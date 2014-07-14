#include <iostream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include "lodepng.h"
#include "bitarray.h"
#include "compression.h"

using namespace std;

//#ifdef M_PI
//#define PI M_PI
//#else
#define PI 3.14159265359
//#endif
#define TWO_PI (2*PI)
#define PI_2 (PI/2)

#define BLOCK_SIZE (16)
#define BLOCK_SIZE2 (BLOCK_SIZE*BLOCK_SIZE)

struct block_t
{
   double luma[BLOCK_SIZE2];
   double chroma_a[BLOCK_SIZE2];
   double chroma_b[BLOCK_SIZE2];

   double luma_res[BLOCK_SIZE2];

   uint8_t luma_pred_mode;
   uint8_t luma_dc[5]; // 8-bit

   uint8_t luma_diff_scale;  // 8-bit
   uint8_t luma_diff_offset; // 8-bit
   uint8_t bit_allocation; // 3-bit

   uint8_t luma_diff[BLOCK_SIZE2]; // 2-8 bit each
   uint8_t luma_comp[BLOCK_SIZE2];
   uint8_t luma_bits[BLOCK_SIZE2*2]; // variable length

   uint16_t chroma_a_c0; // 10-bit
   uint16_t chroma_a_c1; // 10-bit
   uint16_t chroma_b_c0; // 10-bit
   uint16_t chroma_b_c1; // 10-bit

   int   luma_bits_used;
   int   chroma_bits_used;
};

struct image_t {
   unsigned width, height, size;
   block_t* blocks;

   huffman_tree_t* luma_1bit_huffman;
   huffman_tree_t* luma_2bit_huffman;
   huffman_tree_t* luma_3bit_huffman;
   huffman_tree_t* luma_4bit_huffman;
};

float random(int x, int seed)
{
   int n = x * 1619 + seed * 6971;
   n = (n>>8)^n;
   return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float random(int x, int y, int seed)
{
   int n = x * 1619 + y * 31337 + seed * 6971;
   n = (n>>8)^n;
   return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

uint8_t toInt(double v)
{
   return (uint8_t)max(min((int)round(v * 255.0), 255), 0);
}

uint16_t toInt10(double v)
{
   return (uint16_t)max(min((int)round(v * 1023.0), 1023), 0);
}

double fromInt(uint8_t i)
{
   return i / 255.0;
} 

double fromInt10(uint16_t i)
{
   return min((int)i, 1023) / 1023.0;
}

double fromInt12(uint16_t i)
{
   return min((int)i, 4095) / 4095.0;
}

void predict_luma(int mode, block_t& upleft_block, block_t& left_block, block_t& up_block, double* out_luma)
{
   switch (mode)
   {
      case 2:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
               out_luma[x+y*BLOCK_SIZE] = left_block.luma_res[BLOCK_SIZE-1+y*BLOCK_SIZE];
         break;
      case 3:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
               out_luma[x+y*BLOCK_SIZE] = up_block.luma_res[x+(BLOCK_SIZE-1)*BLOCK_SIZE];
         break;
      case 4:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
            {
               int i = x+y*BLOCK_SIZE;
               if (x > y)
                  out_luma[i] = up_block.luma_res[(x-y-1)+(BLOCK_SIZE-1)*BLOCK_SIZE];
               else if (x < y)
                  out_luma[i] = left_block.luma_res[BLOCK_SIZE-1+(y-x-1)*BLOCK_SIZE];
               else
                  out_luma[i] = (up_block.luma_res[(BLOCK_SIZE-1)*BLOCK_SIZE] + left_block.luma_res[BLOCK_SIZE-1]) * 0.5;
            }
         break;
      case 5:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
            {
               int i = x+y*BLOCK_SIZE;
               int j = min(x + y, BLOCK_SIZE-1);
               double up = up_block.luma_res[j+(BLOCK_SIZE-1)*BLOCK_SIZE];
               double left = left_block.luma_res[BLOCK_SIZE-1+j*BLOCK_SIZE];
               out_luma[i] = (left * (x+1.0) + up * (y+1.0)) / (x+y+2.0);
            }
         break;
      default:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
               out_luma[x+y*BLOCK_SIZE] = (x&1) == (y&1) ? 1.0 : 0.0;
   }
}

void decode_dc_luma(int type, double* params, double* out_luma)
{
   if (type == 0)
   {
      for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         out_luma[i] = params[0];
   }
   else
   {
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            unsigned index = x+y*BLOCK_SIZE;
            double luma_p = params[0];
            luma_p += params[1] * x;
            luma_p += params[2] * y;
            out_luma[index] = luma_p;
         }
   }
}

void encode_decode(unsigned width, unsigned height, uint8_t* data)
{
   unsigned size = width * height;
   double*  pixels = new double[size*4];

   unsigned blocks_width = width / BLOCK_SIZE;
   unsigned blocks_height = height / BLOCK_SIZE;
   unsigned blocks_size = blocks_width * blocks_height;

   image_t image;
   image.width = width;
   image.height = height;
   image.blocks = new block_t[blocks_size];

   const double THRESHOLD = 0.04045;
   const double THRESHOLD2 = 0.0031308;
   const double OFFSET = 0.055;
   const double FACTOR1 = 1.055;
   const double FACTOR2 = 12.92;
   const double GAMMA = 2.4;
   const double INV_GAMMA = 1.0 / 2.4;

   for (unsigned i = 0; i < size; i++)
   {
      double x, y, z, alpha;

      {
         double r = data[i*4+0] / 255.0;
         double g = data[i*4+1] / 255.0;
         double b = data[i*4+2] / 255.0;
         double a = data[i*4+3] / 255.0;

         // Convert to linear RGB
         r = (r > THRESHOLD) ? pow((r + OFFSET) / FACTOR1, GAMMA) : r / FACTOR2;
         g = (g > THRESHOLD) ? pow((g + OFFSET) / FACTOR1, GAMMA) : g / FACTOR2;
         b = (b > THRESHOLD) ? pow((b + OFFSET) / FACTOR1, GAMMA) : b / FACTOR2;

         // Convert to XYZ
         x = r * 0.4124 + g * 0.3576 + b * 0.1805;
         y = r * 0.2126 + g * 0.7152 + b * 0.0722;
         z = r * 0.0193 + g * 0.1192 + b * 0.9505;
         alpha = a;
      }

      double xr = x / 0.95047;
      double yr = y / 1.00000;
      double zr = z / 1.08883;

      xr = xr > 0.008856 ? pow(xr, 1.0 / 3.0) : 7.787 * xr + 16.0 / 116.0;
      yr = yr > 0.008856 ? pow(yr, 1.0 / 3.0) : 7.787 * yr + 16.0 / 116.0;
      zr = zr > 0.008856 ? pow(zr, 1.0 / 3.0) : 7.787 * zr + 16.0 / 116.0;

      double L = 1.16 * yr - 0.16;
      double a = 5.00 * (xr - yr);
      double b = 2.00 * (yr - zr);

      pixels[i*4+0] = L;
      pixels[i*4+1] = a;
      pixels[i*4+2] = b;
      pixels[i*4+3] = alpha;
   }



   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // ENCODE //
   ////////////
   int luma_bits_used = 0, chroma_bits_used = 0;
   {
      // Split into blocks
      for (unsigned my = 0; my < blocks_height; my++)
         for (unsigned mx = 0; mx < blocks_width; mx++)
         {
            block_t& block = image.blocks[mx+my*blocks_width];
            for (unsigned by = 0; by < BLOCK_SIZE; by++)
               for (unsigned bx = 0; bx < BLOCK_SIZE; bx++)
               {
                  unsigned x = mx*BLOCK_SIZE + bx, y = my*BLOCK_SIZE + by;
                  unsigned index = x+y*width;
                  unsigned block_index = bx+by*BLOCK_SIZE;
                  block.luma[block_index] = pixels[index*4+0];
                  block.chroma_a[block_index] = pixels[index*4+1];
                  block.chroma_b[block_index] = pixels[index*4+2];
               }
         }

      gsl_vector* luma_vec = gsl_vector_calloc(BLOCK_SIZE2);
      gsl_vector* luma_weight_vec = gsl_vector_calloc(BLOCK_SIZE2);

      gsl_multifit_linear_workspace* fitting_workspace3 = gsl_multifit_linear_alloc(BLOCK_SIZE2, 3);
      gsl_matrix* pos_matrix3 = gsl_matrix_calloc(BLOCK_SIZE2, 3);
      gsl_vector *c3 = gsl_vector_calloc(3);
      gsl_matrix *cov3 = gsl_matrix_calloc(3, 3);

      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            double weight = (x < 2 || y < 2 || x >= BLOCK_SIZE-2 || y >= BLOCK_SIZE-2) ? 5.0 : 0.2;
            gsl_vector_set(luma_weight_vec, x+y*BLOCK_SIZE, weight);
         }

      for (unsigned by = 0; by < blocks_height; by++)
         for (unsigned bx = 0; bx < blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*blocks_width];

            block.chroma_bits_used = 0;
            block.luma_bits_used = 0;

            // Determine roughness
            double vert_diffs = 0.0, hor_diffs = 0.0;
            for (unsigned y = 1; y < BLOCK_SIZE; y++)
               for (unsigned x = 1; x < BLOCK_SIZE; x++)
               {
                  unsigned index = x+y*BLOCK_SIZE;
                  unsigned l_index = (x-1)+y*BLOCK_SIZE;
                  unsigned u_index = x+(y-1)*BLOCK_SIZE;
                  double luma = block.luma[index];
                  double l_luma = block.luma[l_index];
                  double u_luma = block.luma[u_index];
                  hor_diffs += min(abs(luma - l_luma), 0.2);
                  vert_diffs += min(abs(luma - u_luma), 0.2);
               }
            double roughness = max(vert_diffs, hor_diffs) / (BLOCK_SIZE-1)*(BLOCK_SIZE-1);

            // Fit plane to luma
            unsigned mode_bits_used[6];
            unsigned mode_diff_bits_used[6] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
            double   mode_luma_dc_chisq[6] = {1e9, 1e9, 1e9, 1e9, 1e9, 1e9};
            double   mode_luma_dc[5*6];
            for (unsigned mode = 0; mode < 6; mode++)
            {
               if ((bx == 0 || by == 0) && mode > 1) continue;

               double pred_luma[BLOCK_SIZE2] = {0.0};
               
               if (mode > 1)
               {
                  int upleft_index = bx-1+(by-1)*blocks_width;
                  int up_index = bx+(by-1)*blocks_width;
                  int left_index = bx-1+by*blocks_width;
                  assert (upleft_index >= 0 && up_index >= 0 && left_index >= 0);
                  predict_luma(mode, image.blocks[upleft_index], image.blocks[left_index], image.blocks[up_index], pred_luma);
               }

               if (mode == 0)
               {
                  double total = 0.0;
                  for (unsigned i = 0; i < BLOCK_SIZE2; i++) total += block.luma[i];
                  mode_luma_dc[mode*5] = total / BLOCK_SIZE2;
                  mode_bits_used[mode] = 8;
               }
               else
               {
                  for (unsigned y = 0; y < BLOCK_SIZE; y++)
                     for (unsigned x = 0; x < BLOCK_SIZE; x++)
                     {
                        unsigned index = x+y*BLOCK_SIZE;
                        double luma = block.luma[index] - pred_luma[index];
                        gsl_matrix_set(pos_matrix3, index, 0, 1.0);
                        gsl_matrix_set(pos_matrix3, index, 1, x);
                        gsl_matrix_set(pos_matrix3, index, 2, y);
                        gsl_vector_set(luma_vec, index, luma);
                     }

                  gsl_multifit_wlinear(pos_matrix3, luma_weight_vec, luma_vec, c3, cov3, &mode_luma_dc_chisq[mode], fitting_workspace3);
                  for (int i = 0; i < 3; i++) mode_luma_dc[i + mode*5] = gsl_vector_get(c3, i);
                  mode_bits_used[mode] = 3 * 8;
               }

               // Determine luma residual range
               {
                  double luma_dc_fit[BLOCK_SIZE2] = {0.0};
                  decode_dc_luma(mode == 0 ? 0 : 1, &mode_luma_dc[mode*5], luma_dc_fit);

                  double max_luma = -1e6, min_luma = 0; // Min Luma starts with zero to force it to at least zero
                  for (unsigned y = 0; y < BLOCK_SIZE; y++)
                     for (unsigned x = 0; x < BLOCK_SIZE; x++)
                     {
                        unsigned index = x+y*BLOCK_SIZE;
                        double luma_diff = block.luma[index] - max(min(pred_luma[index] + luma_dc_fit[index], 1.0), 0.0);
                        if (max_luma < luma_diff) max_luma = luma_diff;
                        if (min_luma > luma_diff) min_luma = luma_diff;
                     }

                  unsigned q_offset = min(max((int)ceil(-min_luma * 255.0), 0), 255);
                  double   luma_offset = q_offset / 255.0;
                  unsigned q_scale = min(max((int)ceil((max_luma + luma_offset) * 127.0), 0), 255);

                  unsigned luma_bits = 1;
                  while ((1 << luma_bits) < q_scale) luma_bits++;
                  mode_diff_bits_used[mode] = luma_bits;
               }
            }
            
            double luma_dc[5] = {0};
            unsigned luma_dc_score = mode_diff_bits_used[0] * BLOCK_SIZE2 / 2 + mode_bits_used[0];
            for (unsigned i = 0; i < 5; i++) luma_dc[i] = mode_luma_dc[i];

            block.luma_pred_mode = 0;
            for (unsigned i = 1; i < 6; i++)
            {
               if (mode_diff_bits_used[i] >= 0xFF) continue;

               unsigned score = mode_diff_bits_used[i] * BLOCK_SIZE2 / 2 + mode_bits_used[i];
               if (luma_dc_score > score)
               {
                  luma_dc_score = score;
                  for (unsigned j = 0; j < 5; j++) luma_dc[j] = mode_luma_dc[j + i*5];
                  block.luma_pred_mode = i;
               }
            }

            //cout << (int)block.luma_pred_mode << endl;
            block.luma_bits_used += 2 + mode_bits_used[block.luma_pred_mode];

            // Quantize and decode again for reference
            block.luma_dc[0] = toInt(luma_dc[0] * 0.5 + 0.5);
            for (unsigned i = 1; i < 3; i++)
               block.luma_dc[i] = toInt(luma_dc[i] * BLOCK_SIZE + 0.5);

            luma_dc[0] = (fromInt(block.luma_dc[0]) - 0.5) * 2.0;
            for (unsigned i = 1; i < 3; i++)
               luma_dc[i] = (fromInt(block.luma_dc[i]) - 0.5) / BLOCK_SIZE;
         
            block.luma_bits_used += 3 * 8;

            // Generate fit values
            double luma_dc_res[BLOCK_SIZE2] = {0.0};
            {
               unsigned pred_type = block.luma_pred_mode;

               double luma_dc_pred[BLOCK_SIZE2] = {0.0};
               double luma_dc_fit[BLOCK_SIZE2] = {0.0};
               if (pred_type > 1)
                  predict_luma(pred_type, image.blocks[(bx-1)+(by-1)*blocks_width], image.blocks[(bx-1)+by*blocks_width], image.blocks[bx+(by-1)*blocks_width], luma_dc_pred);
               decode_dc_luma(pred_type == 0 ? 0 : 1, luma_dc, luma_dc_fit);

               for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                  block.luma_res[i] = luma_dc_res[i] = max(min(luma_dc_pred[i] + luma_dc_fit[i], 1.0), 0.0);
            }

            // Calculate difference range
            double max_luma = -1e6, min_luma = 0; // Min Luma starts with zero to force it to at least zero
            for (unsigned y = 0; y < BLOCK_SIZE; y++)
               for (unsigned x = 0; x < BLOCK_SIZE; x++)
               {
                  unsigned index = x+y*BLOCK_SIZE;
                  double luma_diff = block.luma[index] - luma_dc_res[index];
                  if (max_luma < luma_diff) max_luma = luma_diff;
                  if (min_luma > luma_diff) min_luma = luma_diff;
               }

            // Lookup threshold per quality level
            unsigned diff_threshold = 2;
            unsigned bit_alloc_offset = 2;
            unsigned lowest_bit_alloc = 1;
            unsigned most_bit_alloc = 4;

            // Calculate luma difference scale and offset
            block.luma_diff_offset = min(max((int)ceil(-min_luma * 255.0), 0), 255);
            double luma_offset = block.luma_diff_offset / 255.0;
            block.luma_diff_scale = min(max((int)ceil((max_luma + luma_offset) * 127.0), 0), 255);
            double luma_scale = block.luma_diff_scale / 127.0;

            if (block.luma_diff_scale < diff_threshold)
            {
               block.luma_diff_scale = 0;
               luma_scale = 0.0;
               block.luma_bits_used += 8;
            }
            else 
               block.luma_bits_used += 16;

            unsigned luma_scale_bits = 1;
            while ((1 << luma_scale_bits) < block.luma_diff_scale) luma_scale_bits++;

            const double inv_luma_scale = 127.0 / block.luma_diff_scale;

            uint16_t factor = 0;
            if (block.luma_diff_scale > 0)
            {
               // Determine luma difference bit allocation
               unsigned bit_allocation = max(luma_scale_bits-bit_alloc_offset, 1U);//min(max((unsigned)luma_scale_bits, 2U), 6U);
               //if (roughness > 15) bit_allocation--;
               if (roughness > 20) bit_allocation--;
               if (bit_allocation == 1 && roughness < 5 && luma_scale_bits >= 3) bit_allocation = 2;
               bit_allocation = min(max((unsigned)bit_allocation, lowest_bit_alloc), most_bit_alloc);

               block.bit_allocation = bit_allocation;
               factor = (1<<block.bit_allocation)-1;
               block.luma_bits_used += 4;

               // Calculate luma difference
               double above_error[BLOCK_SIZE];
               for (unsigned i = 0; i < BLOCK_SIZE; i++) above_error[i] = 0;

               for (unsigned y = 0; y < BLOCK_SIZE; y++)
               {
                  double error = 0.0;
                  for (unsigned x = 0; x < BLOCK_SIZE; x++)
                  {
                     unsigned index = x+y*BLOCK_SIZE;
                     double luma_diff = (block.luma[index] - luma_dc_res[index] + luma_offset) * inv_luma_scale + (error + above_error[x]) * 0.5;
                     int luma_q = (int)round(luma_diff * factor);

                     block.luma_diff[index] = luma_q;

                     // Decode again for reference as well as for error calculation
                     double decode_luma_diff = luma_q / (double)factor;
                     error = luma_diff - decode_luma_diff;
                     above_error[x] = error;
                     block.luma_res[index] = luma_dc_res[index] + decode_luma_diff * luma_scale - luma_offset;
                  }
               }
            }
            else
               for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                  block.luma_res[i] = luma_dc_res[i];

            // Compress luma diffs
            /*uint8_t l = block.luma_comp[0] = block.luma_diff[0];
            for (unsigned i = 1; i < BLOCK_SIZE2; i++)
            {
               block.luma_comp[i] = l ^ block.luma_diff[i];
               l = block.luma_diff[i];
            }*/

            // Determine chroma prediction coefficients
            bool chroma_a_same = true, chroma_b_same = true;
            double chroma_a_c0, chroma_a_c1, chroma_b_c0, chroma_b_c1;
            double cov00, cov01, cov11, sumsq;

            double first_value = block.chroma_a[0];
            for (unsigned i = 1; i < BLOCK_SIZE2; i++)
               if (abs(first_value - block.chroma_a[i]) > 0.001)
               {
                  chroma_a_same = false;
                  break;
               }

            if (chroma_a_same)
            {
               chroma_a_c0 = first_value;
               chroma_a_c1 = 0.0;
            }
            else
            {
               gsl_fit_linear(block.luma, 1, block.chroma_a, 1, BLOCK_SIZE2, 
                  &chroma_a_c0, &chroma_a_c1, 
                  &cov00, &cov01, &cov11, &sumsq);
            }

            first_value = block.chroma_b[0];
            for (unsigned i = 1; i < BLOCK_SIZE2; i++)
               if (abs(first_value - block.chroma_b[i]) > 0.001)
               {
                  chroma_b_same = false;
                  break;
               }

            if (chroma_b_same)
            {
               chroma_b_c0 = first_value;
               chroma_b_c1 = 0.0;
            }
            else
            {
               gsl_fit_linear(block.luma, 1, block.chroma_b, 1, BLOCK_SIZE2, 
                  &chroma_b_c0, &chroma_b_c1, 
                  &cov00, &cov01, &cov11, &sumsq);
            }

            block.chroma_a_c0 = toInt10(chroma_a_c0 * 0.125 + 0.5);
            block.chroma_a_c1 = toInt10(chroma_a_c1 * 0.0625 + 0.5);
            block.chroma_b_c0 = toInt10(chroma_b_c0 * 0.125 + 0.5);
            block.chroma_b_c1 = toInt10(chroma_b_c1 * 0.0625 + 0.5);
            block.chroma_bits_used += 2 * (10 * 2);
         }

      /// Initial transformation and stats collection phase
      unsigned stats_luma_1bit[256] = {0};
      unsigned stats_luma_2bit[4*4*4*4] = {0};
      unsigned stats_luma_3bit[8*8] = {0};
      unsigned stats_luma_4bit[16*16] = {0};

      for (unsigned by = 0; by < blocks_height; by++)
         for (unsigned bx = 0; bx < blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*blocks_width];

            if (block.luma_diff_scale > 0)
            {
               mtf_encode8(block.luma_diff, block.luma_comp, BLOCK_SIZE2);

               if (block.bit_allocation == 1)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 8)
                  {
                     uint8_t index = (block.luma_comp[i] << 7) + (block.luma_comp[i+1] << 6) + (block.luma_comp[i+2] << 5) + (block.luma_comp[i+3] << 4);
                     index += (block.luma_comp[i+4] << 3) + (block.luma_comp[i+5] << 2) + (block.luma_comp[i+6] << 1) + block.luma_comp[i+7];
                     stats_luma_1bit[index]++;
                  }
               }
               else if (block.bit_allocation == 2)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 4)
                     stats_luma_2bit[(block.luma_comp[i] << 6) + (block.luma_comp[i+1] << 4) + (block.luma_comp[i+2] << 2) + block.luma_comp[i+3]]++;
               }
               else if (block.bit_allocation == 3)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
                     stats_luma_3bit[(block.luma_comp[i] << 3) + block.luma_comp[i+1]]++;
               }
               else if (block.bit_allocation == 4)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
                     stats_luma_4bit[(block.luma_comp[i] << 4) + block.luma_comp[i+1]]++;
               }
            }
         }

      /// Huffman tree construction
      {
         unsigned total = 0;
         for (unsigned i = 0; i < 256; i++) total += stats_luma_1bit[i];
         if (total > 0)
         {
            printf("1-bit stats:\n");
            //for (unsigned i = 0; i < 256; i++) printf("%2i: %i\n", i, stats_luma_1bit[i]);
            //printf("\n");

            image.luma_1bit_huffman = build_huffman_tree(stats_luma_1bit, 256);

            unsigned total_bits = 0;
            for (unsigned i = 0; i < 256; i++)
            {
               uint32_t code = 0;
               uint8_t length = get_huffman_code(image.luma_1bit_huffman, i, &code);
               //printf("%2i: ", i);
               //for (int j = length - 1; j >= 0; j--) printf("%i", (code >> j) & 1);
               //printf("\n");
               total_bits += length * stats_luma_1bit[i];
            }
            //printf("\n");

            double bits_per_symbol = (double)total_bits / total;
            printf("bits per symbol: %f\n", bits_per_symbol/8.0);
            printf("compression: %f%%\n\n", bits_per_symbol/0.08);
         }
      }

      {
         unsigned total = 0;
         for (unsigned i = 0; i < 4*4*4*4; i++) total += stats_luma_2bit[i];
         if (total > 0)
         {
            printf("2-bit stats:\n");
            //for (unsigned i = 0; i < 4*4; i++) printf("%2i: %i\n", i, stats_luma_2bit[i]);
            //printf("\n");

            image.luma_2bit_huffman = build_huffman_tree(stats_luma_2bit, 4*4*4*4);

            unsigned total_bits = 0;
            for (unsigned i = 0; i < 4*4*4*4; i++)
            {
               uint32_t code = 0;
               uint8_t length = get_huffman_code(image.luma_2bit_huffman, i, &code);
               //printf("%2i: ", i);
               //for (int j = length - 1; j >= 0; j--) printf("%i", (code >> j) & 1);
               //printf("\n");
               total_bits += length * stats_luma_2bit[i];
            }
            //printf("\n");

            double bits_per_symbol = (double)total_bits / total;
            printf("bits per symbol: %f\n", bits_per_symbol/4.0);
            printf("compression: %f%%\n\n", bits_per_symbol/0.08);
         }
      }

      {
         unsigned total = 0;
         for (unsigned i = 0; i < 8*8; i++) total += stats_luma_3bit[i];

         if (total > 0)
         {
            printf("3-bit stats:\n");
            //for (unsigned i = 0; i < 8*8; i++) printf("%2i: %i\n", i, stats_luma_3bit[i]);
            //printf("\n");

            image.luma_3bit_huffman = build_huffman_tree(stats_luma_3bit, 8*8);

            unsigned total_bits = 0;
            for (unsigned i = 0; i < 8*8; i++)
            {
               uint32_t code = 0;
               uint8_t length = get_huffman_code(image.luma_3bit_huffman, i, &code);
               //printf("%2i: ", i);
               //for (int j = length - 1; j >= 0; j--) printf("%i", (code >> j) & 1);
               //printf("\n");
               total_bits += length * stats_luma_3bit[i];
            }
            //printf("\n");

            double bits_per_symbol = (double)total_bits / total;
            printf("bits per symbol: %f\n", bits_per_symbol/2.0);
            printf("compression: %f%%\n\n", bits_per_symbol/0.06);
         }
      }

      {
         unsigned total = 0;
         for (unsigned i = 0; i < 16*16; i++) total += stats_luma_4bit[i];

         if (total > 0)
         {
            printf("4-bit stats:\n");
            //for (unsigned i = 0; i < 16*16; i++) printf("%2i: %i\n", i, stats_luma_4bit[i]);
            //printf("\n");

            image.luma_4bit_huffman = build_huffman_tree(stats_luma_4bit, 16*16);

            unsigned total_bits = 0;
            for (unsigned i = 0; i < 16*16; i++)
            {
               uint32_t code = 0;
               uint8_t length = get_huffman_code(image.luma_4bit_huffman, i, &code);
               //printf("%2i: ", i);
               //for (int j = length - 1; j >= 0; j--) printf("%i", (code >> j) & 1);
               //printf("\n");
               total_bits += length * stats_luma_4bit[i];
            }
            //printf("\n");

            double bits_per_symbol = (double)total_bits / total;
            printf("bits per symbol: %f\n", bits_per_symbol/2.0);
            printf("compression: %f%%\n\n", bits_per_symbol/0.08);
         }
      }

      /// Block compression and writing
      for (unsigned by = 0; by < blocks_height; by++)
         for (unsigned bx = 0; bx < blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*blocks_width];

            if (block.luma_diff_scale > 0)
            {
               /*
               printf("Off: %.3f - Scale: %.3f - Max Diff: %i - Min Luma: %.3f - Max Luma: %.3f\n", luma_offset, luma_scale, factor, min_luma, max_luma);
               for (unsigned y = 0; y < BLOCK_SIZE; y++)
               {
                  for (unsigned x = 0; x < BLOCK_SIZE; x++)
                     printf("%3i ", block.luma_diff[x+y*BLOCK_SIZE]);
                  printf("\n");
               }
               printf("\n");
               */

               uint8_t mtf_encoded[BLOCK_SIZE2];
               uint8_t temp_encoded[BLOCK_SIZE2*2];

               for (unsigned i = 0; i < BLOCK_SIZE2; i++) mtf_encoded[i] = block.luma_comp[i];
               mtf_encode8(block.luma_diff, mtf_encoded, BLOCK_SIZE2);

               bit_array_writer_t temp_luma_writer(temp_encoded, BLOCK_SIZE2*2);

               bit_array_writer_t luma_writer(block.luma_bits, BLOCK_SIZE2*2);

               if (block.bit_allocation == 1)
               {
                  unsigned orig_size = BLOCK_SIZE2;
                  unsigned huffman_size = 0;
                  unsigned runlength_size = runlength_bit_encode(mtf_encoded, BLOCK_SIZE2, temp_encoded, BLOCK_SIZE2*2, 4);

                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 8)
                  {
                     uint8_t symbol = (mtf_encoded[i] << 7) + (mtf_encoded[i+1] << 6) + (mtf_encoded[i+2] << 5) + (mtf_encoded[i+3] << 4);
                     symbol += (mtf_encoded[i+4] << 3) + (mtf_encoded[i+5] << 2) + (mtf_encoded[i+6] << 1) + mtf_encoded[i+7];
                     write_huffman_symbol(image.luma_1bit_huffman, temp_luma_writer, symbol);
                  }
                  temp_luma_writer.pad_to_byte();
                  huffman_size = temp_luma_writer.get_size();

                  if (runlength_size < huffman_size)
                  {
                     luma_writer.write_bit(1);
                     runlength_bit_encode(mtf_encoded, BLOCK_SIZE2, luma_writer, 4);
                  }
                  else
                  {
                     luma_writer.write_bit(0);
                     for (unsigned i = 0; i < BLOCK_SIZE2; i += 8)
                     {
                        uint8_t symbol = (mtf_encoded[i] << 7) + (mtf_encoded[i+1] << 6) + (mtf_encoded[i+2] << 5) + (mtf_encoded[i+3] << 4);
                        symbol += (mtf_encoded[i+4] << 3) + (mtf_encoded[i+5] << 2) + (mtf_encoded[i+6] << 1) + mtf_encoded[i+7];
                        write_huffman_symbol(image.luma_1bit_huffman, luma_writer, symbol);
                     }
                  }
               }
               else if (block.bit_allocation == 2)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 4)
                  {
                     uint8_t symbol = (mtf_encoded[i] << 6) + (mtf_encoded[i+1] << 4) + (mtf_encoded[i+2] << 2) + mtf_encoded[i+3];
                     write_huffman_symbol(image.luma_2bit_huffman, luma_writer, symbol);
                  }
               }
               else if (block.bit_allocation == 3)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
                  {
                     uint8_t symbol = (mtf_encoded[i] << 3) + mtf_encoded[i+1];
                     write_huffman_symbol(image.luma_3bit_huffman, luma_writer, symbol);
                  }
               }
               else if (block.bit_allocation == 4)
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
                  {
                     uint8_t symbol = (mtf_encoded[i] << 4) + mtf_encoded[i+1];
                     write_huffman_symbol(image.luma_4bit_huffman, luma_writer, symbol);
                  }
               }
               else
               {
                  for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                     luma_writer.write_bits8(mtf_encoded[i], block.bit_allocation);
               }
               luma_writer.pad_to_byte();

               block.luma_bits_used += luma_writer.get_size() * 8;

               if (bx == 16 && by == 9)
               {
                  for (unsigned y = 0; y < BLOCK_SIZE; y++)
                  {
                     for (unsigned x = 0; x < BLOCK_SIZE; x++)
                        printf("%3i ", mtf_encoded[x+y*BLOCK_SIZE]);
                     printf("\n");
                  }
                  printf("\n");

                  printf("Original bytes  : %i\n", block.bit_allocation*BLOCK_SIZE2/8);
                  printf("Compressed bytes: %i\n", luma_writer.get_size());
               }
            }

            luma_bits_used += block.luma_bits_used;
            chroma_bits_used += block.chroma_bits_used;
         }

      gsl_multifit_linear_free(fitting_workspace3);
   }

   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // DECODE //
   ////////////
   {

      for (unsigned by = 0; by < blocks_height; by++)
         for (unsigned bx = 0; bx < blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*blocks_width];

            // Decode luma dc
            int pred_mode = block.luma_pred_mode;

            double luma_dc[5];
            luma_dc[0] = (fromInt(block.luma_dc[0]) - 0.5) * 2.0;
            for (int i = 1; i < 3; i++)
               luma_dc[i] = (fromInt(block.luma_dc[i]) - 0.5) / BLOCK_SIZE;

            double luma_dc_pred[BLOCK_SIZE2] = {0.0};
            double luma_dc_fit[BLOCK_SIZE2] = {0.0};
            double luma_dc_res[BLOCK_SIZE2] = {0.0};

            for (unsigned y = 0; y < BLOCK_SIZE; y++)
               for (unsigned x = 0; x < BLOCK_SIZE; x++)
                  luma_dc_pred[x+y*BLOCK_SIZE] = 0.0;

            if (pred_mode > 1)
               predict_luma(pred_mode, image.blocks[(bx-1)+(by-1)*blocks_width], image.blocks[(bx-1)+by*blocks_width], image.blocks[bx+(by-1)*blocks_width], luma_dc_pred);

            decode_dc_luma(pred_mode == 0 ? 0 : 1, luma_dc, luma_dc_fit);

            for (unsigned i = 0; i < BLOCK_SIZE2; i++)
               luma_dc_res[i] = max(min(luma_dc_pred[i] + luma_dc_fit[i], 1.0), 0.0);

            // Decode luma difference
            if (block.luma_diff_scale > 0)
            {
               double luma_offset = block.luma_diff_offset / 255.0;
               double luma_scale = block.luma_diff_scale / 127.0;
               uint16_t factor = (1<<block.bit_allocation)-1;

               /*
               for (unsigned y = 0; y < BLOCK_SIZE; y+=2)
                  for (unsigned x = 0; x < BLOCK_SIZE; x+=2)
                  {
                     unsigned index = x+y*BLOCK_SIZE;
                     double v = luma_qs[index];
                     double t = v;
                     bool a1 = abs(luma_qs[index+1]-v) < 1.1;
                     bool a2 = abs(luma_qs[index+BLOCK_SIZE]-v) < 1.1;
                     bool a3 = abs(luma_qs[index+1+BLOCK_SIZE]-v) < 1.1;
                     t += a1 ? luma_qs[index+1] : v;
                     t += a2 ? luma_qs[index+BLOCK_SIZE] : v;
                     t += a3 ? luma_qs[index+1+BLOCK_SIZE] : v;
                     t *= 0.25;
                     luma_qs[index] = t;
                     if (a1) luma_qs[index+1] = t;
                     if (a2) luma_qs[index+BLOCK_SIZE] = t;
                     if (a3) luma_qs[index+1+BLOCK_SIZE] = t;
                  }
               */

               bit_array_reader_t luma_reader(block.luma_bits, BLOCK_SIZE2*2);

               if (block.bit_allocation == 1)
               {
                  uint8_t b = luma_reader.read_bit();

                  if (b == 1)
                     runlength_bit_decode(luma_reader, block.luma_comp, BLOCK_SIZE2, 4);
                  else
                  {
                     for (int i = 0; i < BLOCK_SIZE2; i += 8)
                     {
                        uint8_t symbol = read_huffman_symbol(image.luma_1bit_huffman, luma_reader);
                        block.luma_comp[i  ] = (symbol >> 7) & 1;
                        block.luma_comp[i+1] = (symbol >> 6) & 1;
                        block.luma_comp[i+2] = (symbol >> 5) & 1;
                        block.luma_comp[i+3] = (symbol >> 4) & 1;
                        block.luma_comp[i+4] = (symbol >> 3) & 1;
                        block.luma_comp[i+5] = (symbol >> 2) & 1;
                        block.luma_comp[i+6] = (symbol >> 1) & 1;
                        block.luma_comp[i+7] = symbol & 1;
                     }
                  }
               }
               else if (block.bit_allocation == 2)
               {
                  for (int i = 0; i < BLOCK_SIZE2; i += 4)
                  {
                     uint8_t symbol = read_huffman_symbol(image.luma_2bit_huffman, luma_reader);
                     block.luma_comp[i] = (symbol >> 6) & 3;
                     block.luma_comp[i+1] = (symbol >> 4) & 3;
                     block.luma_comp[i+2] = (symbol >> 2) & 3;
                     block.luma_comp[i+3] = symbol & 3;
                  }
               }
               else if (block.bit_allocation == 3)
               {
                  for (int i = 0; i < BLOCK_SIZE2; i += 2)
                  {
                     uint8_t symbol = read_huffman_symbol(image.luma_3bit_huffman, luma_reader);
                     block.luma_comp[i] = (symbol >> 3) & 7;
                     block.luma_comp[i+1] = symbol & 7;
                  }
               }
               else if (block.bit_allocation == 4)
               {
                  for (int i = 0; i < BLOCK_SIZE2; i += 2)
                  {
                     uint8_t symbol = read_huffman_symbol(image.luma_4bit_huffman, luma_reader);
                     block.luma_comp[i] = (symbol >> 4) & 15;
                     block.luma_comp[i+1] = symbol & 15;
                  }
               }
               else
               {
                  for (int i = 0; i < BLOCK_SIZE2; i++)
                     block.luma_comp[i] = luma_reader.read_bits8(block.bit_allocation);
               }

               mtf_decode8(block.luma_comp, block.luma_diff, BLOCK_SIZE2);

               for (unsigned y = 0; y < BLOCK_SIZE; y++)
                  for (unsigned x = 0; x < BLOCK_SIZE; x++)
                  {
                     unsigned index = x+y*BLOCK_SIZE;

                     double luma_q = block.luma_diff[index];
                     double luma_diff = (luma_q / factor) * luma_scale - luma_offset;
                     //double luma = luma_dc_res[index];
                     double luma = luma_dc_res[index] + luma_diff;

                     block.luma[index] = luma;// + (random(x, y, bx+by*8) * 0.5 - 0.5) / factor * luma_scale;
                     block.luma_res[index] = luma;
                     //block.luma[index] = (block.luma[index] - luma) * 2.0 + 0.5;
                  } 
            }
            else
            {
               for (unsigned y = 0; y < BLOCK_SIZE; y++)
                  for (unsigned x = 0; x < BLOCK_SIZE; x++)
                  {
                     unsigned index = x+y*BLOCK_SIZE;

                     double luma = luma_dc_res[index];

                     block.luma[index] = luma;
                     block.luma_res[index] = luma;
                     //block.luma[index] = (block.luma[index] - luma_dc_res[index]) * 2.0 + 0.5;
                  }

               // deblocking
               if (bx > 0)
                  for (unsigned y = 0; y < BLOCK_SIZE; y++)
                  {
                     unsigned index = y*BLOCK_SIZE;
                     unsigned left_index = bx-1+by*blocks_width;
                     block.luma[index] = block.luma[index] * 0.5 + image.blocks[left_index].luma_res[index+BLOCK_SIZE-1] * 0.5;
                     block.luma[index+1] = block.luma[index+1] * 0.67 + image.blocks[left_index].luma_res[index+BLOCK_SIZE-1] * 0.33;
                  }

               if (by > 0)
                  for (unsigned x = 0; x < BLOCK_SIZE; x++)
                  {
                     unsigned up_index = bx+(by-1)*blocks_width;
                     block.luma[x] = block.luma[x] * 0.5 + image.blocks[up_index].luma_res[x+(BLOCK_SIZE-1)*BLOCK_SIZE] * 0.5;
                     block.luma[x+BLOCK_SIZE] = block.luma[x+BLOCK_SIZE] * 0.67 + image.blocks[up_index].luma_res[x+(BLOCK_SIZE-1)*BLOCK_SIZE] * 0.33;
                  }

               // add a bit of noise
               for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                  block.luma[i] += (random(i, bx+by*blocks_width) - 0.5) / 63.0;
            }

            // Predict chroma
            double chroma_a_c0 = fromInt10(block.chroma_a_c0) * 8.0 - 4.0;
            double chroma_a_c1 = fromInt10(block.chroma_a_c1) * 16.0 - 8.0;
            double chroma_b_c0 = fromInt10(block.chroma_b_c0) * 8.0 - 4.0;
            double chroma_b_c1 = fromInt10(block.chroma_b_c1) * 16.0 - 8.0;

            for (unsigned i = 0; i < BLOCK_SIZE2; i++)
            {
               block.chroma_a[i] = 0.0;
               block.chroma_b[i] = 0.0;
            }

            for (unsigned y = 0; y < BLOCK_SIZE; y++)
               for (unsigned x = 0; x < BLOCK_SIZE; x++)
               {
                  unsigned index = x+y*BLOCK_SIZE;
                  /*
                  double L = 0.25*(block.luma[index]+block.luma[index+1]+block.luma[index+1+BLOCK_SIZE]+block.luma[index+BLOCK_SIZE]);

                  block.chroma_a[index] = a;
                  block.chroma_a[index+1] = a;
                  block.chroma_a[index+1+BLOCK_SIZE] = a;
                  block.chroma_a[index+BLOCK_SIZE] = a;
                  block.chroma_b[index] = b;
                  block.chroma_b[index+1] = b;
                  block.chroma_b[index+1+BLOCK_SIZE] = b;
                  block.chroma_b[index+BLOCK_SIZE] = b;
                  */

                  double L = block.luma[index];
                  double a = max(min(L * chroma_a_c1 + chroma_a_c0, 1.0), -1.0);
                  double b = max(min(L * chroma_b_c1 + chroma_b_c0, 1.0), -1.0);
                  block.chroma_a[index] = a;
                  block.chroma_b[index] = b;
               }
         }

         for (unsigned my = 0; my < blocks_height; my++)
            for (unsigned mx = 0; mx < blocks_width; mx++)
            {
               block_t& block = image.blocks[mx+my*blocks_width];
               for (unsigned by = 0; by < BLOCK_SIZE; by++)
                  for (unsigned bx = 0; bx < BLOCK_SIZE; bx++)
                  {
                     unsigned x = mx*BLOCK_SIZE + bx, y = my*BLOCK_SIZE + by;
                     unsigned index = x+y*width;
                     unsigned block_index = bx+by*BLOCK_SIZE;
                     
                     pixels[index*4+0] = block.luma[block_index];
                     pixels[index*4+1] = block.chroma_a[block_index];
                     pixels[index*4+2] = block.chroma_b[block_index];
                  }
            }
   }



   double r_error = 0.0, g_error = 0.0, b_error = 0.0;

   for (unsigned i = 0; i < size; i++)
   {
      double L = pixels[i*4+0];
      double a = pixels[i*4+1];
      double b = pixels[i*4+2];

      double yr = (L + 0.16) / 1.16;
      double xr = a / 5.00 + yr;
      double zr = yr - b / 2.00;

      double xr3 = xr*xr*xr, yr3 = yr*yr*yr, zr3 = zr*zr*zr;
      xr = xr3 > 0.008856 ? xr3 : (xr - 16.0 / 116.0) / 7.87;
      yr = yr3 > 0.008856 ? yr3 : (yr - 16.0 / 116.0) / 7.87;
      zr = zr3 > 0.008856 ? zr3 : (zr - 16.0 / 116.0) / 7.87;

      double x = xr * 0.95047;
      double y = yr * 1.00000;
      double z = zr * 1.08883;

      {
         double r = x *  3.2406 + y * -1.5372 + z * -0.4986;
         double g = x * -0.9689 + y *  1.8758 + z *  0.0415;
         double b = x *  0.0557 + y * -0.2040 + z *  1.0570;

         // Convert to sRGB
         r = (r > THRESHOLD2) ? FACTOR1 * pow(r, INV_GAMMA) - OFFSET : FACTOR2 * r;
         g = (g > THRESHOLD2) ? FACTOR1 * pow(g, INV_GAMMA) - OFFSET : FACTOR2 * g;
         b = (b > THRESHOLD2) ? FACTOR1 * pow(b, INV_GAMMA) - OFFSET : FACTOR2 * b;

         uint8_t ir = toInt(r);
         uint8_t ig = toInt(g);
         uint8_t ib = toInt(b);

         int er = data[i*4+0]-ir;
         int eg = data[i*4+1]-ig;
         int eb = data[i*4+2]-ib;
         r_error += er*er;
         g_error += eg*eg;
         b_error += eb*eb;

         data[i*4+0] = ir;
         data[i*4+1] = ig;
         data[i*4+2] = ib;
         data[i*4+3] = 255;
      }
   }

   int bits_used = luma_bits_used + chroma_bits_used;
   r_error /= size;
   b_error /= size;
   g_error /= size;
   cout << "Luma Bytes used:" << (luma_bits_used/8) << " B" << endl;
   cout << "Chroma Bytes used:" << (chroma_bits_used/8) << " B" << endl;
   cout << "Bytes used:" << (bits_used/8) << " B" << endl;
   cout << "Luma Bits per pixel:" << ((double)luma_bits_used/width/height) << " bps" << endl;
   cout << "Chroma Bits per pixel:" << ((double)chroma_bits_used/width/height) << " bps" << endl;
   cout << "Bits per pixel:" << ((double)bits_used/width/height) << " bps" << endl;
   cout << "Kilobytes used:" << (bits_used/8/1024) << " KB" << endl;
   cout << "R error: " << r_error << endl;
   cout << "G error: " << g_error << endl;
   cout << "B error: " << b_error << endl;
   cout << "W error: " << (r_error * 0.2126 + g_error * 0.7152 + b_error * 0.0722) << endl;

   delete [] pixels;
   delete [] image.blocks;
}

int main()
{
   string   name = "test/coast";
   string   filename = name + ".png";
   string   result_filename = name + "_out.png";

   unsigned width, height;
   uint8_t* image;


   unsigned error;
   error = lodepng_decode32_file(&image, &width, &height, filename.c_str());
   if (error)
   {
      cout << "Could not load file: " << filename << endl;
      return 1;
   }

   encode_decode(width, height, image);

   error = lodepng_encode32_file(result_filename.c_str(), image, width, height);
   if (error)
   {
      cout << "Could not save file: " << result_filename << endl;
      return 2;
   }

   return 0;
}