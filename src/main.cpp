#include <iostream>
#include <cmath>
#include <random>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#include "lodepng.h"

#include "bitio/bitarray.h"
#include "bitio/bitfile.h"
#include "bitio/bitview.h"
#include "bitio/bitdebug.h"

#include "compression/bitmanip.h"
#include "compression/compression.h"
#include "compression/huffman.h"
#include "compression/runlength.h"
#include "compression/arithmetic_coder.h"
#include "compression/rice.h"
#include "compression/elias_gamma.h"

#include "util.h"

using std::min;
using std::max;
using std::cout;
using std::endl;

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
   double contrast[BLOCK_SIZE2];
   double luma[BLOCK_SIZE2];
   double chroma_a[BLOCK_SIZE2];
   double chroma_b[BLOCK_SIZE2];

   double luma_x;
   double luma_x2;

   double luma_res[BLOCK_SIZE2];

   uint8_t luma_pred_mode;
   uint8_t luma_dc[16];

   uint8_t luma_diff_scale; 
   uint8_t luma_diff_offset; 
   uint8_t bit_allocation;

   uint8_t luma_diff[BLOCK_SIZE2]; 
   uint8_t luma_comp[BLOCK_SIZE2];
   uint8_t luma_bits[BLOCK_SIZE2*2];

   uint8_t luma_compress_mode;
   uint8_t luma_length;

   uint8_t chroma_a_c0[2];
   uint8_t chroma_a_c1[2];
   uint8_t chroma_b_c0[2];
   uint8_t chroma_b_c1[2];
};

struct image_t {
   unsigned width, height, size;
   double*  pixels;
   double*  contrast;

   unsigned blocks_width, blocks_height, blocks_size;
   block_t* blocks;

   canonical_huffman_t* luma_mode_coder;
   canonical_huffman_t* luma_1bit_coder;
   canonical_huffman_t* luma_2bit_coder;
   canonical_huffman_t* luma_3bit_coder;
   canonical_huffman_t* luma_4bit_coder;
   canonical_huffman_t* luma_5bit_coder;
   canonical_huffman_t* luma_6bit_coder;
};

gsl_vector* luma_vec = gsl_vector_calloc(BLOCK_SIZE2);
gsl_vector* luma_weight_vec = gsl_vector_calloc(BLOCK_SIZE2);

gsl_multifit_linear_workspace* fitting_workspace3 = gsl_multifit_linear_alloc(BLOCK_SIZE2, 3);
gsl_matrix* pos_matrix3 = gsl_matrix_calloc(BLOCK_SIZE2, 3);
gsl_vector *c3 = gsl_vector_calloc(3);
gsl_matrix *cov3 = gsl_matrix_calloc(3, 3);

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
      case 6:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
            {
               int i = x+y*BLOCK_SIZE;
               double up = up_block.luma_res[i];
               double left = left_block.luma_res[i];
               out_luma[i] = up + left - up * left;
            }
         break;
      default:
         for (int y = 0; y < BLOCK_SIZE; y++)
            for (int x = 0; x < BLOCK_SIZE; x++)
               out_luma[x+y*BLOCK_SIZE] = (x&1) == (y&1) ? 1.0 : 0.0;
   }
}

void decode_dc_luma(int mode, double* params, double* out_luma)
{
   if (mode == 0)
   {
      for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         out_luma[i] = params[0];
   }
   else if (mode == 7)
   {
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            int ix = BLOCK_SIZE - 1 - x;
            int iy = BLOCK_SIZE - 1 - y;

            double v = params[0] * ix * iy + params[1] * x * iy + params[2] * ix * y + params[3] * x * y;
            out_luma[x+y*BLOCK_SIZE] = v / ((BLOCK_SIZE-1) * (BLOCK_SIZE-1));
         }
   }
   else if (mode == 1)
   {
      const float factor = 2.0 / BLOCK_SIZE;
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            unsigned index = x+y*BLOCK_SIZE;
            double luma_p = params[0];
            luma_p += params[1] * x * factor;
            luma_p += params[2] * y * factor;
            out_luma[index] = luma_p;
         }
   }
   else
   {
      for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         out_luma[i] = params[0];
   }
}

void write_variable_length_number(uint32_t value, bit_writer_t& writer)
{
   do
   {
      writer.write_byte((value > 127 ? 0x80 : 0) + (value & 0x7F));
      value >>= 7;
   }
   while (value != 0);
}

uint32_t read_variable_length_number(bit_reader_t& reader)
{
   unsigned b = 0x80;
   unsigned shift = 0;
   uint32_t value = 0;
   while ((b & 0x80) == 0x80)
   {
      b = reader.read_byte();
      value += (b & 0x7F) << shift;
      shift += 7;
   }
   return value;
}

double calculate_score(unsigned bit_cost, double mse, unsigned quality)
{
   double q = sqr(clamp(quality / 7.0, 0.0, 1.0));
   double score = mse * (1.0 - q) * 100.0 + bit_cost * q * 0.001;
   return score;
}

void guassian_blur_filter(double* input, int width, int height, unsigned input_stride, double* output, unsigned output_stride)
{
   double* buffer = new double[width * height];

   const double guassian_filter[7] = {1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0};
   for (unsigned y = 0; y < height; y++)
      for (unsigned x = 0; x < width; x++)
      {
         double t = 0.0;
         for (int fy = -3; fy <= 3; fy++) 
         {
            int sy = clamp(y + fy, 0, height - 1);
            t += input[(x + sy * width) * input_stride] * guassian_filter[fy + 3];
         }

         buffer[x + y * width] = t;
      }

   for (unsigned y = 0; y < height; y++)
      for (unsigned x = 0; x < width; x++)
      {
         double t = 0.0;
         for (int fx = -3; fx <= 3; fx++) 
         {
            int sx = clamp(x + fx, 0, width - 1);
            t += buffer[sx + y * width] * guassian_filter[fx + 3];
         }

         output[(x + y * width) * output_stride] = t / (64.0 * 64.0);
      }

   delete [] buffer;
}

void max_filter(double* input, int width, int height, unsigned input_stride, double* output, unsigned output_stride)
{
   double* buffer = new double[width * height];

   for (unsigned y = 0; y < height; y++)
      for (unsigned x = 0; x < width; x++)
      {
         double t = 0.0;
         for (int fy = -1; fy <= 1; fy++) 
         {
            int sy = y + fy;
            if (sy >= 0 && sy < height)
               t = max(t, input[(x + sy * width) * input_stride]);
         }

         buffer[x + y * width] = t;
      }

   for (unsigned y = 0; y < height; y++)
      for (unsigned x = 0; x < width; x++)
      {
         double t = 0.0;
         for (int fx = -1; fx <= 1; fx++) 
         {
            int sx = x + fx;
            if (sx >= 0 && sx < width)
               t = max(t, buffer[sx + y * width]);
         }
         
         output[(x + y * width) * output_stride] = t;
      } 

   delete [] buffer;
}

struct mode_stats_t
{
   bool     valid = false;
   unsigned bits_used = 0;
   double   luma_range = 0.0;
   double   parameters[16];
};

struct residual_stats_t
{
   bool     valid = false;
   unsigned bits_used = 0;
   double   error = 1e9;
   double   weighted_error = 1e9;
};

double* generate_mode_data(image_t& image, unsigned bx, unsigned by, mode_stats_t& stats, unsigned mode)
{
   block_t&       block = image.blocks[bx+by*image.blocks_width];

   if ((bx == 0 || by == 0) && mode > 1 && mode < 7) return nullptr;

   // Calculate mode predictions
   double pred_luma[BLOCK_SIZE2] = {0.0};
   
   if (mode > 1 && mode < 7)
   {
      int upleft_index = bx-1+(by-1)*image.blocks_width;
      int up_index = bx+(by-1)*image.blocks_width;
      int left_index = bx-1+by*image.blocks_width;
      assert (upleft_index >= 0 && up_index >= 0 && left_index >= 0);
      predict_luma(mode, image.blocks[upleft_index], image.blocks[left_index], image.blocks[up_index], pred_luma);
   }

   if (mode == 0)
   {
      double total = 0.0;
      for (unsigned i = 0; i < BLOCK_SIZE2; i++) total += block.luma[i];
      stats.parameters[0] = total / BLOCK_SIZE2;
      stats.bits_used = 8;
   }
   else if (mode == 7)
   {
      stats.parameters[0] = block.luma[0];
      stats.parameters[1] = block.luma[BLOCK_SIZE-1];
      stats.parameters[2] = block.luma[BLOCK_SIZE2-BLOCK_SIZE];
      stats.parameters[3] = block.luma[BLOCK_SIZE2-1];
      stats.bits_used = 4 * 6;
   }
   else if (mode == 1)
   {
      double factor = 2.0 / BLOCK_SIZE;
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            unsigned index = x+y*BLOCK_SIZE;
            double luma = block.luma[index] - pred_luma[index];
            gsl_matrix_set(pos_matrix3, index, 0, 1.0);
            gsl_matrix_set(pos_matrix3, index, 1, x * factor - 1.0);
            gsl_matrix_set(pos_matrix3, index, 2, y * factor - 1.0);
            gsl_vector_set(luma_vec, index, luma);
         }

      double chisq;
      gsl_multifit_wlinear(pos_matrix3, luma_weight_vec, luma_vec, c3, cov3, &chisq, fitting_workspace3);
      for (int i = 0; i < 3; i++) stats.parameters[i] = gsl_vector_get(c3, i);
      stats.bits_used = 3 * 8;
   }
   else
   {
      double factor = 2.0 / BLOCK_SIZE;
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            unsigned index = x+y*BLOCK_SIZE;
            double luma = block.luma[index] - pred_luma[index];
            gsl_matrix_set(pos_matrix3, index, 0, 1.0);
            gsl_matrix_set(pos_matrix3, index, 1, 0.0);
            gsl_matrix_set(pos_matrix3, index, 2, 0.0);
            gsl_vector_set(luma_vec, index, luma);
         }

      double chisq;
      gsl_multifit_wlinear(pos_matrix3, luma_weight_vec, luma_vec, c3, cov3, &chisq, fitting_workspace3);
      for (int i = 0; i < 3; i++) stats.parameters[i] = gsl_vector_get(c3, i);
      stats.bits_used = 1 * 8;
   }

   // Determine luma error and residual range
   {
      double luma_dc_fit[BLOCK_SIZE2] = {0.0};
      decode_dc_luma(mode, stats.parameters, luma_dc_fit);

      double max_luma = -1e6, min_luma = 0.0; // Min Luma starts with zero to force it to at least zero
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
   }

   // Quantize and decode again for reference
   if (mode == 0)
   {
      block.luma_dc[0] = toInt(stats.parameters[0]);
      stats.parameters[0] = fromInt(block.luma_dc[0]);
   }
   else if (mode == 7)
   {
      for (unsigned i = 0; i < 4; i++)
         block.luma_dc[i] = toInt(stats.parameters[i]);

      for (unsigned i = 0; i < 4; i++)
         stats.parameters[i] = fromInt(block.luma_dc[i]);
   }
   else
   {
      for (unsigned i = 0; i < 3; i++)
         block.luma_dc[i] = toInt(stats.parameters[i] * 0.5 + 0.5);

      for (unsigned i = 0; i < 3; i++)
         stats.parameters[i] = (fromInt(block.luma_dc[i]) - 0.5) * 2.0;
   }

   // Generate fit values
   double* luma_dc_res = new double[BLOCK_SIZE2];
   {
      for (int i = 0; i < BLOCK_SIZE2; i++) luma_dc_res[i] = 0.0;

      double luma_dc_pred[BLOCK_SIZE2] = {0.0};
      double luma_dc_fit[BLOCK_SIZE2] = {0.0};
      if (mode > 1 && mode < 7)
         predict_luma(mode, image.blocks[(bx-1)+(by-1)*image.blocks_width], image.blocks[(bx-1)+by*image.blocks_width], image.blocks[bx+(by-1)*image.blocks_width], luma_dc_pred);

      decode_dc_luma(mode, stats.parameters, luma_dc_fit);

      for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         block.luma_res[i] = luma_dc_res[i] = max(min(luma_dc_pred[i] + luma_dc_fit[i], 1.0), 0.0);
   }

   // Calculate difference range
   double max_luma = -1e6, min_luma = 0.0, sum_x = 0, sum_x2 = 0.0; // Min Luma starts with zero to force it to at least zero
   for (unsigned y = 0; y < BLOCK_SIZE; y++)
      for (unsigned x = 0; x < BLOCK_SIZE; x++)
      {
         unsigned index = x+y*BLOCK_SIZE;
         double luma_diff = block.luma[index] - luma_dc_res[index];
         sum_x += luma_diff;
         sum_x2 += luma_diff * luma_diff;
         if (max_luma < luma_diff) max_luma = luma_diff;
         if (min_luma > luma_diff) min_luma = luma_diff;
      }

   // Calculate luma difference scale and offset
   block.luma_diff_offset = min(max((int)ceil(-min_luma * 255.0), 0), 255);
   double luma_offset = block.luma_diff_offset / 255.0;
   block.luma_diff_scale = min(max((int)ceil((max_luma + luma_offset) * 127.0), 0), 255);
   if (block.luma_diff_scale < 1) block.luma_diff_scale = 1;
   double luma_scale = block.luma_diff_scale / 127.0;
   stats.luma_range = luma_scale;

   block.luma_pred_mode = mode;
   stats.valid = true;

   return luma_dc_res;
}

void generate_bit_allocation_data(image_t& image, unsigned bx, unsigned by, residual_stats_t& stats, double* residual, unsigned bit_allocation)
{
   block_t& block = image.blocks[bx+by*image.blocks_width];

   double luma_offset = block.luma_diff_offset / 255.0;
   double luma_scale = block.luma_diff_scale / 127.0;
   if (bit_allocation == 0)
   {
      block.luma_diff_scale = 0;
      luma_scale = 0.0;
   }

   double luma_scale_log = log2(block.luma_diff_scale);
   const double inv_luma_scale = 127.0 / block.luma_diff_scale;

   uint16_t factor = 0;
   if (block.luma_diff_scale > 0)
   {
      int needed_bit_allocation = max((int)luma_scale_log, 1);
      if (needed_bit_allocation < bit_allocation) return;
      block.bit_allocation = clamp(bit_allocation, 1, 6);
      factor = (1<<block.bit_allocation)-1;

      // Calculate luma difference
      for (unsigned y = 0; y < BLOCK_SIZE; y++)
         for (unsigned x = 0; x < BLOCK_SIZE; x++)
         {
            unsigned index = x+y*BLOCK_SIZE;
            double luma_base = residual[index] - luma_offset;

            double diff = block.luma[index] - residual[index];

            // Positioned dithering
            double dither[64] = {1, 49, 13, 61, 4, 52, 16, 64, 33, 17, 45, 29, 36, 20, 48, 32, 9, 57, 5, 53, 12, 60, 8, 56, 41, 25, 37, 21, 44, 28, 40, 24, 3, 51, 15, 63, 2, 50, 14, 62, 35, 19, 47, 31, 34, 18, 46, 30, 11, 59, 7, 55, 10, 58, 6, 54, 43, 27, 39, 23, 42, 26, 38, 22};

            double luma_diff = (block.luma[index] - luma_base) * inv_luma_scale;
            int luma_q = clamp((int)floor(luma_diff * factor + dither[(x&7)+(y&7)*8] / 65.0), 0, factor);

            block.luma_diff[index] = luma_q;

            // Decode again for reference as well as for error calculation
            double luma_decoded_diff = luma_q / (double)factor;
            double luma_decoded = clamp(luma_decoded_diff * luma_scale + luma_base, 0.0, 1.0);
            block.luma_res[index] = luma_decoded;
         }
   }
   else
      for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         block.luma_res[i] = residual[i];

   // Calculate error
   double sq_error = 0.0;
   for (unsigned i = 0; i < BLOCK_SIZE2; i++)
      sq_error += sqr(block.luma_res[i] - block.luma[i]);
   stats.error = sq_error / BLOCK_SIZE2;

   sq_error = 0.0;
   for (unsigned i = 0; i < BLOCK_SIZE2; i++)
      sq_error += sqr(block.luma_res[i] - block.luma[i]) * clamp(1.0 - block.contrast[i], 0.1, 1.0);
   stats.weighted_error = sq_error / BLOCK_SIZE2;

   stats.valid = true;
   stats.bits_used = bit_allocation;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ENCODE //
////////////
void encode(string output_filename, uint8_t* data, unsigned width, unsigned height, unsigned quality = 2)
{
   image_t image;
   image.width = width;
   image.height = height;
   image.pixels = new double[width*height*4];
   image.contrast = new double[width*height];

   unsigned size = width * height;

   // Transform data to YUV colour space
   for (unsigned i = 0; i < size; i++)
   {
      double r = data[i*4+0] / 255.0;
      double g = data[i*4+1] / 255.0;
      double b = data[i*4+2] / 255.0;
      double a = data[i*4+3] / 255.0;

      image.pixels[i*4+0] = r * 0.299 + g * 0.587 + b * 0.114;
      image.pixels[i*4+1] = 0.5 - r * 0.168736 - g * 0.331264 + b * 0.5;
      image.pixels[i*4+2] = 0.5 + r * 0.5 - g * 0.418688 - b * 0.081312;
      image.pixels[i*4+3] = a;
   }

   // Calculate constrast
   guassian_blur_filter(image.pixels, width, height, 4, image.contrast, 1);

   for (unsigned i = 0; i < size; i++)
      image.contrast[i] = fabs(image.contrast[i] - image.pixels[i*4]);

   max_filter(image.contrast, width, height, 1, image.contrast, 1);

   // for (unsigned i = 0; i < size; i++)
   //    image.pixels[i*4] = min(image.contrast[i], 1.0);

   // for (unsigned i = 0; i < size; i++)
   //    image.contrast[i] = 0.0;

   // Split into blocks
   image.blocks_width = (width + BLOCK_SIZE - 1) / BLOCK_SIZE;
   image.blocks_height = (height + BLOCK_SIZE - 1) / BLOCK_SIZE;
   image.blocks_size = image.blocks_width * image.blocks_height;
   image.blocks = new block_t[image.blocks_size];

   for (unsigned my = 0; my < image.blocks_height; my++)
      for (unsigned mx = 0; mx < image.blocks_width; mx++)
      {
         block_t& block = image.blocks[mx+my*image.blocks_width];

         for (unsigned by = 0; by < BLOCK_SIZE; by++)
            for (unsigned bx = 0; bx < BLOCK_SIZE; bx++)
            {
               unsigned x = min(mx*BLOCK_SIZE + bx, width - 1), y = min(my*BLOCK_SIZE + by, height - 1);
               unsigned index = x+y*width;
               unsigned block_index = bx+by*BLOCK_SIZE;
               block.contrast[block_index] = image.contrast[index];
               block.luma[block_index] = image.pixels[index*4+0];
               block.chroma_a[block_index] = image.pixels[index*4+1];
               block.chroma_b[block_index] = image.pixels[index*4+2];
            }
      }

   for (unsigned y = 0; y < BLOCK_SIZE; y++)
      for (unsigned x = 0; x < BLOCK_SIZE; x++)
      {
         double weight = (x < 2 || y < 2 || x >= BLOCK_SIZE-2 || y >= BLOCK_SIZE-2) ? 5.0 : 0.2;
         gsl_vector_set(luma_weight_vec, x+y*BLOCK_SIZE, weight);
      }

   // Heuristics for compression
   const unsigned block_size_comp[7] = {
      0,
      BLOCK_SIZE2 / 2,
      BLOCK_SIZE2 / 2 + BLOCK_SIZE2 / 8,
      BLOCK_SIZE2 - BLOCK_SIZE2 / 4,
      BLOCK_SIZE2 - BLOCK_SIZE2 / 8,
      BLOCK_SIZE2 - BLOCK_SIZE2 / 8,
      BLOCK_SIZE2 - BLOCK_SIZE2 / 16};

   for (unsigned by = 0; by < image.blocks_height; by++)
      for (unsigned bx = 0; bx < image.blocks_width; bx++)
      {
         unsigned lowest_mode = 0;
         unsigned lowest_bit_alloc = 1;
         double   lowest_score = 1e6;

         for (unsigned mode = 0; mode < 8; mode++)
         {
            mode_stats_t mode_stats;
            double* residual = generate_mode_data(image, bx, by, mode_stats, mode);
            if (!mode_stats.valid || !residual) break;

            for (unsigned bit_allocation = 0; bit_allocation <= 6; bit_allocation++)
            {
               if (bit_allocation == 0 && mode_stats.luma_range > 0.001 && (mode != 7 || mode_stats.luma_range > 0.003)) continue;

               residual_stats_t residual_stats;
               generate_bit_allocation_data(image, bx, by, residual_stats, residual, bit_allocation);
               if (!residual_stats.valid) break;

               unsigned cost  = residual_stats.bits_used * block_size_comp[bit_allocation] + mode_stats.bits_used;
               double   weight = bit_allocation == 0 ? 0.0 : quality / 7.0;
               double   error = residual_stats.weighted_error * weight + residual_stats.error * (1.0 - weight);
               double   score = calculate_score(cost, error, quality);
               if (score < lowest_score)
               {
                  lowest_score = score;
                  lowest_bit_alloc = bit_allocation;
                  lowest_mode = mode;
               }
            }
         }

         mode_stats_t mode_stats;
         double* residual = generate_mode_data(image, bx, by, mode_stats, lowest_mode);
         residual_stats_t residual_stats;
         generate_bit_allocation_data(image, bx, by, residual_stats, residual, lowest_bit_alloc);

         block_t& block = image.blocks[bx+by*image.blocks_width];

         // Determine chroma prediction coefficients
         double chroma_a_c0, chroma_a_c1, chroma_b_c0, chroma_b_c1;
         {
            bool chroma_a_same = true, chroma_b_same = true;
            double cov00, cov01, cov11, sumsq;

            double first_value = block.chroma_a[0];
            for (unsigned i = 1; i < BLOCK_SIZE2; i++)
               if (fabs(first_value - block.chroma_a[i]) > 0.01)
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
            {
               if (fabs(first_value - block.chroma_b[i]) > 0.01)
               {
                  chroma_b_same = false;
                  break;
               }
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
         }

         //printf("%6f %6f %6f %6f\n", chroma_a_c0, chroma_a_c1, chroma_b_c0, chroma_b_c1);

         block.chroma_a_c0[0] = toInt(chroma_a_c0 * 0.25 + 0.5);
         block.chroma_a_c1[0] = toInt(chroma_a_c1 * 0.5 + 0.5);
         block.chroma_b_c0[0] = toInt(chroma_b_c0 * 0.25 + 0.5);
         block.chroma_b_c1[0] = toInt(chroma_b_c1 * 0.5 + 0.5);

         //printf("%3i %3i %3i %3i\n", block.chroma_a_c0[0], block.chroma_a_c1[0], block.chroma_b_c0[0], block.chroma_b_c1[0]);
      }

   /// Initial transformation and stats collection phase
   unsigned long stats_luma_mode[64] = {0};
   unsigned long stats_luma_1bit[256] = {0};
   unsigned long stats_luma_2bit[256] = {0};
   unsigned long stats_luma_3bit[64] = {0};
   unsigned long stats_luma_4bit[256] = {0};
   unsigned long stats_luma_5bit[32] = {0};
   unsigned long stats_luma_6bit[64] = {0};
   unsigned long stats_luma_mode_total = 0;
   unsigned long stats_luma_1bit_total = 0;
   unsigned long stats_luma_2bit_total = 0;
   unsigned long stats_luma_3bit_total = 0;
   unsigned long stats_luma_4bit_total = 0;
   unsigned long stats_luma_5bit_total = 0;
   unsigned long stats_luma_6bit_total = 0;

   stats_luma_mode[63]++; // Hack

   uint8_t* luma_modes = new uint8_t[image.blocks_size/2];
   {
      for (unsigned i = 0; i < image.blocks_size; i += 2)
         luma_modes[i>>1] = (image.blocks[i].luma_pred_mode << 3) + image.blocks[i+1].luma_pred_mode;

      for (unsigned i = 0; i < image.blocks_size/2; i++)
         stats_luma_mode[luma_modes[i]]++;
   }

   for (unsigned by = 0; by < image.blocks_height; by++)
      for (unsigned bx = 0; bx < image.blocks_width; bx++)
      {
         block_t& block = image.blocks[bx+by*image.blocks_width];

         if (block.luma_diff_scale > 0)
         {
            block.luma_compress_mode = 0;

            // if (block.bit_allocation == 4)//(bx == 16 && by == 9)
            // {
            //    printf("%i %i\n", bx, by);
            //    for (unsigned y = 0; y < BLOCK_SIZE; y++)
            //    {
            //       for (unsigned x = 0; x < BLOCK_SIZE; x++)
            //          printf("%3i ", block.luma_diff[x+y*BLOCK_SIZE]);
            //       printf("\n");
            //    }
            //    printf("\n");
            // }

            if (block.bit_allocation == 1)
            {
               for (int y = 0; y < BLOCK_SIZE; y++)
                  for (int x = 0; x < BLOCK_SIZE; x++) {
                     int i = x + y * BLOCK_SIZE;
                     block.luma_comp[i] = block.luma_diff[i] ^ ((x & 1) != (y & 1));
                  }
            }
            else
            {
               mtf_encode8(block.luma_diff, block.luma_comp, BLOCK_SIZE2);
            }

            // if (block.bit_allocation == 4)//(bx == 16 && by == 9)
            // {
            //    for (unsigned y = 0; y < BLOCK_SIZE; y++)
            //    {
            //       for (unsigned x = 0; x < BLOCK_SIZE; x++)
            //          printf("%3i ", block.luma_comp[x+y*BLOCK_SIZE]);
            //       printf("\n");
            //    }
            //    printf("\n");
            // }

            if (block.bit_allocation == 1)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i += 8)
               {
                  unsigned symbol = 0;
                  for (unsigned j = 0; j < 8; j++)
                     symbol = (symbol << 1) + (block.luma_comp[i+j] & 1);
                  stats_luma_1bit[symbol]++;
               }
            }
            else if (block.bit_allocation == 2)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i += 4)
               {
                  unsigned symbol = 0;
                  for (unsigned j = 0; j < 4; j++)
                     symbol = (symbol << 2) + (block.luma_comp[i+j] & 3);
                  stats_luma_2bit[symbol]++;
               }
            }
            else if (block.bit_allocation == 3)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
               {
                  unsigned symbol = 0;
                  for (unsigned j = 0; j < 2; j++)
                     symbol = (symbol << 3) + (block.luma_comp[i+j] & 7);
                  stats_luma_3bit[symbol]++;
               }
            }
            else if (block.bit_allocation == 4)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i += 2)
               {
                  unsigned symbol = 0;
                  for (unsigned j = 0; j < 2; j++)
                     symbol = (symbol << 4) + (block.luma_comp[i+j] & 15);
                  stats_luma_4bit[symbol]++;
               }
            }
            else if (block.bit_allocation == 5)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                  stats_luma_5bit[block.luma_comp[i] & 31]++;
            }
            else if (block.bit_allocation == 6)
            {
               for (unsigned i = 0; i < BLOCK_SIZE2; i++)
                  stats_luma_6bit[block.luma_comp[i] & 63]++;
            }
         }
      }

   /// Huffman tree construction
   {
      for (unsigned i = 0; i < 64; i++) stats_luma_mode_total += stats_luma_mode[i];
      if (stats_luma_mode_total > 0)
         image.luma_mode_coder = build_canonical_huffman(stats_luma_mode, 64);

      for (unsigned i = 0; i < 256; i++) stats_luma_1bit_total += stats_luma_1bit[i];
      if (stats_luma_1bit_total > 0)
         image.luma_1bit_coder = build_canonical_huffman(stats_luma_1bit, 256);

      for (unsigned i = 0; i < 256; i++) stats_luma_2bit_total += stats_luma_2bit[i];
      if (stats_luma_2bit_total > 0)
         image.luma_2bit_coder = build_canonical_huffman(stats_luma_2bit, 256);

      for (unsigned i = 0; i < 64; i++) stats_luma_3bit_total += stats_luma_3bit[i];
      if (stats_luma_3bit_total > 0)
         image.luma_3bit_coder = build_canonical_huffman(stats_luma_3bit, 64);

      for (unsigned i = 0; i < 16; i++) stats_luma_4bit_total += stats_luma_4bit[i];
      if (stats_luma_4bit_total > 0)
         image.luma_4bit_coder = build_canonical_huffman(stats_luma_4bit, 256);

      for (unsigned i = 0; i < 32; i++) stats_luma_5bit_total += stats_luma_5bit[i];
      if (stats_luma_5bit_total > 0)
         image.luma_5bit_coder = build_canonical_huffman(stats_luma_5bit, 32);

      for (unsigned i = 0; i < 64; i++) stats_luma_6bit_total += stats_luma_6bit[i];
      if (stats_luma_6bit_total > 0)
         image.luma_6bit_coder = build_canonical_huffman(stats_luma_6bit, 64);
   }
   
   unsigned long stats_luma_bit_length[6] = {0, 0, 0, 0, 0, 0};
   unsigned long stats_luma_bit_blocks[6] = {0, 0, 0, 0, 0, 0};

   /// Block compression and writing
   bit_file_writer_t writer(output_filename);
   
   {
      // Write magic bytes
      writer.write_byte('P');
      writer.write_byte('A');
      writer.write_byte('C');
      writer.write_byte('%');

      // Write image size
      write_variable_length_number(image.width, writer);
      write_variable_length_number(image.height, writer);

      printf(" block count: %i\n", image.blocks_size);
      printf(" basic header info: %i\n", writer.get_bytes_written());

      // Write block modes   
      write_canonical_huffman(image.luma_mode_coder, writer);
      image.luma_mode_coder->encode(luma_modes, image.blocks_size/2, writer);
      delete [] luma_modes;

      // Write block mode values
      for (size_t i = 0; i < image.blocks_size; i++)
      {
         block_t& block = image.blocks[i];

         uint8_t components = 0;
         if (block.luma_pred_mode == 0)
            components = 1;
         else if (block.luma_pred_mode == 1)
            components = 3;
         else if (block.luma_pred_mode == 7)
            components = 4;
         else
            components = 1;
         for (size_t j = 0; j < components; j++)
            writer.write_byte(block.luma_dc[j]);
      }

      // Write block bit allocation values
      {
         uint8_t* block_bit_allocs = new uint8_t[image.blocks_size];
         uint8_t* block_bit_allocs_encoded = new uint8_t[image.blocks_size];

         for (size_t i = 0; i < image.blocks_size; i++) 
            block_bit_allocs[i] = image.blocks[i].bit_allocation;

         mtf_encode8(block_bit_allocs, block_bit_allocs_encoded, image.blocks_size);

         for (size_t i = 0; i < image.blocks_size; i++) 
            write_rice2_code(block_bit_allocs_encoded[i], writer);

         delete [] block_bit_allocs;
         delete [] block_bit_allocs_encoded;
      }

      printf(" + mode data size: %i\n", writer.get_bytes_written());

      // Write block luma dc and chroma predictors
      for (size_t i = 0; i < image.blocks_size; i++)
      {
         block_t& block = image.blocks[i];
         if (block.bit_allocation > 0)
         {
            writer.write_byte(block.luma_diff_offset);
            writer.write_byte(block.luma_diff_scale);
         }
      }

      printf(" + dc pred data size: %i\n", writer.get_bytes_written());

      uint8_t* buffer1 = new uint8_t[image.blocks_size];
      uint8_t* buffer2 = new uint8_t[image.blocks_size];

      for (size_t i = 0; i < image.blocks_size; i++) buffer1[i] = image.blocks[i].chroma_a_c0[0];
      delta_encode8(buffer1, buffer2, image.blocks_size);

      writer.write_byte(buffer2[0]);
      for (size_t i = 1; i < image.blocks_size; i++) write_elias_gamma_code(buffer2[i], writer);

      for (size_t i = 0; i < image.blocks_size; i++) buffer1[i] = image.blocks[i].chroma_a_c1[0];
      delta_encode8(buffer1, buffer2, image.blocks_size);

      writer.write_byte(buffer2[0]);
      for (size_t i = 1; i < image.blocks_size; i++) write_elias_gamma_code(buffer2[i], writer);

      for (size_t i = 0; i < image.blocks_size; i++) buffer1[i] = image.blocks[i].chroma_b_c0[0];
      delta_encode8(buffer1, buffer2, image.blocks_size);

      writer.write_byte(buffer2[0]);
      for (size_t i = 1; i < image.blocks_size; i++) write_elias_gamma_code(buffer2[i], writer);

      for (size_t i = 0; i < image.blocks_size; i++) buffer1[i] = image.blocks[i].chroma_b_c1[0];
      delta_encode8(buffer1, buffer2, image.blocks_size);

      writer.write_byte(buffer2[0]);
      for (size_t i = 1; i < image.blocks_size; i++) write_elias_gamma_code(buffer2[i], writer);

      writer.pad_to_byte(0);

      delete [] buffer1;
      delete [] buffer2;

      printf(" + chroma pred data size: %i\n", writer.get_bytes_written());

      // Write huffman trees
      uint8_t used_trees = 0;
      used_trees |= (stats_luma_1bit_total > 0 ? 1 : 0) << 0;
      used_trees |= (stats_luma_2bit_total > 0 ? 1 : 0) << 1;
      used_trees |= (stats_luma_3bit_total > 0 ? 1 : 0) << 2;
      used_trees |= (stats_luma_4bit_total > 0 ? 1 : 0) << 3;
      used_trees |= (stats_luma_5bit_total > 0 ? 1 : 0) << 4;
      used_trees |= (stats_luma_6bit_total > 0 ? 1 : 0) << 5;
      writer.write_byte(used_trees);

      if (stats_luma_1bit_total > 0) write_canonical_huffman(image.luma_1bit_coder, writer);
      if (stats_luma_2bit_total > 0) write_canonical_huffman(image.luma_2bit_coder, writer);
      if (stats_luma_3bit_total > 0) write_canonical_huffman(image.luma_3bit_coder, writer);
      if (stats_luma_4bit_total > 0) write_canonical_huffman(image.luma_4bit_coder, writer);
      if (stats_luma_5bit_total > 0) write_canonical_huffman(image.luma_5bit_coder, writer);
      if (stats_luma_6bit_total > 0) write_canonical_huffman(image.luma_6bit_coder, writer);
      writer.pad_to_byte(0);

      printf(" + residual huffman: %i\n", writer.get_bytes_written());
   }

   unsigned block_bit_counts[6] = {0, 0, 0, 0, 0, 0};

   for (unsigned by = 0; by < image.blocks_height; by++)
      for (unsigned bx = 0; bx < image.blocks_width; bx++)
      {
         block_t& block = image.blocks[bx+by*image.blocks_width];

         if (block.luma_diff_scale > 0 && block.bit_allocation <= 6)
            block_bit_counts[block.bit_allocation-1]++;
      }

   for (unsigned b = 0; b < 6; b++)
   {
      if (block_bit_counts[b] == 0) continue;
      unsigned bit_allocation = b + 1;
      unsigned total_length = block_bit_counts[b] * BLOCK_SIZE2;
      unsigned index = 0;
      uint8_t* luma_all = new uint8_t[total_length];

      for (unsigned by = 0; by < image.blocks_height; by++)
         for (unsigned bx = 0; bx < image.blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*image.blocks_width];

            if (block.luma_diff_scale > 0 && block.bit_allocation == bit_allocation)
               for (unsigned i = 0; i < BLOCK_SIZE2; i++) 
                  luma_all[index++] = block.luma_comp[i];
         } 

      unsigned start_bytes = writer.get_bytes_written();

      if (bit_allocation == 1)
      {
         uint8_t* temp = new uint8_t[total_length/8];
         compress_1bit_to_8bit(luma_all, total_length, temp); 
         image.luma_1bit_coder->encode(temp, total_length/8, writer);
         delete [] temp;
      }
      else if (bit_allocation == 2)
      {
         uint8_t* temp = new uint8_t[total_length/4];
         compress_2bit_to_8bit(luma_all, total_length, temp); 
         image.luma_2bit_coder->encode(temp, total_length/4, writer);
         delete [] temp;
      }
      else if (bit_allocation == 3)
      {
         uint8_t* temp = new uint8_t[total_length/2];
         compress_3bit_to_6bit(luma_all, total_length, temp); 
         image.luma_3bit_coder->encode(temp, total_length/2, writer);
         delete [] temp;
      }
      else if (bit_allocation == 4)
      {
         uint8_t* temp = new uint8_t[total_length/2];
         compress_4bit_to_8bit(luma_all, total_length, temp); 
         image.luma_4bit_coder->encode(temp, total_length/2, writer);
         delete [] temp;
      }
      else if (bit_allocation == 5)
         image.luma_5bit_coder->encode(luma_all, total_length, writer);
      else if (bit_allocation == 6)
         image.luma_6bit_coder->encode(luma_all, total_length, writer);

      stats_luma_bit_length[b] = writer.get_bytes_written() - start_bytes;
      stats_luma_bit_blocks[b] = block_bit_counts[b];

      delete [] luma_all;
   }

   writer.pad_to_byte();

   for (int i = 0; i < 6; i++)
      if (stats_luma_bit_blocks[i] > 0) 
      {
         double encoded_per_block = (double)stats_luma_bit_length[i] / stats_luma_bit_blocks[i];
         double original_per_block = (BLOCK_SIZE2 * (i + 1) + 7) / 8;
         printf("%i-bit block:\n", i+1);
         printf(" - Encoded bytes per block:  %f\n", encoded_per_block);
         printf(" - Original bytes per block: %f\n", original_per_block);
         printf(" - Compression: %f%%\n", encoded_per_block * 100.0 / original_per_block);
      }

   int bytes_used = writer.get_bytes_written();
   cout << "Bits per pixel:" << ((double)bytes_used*8/width/height) << " bpp" << endl;
   cout << "Bytes used:" << bytes_used << " bytes" << endl;
   cout << "Kilobytes used:" << (bytes_used/1000.0) << " kB" << endl;

   gsl_multifit_linear_free(fitting_workspace3);
   delete [] image.pixels;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DECODE //
////////////
uint8_t* decode(string input_filename, unsigned& width, unsigned& height)
{
   cout << " ==== DECODING ==== " << endl;

   image_t image;

   unsigned blocks_width;
   unsigned blocks_height;
   unsigned blocks_size;

   unsigned size;

   bit_file_reader_t reader(input_filename);

   if (!reader.is_open())
   {
      cout << "Could not open file" << endl;
      return nullptr;
   }

   uint8_t magic[4];
   magic[0] = reader.read_byte();
   magic[1] = reader.read_byte();
   magic[2] = reader.read_byte();
   magic[3] = reader.read_byte();
   assert(magic[0] == 'P' && magic[1] == 'A' && magic[2] == 'C' && magic[3] == '%');

   // Read image size
   image.width = read_variable_length_number(reader);
   image.height = read_variable_length_number(reader);
   printf("size: %ix%i\n", image.width, image.height);
   
   // Calculate block sizes
   size = image.width * image.height;

   blocks_width = (image.width + BLOCK_SIZE - 1) / BLOCK_SIZE;
   blocks_height = (image.height + BLOCK_SIZE - 1) / BLOCK_SIZE;
   blocks_size = blocks_width * blocks_height;
   image.blocks = new block_t[blocks_size];

   // Read in block modes
   printf("read: %i\n", reader.get_bytes_read());
   uint8_t* block_modes = new uint8_t[blocks_size/2];
   image.luma_mode_coder = read_canonical_huffman(reader);
   image.luma_mode_coder->decode(reader, block_modes, blocks_size/2);

   for (size_t i = 0; i < blocks_size; i += 2)
   {
      image.blocks[i  ].luma_pred_mode = block_modes[i>>1] >> 3;
      image.blocks[i+1].luma_pred_mode = block_modes[i>>1] & 7;
   }
   delete [] block_modes;

   // Read in block values
   for (size_t i = 0; i < blocks_size; i++)
   {
      block_t& block = image.blocks[i];

      uint8_t components = 0;
      if (block.luma_pred_mode == 0)
         components = 1;
      else if (block.luma_pred_mode == 1)
         components = 3;
      else if (block.luma_pred_mode == 7)
         components = 4;
      else
         components = 1;
      for (size_t j = 0; j < components; j++)
         block.luma_dc[j] = reader.read_byte();
   }

   {
      uint8_t* block_bit_allocs_encoded = new uint8_t[blocks_size];
      uint8_t* block_bit_allocs = new uint8_t[blocks_size];

      for (size_t i = 0; i < blocks_size; i++) 
         block_bit_allocs_encoded[i] = read_rice2_code(reader);

      mtf_decode8(block_bit_allocs_encoded, block_bit_allocs, blocks_size);

      for (size_t i = 0; i < blocks_size; i++) 
         image.blocks[i].bit_allocation = block_bit_allocs[i];

      delete [] block_bit_allocs;
      delete [] block_bit_allocs_encoded;
   }

   for (size_t i = 0; i < blocks_size; i++)
   {
      block_t& block = image.blocks[i];

      if (block.bit_allocation > 0)
      {
         block.luma_diff_offset = reader.read_byte();
         block.luma_diff_scale = reader.read_byte();
      }
   }

   uint8_t* buffer1 = new uint8_t[blocks_size];
   uint8_t* buffer2 = new uint8_t[blocks_size];

   buffer1[0] = reader.read_byte();
   for (unsigned i = 1; i < blocks_size; i++)
      buffer1[i] = read_elias_gamma_code(reader);

   delta_decode8(buffer1, buffer2, blocks_size);

   for (unsigned i = 0; i < blocks_size; i++) image.blocks[i].chroma_a_c0[0] = buffer2[i];

   buffer1[0] = reader.read_byte();
   for (unsigned i = 1; i < blocks_size; i++)
      buffer1[i] = read_elias_gamma_code(reader);

   delta_decode8(buffer1, buffer2, blocks_size);

   for (unsigned i = 0; i < blocks_size; i++) image.blocks[i].chroma_a_c1[0] = buffer2[i];

   buffer1[0] = reader.read_byte();
   for (unsigned i = 1; i < blocks_size; i++)
      buffer1[i] = read_elias_gamma_code(reader);

   delta_decode8(buffer1, buffer2, blocks_size);

   for (unsigned i = 0; i < blocks_size; i++) image.blocks[i].chroma_b_c0[0] = buffer2[i];

   buffer1[0] = reader.read_byte();
   for (unsigned i = 1; i < blocks_size; i++)
      buffer1[i] = read_elias_gamma_code(reader);

   delta_decode8(buffer1, buffer2, blocks_size);

   for (unsigned i = 0; i < blocks_size; i++) image.blocks[i].chroma_b_c1[0] = buffer2[i];

   reader.skip_to_next_byte();

   delete [] buffer1;
   delete [] buffer2;

   // Read in huffman trees
   uint8_t used_trees = reader.read_byte();
   printf("used: %4x\n", used_trees);
   printf("read: %i\n", reader.get_bytes_read());
   if (((used_trees >> 0) & 1) == 1) image.luma_1bit_coder = read_canonical_huffman(reader);
   if (((used_trees >> 1) & 1) == 1) image.luma_2bit_coder = read_canonical_huffman(reader);
   if (((used_trees >> 2) & 1) == 1) image.luma_3bit_coder = read_canonical_huffman(reader);
   if (((used_trees >> 3) & 1) == 1) image.luma_4bit_coder = read_canonical_huffman(reader);
   if (((used_trees >> 4) & 1) == 1) image.luma_5bit_coder = read_canonical_huffman(reader);
   if (((used_trees >> 5) & 1) == 1) image.luma_6bit_coder = read_canonical_huffman(reader);
   reader.skip_to_next_byte();

   // Decode residual
   unsigned block_bit_counts[6] = {0, 0, 0, 0, 0, 0};

   for (unsigned by = 0; by < blocks_height; by++)
      for (unsigned bx = 0; bx < blocks_width; bx++)
      {
         block_t& block = image.blocks[bx+by*blocks_width];

         if (block.luma_diff_scale > 0 && block.bit_allocation <= 6)
            block_bit_counts[block.bit_allocation-1]++;
      }

   for (unsigned b = 0; b < 6; b++)
   {
      if (block_bit_counts[b] == 0) continue;
      unsigned bit_allocation = b + 1;
      unsigned total_length = block_bit_counts[b] * BLOCK_SIZE2;
      unsigned index = 0;
      uint8_t* luma_all = new uint8_t[total_length];

      if (bit_allocation == 1)
      {
         uint8_t* temp = new uint8_t[total_length/8];
         image.luma_1bit_coder->decode(reader, temp, total_length/8);
         expand_8bit_to_1bit(temp, total_length/8, luma_all);
         delete [] temp;
      }
      else if (bit_allocation == 2)
      {
         uint8_t* temp = new uint8_t[total_length/4];
         image.luma_2bit_coder->decode(reader, temp, total_length/4);
         expand_8bit_to_2bit(temp, total_length/4, luma_all);
         delete [] temp;
      }
      else if (bit_allocation == 3)
      {
         uint8_t* temp = new uint8_t[total_length/2];
         image.luma_3bit_coder->decode(reader, temp, total_length/2);
         expand_6bit_to_3bit(temp, total_length/2, luma_all);
         delete [] temp;
      }
      else if (bit_allocation == 4)
      {
         uint8_t* temp = new uint8_t[total_length/2];
         image.luma_4bit_coder->decode(reader, temp, total_length/2);
         expand_8bit_to_4bit(temp, total_length/2, luma_all);
         delete [] temp;
      }
      else if (bit_allocation == 5)
         image.luma_5bit_coder->decode(reader, luma_all, total_length);
      else if (bit_allocation == 6)
         image.luma_6bit_coder->decode(reader, luma_all, total_length);

      for (unsigned by = 0; by < blocks_height; by++)
         for (unsigned bx = 0; bx < blocks_width; bx++)
         {
            block_t& block = image.blocks[bx+by*blocks_width];

            // Decode luma difference
            if (block.luma_diff_scale > 0 && bit_allocation == block.bit_allocation)
               for (unsigned i = 0; i < BLOCK_SIZE2; i++) 
                  block.luma_comp[i] = luma_all[index++];
         }
   }

   // Calculate final values
   for (unsigned by = 0; by < blocks_height; by++)
      for (unsigned bx = 0; bx < blocks_width; bx++)
      {     
         block_t& block = image.blocks[bx+by*blocks_width];

         int mode = block.luma_pred_mode;

         double luma_dc[16];
         if (mode == 0)
         {
            luma_dc[0] = fromInt(block.luma_dc[0]);
         }
         else if (mode == 7)
         {
            for (int i = 0; i < 4; i++)
               luma_dc[i] = fromInt(block.luma_dc[i]);
         }
         else
         {
            for (int i = 0; i < 3; i++)
               luma_dc[i] = (fromInt(block.luma_dc[i]) - 0.5) * 2.0;
         }

         double luma_dc_pred[BLOCK_SIZE2] = {0.0};
         double luma_dc_fit[BLOCK_SIZE2] = {0.0};
         double luma_dc_res[BLOCK_SIZE2] = {0.0};

         for (unsigned y = 0; y < BLOCK_SIZE; y++)
            for (unsigned x = 0; x < BLOCK_SIZE; x++)
               luma_dc_pred[x+y*BLOCK_SIZE] = 0.0;

         if (mode > 1 && mode < 7)
            predict_luma(mode, image.blocks[(bx-1)+(by-1)*blocks_width], image.blocks[(bx-1)+by*blocks_width], image.blocks[bx+(by-1)*blocks_width], luma_dc_pred);

         decode_dc_luma(mode, luma_dc, luma_dc_fit);

         for (unsigned i = 0; i < BLOCK_SIZE2; i++)
            luma_dc_res[i] = max(min(luma_dc_pred[i] + luma_dc_fit[i], 1.0), 0.0);

         if (block.luma_diff_scale > 0)
         {
            double luma_offset = block.luma_diff_offset / 255.0;
            double luma_scale = block.luma_diff_scale / 127.0;
            uint16_t factor = (1<<block.bit_allocation)-1;

            if (block.bit_allocation == 1) 
            {
               for (int y = 0; y < BLOCK_SIZE; y++)
                  for (int x = 0; x < BLOCK_SIZE; x++) {
                     int i = x + y * BLOCK_SIZE;
                     block.luma_diff[i] = block.luma_comp[i] ^ ((x & 1) != (y & 1));
                  }
            }
            else 
            {
               mtf_decode8(block.luma_comp, block.luma_diff, BLOCK_SIZE2);
            }

            for (unsigned y = 0; y < BLOCK_SIZE; y++)
               for (unsigned x = 0; x < BLOCK_SIZE; x++)
               {
                  unsigned index = x+y*BLOCK_SIZE;
                  double luma_q = block.luma_diff[index];
                  double luma_diff = (luma_q / factor) * luma_scale - luma_offset;
                  double luma = luma_diff + luma_dc_res[index];
                  block.luma[index] = luma;
                  block.luma_res[index] = luma;
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
               }
         }

         // Predict chroma
         double chroma_a_c0 = fromInt(block.chroma_a_c0[0]) * 4.0 - 2.0;
         double chroma_a_c1 = fromInt(block.chroma_a_c1[0]) * 2.0 - 1.0;
         double chroma_b_c0 = fromInt(block.chroma_b_c0[0]) * 4.0 - 2.0;
         double chroma_b_c1 = fromInt(block.chroma_b_c1[0]) * 2.0 - 1.0;

         for (unsigned i = 0; i < BLOCK_SIZE2; i++)
         {
            block.chroma_a[i] = 0.5;
            block.chroma_b[i] = 0.5;
         }

         for (unsigned y = 0; y < BLOCK_SIZE; y++)
            for (unsigned x = 0; x < BLOCK_SIZE; x++)
            {
               unsigned index = x+y*BLOCK_SIZE;
               double L = block.luma[index];
               double a = max(min(L * chroma_a_c1 + chroma_a_c0, 1.0), 0.0);
               double b = max(min(L * chroma_b_c1 + chroma_b_c0, 1.0), 0.0);
               //block.chroma_a[index] = a;
               //block.chroma_b[index] = b;
            }
      }

   image.pixels = new double[size*4];

   for (unsigned my = 0; my < blocks_height; my++)
      for (unsigned mx = 0; mx < blocks_width; mx++)
      {
         block_t& block = image.blocks[mx+my*blocks_width];
         for (unsigned by = 0; by < BLOCK_SIZE; by++)
            for (unsigned bx = 0; bx < BLOCK_SIZE; bx++)
            {
               unsigned x = mx*BLOCK_SIZE + bx, y = my*BLOCK_SIZE + by;
               unsigned index = x+y*image.width;
               unsigned block_index = bx+by*BLOCK_SIZE;

               if (x >= image.width || y >= image.height) continue;
               
               image.pixels[index*4+0] = block.luma[block_index];
               image.pixels[index*4+1] = block.chroma_a[block_index];
               image.pixels[index*4+2] = block.chroma_b[block_index];
            }
      }

   delete [] image.blocks;

   uint8_t* data = new uint8_t[image.width*image.height*4];
   for (unsigned i = 0; i < size; i++)
   {
      double y = image.pixels[i*4+0];
      double cb = image.pixels[i*4+1];
      double cr = image.pixels[i*4+2];

      double r = y + (cr - 0.5) * 1.402;
      double g = y - (cb - 0.5) * 0.34414 - (cr - 0.5) * 0.71414;
      double b = y + (cb - 0.5) * 1.722;

      data[i*4+0] = toInt(r);
      data[i*4+1] = toInt(g);
      data[i*4+2] = toInt(b);
      data[i*4+3] = 255;
   }

   delete [] image.pixels;

   width = image.width;
   height = image.height;
   return data;
}

int main(int argc, char* args[])
{
   if (argc < 2)
   {
      printf("Usage: %s <filename> [<quality>]", args[0]);
      return 1;
   }

   unsigned quality = (unsigned)(argc >= 3 ? atoi(args[2]) : 4);
   string   name = string(args[1]);
   string   filename = name + ".png";
   string   result_filename = name + "_out.png";

   {
      {
         unsigned width, height;
         uint8_t* image;
         unsigned error;

         error = lodepng_decode32_file(&image, &width, &height, filename.c_str());
         if (error)
         {
            cout << "Could not load file: " << filename << endl;
            return 1;
         }

         encode("test/a.pac", image, width, height, quality);
      }

      {
         unsigned width, height;
         uint8_t* image;
         unsigned error;

         image = decode("test/a.pac", width, height);

         error = lodepng_encode32_file(result_filename.c_str(), image, width, height);
         if (error)
         {
            cout << "Could not save file: " << result_filename << endl;
            return 2;
         }

         cout << "Output: " << result_filename << endl;
      }
   }

   return 0;
}