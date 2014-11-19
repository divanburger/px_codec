#ifndef COMPRESSION_H_
#define COMPRESSION_H_

#include "../bitio/bit.h"

#include "math.h"

// ------------------------
// Move-to-front transform
// ------------------------

void mtf_encode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t symbols[256];
   for (unsigned i = 0; i < 256; i++) symbols[i] = i;

   for (unsigned i = 0; i < n; i++)
   {
      uint8_t in = input[i];

      // Look-up symbol
      uint8_t sym = 0;
      for (unsigned j = 0; j < 256; j++)
         if (symbols[j] == in)
         {
            sym = j;
            break;
         }

      // Output symbol
      *(output++) = sym;

      // Move up symbol table
      for (unsigned j = sym; j > 0; j--) symbols[j] = symbols[j-1];
      symbols[0] = in;
   }
}

void mtf_decode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t symbols[256];
   for (unsigned i = 0; i < 256; i++) symbols[i] = i;

   for (unsigned i = 0; i < n; i++)
   {
      uint8_t in = input[i];
      uint8_t sym = symbols[in];

      // Output symbol
      *(output++) = sym;

      // Move up symbol table
      for (unsigned j = in; j > 0; j--) symbols[j] = symbols[j-1];
      symbols[0] = sym;
   }
}

// --------------
// Delta-encoding
// --------------

void delta_encode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n == 0) return;
   output[0] = input[0];

   for (unsigned i = 1; i < n; i++)
   {
      int b = (int)input[i-1];
      int m = b <= 127 ? b : 256 - b;
      int d = (int)input[i] - b;
      if (d >= 0)
         output[i] = imin(d, m) * 2 + imax(d - m, 0);
      else
         output[i] = imin(-d, m) * 2 - 1 + imax(-d - m, 0);
   }
}

void delta_decode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n == 0) return;
   output[0] = input[0];

   int b = output[0];
   for (unsigned i = 1; i < n; i++)
   {
      int m = b <= 127 ? b : 256 - b;
      int d = (int)input[i];
      int nd = 0;
      if ((d <= m*2-1) ? (d & 1) == 0 : b <= 127)
         nd = imin(d/2, m) + imax(d - m * 2, 0);
      else
         nd = -imin((d+1)/2, m) - imax(d - m * 2 + 1, 0);

      b += nd;
      output[i] = b;
   }
}

void delta_encode89(uint8_t* input, uint16_t* output, unsigned n)
{
   if (n == 0) return;
   output[0] = input[0];

   for (unsigned i = 1; i < n; i++)
   {
      int d = (int)input[i] - (int)input[i-1];
      output[i] = (d >= 0) ? d << 1 : ((-d) << 2) + 1;
   }
}

void delta_decode89(uint16_t* input, uint8_t* output, unsigned n)
{
   if (n == 0) return;
   output[0] = input[0];

   int b = output[0];
   for (unsigned i = 1; i < n; i++)
   {
      int d = input[i];
      b += (d & 1) ? d >> 1 : -((d-1) >> 1);
      output[i] = b;
   }
}

// -------------
// Miscellaneous
// -------------

void xor_code8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t last = input[0];
   for (unsigned i = 1; i < n; i++)
   {
      output[i] = last ^ input[i];
      last = input[i];
   }
}

#endif