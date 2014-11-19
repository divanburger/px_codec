#ifndef BITMANIP_H_
#define BITMANIP_H_

inline unsigned compress_4bit_to_8bit(uint8_t* input, unsigned n, uint8_t* output)
{
   unsigned n2 = n >> 1;
   for (unsigned i = 0; i < n2; i++)
   {
      unsigned j = i << 1;
      output[i] = (input[j] << 4) + input[j+1];
   }

   if (n & 1)
   {
      output[n2] = input[n2 << 1] << 4;
      return n2 + 1;
   }

   return n2;
}

inline unsigned compress_3bit_to_6bit(uint8_t* input, unsigned n, uint8_t* output)
{
   unsigned n2 = n >> 1;
   for (unsigned i = 0; i < n2; i++)
   {
      unsigned j = i << 1;
      output[i] = (input[j] << 3) + input[j+1];
   }

   if (n & 1)
   {
      output[n2] = input[n2 << 1] << 3;
      return n2 + 1;
   }

   return n2;
}

inline unsigned compress_2bit_to_8bit(uint8_t* input, unsigned n, uint8_t* output)
{
   unsigned n4 = n >> 2;
   for (unsigned i = 0; i < n4; i++)
   {
      unsigned j = i << 2;
      output[i] = (input[j] << 6) + (input[j+1] << 4) + (input[j+2] << 2) + input[j+3];
   }

   if (n & 3)
   {
      unsigned l = n4 << 2, d = n - l;
      uint8_t v = input[l] << 3;
      v += (d >= 1) ? (input[l+1] << 4) : 0;
      v += (d >= 2) ? (input[l+2] << 2) : 0;
      v += (d >= 3) ?  input[l+3]       : 0;
      output[n4] = v;
      return n4 + 1;
   }

   return n4;
}

inline unsigned compress_1bit_to_4bit(uint8_t* input, unsigned n, uint8_t* output)
{
   unsigned n4 = n >> 2;
   for (unsigned i = 0; i < n4; i++)
   {
      unsigned j = i << 2;
      output[i] = (input[j] << 3) + (input[j+1] << 2) + (input[j+2] << 1) + input[j+3];
   }

   if (n & 3)
   {
      unsigned l = n4 << 2, d = n - l;
      uint8_t v = input[l] << 3;
      v += (d >= 1) ? (input[l+1] << 2) : 0;
      v += (d >= 2) ? (input[l+2] << 1) : 0;
      v += (d >= 3) ?  input[l+3]       : 0;
      output[n4] = v;
      return n4 + 1;
   }

   return n4;
}

inline unsigned compress_1bit_to_8bit(uint8_t* input, unsigned n, uint8_t* output)
{
   for (unsigned i = 0; i < n; i += 8)
      output[i>>3] = (input[i] << 7) + (input[i+1] << 6) + (input[i+2] << 5) + (input[i+3] << 4) + (input[i+4] << 3) + (input[i+5] << 2) + (input[i+6] << 1) + input[i+7];
   
   unsigned n8 = n >> 3;
   if (n & 7)
   {
      unsigned l = n8 << 3, d = n - l;
      uint8_t v = input[l] << 7;
      v += (d >= 1) ? (input[l+1] << 6) : 0;
      v += (d >= 2) ? (input[l+2] << 5) : 0;
      v += (d >= 3) ? (input[l+3] << 4) : 0;
      v += (d >= 4) ? (input[l+4] << 3) : 0;
      v += (d >= 5) ? (input[l+5] << 2) : 0;
      v += (d >= 6) ? (input[l+6] << 1) : 0;
      v += (d >= 7) ?  input[l+7]       : 0;
      output[n8] = v;
      return n8 + 1;
   }

   return n8;
}

inline unsigned expand_8bit_to_1bit(uint8_t* input, unsigned n, uint8_t* output)
{
   for (unsigned i = 0; i < n; i++)
   {
      unsigned j = i<<3;
      uint8_t b = input[i];
      output[j+0] = (b >> 7) & 1;
      output[j+1] = (b >> 6) & 1;
      output[j+2] = (b >> 5) & 1;
      output[j+3] = (b >> 4) & 1;
      output[j+4] = (b >> 3) & 1;
      output[j+5] = (b >> 2) & 1;
      output[j+6] = (b >> 1) & 1;
      output[j+7] = (b >> 0) & 1;
   }

   return n*8;
}

inline unsigned expand_8bit_to_2bit(uint8_t* input, unsigned n, uint8_t* output)
{
   for (unsigned i = 0; i < n; i++)
   {
      unsigned j = i<<2;
      uint8_t b = input[i];
      output[j+0] = (b >> 6) & 3;
      output[j+1] = (b >> 4) & 3;
      output[j+2] = (b >> 2) & 3;
      output[j+3] =  b       & 3;
   }

   return n*4;
}

inline unsigned expand_6bit_to_3bit(uint8_t* input, unsigned n, uint8_t* output)
{
   for (unsigned i = 0; i < n; i++)
   {
      unsigned j = i<<1;
      uint8_t b = input[i];
      output[j+0] = (b >> 3) & 7;
      output[j+1] =  b       & 7;
   }

   return n*2;
}

inline unsigned expand_8bit_to_4bit(uint8_t* input, unsigned n, uint8_t* output)
{
   for (unsigned i = 0; i < n; i++)
   {
      unsigned j = i<<1;
      uint8_t b = input[i];
      output[j+0] = (b >> 4) & 15;
      output[j+1] =  b       & 15;
   }

   return n*2;
}

// packs 2 5-bit indices into a 10-bit morton code
inline unsigned morton_2d_encode_5bit(unsigned x, unsigned y)
{ 
   x &= 0x0000001f;
   y &= 0x0000001f;
   x *= 0x01041041;
   y *= 0x01041041;
   x &= 0x10204081;
   y &= 0x10204081;
   x *= 0x00108421;
   y *= 0x00108421;
   x &= 0x15500000;
   y &= 0x15500000;
   return (x >> 20) | (y >> 19);
}

#endif