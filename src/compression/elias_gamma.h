#ifndef COMPRESSION_ELIAS_GAMMA_H_
#define COMPRESSION_ELIAS_GAMMA_H_

#include "../bitio/bit.h"

void write_elias_gamma_code(unsigned symbol, bit_writer_t& writer)
{
   symbol++;
   unsigned s = symbol;
   unsigned len = 0;
   while (s > 1)
   {
      s >>= 1;
      len++;
   }
   for (unsigned i = 0; i < len; i++) writer.write_bit(0);
   writer.write_bit(1);
   for (unsigned i = 0; i < len; i++) writer.write_bit((symbol>>(len-1-i))&1);
}

unsigned read_elias_gamma_code(bit_reader_t& reader)
{
   unsigned len = 0;
   while (true)
   {
      uint8_t b = reader.read_bit();
      if (b == 0) 
         len++;
      else
         break;
   }
   unsigned symbol = 1;
   for (unsigned i = 0; i < len; i++)
      symbol = (symbol << 1) + reader.read_bit();
   return symbol - 1;
}

#endif