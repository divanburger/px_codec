#ifndef COMPRESSION_RICE_H_
#define COMPRESSION_RICE_H_

#include "../bitio/bit.h"

void write_rice2_code(unsigned symbol, bit_writer_t& writer)
{
   unsigned q = symbol >> 2;
   unsigned r = symbol & 0x3;
   for (unsigned i = 0; i < q; i++) writer.write_bit(1);
   writer.write_bit(0);
   writer.write_bits8(r, 2);
}

void write_rice3_code(unsigned symbol, bit_writer_t& writer)
{
   unsigned q = symbol >> 3;
   unsigned r = symbol & 0x7;
   for (unsigned i = 0; i < q; i++) writer.write_bit(1);
   writer.write_bit(0);
   writer.write_bits8(r, 3);
}

void write_rice4_code(unsigned symbol, bit_writer_t& writer)
{
   unsigned q = symbol >> 4;
   unsigned r = symbol & 0xF;
   for (unsigned i = 0; i < q; i++) writer.write_bit(1);
   writer.write_bit(0);
   writer.write_bits8(r, 4);
}

unsigned read_rice2_code(bit_reader_t& reader)
{
   unsigned q = 0;
   while (true)
   {
      uint8_t b = reader.read_bit();
      if (b == 1) 
         q++;
      else
         break;
   }
   return (q << 2) + reader.read_bits8(2);
}

unsigned read_rice3_code(bit_reader_t& reader)
{
   unsigned q = 0;
   while (true)
   {
      uint8_t b = reader.read_bit();
      if (b == 1) 
         q++;
      else
         break;
   }
   return (q << 3) + reader.read_bits8(3);
}

unsigned read_rice4_code(bit_reader_t& reader)
{
   unsigned q = 0;
   while (true)
   {
      uint8_t b = reader.read_bit();
      if (b == 1) 
         q++;
      else
         break;
   }
   return (q << 4) + reader.read_bits8(4);
}

#endif