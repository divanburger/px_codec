#ifndef BITDEBUG_H_
#define BITDEBUG_H_

#include <cassert>

#include "bit.h"

class bit_debug_writer_t : public bit_writer_t
{
public:
   void write_bytes(uint8_t* bytes, unsigned length)
   {
      for (unsigned i = 0; i < length; i++)
         write_byte(bytes[i]);
   }

   void write_byte(uint8_t byte)
   {
      write_bits32(byte, 8);
   }

   void write_bits8(uint8_t bits, uint8_t bit_count)
   {
      write_bits32(bits, bit_count);
   }

   void write_bits16(uint16_t bits, uint8_t bit_count)
   {
      write_bits32(bits, bit_count);
   }

   void write_bits32(uint32_t bits, uint8_t bit_count)
   {
      for (uint8_t i = 0; i < bit_count; i++)
         write_bit((bits >> (bit_count-i-1)) & 1);
   }

   void write_bit(uint8_t bit)
   {
      printf("%c", (bit == 0) ? '0' : '1');
   }

   void pad_to_byte(uint8_t pad_bit = 0)
   {
      printf("P");
   }
};

#endif