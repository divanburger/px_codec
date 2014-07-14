#ifndef BITARRAY_H_
#define BITARRAY_H_

#include <cassert>

class bit_array_writer_t
{
public:
   bit_array_writer_t(uint8_t* bytes, unsigned max_size) : buffer(0), bit_size(0), bytes(bytes), max_size(max_size)
   {
   }

   unsigned get_size()
   {
      return bytes_index;
   }

   unsigned get_max_size()
   {
      return max_size;
   }

   uint8_t* get_data()
   {
      return bytes;
   }

   void write_bits8(uint8_t bits, uint8_t bit_count)
   {
      for (uint8_t i = 0; i < bit_count; i++)
      {
         write_bit((bits >> (bit_count-i-1)) & 1);
      }
   }

   void write_bits16(uint16_t bits, uint8_t bit_count)
   {
      for (uint8_t i = 0; i < bit_count; i++)
      {
         write_bit((bits >> (bit_count-i-1)) & 1);
      }
   }

   void write_bits32(uint32_t bits, uint8_t bit_count)
   {
      for (uint8_t i = 0; i < bit_count; i++)
      {
         write_bit((bits >> (bit_count-i-1)) & 1);
      }
   }

   void write_bit(uint8_t bit)
   {
      buffer |= (bit&1) << (7-bit_size);
      bit_size++;
      if (bit_size == 8)
      {
         bytes[bytes_index++] = buffer;
         assert(bytes_index <= max_size);
         bit_size = 0;
         buffer = 0;
      }
   }

   void pad_to_byte(uint8_t pad_bit = 0)
   {
      if (bit_size == 0) return;
      if (pad_bit == 1) buffer |= (1 << (7-bit_size+1)) - 1;
      bytes[bytes_index++] = buffer;
      assert(bytes_index <= max_size);
      bit_size = 0;
      buffer = 0;
   }

public:
   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   uint8_t* bytes = nullptr;
   unsigned bytes_index = 0;
   unsigned max_size = 0;
};

class bit_array_reader_t
{
public:
   bit_array_reader_t(uint8_t* bytes, unsigned size) : bytes(bytes), size(size)
   {
   }

   uint8_t read_bits8(uint8_t bit_count)
   {
      uint8_t bits = 0;
      for (uint8_t i = 0; i < bit_count; i++)
         bits |= (read_bit() << (bit_count-i-1));
      return bits;
   }

   uint16_t read_bits16(uint8_t bit_count)
   {
      uint16_t bits = 0;
      for (uint8_t i = 0; i < bit_count; i++)
         bits |= (read_bit() << (bit_count-i-1));
      return bits;
   }

   uint32_t read_bits32(uint8_t bit_count)
   {
      uint32_t bits = 0;
      for (uint8_t i = 0; i < bit_count; i++)
         bits |= (read_bit() << (bit_count-i-1));
      return bits;
   }

   uint8_t read_bit()
   {
      if (bit_size == 0)
      {
         buffer = bytes[bytes_index++];
         assert(bytes_index <= size);
         bit_size = 8;
      }
      uint8_t b = (buffer >> 7) & 1;
      buffer <<= 1;
      bit_size--;
      return b;
   }

private:
   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   uint8_t* bytes = nullptr;
   unsigned bytes_index = 0;
   unsigned size = 0;
};

#endif