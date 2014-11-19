#ifndef BITARRAY_H_
#define BITARRAY_H_

#include <cassert>

#include "bit.h"

class bit_array_writer_t : public bit_writer_t
{
public:
   bit_array_writer_t(uint8_t* bytes, unsigned max_size, unsigned offset = 0) : buffer(0), bit_size(0), bytes(bytes), max_size(max_size)
   {
      if (offset > 0) bytes = &bytes[offset];
   }

   unsigned get_bytes_written()
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

   void write_bytes(uint8_t* bytes, unsigned length)
   {
      for (unsigned i = 0; i < length; i++)
         write_byte(bytes[i]);
   }

   void write_byte(uint8_t byte)
   {
      if (bit_size == 0)
      {
         assert(bytes_index < max_size);
         bytes[bytes_index++] = byte;
      }
      else
         write_bits8(byte, 8);
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
         assert(bytes_index < max_size);
         bytes[bytes_index++] = buffer;
         bit_size = 0;
         buffer = 0;
      }
   }

   void pad_to_byte(uint8_t pad_bit = 0)
   {
      if (bit_size == 0) return;
      if (pad_bit == 1) buffer |= (1 << (7-bit_size+1)) - 1;
      assert(bytes_index < max_size);
      bytes[bytes_index++] = buffer;
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

class bit_array_reader_t : public bit_reader_t
{
public:
   bit_array_reader_t(uint8_t* bytes, unsigned size, unsigned offset = 0) : bytes(bytes), size(size)
   {
      if (offset > 0) bytes = &bytes[offset];
   }

   uint8_t read_byte()
   {
      if (bit_size == 0)
      {
         if (bytes_index >= size) throw end_of_file_exception_t();
         buffer = bytes[bytes_index++];
         return buffer;
      }
      return read_bits8(8);
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
         if (bytes_index >= size) throw end_of_file_exception_t();
         buffer = bytes[bytes_index++];
         bit_size = 8;
      }
      uint8_t b = (buffer >> 7) & 1;
      buffer <<= 1;
      bit_size--;
      return b;
   }

   void skip_to_next_byte()
   {
      bit_size = 0;
   }

private:
   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   uint8_t* bytes = nullptr;
   unsigned bytes_index = 0;
   unsigned size = 0;
};

#endif