#ifndef BITVIEW_H_
#define BITVIEW_H_

#include <cassert>

#include "bit.h"

class bit_view_writer_t : public bit_writer_t
{
public:
   bit_view_writer_t(bit_writer_t* sink, unsigned max_size) : sink(sink), buffer(0), bit_size(0), max_size(max_size)
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
         bytes_index++;
         sink->write_byte(byte);
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
         bytes_index++;
         sink->write_byte(buffer);
         bit_size = 0;
         buffer = 0;
      }
   }

   void pad_to_byte(uint8_t pad_bit = 0)
   {
      if (bit_size == 0) return;
      if (pad_bit == 1) buffer |= (1 << (7-bit_size+1)) - 1;
      assert(bytes_index < max_size);
      bytes_index++;
      sink->write_byte(buffer);
      bit_size = 0;
      buffer = 0;
   }

public:
   bit_writer_t* sink;

   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   unsigned bytes_index = 0;
   unsigned max_size = 0;
};

class bit_view_reader_t : public bit_reader_t
{
public:
   bit_view_reader_t(bit_reader_t* source, unsigned size) : source(source), size(size)
   {
   }

   uint8_t read_byte()
   {
      if (bit_size == 0)
      {
         if (bytes_index < size)
         {
            bytes_index++;
            buffer = source->read_byte();
         }
         else
            buffer = 0;
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
         if (bytes_index < size)
         {
            bytes_index++;
            buffer = source->read_byte();
         }
         else
            buffer = 0;
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
   bit_reader_t* source;

   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   unsigned bytes_index = 0;
   unsigned size = 0;
};

#endif