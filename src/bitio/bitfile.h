#ifndef BITFILE_H_
#define BITFILE_H_

#include <cassert>
#include <cstdio>
#include <string>

#include "bit.h"

using std::string;

class bit_file_writer_t : public bit_writer_t
{
public:
   bit_file_writer_t(string filename) : buffer(0), bit_size(0)
   {
      file = fopen(filename.c_str(), "wb");
   }

   ~bit_file_writer_t()
   {
      fclose(file);
   }

   unsigned get_bytes_written()
   {
      return bytes_written;
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
         bytes_written++;
         size_t res = fwrite(&byte, 1, 1, file);
         assert(res == 1);
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
         bytes_written++;
         size_t res = fwrite(&buffer, 1, 1, file);
         assert(res == 1);
         bit_size = 0;
         buffer = 0;
      }
   }

   void pad_to_byte(uint8_t pad_bit = 0)
   {
      if (bit_size == 0) return;
      if (pad_bit == 1) buffer |= (1 << (7-bit_size+1)) - 1;
      bytes_written++;
      size_t res = fwrite(&buffer, 1, 1, file);
      assert(res == 1);
      bit_size = 0;
      buffer = 0;
   }

public:
   uint8_t  buffer = 0;
   uint8_t  bit_size = 0;

   FILE*    file = nullptr;
   unsigned bytes_written = 0;
};

class bit_file_reader_t : public bit_reader_t
{
public:
   bit_file_reader_t(string filename)
   {
      file = fopen(filename.c_str(), "rb");
   }

   ~bit_file_reader_t()
   {
      fclose(file);
   }

   bool is_open()
   {
      return file != nullptr;
   }

   unsigned get_bytes_read()
   {
      return bytes_read;
   }

   uint8_t read_byte()
   {
      if (bit_size == 0)
      {
         size_t res = fread(&buffer, 1, 1, file);
         bytes_read++;
         assert(res == 1);
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
         size_t res = fread(&buffer, 1, 1, file);
         assert(res == 1);
         bytes_read++;
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

   FILE*    file = nullptr;
   unsigned bytes_read = 0;
};

#endif