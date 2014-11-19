#ifndef COMPRESSION_RUNLENGTH_H_
#define COMPRESSION_RUNLENGTH_H_

#include "../bitio/bit.h"

// -----------------------
// Bit Run-length encoder
// -----------------------
// Assumes 0-bit is much more common than 1-bit
// Only codes runs of 0-bits
//

void runlength_bit_encode(uint8_t* input, unsigned n, bit_writer_t& writer, unsigned run_bits = 3)
{
   const unsigned max_run_length = (1 << run_bits);
   unsigned run = 0;
   for (unsigned i = 0; i < n; i++)
   {
      uint8_t bit = input[i];

      if (bit == 1)
      {
         if (run > 0)
         {
            writer.write_bit(0);
            writer.write_bits8(run-1, run_bits);
            run = 0;
         }
         writer.write_bit(1);
      }
      else
      {
         if (++run == max_run_length)
         {
            writer.write_bit(0);
            writer.write_bits8(run-1, run_bits);
            run = 0;
         }
      }
   }

   if (run > 0)
   {
      writer.write_bit(0);
      writer.write_bits8(run-1, run_bits);
   }
}

// Input must have one-bit values (ie. 0 or 1)
// Returns the length in bytes of the encoded result
unsigned runlength_bit_encode(uint8_t* input, unsigned n, uint8_t* output, unsigned max_out, unsigned run_bits = 3)
{
   bit_array_writer_t writer(output, max_out);
   runlength_bit_encode(input, n, writer, run_bits);
   writer.pad_to_byte(0);
   return writer.get_bytes_written();
}

// Input must have one-bit values (ie. 0 or 1)
void runlength_bit_decode(bit_reader_t& reader, uint8_t* output, unsigned n, unsigned run_bits = 3)
{
   for (unsigned i = 0; i < n;)
   {
      uint8_t b = reader.read_bit();
      if (b == 1)
         output[i++] = 1;
      else
      {
         uint8_t run = reader.read_bits8(run_bits) + 1;
         for (unsigned j = 0; j < run; j++, i++)
            output[i] = 0;
      }
   }
}

// ----------------------------
// Run-length bit-pair encoder
// ----------------------------
unsigned runlength_bitpair_encoder(uint8_t* input, unsigned n, uint8_t* output)
{
   uint8_t last_value = 0;
   unsigned zero_run = 0;
   unsigned one_run = 0;
   unsigned length = 0;
   for (unsigned i = 0; i < n; i++)
   {
      uint8_t b = input[i];

      if (zero_run == 31 || one_run == 7 || (last_value == 1 && b == 0))
      {
         output[length++] = (one_run << 5) + zero_run;
         zero_run = 0;
         one_run = 0;
         last_value = 0;
      }

      if (b == 1)
         one_run++;
      else
         zero_run++;

      last_value = b;
   }

   if (zero_run > 0 || one_run > 0) output[length++] = (one_run << 5) + zero_run;

   return length;
}

// ----------------------
// Run-length Modified AB 
// ----------------------
// Based on bzip's second run-length encoder with RUNA, RUNB
// Input must be 4-bit values, output contains an extra two values
// An entropy encoder is recommended to be applied to the output
//

unsigned runlength_ab_encode4(uint8_t* input, unsigned n, uint8_t* output)
{
   const int RUNA = 16;
   const int RUNB = 17;

   unsigned out_length = 0;
   unsigned run = 0;
   uint8_t run_value = 0xFF;
   for (unsigned i = 0; i < n; i++)
   {
      uint8_t value = input[i];

      if (value == run_value)
      {
         run++;
         continue;
      }

      while (run != 0)
      {
         unsigned b = run & 1;
         run >>= 1;
         if (b == 0)
         {
            output[out_length++] = RUNB;
            run--;
         }
         else
            output[out_length++] = RUNA;
      }

      output[out_length++] = value;
      run_value = value;
      run = 0;
   }

   while (run != 0)
   {
      unsigned b = run & 1;
      run >>= 1;
      if (b == 0)
      {
         output[out_length++] = RUNB;
         run--;
      }
      else
         output[out_length++] = RUNA;
   }

   return out_length;
}

#endif