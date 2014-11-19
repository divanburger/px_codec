#ifndef COMPRESSION_ARITHCODER_H_
#define COMPRESSION_ARITHCODER_H_

#include "../bitio/bit.h"

#include "rice.h"

struct arithmetic_coder_t
{
   arithmetic_coder_t(unsigned size) : size(size)
   {
      ranges = new uint32_t[size+1];
   }

   ~arithmetic_coder_t()
   {
      delete [] ranges;
   }

   void encode(uint8_t* data, unsigned length, bit_writer_t& writer) const
   {
      uint32_t low = 0;
      uint32_t high = 0xFFFF;
      uint32_t underflow = 0;
      uint32_t wrote_bits = 0;

      uint8_t highest_prob = ranges[size];
      for (unsigned i = 0; i < length; i++)
      {
         uint8_t symbol = data[i];

         // Apply range
         uint32_t range = high - low + 1;
         high = low + (range * ranges[symbol+1]) / highest_prob - 1;
         low = low + (range * ranges[symbol]) / highest_prob;
         assert(high >= low);

         // Write out bits
         while (true)
         {
            if ((low & 0x8000) == (high & 0x8000))
            {
               writer.write_bit((low & 0x8000) >> 15);
               wrote_bits++;

               // Write underflow
               while (underflow > 0)
               {
                  writer.write_bit(((low & 0x8000) >> 15) ^ 0x1);
                  wrote_bits++;
                  underflow--;
               }
               underflow = 0;
            }
            else if (((low & 0x4000) == 0x4000) && ((high & 0x4000) != 0x4000))
            {
               // Underflow
               underflow++;
               low &= 0x3FFF;
               high |= 0x4000;
            }
            else
               break;

            low = (low << 1) & 0xFFFE;
            high = ((high << 1) & 0xFFFE) | 1;
         }
      }

      writer.write_bit((low & 0x4000) >> 15);
      wrote_bits++;

      if ((((low & 0x8000) >> 15) ^ 0x1) == 0)
      {
         underflow++;
         while (underflow > 0)
         {
            writer.write_bit(((low & 0x8000) >> 15) ^ 0x1);
            wrote_bits++;
            underflow--;
         }
      }
   }

   void decode(bit_reader_t& reader, uint8_t* data, unsigned length) const
   {
      uint32_t low = 0;
      uint32_t high = 0xFFFF;
      uint32_t index = 0;
      uint32_t code = 0;
      uint32_t read_bits = 0;

      // Read in initial bits
      for (unsigned i = 0; i < 16; i++)
      {
         code <<= 1;
         code |= reader.read_bit() & 1;
         read_bits++;
      }

      uint8_t highest_prob = ranges[size];
      for (unsigned i = 0; i < length;)
      {
         uint32_t range = high - low + 1;
         uint32_t prob = code - low + 1;
         prob = (prob * highest_prob - 1) / range;

         // Find symbol
         uint32_t symbol = ~0U;
         for (unsigned j = 0; j < size; j++)
            if (prob < ranges[j+1])
            {
               symbol = j;
               break;
            }

         if (symbol == ~0U)
         {
            printf("Error!!!\n");
            break;
         }
         //printf("sym: %3i\n", symbol);

         // Output symbol
         data[index++] = symbol;
         if (index >= length) break;

         // Apply range
         high = low + (range * ranges[symbol+1]) / highest_prob - 1;
         low = low + (range * ranges[symbol]) / highest_prob;
         assert(high >= low);

         // Read more bits
         while (true)
         {
            if ((low & 0x8000) == (high & 0x8000))
            {
               // Ignore MSB
            }
            else if (((low & 0x4000) == 0x4000) && ((high & 0x4000) != 0x4000))
            {
               // Underflow
               low &= 0x3FFF;
               high |= 0x4000;
               code ^= 0x4000;
            }
            else
               break;

            low = (low << 1) & 0xFFFE;
            high = ((high << 1) & 0xFFFE) | 1;
            code = (code << 1) & 0xFFFE;

            // Read next bit
            code |= reader.read_bit();
            read_bits++;
         }
      }
   }

   unsigned    size;
   uint32_t*   ranges;
};

void write_arithmetic_coder(arithmetic_coder_t* coder, bit_writer_t& writer)
{
   writer.write_byte(coder->size-1);

   writer.write_byte(coder->ranges[1]);
   unsigned last = coder->ranges[1];
   for (unsigned i = 1; i < coder->size; i++)
   {
      unsigned range = coder->ranges[i+1];
      write_rice2_code(range - last, writer);
      last = range;
   }
}

arithmetic_coder_t* read_arithmetic_coder(bit_reader_t& reader)
{
   unsigned size = reader.read_byte() + 1;
   arithmetic_coder_t* coder = new arithmetic_coder_t(size);

   coder->ranges[0] = 0;
   coder->ranges[1] = reader.read_byte();
   unsigned last = coder->ranges[1];
   for (unsigned i = 1; i < coder->size; i++)
   {
      unsigned range = last + read_rice2_code(reader);
      coder->ranges[i+1] = range;
      last = range;
   }

   return coder;
}

arithmetic_coder_t* build_arithmetic_coder(uint64_t* counts, unsigned size)
{
   const unsigned max_total = (1 << (8 - 2));

   arithmetic_coder_t* coder = new arithmetic_coder_t(size);

   uint64_t total = 0;
   for (unsigned i = 0; i < size; i++)
      total += counts[i];

   uint64_t scale_factor = total < max_total ? 1 : total / max_total + 1;
   unsigned range = 0;
   coder->ranges[0] = 0;
   for (unsigned i = 0; i < size; i++)
   {
      if (counts[i] != 0)
         range += std::max(counts[i] / scale_factor, (uint64_t)1U);
      coder->ranges[i+1] = range;
   }

   return coder;
}

#endif