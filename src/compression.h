#ifndef COMPRESSION_H_
#define COMPRESSION_H_

#include "bitarray.h"

void mtf_encode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t symbols[256];
   for (unsigned i = 0; i < 256; i++) symbols[i] = i;

   for (unsigned i = 0; i < n; i++)
   {
      uint8_t in = input[i];

      // Look-up symbol
      uint8_t sym = 0;
      for (unsigned j = 0; j < 256; j++)
         if (symbols[j] == in)
         {
            sym = j;
            break;
         }

      // Output symbol
      *(output++) = sym;

      // Move up symbol table
      for (unsigned j = sym; j > 0; j--) symbols[j] = symbols[j-1];
      symbols[0] = in;
   }
}

void mtf_decode8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t symbols[256];
   for (unsigned i = 0; i < 256; i++) symbols[i] = i;

   for (unsigned i = 0; i < n; i++)
   {
      uint8_t in = input[i];
      uint8_t sym = symbols[in];

      // Output symbol
      *(output++) = sym;

      // Move up symbol table
      for (unsigned j = in; j > 0; j--) symbols[j] = symbols[j-1];
      symbols[0] = sym;
   }
}

void runlength_bit_encode(uint8_t* input, unsigned n, bit_array_writer_t& writer, unsigned run_bits = 3)
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

// Returns the length in bytes of the encoded result
unsigned runlength_bit_encode(uint8_t* input, unsigned n, uint8_t* output, unsigned max_out, unsigned run_bits = 3)
{
   bit_array_writer_t writer(output, max_out);
   runlength_bit_encode(input, n, writer, run_bits);
   writer.pad_to_byte(0);
   return writer.get_size();
}

void runlength_bit_decode(bit_array_reader_t& reader, uint8_t* output, unsigned n, unsigned run_bits = 3)
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

void xor_code8(uint8_t* input, uint8_t* output, unsigned n)
{
   if (n <= 0) return;

   uint8_t last = input[0];
   for (unsigned i = 1; i < n; i++)
   {
      output[i] = last ^ input[i];
      last = input[i];
   }
}

struct huffman_tree_t
{
   huffman_tree_t(unsigned symbols) : size(0), symbols(symbols)
   {
      bits = new uint8_t[symbols*2];
      parents = new unsigned[symbols*2];
   }

   ~huffman_tree_t()
   {
      delete [] bits;
      delete [] parents;
   }

   unsigned  size;
   unsigned  symbols;
   uint8_t*  bits;
   unsigned* parents;  
};

// returns the length of the code
uint8_t get_huffman_code(huffman_tree_t* huffman_tree, uint8_t symbol, uint32_t* code)
{
   uint32_t c = 0;

   uint8_t length = 0;
   unsigned index = symbol;
   while (index != huffman_tree->size - 1)
   {
      uint8_t bit = huffman_tree->bits[index];
      unsigned parent = huffman_tree->parents[index];

      c |= (uint32_t)(bit & 1) << length;
      length++;

      index = parent;
   }

   *code = c;
   return length;
}

void write_huffman_symbol(huffman_tree_t* huffman_tree, bit_array_writer_t& writer, uint8_t symbol)
{
   uint32_t code;
   uint8_t length = get_huffman_code(huffman_tree, symbol, &code);
   writer.write_bits32(code, length);
}

uint8_t read_huffman_symbol(huffman_tree_t* huffman_tree, bit_array_reader_t& reader)
{
   unsigned index = huffman_tree->size - 1;
   while (index >= huffman_tree->symbols)
   {
      uint8_t b = reader.read_bit();

      for (unsigned i = 0; i < huffman_tree->size; i++)
         if (huffman_tree->parents[i] == index && huffman_tree->bits[i] == b)
         {
            index = i;
            break;
         }
   }

   return index;
}

huffman_tree_t* build_huffman_tree(unsigned* frequencies, unsigned n)
{
   // Initialize tree
   const unsigned unassigned = 0xFFFFFFFF;

   huffman_tree_t* huffman_tree = new huffman_tree_t(n);
   unsigned*       freqs = new unsigned[n*2];

   for (unsigned i = 0; i < n; i++)
   {
      huffman_tree->parents[i] = unassigned;
      huffman_tree->bits[i] = 0;
      freqs[i] = frequencies[i];
   }
   huffman_tree->size = n;

   // Create branch nodes
   while (true)
   {
      unsigned lowest = unassigned,      second_lowest = unassigned;
      unsigned lowest_freq = unassigned, second_lowest_freq = unassigned;
      for (unsigned i = 0; i < huffman_tree->size; i++)
      {
         if (huffman_tree->parents[i] != unassigned) continue;
         unsigned freq = freqs[i];
         if (freq < lowest_freq)
         {
            second_lowest_freq = lowest_freq;
            lowest_freq = freq;
            second_lowest = lowest;
            lowest = i;
         }
         else if (freq < second_lowest_freq)
         {
            second_lowest_freq = freq;
            second_lowest = i;
         }
      }

      if (second_lowest == unassigned) break;

      unsigned new_index = huffman_tree->size++;

      huffman_tree->parents[lowest] = new_index;
      huffman_tree->bits[lowest] = 0;
      huffman_tree->parents[second_lowest] = new_index;
      huffman_tree->bits[second_lowest] = 1;

      huffman_tree->parents[new_index] = unassigned;
      huffman_tree->bits[new_index] = 0;
      freqs[new_index] = lowest_freq + second_lowest_freq;
   }

   delete [] freqs;
   return huffman_tree;
}

#endif