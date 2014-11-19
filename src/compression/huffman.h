#ifndef COMPRESSION_HUFFMAN_H_
#define COMPRESSION_HUFFMAN_H_

#include <algorithm>

#include "../bitio/bit.h"

#include "rice.h"
#include "elias_gamma.h"
#include "compression.h"

using std::sort;

// -------------
// Huffman tree
// -------------

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

   uint8_t get_code(uint8_t symbol, unsigned* code)
   {
      const unsigned unassigned = 0xFFFFFFFF;

      unsigned c = 0;

      uint8_t length = 0;
      unsigned index = symbol;
      while (parents[index] != unassigned)
      {
         uint8_t bit = bits[index];
         unsigned parent = parents[index];

         c |= (unsigned)(bit & 1) << length;
         length++;

         index = parent;
      }

      *code = c;
      return length;
   }

   unsigned  size;
   unsigned  symbols;
   uint8_t*  bits;
   unsigned* parents;  
};

// returns the length of the code
uint8_t get_huffman_code(huffman_tree_t* huffman_tree, uint8_t symbol, unsigned* code)
{
   return huffman_tree->get_code(symbol, code);
}

void write_huffman_symbol(huffman_tree_t* huffman_tree, bit_writer_t& writer, uint8_t symbol)
{
   unsigned code;
   uint8_t length = huffman_tree->get_code(symbol, &code);
   writer.write_bits32(code, length);
}

uint8_t read_huffman_symbol(huffman_tree_t* huffman_tree, bit_reader_t& reader)
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

huffman_tree_t* build_huffman_tree(unsigned long* frequencies, unsigned n)
{
   // Initialize tree
   const unsigned unassigned = 0xFFFFFFFF;

   huffman_tree_t* huffman_tree = new huffman_tree_t(n);
   unsigned long*  freqs = new unsigned long[n*2];

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
      unsigned long lowest = unassigned,      second_lowest = unassigned;
      unsigned long lowest_freq = unassigned, second_lowest_freq = unassigned;
      for (unsigned i = 0; i < huffman_tree->size; i++)
      {
         if (huffman_tree->parents[i] != unassigned) continue;
         unsigned freq = freqs[i];
         if (freq == 0) continue;
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

// -----------------
// Canonical Huffman
// -----------------

struct canonical_huffman_t
{
   struct symbol_t
   {
      unsigned code;
      uint8_t  length;
   };

   canonical_huffman_t(unsigned size) : size(size)
   {
      symbols = new symbol_t[size];
   }

   ~canonical_huffman_t()
   {
      delete [] symbols;
   }

   uint8_t get_code(uint8_t symbol, unsigned* code) const
   {
      *code = symbols[symbol].code;
      return symbols[symbol].length;
   }

   void encode(uint8_t* data, unsigned length, bit_writer_t& writer) const
   {
      for (unsigned i = 0; i < length; i++) write_symbol(writer, data[i]);
   }

   void decode(bit_reader_t& reader, uint8_t* data, unsigned length) const
   {
      for (unsigned i = 0; i < length; i++) data[i] = read_symbol(reader);
   }

   void write_symbol(bit_writer_t& writer, uint8_t symbol) const
   {
      unsigned code;
      uint8_t length = get_code(symbol, &code);
      writer.write_bits32(code, length);
   }

   uint8_t read_symbol(bit_reader_t& reader) const
   {
      unsigned code = 0;
      for (size_t length = 1; length <= 256; length++)
      {
         unsigned b = reader.read_bit();
         code = (code << 1) | b;

         for (unsigned i = 0; i < size; i++)
            if (symbols[i].code == code && symbols[i].length == length)
               return i;
      }

      assert(1 == 2);
      return 0;
   }

   void reconstruct_codes()
   {
      struct entry_t
      {
         uint8_t  symbol;
         unsigned code;
         uint8_t  length;

         bool operator<(const entry_t& other) const
         {
            if (length == other.length) return symbol < other.symbol;
            return length < other.length;
         }
      };

      entry_t* entries = new entry_t[size];

      for (size_t i = 0; i < size; i++)
      {
         entries[i].symbol = i;
         entries[i].length = symbols[i].length;
      }

      sort(&entries[0], &entries[size]);

      unsigned length = 0;
      unsigned code = 0;

      for (size_t i = 0; i < size; i++)
      {
         if (length > 0) code = (code + 1) << (entries[i].length - length);

         uint8_t symbol = entries[i].symbol;
         symbols[symbol].code = code;
         symbols[symbol].length = entries[i].length;
         length = entries[i].length;
      }

      delete [] entries;
   }

   symbol_t* symbols;
   unsigned  size;
};

canonical_huffman_t* build_canonical_huffman_from_tree(huffman_tree_t* huffman_tree)
{
   if (huffman_tree == nullptr || huffman_tree->symbols == 0) return nullptr;
   
   canonical_huffman_t* ch = new canonical_huffman_t(huffman_tree->symbols);

   for (size_t i = 0; i < ch->size; i++)
   {
      unsigned code;
      ch->symbols[i].length = huffman_tree->get_code(i, &code);
   }

   ch->reconstruct_codes();   
   return ch;
}

canonical_huffman_t* build_canonical_huffman(unsigned long* frequencies, unsigned n)
{
   huffman_tree_t* huffman_tree = build_huffman_tree(frequencies, n);
   canonical_huffman_t* ch = new canonical_huffman_t(huffman_tree->symbols);

   for (size_t i = 0; i < ch->size; i++)
   {
      unsigned code;
      ch->symbols[i].length = huffman_tree->get_code(i, &code);
   }

   delete huffman_tree;

   ch->reconstruct_codes();   
   return ch;
}

void write_canonical_huffman(canonical_huffman_t* ch, bit_writer_t& writer)
{
   assert(ch->size <= 256);
   writer.write_byte(ch->size-1);

   uint8_t lengths[256];
   uint8_t encoded[256];
   for (unsigned i = 0; i < ch->size; i++) lengths[i] = ch->symbols[i].length;
   mtf_encode8(lengths, encoded, ch->size);

   for (unsigned i = 0; i < ch->size; i++)
      write_elias_gamma_code(encoded[i], writer);
}

canonical_huffman_t* read_canonical_huffman(bit_reader_t& reader)
{
   unsigned size = reader.read_byte()+1;

   canonical_huffman_t* ch = new canonical_huffman_t(size);

   uint8_t encoded[256];
   for (unsigned i = 0; i < ch->size; i++)
      encoded[i] = read_elias_gamma_code(reader);

   uint8_t lengths[256];
   mtf_decode8(encoded, lengths, ch->size);
   for (unsigned i = 0; i < ch->size; i++) ch->symbols[i].length = lengths[i];

   ch->reconstruct_codes();
   return ch;
}

#endif