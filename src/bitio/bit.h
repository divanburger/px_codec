#ifndef BIT_H_
#define BIT_H_

class end_of_file_exception_t : std::exception
{
   virtual const char* what() const throw()
   {
      return "end of file";
   }
};

class bit_writer_t
{
public:
   virtual void write_bytes(uint8_t* bytes, unsigned length) = 0;
   virtual void write_byte(uint8_t byte) = 0;
   virtual void write_bits8(uint8_t bits, uint8_t bit_count) = 0;
   virtual void write_bits16(uint16_t bits, uint8_t bit_count) = 0;
   virtual void write_bits32(uint32_t bits, uint8_t bit_count) = 0;
   virtual void write_bit(uint8_t bit) = 0;
   virtual void pad_to_byte(uint8_t pad_bit = 0) = 0;
};

class bit_reader_t
{
public:
   virtual uint8_t read_byte() = 0;
   virtual uint8_t read_bits8(uint8_t bit_count) = 0;
   virtual uint16_t read_bits16(uint8_t bit_count) = 0;
   virtual uint32_t read_bits32(uint8_t bit_count) = 0;
   virtual void skip_to_next_byte() = 0;
   virtual uint8_t read_bit() = 0;
};

#endif