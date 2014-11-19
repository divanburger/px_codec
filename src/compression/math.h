#ifndef COMPRESSION_MATH_H_
#define COMPRESSION_MATH_H_

uint16_t imin(uint16_t a, uint16_t b)
{
   return a < b ? a : b;
}

uint16_t imax(uint16_t a, uint16_t b)
{
   return a > b ? a : b;
}

int imin(int a, int b)
{
   return a < b ? a : b;
}

int imax(int a, int b)
{
   return a > b ? a : b;
}

#endif
