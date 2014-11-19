#ifndef UTIL_H_
#define UTIL_H_

#include <cmath>

using std::min;
using std::max;

template <typename T, typename U, typename V>
auto clamp(T value, U minimum, V maximum) -> decltype(value)
{
   return value < maximum ? (value > minimum ? value : minimum) : maximum;
}

template <typename T>
T sqr(T v)
{
   return v * v;
}

float random(int x, int seed)
{
   int n = x * 1619 + seed * 6971;
   n = (n>>8)^n;
   return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

float random(int x, int y, int seed)
{
   int n = x * 1619 + y * 31337 + seed * 6971;
   n = (n>>8)^n;
   return ((n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

uint8_t toInt(double v)
{
   return (uint8_t)max(min((int)(v * 255.0 + 0.5), 255), 0);
}

uint16_t toInt10(double v)
{
   return (uint16_t)max(min((int)(v * 1023.0 + 0.5), 1023), 0);
}

double fromInt(uint8_t i)
{
   return i / 255.0;
} 

double fromInt10(uint16_t i)
{
   return min((int)i, 1023) / 1023.0;
}

double fromInt12(uint16_t i)
{
   return min((int)i, 4095) / 4095.0;
}

#endif