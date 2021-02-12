#pragma once

#include <vector>

// As real numbers
using real = long double;

// As vector or array
using vec = std::vector<real>;

// As matrix
using mat = std::vector<vec>;

// Instance class
struct inst {
  size_t n;
  vec x, y;
};
