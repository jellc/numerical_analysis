#pragma once
#include <queue>

#include "poly.h"

// Naive implementation of Lagrange interpolation.
// Time complexity: O(n^3).
vec lagrange_naive(const inst &in) {
  const auto &[n, x, y] = in;
  vec f(n - 1);
  for (size_t i = 0; i != n; ++i) {
    vec l{1};
    for (size_t j = 0; j != n; ++j) {
      if (i == j) continue;
      l *= {-x[j], 1};
      l *= {1 / (x[i] - x[j])};
    }
    l *= {y[i]};
    f += l;
  }
  return f;
}

// Faster than above one
// Time complexity: O(n^2 log n).
vec lagrange_fft(const inst &in) {
  const auto &[n, x, y] = in;
  vec f;
  for (size_t i = 0; i != n; ++i) {
    std::queue<vec> Q;
    Q.push({1});
    for (size_t j = 0; j != n; ++j) {
      if (i == j) continue;
      vec p{-x[j], 1};
      p *= {1 / (x[i] - x[j])};
      Q.push(p);
    }
    while (Q.size() > 1) {
      vec p = Q.front();
      Q.pop();
      p *= Q.front();
      Q.pop();
      Q.push(p);
    }
    f += Q.front() * vec{y[i]};
  }
  return f;
}

// Time complexity: O(n^2).
vec lagrange_div(const inst &in) {
  const auto &[n, x, y] = in;
  vec p{1};
  for (size_t i = 0; i != n; ++i) p *= {-x[i], 1};
  vec f;
  for (size_t i = 0; i != n; ++i) {
    real coef = y[i];
    for (size_t j = 0; j != n; ++j)
      if (i != j) coef /= x[i] - x[j];
    f += p / vec{-x[i], 1} * vec{coef};
  }
  return f;
}
