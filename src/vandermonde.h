#pragma once
#include "poly.h"

// Based on LU-decomposition of Vandermonde matrix.
// Time complexity: O(n^2).
vec vandermonde(const inst& in) {
  auto [n, x, y] = in;

  // Backward substitution with L
  vec c(n, 1);
  for (size_t i = 0; i != n; ++i) {
    y[i] /= c[i];
    for (size_t j = i + 1; j != n; ++j) {
      y[j] -= y[i] * c[j];
      c[j] *= x[j] - x[i];
    }
  }

  // Calculate U
  mat u(n, vec(n));
  u[0][0] = 1;
  for (size_t i = 1; i != n; ++i) u[0][i] = u[0][i - 1] * x[0];
  for (size_t i = 1; i != n; ++i) {
    u[i][i] = 1;
    for (size_t j = i + 1; j != n; ++j)
      u[i][j] = u[i - 1][j - 1] + u[i][j - 1] * x[i];
  }

  // Backward substitution with U
  for (size_t i = n; i--;) {
    y[i] /= u[i][i];
    for (size_t j = i; j--;) y[j] -= y[i] * u[j][i];
  }

  return y;
}
