#pragma once
#include "poly.h"

// Multipoint evaluation.
// Time complexity: O(n log(n)^2).
vec eval(const vec& a, const vec& x) {
  size_t n = a.size();
  std::vector<vec> tr(n * 2);
  for (size_t i = 0; i != n; ++i) tr[i + n] = {-x[i], 1};
  for (size_t i = n; --i;) tr[i] = tr[i * 2] * tr[i * 2 + 1];
  tr[1] = a % tr[1];
  for (size_t i = 2; i < n * 2; ++i) tr[i] = tr[i / 2] % tr[i];
  vec y(n);
  for (size_t i = 0; i < n; ++i) y[i] = tr[i + n][0];
  return y;
}

// Fast interpolation.
// Time complexity: O(n log(n)^2).
vec fast(const inst& in) {
  const auto& [n, x, y] = in;
  std::vector<std::pair<vec, vec>> tr(n * 2);
  for (size_t i = 0; i != n; ++i) tr[i + n] = {{-x[i], 1}, {1}};
  for (size_t i = n; --i;) {
    auto& [pl, ql] = tr[i * 2];
    auto& [pr, qr] = tr[i * 2 + 1];
    tr[i] = {pl * pr, pl * qr + pr * ql};
  }
  auto coef = eval(tr[1].second, x);
  std::vector<vec> tr2(n * 2);
  for (size_t i = 0; i != n; ++i) tr2[i + n] = {y[i] / coef[i]};
  for (size_t i = n; --i;) {
    tr2[i] =
        tr[i * 2].first * tr2[i * 2 + 1] + tr[i * 2 + 1].first * tr2[i * 2];
  }
  return tr2[1];
}
