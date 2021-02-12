#pragma once

#include "fft.h"

#define FFT_THRESHOLD 100

vec& operator+=(vec& p, const vec& q) {
  if (p.size() < q.size()) p.resize(q.size());
  for (size_t i = 0; i != q.size(); ++i) p[i] += q[i];
  return p;
}

vec operator+(vec p, const vec& q) { return p += q; }

vec operator-(vec p) {
  for (auto&& x : p) {
    x = -x;
  }
  return p;
}

vec& operator-=(vec& p, const vec& q) { return p += -q; }

vec operator-(vec p, const vec& q) { return p -= q; }

vec& operator*=(vec& p, vec q) {
  if (p.size() < q.size()) p.swap(q);
  if (q.empty()) return p = {};

  if (q.size() < FFT_THRESHOLD) {  // naive convolution
    p.resize(p.size() + q.size() - 1);
    for (size_t i = p.size(); i--;) {
      p[i] *= q[0];
      for (size_t j = 1; j <= i && j != q.size(); ++j) p[i] += p[i - j] * q[j];
    }
    return p;
  }

  return p = conv(p, q);
}

vec operator*(vec p, const vec& q) { return p *= q; }

vec pref(const vec& p, size_t n) { return vec(p.begin(), p.begin() + n); }

vec rev(vec p, size_t n) {
  p.resize(n);
  std::reverse(p.begin(), p.end());
  return p;
}

vec inv(const vec& p, size_t n) {
  vec ip{1 / p[0]};
  for (size_t i = 1; i < n;) {
    i <<= 1;
    ip *= pref(-ip * pref(p, i) + vec{2}, i);
  }
  return pref(p, n);
}

vec& operator/=(vec& p, const vec& q) {
  if (p.size() < q.size()) {
    p.clear();
    return p;
  }

  size_t n = p.size() - q.size() + 1;

  if (q.size() < FFT_THRESHOLD) {
    for (size_t i = p.size(); i-- >= q.size();) {
      p[i] /= q.back();
      for (size_t j = 1; j != q.size(); ++j) {
        p[i - j] -= p[i] * q[q.size() - 1 - j];
      }
    }
    p.erase(p.begin(), p.end() - n);
    return p;
  }

  return p = rev(pref(pref(rev(p, n), n) * inv(rev(p, n), n), n), n);
}

vec operator/(vec p, const vec& q) { return p /= q; }

vec& operator%=(vec& p, const vec& q) { return p -= p / q * q; }

vec operator%(vec p, const vec& q) { return p %= q; }
