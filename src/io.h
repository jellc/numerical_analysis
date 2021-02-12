#pragma once
#include <iomanip>
#include <iostream>

#include "def.h"

inst read() {
  size_t n;
  std::cin >> n;
  vec x(n), y(n);
  for (auto&& a : x) std::cin >> a;
  for (auto&& a : y) std::cin >> a;
  return {n, x, y};
}

template <class _Vec> void print(const _Vec& __v) {
  for (const auto& __x : __v) std::cout << __x << ' ';
  std::cout << '\n';
}

struct io_setup {
  io_setup() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout << std::fixed << std::setprecision(17);
  }
} _io_setup;
