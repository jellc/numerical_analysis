#pragma once

#include <algorithm>
#include <cmath>

#include "def.h"

struct cplx {
  real re, im;

  constexpr cplx() : re{real{}}, im{real{}} {}
  constexpr cplx(real _re) : re{_re}, im{real{}} {}
  constexpr cplx(real _re, real _im) : re{_re}, im{_im} {}

  constexpr real norm() const { return re * re + im * im; }
  constexpr cplx conj() const { return cplx(re, -im); }
  constexpr cplx operator-() const { return cplx(-re, -im); }

  constexpr cplx &operator+=(const cplx &x) {
    return re += x.re, im += x.im, *this;
  }
  constexpr cplx &operator-=(const cplx &x) { return *this += -x; }
  constexpr cplx &operator*=(const cplx &x) {
    real _re{re * x.re - im * x.im};
    return im = im * x.re + x.im * re, re = _re, *this;
  }
  constexpr cplx &operator*=(real x) { return re *= x, im *= x, *this; }
  constexpr cplx &operator/=(const cplx &x) {
    return (*this *= x.conj()) /= x.re * x.re + x.im * x.im;
  }
  constexpr cplx &operator/=(real x) { return re /= x, im /= x, *this; }
  constexpr cplx operator+(const cplx &x) const { return cplx(*this) += x; }
  constexpr cplx operator-(const cplx &x) const { return cplx(*this) -= x; }
  constexpr cplx operator*(const cplx &x) const { return cplx(*this) *= x; }
  constexpr cplx operator*(real x) const { return cplx(*this) *= x; }
  constexpr cplx operator/(const cplx &x) const { return cplx(*this) /= x; }
  constexpr cplx operator/(real x) const { return cplx(*this) /= x; }
};

using vec_cplx = std::vector<cplx>;

void fft(vec_cplx &f) {
  constexpr cplx zeta[31] = {
      {1, 0},
      {-1, 0},
      {0, 1},
      {0.70710678118654752438189403651, 0.70710678118654752443610414514},
      {0.92387953251128675610142140795, 0.38268343236508977172325753068},
      {0.98078528040323044911909938781, 0.19509032201612826785692544201},
      {0.99518472667219688623102546998, 0.09801714032956060199569840382},
      {0.99879545620517239270077028412, 0.04906767432741801425693899119},
      {0.99969881869620422009748220149, 0.02454122852291228803212346128},
      {0.99992470183914454093764001552, 0.01227153828571992607945510345},
      {0.99998117528260114264494415325, 0.00613588464915447535972750246},
      {0.99999529380957617150137498041, 0.00306795676296597627029751672},
      {0.99999882345170190993313003025, 0.00153398018628476561237225788},
      {0.99999970586288221914474799723, 0.00076699031874270452695124765},
      {0.99999992646571785113833452651, 0.00038349518757139558906815188},
      {0.99999998161642929381167504976, 0.00019174759731070330743679009},
      {0.99999999540410731290905263501, 0.00009587379909597734587360460},
      {0.99999999885102682753608427379, 0.00004793689960306688454884772},
      {0.99999999971275670682981095982, 0.00002396844980841821872882467},
      {0.99999999992818917670745273995, 0.00001198422490506970642183282},
      {0.99999999998204729416331065783, 0.00000599211245264242784278378},
      {0.99999999999551182356793271877, 0.00000299605622633466075058210},
      {0.99999999999887795586487812538, 0.00000149802811316901122883643},
      {0.99999999999971948897977205850, 0.00000074901405658471572113723},
      {0.99999999999992987223139048746, 0.00000037450702829238412391495},
      {0.99999999999998246807140014902, 0.00000018725351414619534486931},
      {0.99999999999999561700429751010, 0.00000009362675707309808280024},
      {0.99999999999999890425107437752, 0.00000004681337853654909269501},
      {0.99999999999999972607632112153, 0.00000002340668926827455275977},
      {0.99999999999999993153263280754, 0.00000001170334463413727718121},
      {0.99999999999999998286960567472, 0.00000000585167231706863869077}};

  const size_t n = f.size(), mask{n - 1};
  vec_cplx g(n);

  for (size_t i{n >> 1}, ii{1}; i; i >>= 1, ++ii, f.swap(g)) {
    cplx powzeta{1};
    for (size_t j{}; j != n; powzeta *= zeta[ii])
      for (size_t k{}, x{mask & j << 1}, y{mask & (i + (j << 1))}; k != i;
           ++k, ++j, ++x, ++y)
        g[j] = powzeta * f[y] + f[x];
  }
}

void ifft(vec_cplx &f) {
  fft(f);
  std::reverse(std::next(f.begin()), f.end());
  for (cplx &e : f) e /= f.size();
}

vec conv(const vec &f, const vec &g) {
  size_t len = 1;
  while (len + 1 < f.size() + g.size()) len <<= 1;

  vec_cplx p(len);
  for (size_t i = 0; i != f.size(); ++i) p[i].re = f[i];
  for (size_t i = 0; i != g.size(); ++i) p[i].im = g[i];
  fft(p);

  vec_cplx q(len);
  q[0] = p[0].norm() * 4;
  for (size_t i = 1; i != len; ++i)
    q[i] = (p[i] + q[len - i].conj()) * (p[i] - q[len - i].conj());
  ifft(q);

  vec fg(f.size() + g.size() - 1);
  for (size_t i = 0; i != fg.size(); ++i) {
    fg[i] = round(q[i].im / 4);
  }
  return fg;
}
