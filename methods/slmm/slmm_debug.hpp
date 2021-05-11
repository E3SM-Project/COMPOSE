#ifndef INCLUDE_SLMM_DEBUG_HPP
#define INCLUDE_SLMM_DEBUG_HPP

#include <cstdio>
#include <string>

namespace slmm {

template <typename CV3s>
void write_matlab (const std::string& name, const CV3s& p) {
  std::cout << "mat=1; " << name << " = [";
  for (Int ip = 0; ip < nslices(p); ++ip) {
    for (Int k = 0; k < szslice(p); ++k)
      std::cout << " " << p(ip,k);
    std::cout << ";";
  }
  std::cout << "].';\n";
}

template <typename CV3s, typename CIs>
void write_matlab (const std::string& name, const CV3s& p, const CIs& e) {
  printf("mat=1; %s.p = [", name.c_str());
  for (Int ip = 0; ip < nslices(p); ++ip)
    printf(" %1.15e %1.15e %1.15e;", p(ip,0), p(ip,1), p(ip,2));
  printf("].';\n");
  printf("mat=1; %s.n = [", name.c_str());
  for (Int ie = 0; ie < nslices(e); ++ie) {
    for (Int k = 0; k < szslice(e); ++k)
      printf(" %d", e(ie,k)+1);
    printf(";");
  }
  printf("].';\n");
}

} // namespace slmm

#endif
