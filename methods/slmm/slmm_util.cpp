#include "slmm_util.hpp"

#include <cctype>
#include <sys/time.h>
#include <sstream>

namespace slmm {

double wall_time () {
  static const double us = 1.0e6;
  timeval t;
  gettimeofday(&t, 0);
  return (t.tv_sec*us + t.tv_usec)/us;
}

std::string& tolower (std::string& s) {
  for (auto& c: s)
    c = std::tolower(c);
  return s;
}

std::string format_strings_as_list (const char** strings, const Size n) {
  std::stringstream ss;
  ss << "{";
  for (Size i = 0; i < n-1; ++i) ss << strings[i] << ", ";
  ss << strings[n-1] << "}";
  return ss.str();
}

void fill_normals (const AVec3s& p, const AIdxs& e, AVec3s& nml, AIdxs& en) {
  using geo = siqk::SphereGeometry;
  // Count number of edges.
  Int ne = 0;
  for (Int ip = 0; ip < nslices(e); ++ip)
    for (Int iv = 0; iv < szslice(e); ++iv)
      if (e(ip,iv) == -1) break; else ++ne;
  // Fill.
  resize(en, nslices(e), szslice(e));
  copy(en, -1);
  resize(nml, ne);
  Int ie = 0;
  for (Int ip = 0; ip < nslices(e); ++ip)
    for (Int iv = 0; iv < szslice(e); ++iv)
      if (e(ip,iv) == -1)
        break;
      else {
        // Somewhat complicated next node index.
        const Int iv_next = (iv+1 == szslice(e) ? 0 :
                             (e(ip,iv+1) == -1 ? 0 : iv+1));
        geo::edge_normal(slice(p, e(ip, iv)), slice(p, e(ip, iv_next)),
                         slice(nml, ie));
        en(ip,iv) = ie;
        ++ie;
      }
}

void form_rotation(const Real axis[3], const Real angle, Real r[9]) {
  const Real nrm = std::sqrt(square(axis[0]) + square(axis[1]) +
                             square(axis[2]));
  const Real x = axis[0]/nrm, y = axis[1]/nrm, z = axis[2]/nrm, th = angle;
  const Real cth = std::cos(th), sth = std::sin(th), omcth = 1 - cth;
  r[0] = cth + x*x*omcth;
  r[3] = y*x*omcth + z*sth;
  r[6] = z*x*omcth - y*sth;
  r[1] = x*y*omcth - z*sth;
  r[4] = cth + y*y*omcth;
  r[7] = z*y*omcth + x*sth;
  r[2] = x*z*omcth + y*sth;
  r[5] = y*z*omcth - x*sth;
  r[8] = cth + z*z*omcth;
}

} // namespace slmm
