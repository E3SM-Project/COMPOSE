#ifndef INCLUDE_SLMM_UTIL_HPP
#define INCLUDE_SLMM_UTIL_HPP

#include "slmm_defs.hpp"

#include <stdexcept>

namespace slmm {
using siqk::square;
template<typename T> inline constexpr T cube (const T& x) { return x*x*x; }

struct consts {
  static constexpr Real earth_radius_m = 6.37122e6;
};

template<typename T> inline T sign (const T& a) { return a >= 0 ? 1 : -1; }

inline Real sec2day (const Real sec) { return sec/(24*3600); }
inline Real day2sec (const Real day) { return day*(24*3600); }

// Output is in radians.
//todo Make a version that lets you pass R = mag(x,y,z).
inline void xyz2ll (const Real x, const Real y, const Real z,
                    Real& lat, Real& lon) {
  const Real r = std::sqrt(square(x) + square(y) + square(z));
  lat = std::asin(z/r);
  lon = std::atan2(y, x);
}

// Input is in radians.
inline void ll2xyz (const Real lat, const Real lon, Real& x, Real& y, Real& z,
                    const Real radius = 1) {
  const Real sinl = std::sin(lat), cosl = std::cos(lat);
  x = radius*std::cos(lon)*cosl;
  y = radius*std::sin(lon)*cosl;
  z = radius*sinl;
}

// Eq after eq 10 in Lauritzen et al test cases paper.
inline Real great_circle_dist (
  const Real lat1, const Real lon1, const Real lat2, const Real lon2,
  const Real R = 1)
{
  return R*std::acos(std::sin(lat1)*std::sin(lat2) +
                     std::cos(lat1)*std::cos(lat2)*std::cos(lon1 - lon2));
}

inline constexpr Real m2radlat (const Real m)
{ return m/consts::earth_radius_m; }

inline Real m2radlon(const Real lat, const Real m)
{ return m2radlat(m)/std::abs(std::cos(lat)); }

inline constexpr Real deg2rad (const Real v) { return v * (M_PI/180); }
inline constexpr Real rad2deg (const Real v) { return v * (180/M_PI); }

inline Real reldif (const Real a, const Real b, const Real abstol = 0)
{ return std::abs(b - a)/(abstol + std::abs(a)); }

// Row-major R.
inline void form_rotation (const Real axis[3], const Real angle, Real r[9]) {
  const Real nrm = std::sqrt(square(axis[0]) + square(axis[1]) +
                             square(axis[2]));
  const Real& x = axis[0] / nrm, & y = axis[1] / nrm, & z = axis[2] / nrm,
    & th = angle;
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

/*! \brief RAII std stream state saver.
 *
 * Example: Preserve std::cout's state so manipulations don't affect others' use
 * of cout.
 */
template<typename Stream> class IosSaver {
  Stream& s_;
  std::ios state_;
public:
  IosSaver (Stream& s) : s_(s), state_(nullptr) { state_.copyfmt(s); }
  IosSaver (const IosSaver& ios) : s_(ios.s_), state_(nullptr)
  { state_.copyfmt(ios.state_); }
  IosSaver operator= (const IosSaver&) = delete;
  ~IosSaver () { s_.copyfmt(state_); }
};
template<typename Stream> inline IosSaver<Stream> save_ios (Stream& s)
{ return IosSaver<Stream>(s); }

template<typename T>
inline T* tin (T* const p, const char* const msg="") {
  if ( ! p)
    throw std::runtime_error(std::string(std::string("Null pointer: ") + msg));
  return p;
}

inline bool
eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

std::string& tolower(std::string& s);

std::string format_strings_as_list(const char** strings, const Size n);

double wall_time();

template<typename V> inline Int len (const V& v)
{ return static_cast<Int>(v.dimension_0()); }

template<typename T> inline Int len (const std::vector<T>& v)
{ return static_cast<Int>(v.size()); }

class ProgressBar {
  std::string name_;
  const Int nits_; // total # iterations
  const Real wf_;  // write frequency in percentage points
  Int it_;
  Real next_;
  std::ostream& os_;

public:
  ProgressBar (const std::string& name, const Int niterations,
               const Real write_freq = 1.0, std::ostream& os = std::cout)
    : name_(name), nits_(niterations), wf_(write_freq), it_(0), next_(0),
      os_(os)
  {
    os_ << name_ << ":";
    os_.flush();
  }

  void update () {
    ++it_;
    const Real p = 100 * it_ / nits_;
    if (p >= next_ || it_ == nits_) {
      os_ << " " << p;
      if (it_ == nits_) os_ << "\n";
      os_.flush();
      next_ += wf_;
    }
  }
};

} // namespace slmm

#endif
