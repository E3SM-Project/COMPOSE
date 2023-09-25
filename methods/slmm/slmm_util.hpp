#ifndef INCLUDE_SLMM_UTIL_HPP
#define INCLUDE_SLMM_UTIL_HPP

#include "slmm_defs.hpp"

#include <stdexcept>
#include <fstream>

#define require(condition) do {                                         \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

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
  Real xA, yA, zA;
  Real xB, yB, zB;
  ll2xyz(lat1, lon1, xA, yA, zA);
  ll2xyz(lat2, lon2, xB, yB, zB);
  Real cp1, cp2, cp3, cpnorm, dotprod;
  cp1 = yA*zB - yB*zA;
	cp2 = xB*zA - xA*zB;
	cp3 = xA*yB - xB*yA;
	cpnorm = std::sqrt(cp1*cp1 + cp2*cp2 + cp3*cp3);
	dotprod = xA*xB + yA*yB + zA*zB;
	
	return R * std::atan2(cpnorm, dotprod);
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
void form_rotation(const Real axis[3], const Real angle, Real r[9]);

void fill_normals(const AVec3s& p, const AIdxs& e, AVec3s& nml, AIdxs& en);

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
  return ((a == std::string(b1)) ||
          (b2 && a == std::string(b2)) ||
          (a == std::string("-") + std::string(b1)));
}

std::string& tolower(std::string& s);

std::string format_strings_as_list(const char** strings, const Size n);

double wall_time();

template<typename V> inline Int len (const V& v)
{ return static_cast<Int>(v.extent_int(0)); }

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

inline double urand () { return rand() / ((double) RAND_MAX + 1.0); }

template <typename T>
bool write (std::ofstream& os, const T& s) {
  return ! os.write((const char*) &s, sizeof(T)).bad();
}

template <typename T>
bool write (std::ofstream& os, const Int n, const T* const d) {
  return (write(os, n) &&
          ! os.write((const char*) d, n*sizeof(T)).bad());
}

} // namespace slmm

#endif
