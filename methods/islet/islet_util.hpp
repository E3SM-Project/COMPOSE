#ifndef INCLUDE_ISLET_UTIL_HPP
#define INCLUDE_ISLET_UTIL_HPP

#include <cassert>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <memory>

#include "islet_types.hpp"

namespace islet {
#define throw_if(condition, message) do {                               \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

#define require(condition) do {                                         \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)
#define require_msg(condition, message) do {                            \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
           << "\nmessage:\n" << message << "\n";                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

template<typename T> inline T sign (const T& a) { return a >= 0 ? 1 : -1; }
template<typename T> inline constexpr T square (const T& x) { return x*x; }
inline Real reldif (const Real a, const Real b, const Real abstol = 0)
{ return std::abs(b - a)/(abstol + std::abs(a)); }

#define pr(m) do {                              \
    std::stringstream _ss_;                     \
    _ss_ << m << std::endl;                     \
    std::cerr << _ss_.str();                    \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)

template <typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss << name << " = [";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  ss << "]";
  pr(ss.str());
}
template <typename Array>
static void prarr (const std::string& name, const Array& a) {
  prarr(name, a.data(), a.size());
}

#define mprarr(a) prarr(#a, a)

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

inline double urand () { return rand() / ((double) RAND_MAX + 1.0); }

template <typename T>
bool write (std::ofstream& os, const T s) {
  return ! os.write((const char*) &s, sizeof(T)).bad();
}

template <typename T>
bool write (std::ofstream& os, const Int n, const T* const d) {
  return (write(os, n) &&
          ! os.write((const char*) d, n*sizeof(T)).bad());
}

template <typename T>
bool read (std::ifstream& os, T& s) {
  return ! os.read((char*) &s, sizeof(T)).bad();
}

template <typename T>
bool read (std::ifstream& os, Int& n, T* const d) {
  return (read(os, n) &&
          ! os.read((char*) d, n*sizeof(T)).bad());
}

template <typename T> class Array {
  T* p_;
  std::size_t n_, cap_;
public:
  Array () { init(); }
  Array(std::size_t n);
  Array(std::size_t n, const T& init);
  Array(const Array<T>& a);
  ~Array () { clear(); }
  // Initialize the object with the assumption that all variables are uninit'ed
  // prior to calling.
  void init();
  void clear();
  // optclear means optionally clear. The function has the semantics of
  // clearing, but it may not actually release the memory.
  void optclear_and_resize(std::size_t n);
  // _ft indicates first touch.
  void optclear_and_resize_ft(std::size_t n);
  void optclear_and_resize(std::size_t n, const T& i);
  void optclear_and_reserve(std::size_t n);
  void optclear_and_reserve_ft(std::size_t n);
  T& operator[] (std::size_t i) { return p_[i]; }
  const T& operator[] (std::size_t i) const { return p_[i]; }
  T& back () { return p_[n_-1]; }
  const T& back () const { return p_[n_-1]; }
  std::size_t size () const { return n_; }
  bool empty () const { return size() == 0; }
  T* data () const { return p_; }
  // This does not realloc; reserve must provide the necessary memory. It does
  // not throw, either. It asserts.
  void unsafe_push_back(const T& e);
  T* begin () { return p_; }
  T* end () { return p_ + n_; }
  const T* begin () const { return p_; }
  const T* end () const { return p_ + n_; }
  void set (const T& v) { for (std::size_t i = 0; i < n_; ++i) p_[i] = v; }
};

template<typename T> inline int len (const Array<T>& v)
{ return static_cast<int>(v.size()); }

template<typename T> inline void touch (T* const p, const size_t n,
                                        const T& init = T()) {
  // 1 KB should be a safe lower bound on page size. Touch enough to touch every
  // page; I don't think there's any need to touch more memory than that.
  for (size_t i = 0; i < n; i += 1024 / sizeof(T))
    p[i] = init;
  // Make sure the last part is touched.
  if (n) p[n-1] = init;
}
template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) {
  if ( ! n) return 0;
  T* p = new T[n];
  if (first_touch) touch(p, n);
  return p;
}
template<typename T> inline void deln (T*& p) {
  if (p) delete[] p;
  p = 0;
}
template<typename T> inline void deln_const (const T* p) {
  if (p) delete[] p;
}
template<typename T> inline void del (T*& p) {
  if (p) delete p;
  p = 0;
}

template<typename T>
inline void Array<T>::init () {
  n_ = cap_ = 0;
  p_ = 0;
}

template<typename T>
inline Array<T>::Array (std::size_t n)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n); }

template<typename T>
inline Array<T>::Array (std::size_t n, const T& init)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n, init); }

template<typename T>
inline Array<T>::Array (const Array<T>& a) {
  init();
  optclear_and_resize(a.size());
  std::copy(a.begin(), a.end(), begin());
}

template<typename T>
inline void Array<T>::clear () {
  n_ = cap_ = 0;
  deln(p_);
}

template<typename T>
inline void Array<T>::optclear_and_reserve (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_reserve_ft (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n, true);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize_ft (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve_ft(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n, const T& init) {
  optclear_and_resize(n);
  for (std::size_t i = 0; i < n_; ++i)
    memcpy(p_ + i, &init, sizeof(init));
}

template<typename T>
inline void Array<T>::unsafe_push_back (const T& e) {
  assert(n_ < cap_);
  p_[n_++] = e;
}

} // namespace islet

#endif
