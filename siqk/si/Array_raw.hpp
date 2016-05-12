#ifndef INCLUDE_ARRAY_RAW_HPP
#define INCLUDE_ARRAY_RAW_HPP

#include <memory>
#include <type_traits>
#include <stdexcept>
#include <sstream>

static inline void error(const std::string& msg)
{ throw std::runtime_error(msg.c_str()); }

template <typename T> static inline void share_nodelete_delete (T* p) {}
template <typename T> inline std::shared_ptr<T> share_nodelete (T* o)
{ return std::shared_ptr<T>(o, share_nodelete_delete<T>); }

template <typename T>
class Array1D {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array1D<const T_nonconst>;
  int n_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array1D () : n_(0) {}
  Array1D (const int n) { reset(n); }
  Array1D (const int n, T* const a) { reset(n, a); }
  Array1D (const Array1D<T_nonconst>& v)
    : n_(v.n_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int n) {
    n_ = n;
    a_p_ = std::shared_ptr<T>(new T[n], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int n, T* const a) { n_ = n; a_p_ = nullptr; a_ = a; }
  const int& n () const { return n_; }
  T* data () { return a_; }
  const T* data () const { return a_; }
  T& operator[] (const int i) { debug(i); return a_[i]; }
  const T& operator[] (const int i) const { debug(i); return a_[i]; }
  void set (const T& init) { for (int i = 0; i < n_; ++i) a_[i] = init; }
  Array1D<T>& device () { return *this; }
  const Array1D<T>& device () const { return *this; }
  void sync () {}
  void modify () {}
private:
#ifdef SIQK_DEBUG
  void debug (const int& i) const {
    if (i < 0 || i >= m_) {
      std::stringstream ss;
      ss << "Array1D: i is " << i << " but n_ is " << n_ << "\n";
      error(ss.str().c_str());
    }
  }
#else
  static void debug (const int& i) {}
#endif
};

template <typename T>
class Array2D {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array2D<const T_nonconst>;
  int m_, n_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array2D () : m_(0), n_(0) {}
  Array2D (const int m, const int n) { reset(m, n); }
  Array2D (const int m, const int n, T* const a) { reset(m, n, a); }
  Array2D (const Array2D<T_nonconst>& v)
    : m_(v.m_), n_(v.n_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int m, const int n) {
    m_ = m; n_ = n;
    a_p_ = std::shared_ptr<T>(new T[m*n], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int m, const int n, T* const a) { m_ = m; n_ = n; a_p_ = nullptr; a_ = a; }
  const int& m () const { return m_; }
  const int& n () const { return n_; }
  T* data () { return a_; }
  const T* data () const { return a_; }
  T& operator() (const int r, const int c) { debug(r, c); return a_[c*m_ + r]; }
  const T& operator() (const int r, const int c) const { debug(r, c); return a_[c*m_ + r]; }
  T* operator() (const int c) { debug(0, c); return a_ + m_*c; }
  const T* operator() (const int c) const { debug(0, c); return a_ + m_*c; }
  void set (const T& init) { for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
  Array2D<T>& device () { return *this; }
  const Array2D<T>& device () const { return *this; }
  void sync () {}
  void modify () {}
private:
#ifdef SIQK_DEBUG
  void debug (const int& r, const int& c) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array2D: r is " << r << " but m_ is " << m_ << "\n";
      error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array2D: c is " << c << " but n_ is " << n_ << "\n";
      error(ss.str().c_str());
    }
  }
#else
  static void debug (const int& r, const int& c) {}
#endif
};

// Define a few things to minimize KOKKOS guards.
# ifndef KOKKOS_FUNCTION
#  define KOKKOS_FUNCTION
# endif
# ifndef KOKKOS_INLINE_FUNCTION
#  define KOKKOS_INLINE_FUNCTION inline
# endif
# ifndef KOKKOS_FORCEINLINE_FUNCTION
#  define KOKKOS_FORCEINLINE_FUNCTION inline
# endif

namespace Kokkos {
typedef void DefaultExecutionSpace;
inline void fence() {}
}

namespace ko {
using std::min;
using std::max;

template <typename Functor, typename Scalar>
void parallel_reduce (const int n, Functor f, Scalar& r) {
  for (int i = 0; i < n; ++i)
    f(i, r);
}
}

#endif // INCLUDE_ARRAY_RAW_HPP
