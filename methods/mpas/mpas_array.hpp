#ifndef MPAS_INCLUDE_ARRAY_HPP
#define MPAS_INCLUDE_ARRAY_HPP

#include "mpas_macros.hpp"
#include <sstream>
#include <memory>

template <typename T>
class Array1D {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array1D<const T_nonconst>;
  size_t m_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array1D () : m_(0) {}
  Array1D (const size_t m) { reset(m); }
  Array1D (const size_t m, T* const a) { reset(m, a); }
  Array1D (const Array1D<T_nonconst>& v)
    : m_(v.m_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const size_t m) {
    m_ = m;
    assert(size() >= 0);
    a_p_ = std::shared_ptr<T>(new T[m], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const size_t m, T* const a) {
    m_ = m; a_p_ = nullptr; a_ = a;
    assert(size() >= 0);
  }
  int m () const { return static_cast<int>(m_); }
  size_t size () const { return m_; }
  bool empty () const { return size() == 0; }
  T* data () const { return a_; }
  T& operator[] (const int i) const { debug(i); return a_[i]; }
  T& operator() (const int i) const { debug(i); return a_[i]; }
  void set (const T& init) const { for (size_t i = 0; i < m_; ++i) a_[i] = init; }
  void set_ompfor (const T& init) const
  { ompfor for (size_t i = 0; i < m_; ++i) a_[i] = init; }
  void set_ompparfor (const T& init) const
  { ompparfor for (size_t i = 0; i < m_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const size_t i) const {
    if (i < 0 || i >= m_) {
      std::stringstream ss;
      ss << "Array1D: i is " << i << " but m_ is " << m_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int) {}
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
    assert(size() >= 0);
    a_p_ = std::shared_ptr<T>(new T[m*n], std::default_delete<T[]>());
    a_ = a_p_.get();
    assert(size() >= 0);
  }
  void reset (const int m, const int n, T* const a) {
    m_ = m; n_ = n; a_p_ = nullptr; a_ = a;
    assert(size() >= 0);
  }
  int m () const { return m_; }
  int n () const { return n_; }
  size_t size () const { return size_t(n())*m(); }
  bool empty () const { return size() == 0; }
  T* data () const { return a_; }
  T& operator() (const int r, const int c) const { debug(r, c); return a_[c*m_ + r]; }
  T* operator() (const int c) const { debug(0, c); return a_ + m_*c; }
  void set (const T& init) const { for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
  void set_ompfor (const T& init) const
  { ompfor for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
  void set_ompparfor (const T& init) const
  { ompparfor for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int r, const int c) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array2D: r is " << r << " but m_ is " << m_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array2D: c is " << c << " but n_ is " << n_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int, const int) {}
#endif
};

template <typename T>
class Array3D {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array3D<const T_nonconst>;
  int m_, n_, o_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array3D () : m_(0), n_(0), o_(0) {}
  Array3D (const int m, const int n, const int o) { reset(m, n, o); }
  Array3D (const int m, const int n, const int o, T* const a) { reset(m, n, o, a); }
  Array3D (const Array3D<T_nonconst>& v)
    : m_(v.m_), n_(v.n_), o_(v.o_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int m, const int n, const int o) {
    m_ = m; n_ = n; o_ = o;
    assert(size() >= 0);
    a_p_ = std::shared_ptr<T>(new T[size_t(m)*n*o], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int m, const int n, const int o, T* const a) {
    m_ = m; n_ = n; o_ = o; a_p_ = nullptr; a_ = a;
    assert(size() >= 0);
  }
  int m () const { return m_; }
  int n () const { return n_; }
  int o () const { return o_; }
  size_t size () const { return size_t(n())*m()*o(); }
  bool empty () const { return size() == 0; }
  T* data () const { return a_; }
  T& operator() (const int r, const int c, const int s) const
  { debug(r, c, s); return a_[(s*n_ + c)*m_ + r]; }
  T* operator() (const int c, const int s) const
  { debug(0, c, s); return a_ + (s*n_ + c)*m_; }
  void set (const T& init) const
  { for (size_t i = 0; i < size_t(m_)*n_*o_; ++i) a_[i] = init; }
  void set_ompfor (const T& init) const
  { ompfor for (size_t i = 0; i < size_t(m_)*n_*o_; ++i) a_[i] = init; }
  void set_ompparfor (const T& init) const
  { ompparfor for (size_t i = 0; i < size_t(m_)*n_*o_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int r, const int c, const int s) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array3D: r is " << r << " but m_ is " << m_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array3D: c is " << c << " but n_ is " << n_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (s < 0 || s >= o_) {
      std::stringstream ss;
      ss << "Array3D: s is " << s << " but o_ is " << o_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int, const int, const int) {}
#endif
};

template <typename T>
class Array4D {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array4D<const T_nonconst>;
  int m_, n_, o_, p_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array4D () : m_(0), n_(0), o_(0), p_(0) {}
  Array4D (const int m, const int n, const int o, const int p)
  { reset(m, n, o, p); }
  Array4D (const int m, const int n, const int o, const int  p, T* const a)
  { reset(m, n, o, p, a); }
  Array4D (const Array4D<T_nonconst>& v)
    : m_(v.m_), n_(v.n_), o_(v.o_), p_(v.p_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int m, const int n, const int o, const int p) {
    m_ = m; n_ = n; o_ = o; p_ = p;
    assert(size() >= 0);
    a_p_ = std::shared_ptr<T>(new T[size_t(m)*n*o*p], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int m, const int n, const int o, const int p, T* const a) {
    m_ = m; n_ = n; o_ = o; p_ = p; a_p_ = nullptr; a_ = a;
    assert(size() >= 0);
  }
  int m () const { return m_; }
  int n () const { return n_; }
  int o () const { return o_; }
  int p () const { return p_; }
  size_t size () const { return size_t(n())*m()*o()*p(); }
  bool empty () const { return size() == 0; }
  T* data () const { return a_; }
  T& operator() (const int r, const int c, const int s, const int t) const
  { debug(r, c, s, t); return a_[((t*o_ + s)*n_ + c)*m_ + r]; }
  T* operator() (const int c, const int s, const int t) const
  { debug(0, c, s, t); return a_ + ((t*o_ + s)*n_ + c)*m_; }
  void set (const T& init) const
  { for (size_t i = 0; i < size_t(m_)*n_*o_*p_; ++i) a_[i] = init; }
  void set_ompfor (const T& init) const
  { ompfor for (size_t i = 0; i < size_t(m_)*n_*o_*p_; ++i) a_[i] = init; }
  void set_ompparfor (const T& init) const
  { ompparfor for (size_t i = 0; i < size_t(m_)*n_*o_*p_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int r, const int c, const int s, const int t) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array4D: r is " << r << " but m_ is " << m_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array4D: c is " << c << " but n_ is " << n_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (s < 0 || s >= o_) {
      std::stringstream ss;
      ss << "Array4D: s is " << s << " but o_ is " << o_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
    if (t < 0 || t >= p_) {
      std::stringstream ss;
      ss << "Array4D: t is " << t << " but p_ is " << p_ << "\n";
      throw std::runtime_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int, const int, const int, const int) {}
#endif
};

template <typename T> inline T* slice (Array2D<T>& a, const int c) { return a(c); }
template <typename T> inline const T* slice (const Array2D<T>& a, const int c) { return a(c); }
template <typename T> inline int nslices (const Array2D<T>& a) { return a.n(); }
template <typename T> inline int szslice (const Array2D<T>& a) { return a.m(); }

template <typename T> inline Array1D<T> offset (Array1D<T>& a, const int os)
{ return Array1D<T>(a.m() - os, a.data() + os); }
template <typename T> inline Array2D<T> offset (Array2D<T>& a, const int os)
{ return Array2D<T>(a.m(), a.n() - os, a.data() + a.m()*os); }

template <typename T> inline T* slice (Array3D<T>& a, const int c, const int s)
{ return a(c, s); }
template <typename T> inline const T* slice (const Array3D<T>& a, const int c, const int s)
{ return a(c, s); }
template <typename T> inline int szslice (const Array3D<T>& a) { return a.m(); }

template <typename T> inline T* slice (Array4D<T>& a, const int c, const int s, const int t)
{ return a(c, s, t); }
template <typename T> inline
const T* slice (const Array4D<T>& a, const int c, const int s, const int t)
{ return a(c, s, t); }
template <typename T> inline int szslice (const Array4D<T>& a) { return a.m(); }

// Unmanaged views.
template <typename T> inline Array1D<T> unmgd (const Array1D<T>& a)
{ return Array1D<T>(a.m(), a.data()); }
template <typename T> inline Array2D<T> unmgd (const Array2D<T>& a)
{ return Array2D<T>(a.m(), a.n(), a.data()); }
template <typename T> inline Array3D<T> unmgd (const Array3D<T>& a)
{ return Array3D<T>(a.m(), a.n(), a.o(), a.data()); }
template <typename T> inline Array4D<T> unmgd (const Array4D<T>& a)
{ return Array4D<T>(a.m(), a.n(), a.o(), a.p(), a.data()); }
template <typename T> inline Array1D<const T> const_unmgd (const Array1D<T>& a)
{ return Array1D<const T>(a.m(), a.data()); }
template <typename T> inline Array2D<const T> const_unmgd (const Array2D<T>& a)
{ return Array2D<const T>(a.m(), a.n(), a.data()); }
template <typename T> inline Array3D<const T> const_unmgd (const Array3D<T>& a)
{ return Array3D<const T>(a.m(), a.n(), a.o(), a.data()); }
template <typename T> inline Array4D<const T> const_unmgd (const Array4D<T>& a)
{ return Array4D<const T>(a.m(), a.n(), a.o(), a.p(), a.data()); }

#endif
