#ifndef INCLUDE_SLMM_ARRAY_HPP
#define INCLUDE_SLMM_ARRAY_HPP

#include <memory>
#include <exception>
#include <cassert>
#include <sstream>
#include <cstring>

#define SLMM_FORCEINLINE inline __attribute__((always_inline))

namespace slmm {

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
  Array1D (T* const a, const int n) { reset(n, a); }
  Array1D (const Array1D<T_nonconst>& v)
    : n_(v.n_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int n) {
    n_ = n;
    a_p_ = std::shared_ptr<T>(new T[n], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int n, T* const a) { n_ = n; a_p_ = nullptr; a_ = a; }
  SLMM_FORCEINLINE int n () const { return n_; }
  SLMM_FORCEINLINE int size () const { return n_; }
  SLMM_FORCEINLINE T* data () { return a_; }
  const T* data () const { return a_; }
  SLMM_FORCEINLINE T& operator[] (const int i) {
#ifndef NDEBUG
    debug(i);
#endif
    return a_[i];
  }
  SLMM_FORCEINLINE const T& operator[] (const int i) const {
#ifndef NDEBUG
    debug(i);
#endif
    return a_[i];
  }
  SLMM_FORCEINLINE T& operator() (const int i) {
#ifndef NDEBUG
    debug(i);
#endif
    return a_[i];
  }
  SLMM_FORCEINLINE const T& operator() (const int i) const {
#ifndef NDEBUG
    debug(i);
#endif
    return a_[i];
  }
  void set (const T& init) { for (int i = 0; i < n_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int& i) const {
    if (i < 0 || i >= n_) {
      std::stringstream ss;
      ss << "Array1D: i is " << i << " but n_ is " << n_ << "\n";
      throw std::logic_error(ss.str().c_str());
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
  Array2D (T* const a, const int m, const int n) { reset(m, n, a); }
  Array2D (const Array2D<T_nonconst>& v)
    : m_(v.m_), n_(v.n_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int m, const int n) {
    m_ = m; n_ = n;
    a_p_ = std::shared_ptr<T>(new T[m*n], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int m, const int n, T* const a) {
    m_ = m; n_ = n; a_p_ = nullptr; a_ = a;
  }
  SLMM_FORCEINLINE int m () const { return m_; }
  SLMM_FORCEINLINE int n () const { return n_; }
  SLMM_FORCEINLINE size_type size () const { return m_*n_; }
  SLMM_FORCEINLINE T* data () { return a_; }
  SLMM_FORCEINLINE const T* data () const { return a_; }
  SLMM_FORCEINLINE T& operator() (const int r, const int c) {
#ifndef NDEBUG
    debug(r, c);
#endif
    return a_[r*n_ + c];
  }
  SLMM_FORCEINLINE const T& operator() (const int r, const int c) const {
#ifndef NDEBUG
    debug(r, c);
#endif
    return a_[r*n_ + c];
  }
  SLMM_FORCEINLINE T* operator() (const int r) {
#ifndef NDEBUG
    debug(r, 0);
#endif
    return a_ + n_*r;
  }
  SLMM_FORCEINLINE const T* operator() (const int r) const {
#ifndef NDEBUG
    debug(r, 0);
#endif
    return a_ + n_*r;
  }
  void set (const T& init) { for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int& r, const int& c) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array2D: r is " << r << " but m_ is " << m_ << "\n";
      throw std::logic_error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array2D: c is " << c << " but n_ is " << n_ << "\n";
      throw std::logic_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int& r, const int& c) {}
#endif
};

template <typename T, int dim>
class Array2Dim {
  typedef typename std::remove_const<T>::type T_nonconst;
  friend class Array2Dim<const T_nonconst,dim>;
  enum : int { n_ = dim };
  int m_;
  std::shared_ptr<T> a_p_;
  T* a_;
public:
  typedef int size_type;
  Array2Dim () : m_(0) {}
  Array2Dim (const int m) { reset(m); }
  Array2Dim (T* const a, const int m) { reset(m, a); }
  Array2Dim (const Array2D<T_nonconst>& v)
    : m_(v.m_), a_p_(v.a_p_), a_(v.a_)
  {}
  void reset (const int m) {
    m_ = m;
    a_p_ = std::shared_ptr<T>(new T[m*n_], std::default_delete<T[]>());
    a_ = a_p_.get();
  }
  void reset (const int m, T* const a) {
    m_ = m; a_p_ = nullptr; a_ = a;
  }
  SLMM_FORCEINLINE int m () const { return m_; }
  SLMM_FORCEINLINE int n () const { return n_; }
  SLMM_FORCEINLINE size_type size () const { return m_*n_; }
  SLMM_FORCEINLINE T* data () { return a_; }
  SLMM_FORCEINLINE const T* data () const { return a_; }
  SLMM_FORCEINLINE T& operator() (const int r, const int c) {
#ifndef NDEBUG
    debug(r, c);
#endif
    return a_[r*n_ + c];
  }
  SLMM_FORCEINLINE const T& operator() (const int r, const int c) const {
#ifndef NDEBUG
    debug(r, c);
#endif
    return a_[r*n_ + c];
  }
  SLMM_FORCEINLINE T* operator() (const int r) {
#ifndef NDEBUG
    debug(r, 0);
#endif
    return a_ + n_*r;
  }
  SLMM_FORCEINLINE const T* operator() (const int r) const {
#ifndef NDEBUG
    debug(r, 0);
#endif
    return a_ + n_*r;
  }
  void set (const T& init) { for (int i = 0; i < m_*n_; ++i) a_[i] = init; }
private:
#ifndef NDEBUG
  void debug (const int& r, const int& c) const {
    if (r < 0 || r >= m_) {
      std::stringstream ss;
      ss << "Array2Dim: r is " << r << " but m_ is " << m_ << "\n";
      throw std::logic_error(ss.str().c_str());
    }
    if (c < 0 || c >= n_) {
      std::stringstream ss;
      ss << "Array2Dim: c is " << c << " but n_ is " << n_ << "\n";
      throw std::logic_error(ss.str().c_str());
    }
  }
#else
  static void debug (const int& r, const int& c) {}
#endif
};

template <typename T> inline T* slice (Array2D<T>& a, const int r) { return a(r); }
template <typename T> inline const T* slice (const Array2D<T>& a, const int r) { return a(r); }
template <typename T, int d> inline T* slice (Array2Dim<T,d>& a, const int r) { return a(r); }
template <typename T, int d> inline const T* slice (const Array2Dim<T,d>& a, const int r) { return a(r); }
template <typename T> inline int nslices (const Array1D<T>& a) { return a.size(); }
template <typename T> inline int len     (const Array1D<T>& a) { return a.size(); }
template <typename T> inline int nslices (const Array2D<T>& a) { return a.m(); }
template <typename T> inline int szslice (const Array2D<T>& a) { return a.n(); }
template <typename T, int d> inline int nslices (const Array2Dim<T,d>& a) { return a.m(); }
template <typename T, int d> inline int szslice (const Array2Dim<T,d>& a) { return a.n(); }

template <typename T> void resize (Array1D<T>& a, const int n) { a.reset(n); }
template <typename T> void resize (Array2D<T>& a, const int m, const int n) { a.reset(m, n); }
template <typename T, int d> void resize (Array2Dim<T,d>& a, const int m) { a.reset(m); }
template <typename T> void copy (Array1D<T>& d, const Array1D<T>& s)
{ std::copy(s.data(), s.data() + s.size(), d.data()); }
template <typename T> void copy (Array2D<T>& d, const Array2D<T>& s)
{ std::copy(s.data(), s.data() + s.size(), d.data()); }
template <typename T, int D> void copy (Array2Dim<T,D>& d, const Array2Dim<T,D>& s)
{ std::copy(s.data(), s.data() + s.size(), d.data()); }
template <typename T, typename S> void copy (Array1D<T>& d, const S& s) { d.set(s); }
template <typename T, typename S> void copy (Array2D<T>& d, const S& s) { d.set(s); }
template <typename T, typename S, int D> void copy (Array2Dim<T,D>& d, const S& s) { d.set(s); }

} // namespace slmm

#endif
