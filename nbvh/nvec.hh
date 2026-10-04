#ifndef NBVH_VECTOR_HH
#define NBVH_VECTOR_HH

#if __cplusplus >= 201703L
#define ENABLED_CPP_STD_17
#endif

#include <cstddef>
#include <utility> // std::index_sequence
#include <array>
#include <cmath>
#include <iostream>

////////////////////////////////////////////////////////////////
/// ND vector
////////////////////////////////////////////////////////////////

template <typename T, size_t N> struct VectorN
{
  typedef T value_type;

  constexpr VectorN(): v_() {}

  template <typename U, typename ... Us>
  constexpr VectorN(const U &_v, Us... _vs):
  v_ { { static_cast<T>(_v), static_cast<T>(_vs)... } } {}

  static constexpr size_t dimension() noexcept { return N; }

  inline const T &operator[](size_t _i) const { return v_[_i]; }
  inline T &operator[](size_t _i) { return v_[_i]; }

protected:

  std::array<T, N> v_;
};

////////////////////////////////////////////////////////////////
/// ND vector implementation
////////////////////////////////////////////////////////////////

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_neg(const VectorN<T, N> &p, std::index_sequence<I...>)
{ return { (-p[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_add(const VectorN<T, N> &a, const VectorN<T, N> &b, std::index_sequence<I...>)
{ return { (a[I] + b[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_add(const VectorN<T, N> &p, const T &s, std::index_sequence<I...>)
{ return { (p[I] + s)... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_sub(const VectorN<T, N> &a, const VectorN<T, N> &b, std::index_sequence<I...>)
{ return { (a[I] - b[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_sub(const VectorN<T, N> &p, const T &s, std::index_sequence<I...>)
{ return { (p[I] - s)... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_mul(const VectorN<T, N> &a, const VectorN<T, N> &b, std::index_sequence<I...>)
{ return { (a[I] * b[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_mul(const VectorN<T, N> &p, const T &s, std::index_sequence<I...>)
{ return { (p[I] * s)... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_div(const VectorN<T, N> &a, const VectorN<T, N> &b, std::index_sequence<I...>)
{ return { (a[I] / b[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_div(const VectorN<T, N> &p, const T &s, std::index_sequence<I...>)
{ return { (p[I] / s)... }; } // to optimize for large size

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_set(const T &s, std::index_sequence<I...>)
{ return { (s + (T)0 * I)... }; }

template <typename T, size_t N, size_t M, size_t... I>
inline VectorN<T, N> op_impl_set(const VectorN<T, M> &p, std::index_sequence<I...>)
{ return { (I < M ? p[I] : (T)0)... }; }

template <typename T, size_t N, size_t... I>
inline T op_impl_rdc(const VectorN<T, N> &p, std::index_sequence<I...>)
#ifdef ENABLED_CPP_STD_17
{ return ((p[I]) + ...); }
#else
{ using _ = int[]; T r {}; (void)_{ (r += p[I], 0)... }; return r; }
#endif

template <typename T, size_t N, size_t... I>
inline T op_impl_dot(const VectorN<T, N> &a, const VectorN<T, N> &b, std::index_sequence<I...>)
#ifdef ENABLED_CPP_STD_17
{ return ((a[I] * b[I]) + ...); }
#else
{ using _ = int[]; T r {}; (void)_{ (r += a[I] * b[I], 0)... }; return r; }
#endif

template <typename T, size_t N, class Stream, size_t... I>
inline Stream &op_impl_out(Stream &os, const VectorN<T, N>& p, const char *sep, std::index_sequence<I...>)
#ifdef ENABLED_CPP_STD_17
{ ((os << (I == 0 ? "" : sep) << p[I]), ...); return os; }
#else
{ using _ = int[]; (void)_{ (os << (I == 0 ? "" : sep) << p[I], 0)... }; return os; }
#endif

//template <typename T, size_t N, size_t... I>
//inline T op_impl_max(const VectorN<T, N> &p, std::index_sequence<I...>)
//{ using _ = int[]; T r = p[0]; (void)_{ (r = std::max(r, p[I]), 0)... }; return r; }

//template <typename T, size_t N, size_t... I>
//inline size_t op_impl_argmax(const VectorN<T, N> &p, std::index_sequence<I...>)
//{ using _ = int[]; size_t k {}; (void)_{ (k = p[k] < p[I] ? I : k, 0)... }; return k; }

//template <typename T, size_t N, size_t... I>
//inline T op_impl_min(const VectorN<T, N> &p, std::index_sequence<I...>)
//{ using _ = int[]; T r = p[0]; (void)_{ (r = std::min(r, p[I]), 0)... }; return r; }

//template <typename T, size_t N, size_t... I>
//inline size_t op_impl_argmin(const VectorN<T, N> &p, std::index_sequence<I...>)
//{ using _ = int[]; size_t k {}; (void)_{ (k = p[k] > p[I] ? I : k, 0)... }; return k; }

////////////////////////////////////////////////////////////////
/// ND vector implementation extended
////////////////////////////////////////////////////////////////

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_uop(const VectorN<T, N> &p, T (*fun)(T), std::index_sequence<I...>)
{ return { fun(p[I])... }; }

template <typename T, size_t N, typename R, size_t... I>
inline VectorN<R, N> op_impl_uop(const VectorN<T, N> &p, R (*fun)(T), std::index_sequence<I...>)
{ return { fun(p[I])... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_bop(const VectorN<T, N> &p, const T &s, T (*fun)(T, T), std::index_sequence<I...>)
{ return { fun(p[I], s)... }; }

template <typename T, size_t N, typename R, size_t... I>
inline VectorN<R, N> op_impl_bop(const VectorN<T, N> &p, const T &s, R (*fun)(T, T), std::index_sequence<I...>)
{ return { fun(p[I], s)... }; }

template <typename T, size_t N, size_t... I>
inline VectorN<T, N> op_impl_bop(const VectorN<T, N> &a, const VectorN<T, N> &b, T (*fun)(T, T), std::index_sequence<I...>)
{ return { fun(a[I], b[I])... }; }

template <typename T, size_t N, typename R, size_t... I>
inline VectorN<R, N> op_impl_bop(const VectorN<T, N> &a, const VectorN<T, N> &b, R (*fun)(T, T), std::index_sequence<I...>)
{ return { fun(a[I], b[I])... }; }

template <typename T, size_t N, size_t... I>
inline T op_impl_rdc(const VectorN<T, N> &p, T (*fun)(T), std::index_sequence<I...>)
{ using _ = int[]; T r {}; (void)_{ (r += fun(p[I]), 0)... }; return r; }

template <typename T, size_t N, size_t... I>
inline T op_impl_rdc(const VectorN<T, N> &p, T (*fun)(T, T), T init, std::index_sequence<I...>)
{ using _ = int[]; T r { init }; (void)_{ (r = fun(r, p[I]), 0)... }; return r; }

template <typename T, size_t N, typename R, size_t... I>
inline R op_impl_rdc(const VectorN<T, N> &p, R (*fun)(R, T), R init, std::index_sequence<I...>)
{ using _ = int[]; R r { init }; (void)_{ (r = fun(r, p[I]), 0)... }; return r; }

template <typename T, size_t N, size_t... I>
inline T op_impl_dot(const VectorN<T, N> &a, const VectorN<T, N> &b, T (*fun)(T, T), std::index_sequence<I...>)
{ using _ = int[]; T r {}; (void)_{ (r += fun(a[I], b[I]), 0)... }; return r; }

template <typename T, size_t N, size_t... I>
inline T op_impl_dot(const VectorN<T, N> &a, const VectorN<T, N> &b, T (*fun)(T, T, T), T init, std::index_sequence<I...>)
{ using _ = int[]; T r { init }; (void)_{ (r = fun(r, a[I], b[I]), 0)... }; return r; }

template <typename T, size_t N, typename R, size_t... I>
inline R op_impl_dot(const VectorN<T, N> &a, const VectorN<T, N> &b, R (*fun)(R, T, T), R init, std::index_sequence<I...>)
{ using _ = int[]; R r { init }; (void)_{ (r = fun(r, a[I], b[I]), 0)... }; return r; }

template <typename T, size_t N, size_t... I>
inline size_t op_impl_arg(const VectorN<T, N> &p, bool (*cmp)(T, T), std::index_sequence<I...>)
{ using _ = int[]; size_t k {}; (void)_{ (k = cmp(p[k], p[I]) ? I : k, 0)... }; return k; }

////////////////////////////////////////////////////////////////
/// ND vector generators
////////////////////////////////////////////////////////////////

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> make_vector(const T &s)
{ return op_impl_set<T, N>(s, Indices{}); }

template <typename T, size_t N, size_t M, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> make_vector(const VectorN<T, M> &p)
{ return op_impl_set<T, N, M>(p, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> make_vector(const VectorN<T, N> &p, T (*fun)(T))
{ return op_impl_uop<T, N>(p, fun, Indices{}); }

template <typename T, size_t N, typename R, typename Indices = std::make_index_sequence<N>>
inline VectorN<R, N> make_vector(const VectorN<T, N> &p, R (*fun)(T))
{ return op_impl_uop<T, N, R>(p, fun, Indices{}); }

////////////////////////////////////////////////////////////////
/// ND vector operations
////////////////////////////////////////////////////////////////

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator+(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_add(a, b, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator+(const VectorN<T, N> &p, const T &s)
{ return op_impl_add(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator+(const T &s, const VectorN<T, N> &p)
{ return op_impl_add(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator-(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_sub(a, b, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator-(const VectorN<T, N> &p, const T &s)
{ return op_impl_sub(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator-(const T &s, const VectorN<T, N> &p)
{ return op_impl_neg(op_impl_sub(p, s, Indices{}), Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator*(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_mul(a, b, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator*(const VectorN<T, N> &p, const T &s)
{ return op_impl_mul(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator*(const T &s, const VectorN<T, N> &p)
{ return op_impl_mul(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator/(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_div(a, b, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator/(const VectorN<T, N> &p, const T &s)
{ return op_impl_div(p, s, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> operator-(const VectorN<T, N> &p)
{ return op_impl_neg(p, Indices{}); }

template <typename T, size_t N>
inline VectorN<T, N> &operator+=(VectorN<T, N> &a, const VectorN<T, N> &b)
{ a = a + b; return a; }

template <typename T, size_t N>
inline VectorN<T, N> &operator+=(VectorN<T, N> &p, const T &s)
{ p = p + s; return p; }

template <typename T, size_t N>
inline VectorN<T, N> &operator-=(VectorN<T, N> &a, const VectorN<T, N> &b)
{ a = a - b; return a; }

template <typename T, size_t N>
inline VectorN<T, N> &operator-=(VectorN<T, N> &p, const T &s)
{ p = p - s; return p; }

template <typename T, size_t N>
inline VectorN<T, N> &operator*=(VectorN<T, N> &p, const T &s)
{ p = p * s; return p; }

template <typename T, size_t N>
inline VectorN<T, N> &operator/=(VectorN<T, N> &p, const T &s)
{ p = p / s; return p; }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline T sum(const VectorN<T, N> &p)
{ return op_impl_rdc(p, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline T dot(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_dot(a, b, Indices{}); }

////////////////////////////////////////////////////////////////
/// ND vector operations extended
////////////////////////////////////////////////////////////////

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline T max(const VectorN<T, N> &p)
{ return op_impl_rdc<T, N>(p, [] (T x, T y) { return x<y ? y : x; }, p[0], Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline size_t argmax(const VectorN<T, N> &p)
{ return op_impl_arg<T, N>(p, [] (T x, T y) { return x < y; }, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline T min(const VectorN<T, N> &p)
{ return op_impl_rdc<T, N>(p, [] (T x, T y) { return x<y ? x : y; }, p[0], Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline size_t argmin(const VectorN<T, N> &p)
{ return op_impl_arg<T, N>(p, [] (T x, T y) { return x > y; }, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> max(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_bop<T, N>(a, b, [] (T x, T y) { return x<y ? y : x; }, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> min(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_bop<T, N>(a, b, [] (T x, T y) { return x<y ? x : y; }, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline VectorN<T, N> abs(const VectorN<T, N> &p)
{ return op_impl_uop<T, N>(p, [] (T x) { return x<0 ? -x : x; }, Indices{}); }

//template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
//inline VectorN<T, N> pow(const VectorN<T, N> &p, const T &s)
//{ return op_impl_bop<T, N>(p, s, [] (T x, T y) { return std::pow(x, y); }, Indices{}); }

//template <typename T, size_t N>
//inline T norm(const VectorN<T, N> &p, const T &s)
//{ return std::pow(sum(pow(abs(p), s)), (T)1 / s); }

template <typename T, size_t N>
inline T norm1(const VectorN<T, N> &p)
{ return sum(abs(p)); }

template <typename T, size_t N>
inline T norm2(const VectorN<T, N> &p)
{ return sqrt(dot(p, p)); }

template <typename T, size_t N>
inline T norm8(const VectorN<T, N> &p)
{ return max(abs(p)); }

template <typename T, size_t N>
inline VectorN<T, N> normalize(const VectorN<T, N> &p)
{ return p / norm2(p); }

////////////////////////////////////////////////////////////////
/// ND vector IO
////////////////////////////////////////////////////////////////

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline std::ostream &operator<<(std::ostream &os, const VectorN<T, N>& p)
{ os << "("; op_impl_out(os, p, ", ", Indices{}) << ")"; return os; }

////////////////////////////////////////////////////////////////
/// 3D vector specializations
////////////////////////////////////////////////////////////////

template <typename T>
inline VectorN<T, 3> cross(const VectorN<T, 3> &a, const VectorN<T, 3> &b)
{ return { a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] }; }

////////////////////////////////////////////////////////////////
/// 2D vector specializations
////////////////////////////////////////////////////////////////

template <typename T>
inline T cross(const VectorN<T, 2> &a, const VectorN<T, 2> &b)
{ return a[0]*b[1] - a[1]*b[0]; }

////////////////////////////////////////////////////////////////
/// ND vector specializations
////////////////////////////////////////////////////////////////

using Vec2d = VectorN<double, 2>;
using Vec2f = VectorN<float, 2>;
using Vec2i = VectorN<int, 2>;
using Vec2b = VectorN<bool, 2>;
using Vec2c = VectorN<char, 2>;
using Vec2uc= VectorN<unsigned char, 2>;

using Vec3d = VectorN<double, 3>;
using Vec3f = VectorN<float, 3>;
using Vec3i = VectorN<int, 3>;
using Vec3b = VectorN<bool, 3>;
using Vec3c = VectorN<char, 3>;
using Vec3uc= VectorN<unsigned char, 3>;

using Vec4d = VectorN<double, 4>;
using Vec4f = VectorN<float, 4>;
using Vec4i = VectorN<int, 4>;
using Vec4b = VectorN<bool, 4>;
using Vec4c = VectorN<char, 4>;
using Vec4uc= VectorN<unsigned char, 4>;

using Vec5d = VectorN<double, 5>;
using Vec5f = VectorN<float, 5>;
using Vec5i = VectorN<int, 5>;
using Vec5b = VectorN<bool, 5>;
using Vec5c = VectorN<char, 5>;
using Vec5uc= VectorN<unsigned char, 5>;

using Vec6d = VectorN<double, 6>;
using Vec6f = VectorN<float, 6>;
using Vec6i = VectorN<int, 6>;
using Vec6b = VectorN<bool, 6>;
using Vec6c = VectorN<char, 6>;
using Vec6uc= VectorN<unsigned char, 6>;

using Vec7d = VectorN<double, 7>;
using Vec7f = VectorN<float, 7>;
using Vec7i = VectorN<int, 7>;
using Vec7b = VectorN<bool, 7>;
using Vec7c = VectorN<char, 7>;
using Vec7uc= VectorN<unsigned char, 7>;

using Vec8d = VectorN<double, 8>;
using Vec8f = VectorN<float, 8>;
using Vec8i = VectorN<int, 8>;
using Vec8b = VectorN<bool, 8>;
using Vec8c = VectorN<char, 8>;
using Vec8uc= VectorN<unsigned char, 8>;

#endif // NBVH_VECTOR_HH