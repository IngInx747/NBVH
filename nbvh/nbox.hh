#ifndef NBVH_AABB_HH
#define NBVH_AABB_HH

#include <limits>
#include "nvec.hh"

////////////////////////////////////////////////////////////////
/// Axis-aligned Bounding Box
////////////////////////////////////////////////////////////////

template <typename T, size_t N> struct Aabb
{
  typedef VectorN<T, N> vector_type;
  typedef T value_type;

  static constexpr size_t dimension() noexcept { return N; }

  inline const vector_type &operator[](size_t i) const { return v_[i]; }
  inline vector_type &operator[](size_t i) { return v_[i]; }

  vector_type v_[2];
};

////////////////////////////////////////////////////////////////
/// AABB impl utilities
////////////////////////////////////////////////////////////////

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline bool all_less(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_dot<T, N, bool>(a, b, [] (bool r, T x, T y) { return r && x < y; }, true, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline bool all_leq(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_dot<T, N, bool>(a, b, [] (bool r, T x, T y) { return r && x <= y; }, true, Indices{}); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline bool any_less(const VectorN<T, N> &a, const VectorN<T, N> &b)
{ return op_impl_dot<T, N, bool>(a, b, [] (bool r, T x, T y) { return r || x < y; }, false, Indices{}); }

template <typename T, size_t N, typename... R>
inline VectorN<T, N> max(const VectorN<T, N> &a, const VectorN<T, N> &b, const VectorN<R, N> &... rest)
{ return max(max(a, b), rest...); }

template <typename T, size_t N, typename... R>
inline VectorN<T, N> min(const VectorN<T, N> &a, const VectorN<T, N> &b, const VectorN<R, N> &... rest)
{ return min(min(a, b), rest...); }

////////////////////////////////////////////////////////////////
/// AABB boolean impls
////////////////////////////////////////////////////////////////

template <typename T, size_t N>
inline bool valid(const Aabb<T, N> &b)
{ return all_leq(b[0], b[1]); }

template <typename T, size_t N>
inline bool alldim(const Aabb<T, N> &b) // all dimensions are positive
{ return all_less(b[0], b[1]); }

template <typename T, size_t N>
inline bool anydim(const Aabb<T, N> &b) // any dimension is positive
{ return any_less(b[0], b[1]); }

template <typename T, size_t N>
inline bool inside(const Aabb<T, N> &b, const VectorN<T, N> &v)
{ return all_leq(b[0], v) && all_leq(v, b[1]); }

template <typename T, size_t N>
inline bool inside(const Aabb<T, N> &b, const VectorN<T, N> &v, bool)
{ return all_less(b[0], v) && all_less(v, b[1]); }

template <typename T, size_t N>
inline bool inside(const Aabb<T, N> &B, const Aabb<T, N> &b)
{ return all_leq(B[0], b[0]) && all_leq(b[1], B[1]); }

template <typename T, size_t N>
inline bool inside(const Aabb<T, N> &B, const Aabb<T, N> &b, bool)
{ return all_less(B[0], b[0]) && all_less(b[1], B[1]); }

template <typename T, size_t N>
inline bool intersecting(const Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ return all_leq(b0[0], b1[1]) && all_leq(b1[0], b0[1]); }

template <typename T, size_t N>
inline bool intersecting(const Aabb<T, N> &b0, const Aabb<T, N> &b1, bool)
{ return all_less(b0[0], b1[1]) && all_less(b1[0], b0[1]); }

// We can rely on the IEEE 754 floating-point properties which also implicitly
//   handles the corner case where a component of the direction is zero - the
//   tx0 and tx1 values will be infinities of opposite signs if the ray is in
//   the slabs and leave t0 and t1 unchanged.
// If the ray is out of the slabs, tx0 and tx1 will be infinities of the same
//   sign, turning t0 to +inf or t1 to -inf and fails the intersection test.
// AABB is not an entity in the space, hence the distance does not need to be
//   updated in an intersection test between an AABB and a ray.

template <typename T, size_t N>
inline bool intersecting(
  const Aabb<T, N> &b,
  const VectorN<T, N> &org,
  const VectorN<T, N> &dir,
  const T &dist)
{
  const auto k0 = (b[0] - org)/dir;
  const auto k1 = (b[1] - org)/dir;
  const auto t0 = max(min(k0, k1));
  const auto t1 = min(max(k0, k1));
  return t1 > 0 && t1 >= t0 && dist > t0;
}

template <typename T, size_t N>
inline bool intersecting(
  const Aabb<T, N> &b,
  const VectorN<T, N> &org,
  const VectorN<T, N> &inv,
  const T &dist, bool)
{
  const auto k0 = (b[0] - org)*inv;
  const auto k1 = (b[1] - org)*inv;
  const auto t0 = max(min(k0, k1));
  const auto t1 = min(max(k0, k1));
  return t1 > 0 && t1 >= t0 && dist > t0;
}

////////////////////////////////////////////////////////////////
/// AABB property impls
////////////////////////////////////////////////////////////////

template <typename T, size_t N>
inline VectorN<T, N> centroid(const Aabb<T, N> &b)
{ return (b[0] + b[1])/(T)(2); }

template <typename T, size_t N>
inline VectorN<T, N> diagonal(const Aabb<T, N> &b)
{ return b[1] - b[0]; }

template <typename T, size_t N>
inline T length(const Aabb<T, N> &b, size_t dim)
{ return diagonal(b)[dim]; }

template <typename T, size_t N>
inline T lenmax(const Aabb<T, N> &b)
{ return max(diagonal(b)); }

template <typename T, size_t N>
inline size_t longest(const Aabb<T, N> &b)
{ return argmax(diagonal(b)); }

template <typename T, size_t N, typename Indices = std::make_index_sequence<N>>
inline T volume(const Aabb<T, N> &b)
{ return op_impl_rdc<T, N>(diagonal(b), [] (T x, T y) { return x*y; }, (T)1, Indices{}); }

template <typename T, size_t N>
inline VectorN<T, N> nearest(const Aabb<T, N> &b, const VectorN<T, N> &v)
{ return max(b[0], min(b[1], v)); }

////////////////////////////////////////////////////////////////
/// AABB operation impls
////////////////////////////////////////////////////////////////

template <typename T, size_t N>
inline Aabb<T, N> merge(const Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ return { min(b0[0], b1[0]), max(b0[1], b1[1]) }; }

template <typename T, size_t N>
inline Aabb<T, N> intersect(const Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ return { max(b0[0], b1[0]), min(b0[1], b1[1]) }; }

template <typename T, size_t N>
inline Aabb<T, N> operator|(const Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ return merge(b0, b1); }

template <typename T, size_t N>
inline Aabb<T, N> operator&(const Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ return intersect(b0, b1); }

template <typename T, size_t N>
inline Aabb<T, N> &operator|=(Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ b0 = merge(b0, b1); return b0; }

template <typename T, size_t N>
inline Aabb<T, N> &operator&=(Aabb<T, N> &b0, const Aabb<T, N> &b1)
{ b0 = intersect(b0, b1); return b0; }

////////////////////////////////////////////////////////////////
/// AABB ctors
////////////////////////////////////////////////////////////////

template <typename T, size_t N>
inline Aabb<T, N> make_aabb()
{
  return {
    make_vector<T, N>(+std::numeric_limits<T>::max()),
    make_vector<T, N>(-std::numeric_limits<T>::max())
  };
}

template <typename T, size_t N>
inline Aabb<T, N> make_aabb(const VectorN<T, N> &v)
{ return { v, v }; }

template <typename T, size_t N, typename... R>
inline Aabb<T, N> make_aabb(const VectorN<T, N> &v, const VectorN<R, N> &... vs)
{ return { min(v, vs...), max(v, vs...) }; }

template <class BoxT>
inline BoxT make_aabb()
{ return make_aabb<typename BoxT::value_type, BoxT::dimension()>(); }

template <class BoxT, class VecT>
inline BoxT make_aabb(const VecT &v)
{ return make_aabb<typename BoxT::value_type, BoxT::dimension()>(v); }

template <class BoxT, class VecT, class ...Vecs>
inline BoxT make_aabb(const VecT &v, const Vecs &... vs)
{ return make_aabb<typename BoxT::value_type, BoxT::dimension()>(v, vs...); }

////////////////////////////////////////////////////////////////
/// 3D AABB property impls
////////////////////////////////////////////////////////////////

template <typename T>
inline T area(const Aabb<T, 3> &b)
{
  const auto d = diagonal(b);
  return (d[0]*d[1] + d[0]*d[2] + d[1]*d[2])*(T)2;
}

////////////////////////////////////////////////////////////////
/// 2D AABB property impls
////////////////////////////////////////////////////////////////

template <typename T>
inline T area(const Aabb<T, 2> &b)
{
  const auto d = diagonal(b);
  return (d[0] + d[1])*(T)2;
}

#endif // !NBVH_AABB_HH