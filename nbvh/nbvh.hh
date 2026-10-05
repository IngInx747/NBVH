// ======================================================================== //
// Copyright (c) 2022 Ingram Inxent                                         //
//                                                                          //
// Permission is hereby granted, free of charge, to any person obtaining    //
// a copy of this software and associated documentation files (the          //
// "Software"), to deal in the Software without restriction, including      //
// without limitation the rights to use, copy, modify, merge, publish,      //
// distribute, sublicense, and/or sell copies of the Software, and to       //
// permit persons to whom the Software is furnished to do so, subject to    //
// the following conditions:                                                //
//                                                                          //
// The above copyright notice and this permission notice shall be           //
// included in all copies or substantial portions of the Software.          //
//                                                                          //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,          //
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                    //
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE   //
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   //
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION    //
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.          //
// ======================================================================== //

#ifndef NBVH_HH
#define NBVH_HH

#include <stack>
#include <vector>
#include <algorithm>
#include "nbox.hh"

////////////////////////////////////////////////////////////////
/// Bounding volume hierarchy
////////////////////////////////////////////////////////////////

template <class BoxT, typename IndexT = int> class Bvh
{
public:

  typedef IndexT index_type;

public:

  /// For inner nodes, i_{0,1} = index of left and right children;
  /// for leaf nodes, i_0 = offset to the base, i_1 = -|indices|.
  struct Node
  {
    BoxT b {};
    IndexT i[2] {};

    inline const auto &left()  const { return i[0]; }
    inline const auto &right() const { return i[1]; }
    inline auto &left()  { return i[0]; }
    inline auto &right() { return i[1]; }
    inline auto offset() const { return i[0]; }
    inline auto length() const { return-i[1]; }
    inline bool leaf() const { return i[1] < 0; }
    inline void leaf(IndexT _i, IndexT _n) { i[0] = _i; i[1] = -_n; }
  };

public:

  template <class BoundT, class SplitT, class Iter>
  inline void build(
    const BoundT &bound,
    const SplitT &split,
    const Iter &begin,
    const Iter &end,
    const IndexT threshold = 1);

public:

  inline bool empty() const
  { return nodes_.empty(); }

  inline const auto &nodes() const
  { return nodes_; }

  inline void clear()
  { nodes_.clear(); }

protected:

  std::vector<Node> nodes_;
};

////////////////////////////////////////////////////////////////
/// Bvh build
////////////////////////////////////////////////////////////////

template <class BoxT, typename IndexT>
template <class BoundT, class SplitT, class Iter>
inline void Bvh<BoxT, IndexT>::build(
  const BoundT &bound,
  const SplitT &split,
  const Iter &base,
  const Iter &_end,
  const IndexT threshold)
{
  if (base == _end) return;

  struct SE { Iter begin, end; IndexT id; };
  std::stack<SE> se({SE { base, _end, 0 }});
  nodes_.assign(1, Node {});

  while (!se.empty())
  {
    const auto entry = se.top(); se.pop();
    const auto &begin = entry.begin;
    const auto &end = entry.end;
    auto &node = nodes_[entry.id];

    // split current set into two subsets
    const auto len = (IndexT)std::distance(begin, end);
    const auto off = (IndexT)std::distance(base, begin);
    auto pivot = end;
    if (len > threshold)
      pivot = split(begin, end);

    // The node is leaf iff:
    // 1. meets granularity;
    // 2. trivial splitting.
    if (pivot == begin || pivot == end)
    {
      node.leaf(off, len);
      node.b = bound(*begin);
      for (auto it = begin + 1; it != end; ++it)
        node.b |= bound(*it);
    }
    else // the node is inner
    {
      const auto left = (IndexT)nodes_.size();
      const auto right = left + 1;
      node.left()  = left;
      node.right() = right;
      nodes_.emplace_back();
      nodes_.emplace_back();
      se.push({ pivot, end,  right });
      se.push({ begin, pivot, left });
    }
  }

  // Build bounding box of inner nodes.
  // Assume parent_index < child_index:
  // node.left > i and node.right > i,
  // so their boxes have been computed.
  for (IndexT i = (IndexT)nodes_.size() - 1; i >= 0; --i)
  {
    auto &node = nodes_[i]; if (node.leaf()) continue;
    node.b = nodes_[node.left()].b | nodes_[node.right()].b;
  }
}

////////////////////////////////////////////////////////////////
/// Bvh Split Methods
////////////////////////////////////////////////////////////////

/// Split Method: EqualCounts
/// Partition primitives into equally-sized subsets
template <class BoundT, class BoxT, class Iter> struct EqualCountSplit
{
  EqualCountSplit(const BoundT &_bound): bound_(_bound) {}

  inline Iter operator()(const Iter&, const Iter&) const;

  const BoundT &bound_;
};

template <class BoundT, class BoxT, class Iter>
inline Iter EqualCountSplit<BoundT, BoxT, Iter>::operator()(const Iter &begin, const Iter &end) const
{
  auto bc = bound_(*begin); // centroid box
  for (auto it = begin + 1; it != end; ++it)
    bc |= bound_(*it);

  const auto axis = longest(bc);
  const auto pivot = begin + std::distance(begin, end)/2;

  std::nth_element(begin, pivot, end, [&](const auto &a, const auto &b)
  { return centroid(bound_(a))[axis] < centroid(bound_(b))[axis]; });

  return pivot;
}

/// Split Method: MiddlePoint
/// Partition primitives through node's midpoint
template <class BoundT, class BoxT, class Iter> struct MiddlePointSplit
{
  MiddlePointSplit(const BoundT &_bound): bound_(_bound) {}

  inline Iter operator()(const Iter&, const Iter&) const;

  const BoundT &bound_;
};

template <class BoundT, class BoxT, class Iter>
inline Iter MiddlePointSplit<BoundT, BoxT, Iter>::operator()(const Iter &begin, const Iter &end) const
{
  auto bc = bound_(*begin); // centroid box
  for (auto it = begin + 1; it != end; ++it)
    bc |= bound_(*it);

  const auto axis = longest(bc);
  const auto mv = centroid(bc)[axis];

  auto pivot = std::partition(begin, end, [&](const auto &val)
  { return centroid(bound_(val))[axis] < mv; });

  if (pivot == begin || pivot == end) // fallback to EqualCount
  {
    pivot = begin + std::distance(begin, end)/2;
    std::nth_element(begin, pivot, end, [&](const auto &a, const auto &b)
    { return centroid(bound_(a))[axis] < centroid(bound_(b))[axis]; });
  }

  return pivot;
}

/// Split Method: SAH
/// Partition primitives via surface area heuristic
template <class BoundT, class BoxT, class Iter> struct SAHSplit
{
  typedef typename std::iterator_traits<Iter>::value_type value_type;

  SAHSplit(const BoundT &_bound): bound_(_bound) {}

  inline Iter operator()(const Iter &_begin, const Iter &_end) const;

  const BoundT &bound_;
  size_t bsize_ { 16 };
};

template <class BoundT, class BoxT, class Iter>
inline Iter SAHSplit<BoundT, BoxT, Iter>::operator()(const Iter &begin, const Iter &end) const
{
  auto bc = bound_(*begin); // centroid box
  for (auto it = begin + 1; it != end; ++it)
    bc |= bound_(*it);

  // all-zero box, stop splitting
  if (!anydim(bc)) return begin;

  const auto axis = longest(bc);
  const auto inv = 1/diagonal(bc)[axis];

  std::vector<BoxT> bs(bsize_, make_aabb<BoxT>());
  std::vector<size_t> counts(bsize_, 0);

  for (auto it = begin; it != end; ++it)
  {
    const auto d = centroid(bound_(*it)) - bc[0];
    size_t k = (size_t)(bsize_*(d[axis]*inv));
    if (k >= bsize_) k = bsize_ - 1;
    bs[k] |= bound_(*it);
    ++counts[k];
  }

  // the cost of splitting buckets into [0, b] and [b+1, :]
  auto minc = std::numeric_limits<decltype(inv)>::max();
  size_t argminc = 0; // the bucket id to split

  // find bucket id that minimizes SAH metric
  for (size_t k = 0; k < bsize_ - 1; ++k)
  {
    auto b0 = bs[0], b1 = bs[k + 1];
    size_t count0 {}, count1 {};

    for (size_t i = 0; i <= k; ++i)
    {
      b0 |= bs[i];
      count0 += counts[i];
    }
    for (size_t i = k + 1; i < bsize_; ++i)
    {
      b1 |= bs[i];
      count1 += counts[i];
    }

    const auto cost = area(b0)*count0 + area(b1)*count1;
    if (minc > cost) { minc = cost; argminc = k; }
  }

  // split at the position of minimum cost
  auto pivot = std::partition(begin, end, [&](const auto &val)
  {
    const auto d = centroid(bound_(val)) - bc[0];
    size_t k = (size_t)(bsize_*(d[axis]*inv));
    if (k >= bsize_) k = bsize_ - 1;
    return k <= argminc;
  });

  if (pivot == begin || pivot == end) // fallback to EqualCount
  {
    pivot = begin + std::distance(begin, end)/2;
    std::nth_element(begin, pivot, end, [&](const auto &a, const auto &b)
    { return centroid(bound_(a))[axis] < centroid(bound_(b))[axis]; });
  }

  return pivot;
}

////////////////////////////////////////////////////////////////
/// Bvh queries
////////////////////////////////////////////////////////////////

template <class BvhT, class PredicateT, class Iter>
inline bool query(
  const BvhT &bvh,
  PredicateT &pred,
  const Iter &base)
{
  const auto &nodes = bvh.nodes();
  if (nodes.empty()) return false;

  bool hit {};
  std::stack<typename BvhT::index_type> si({ 0 });

  while (!si.empty())
  {
    const auto curr = si.top(); si.pop();
    const auto &node = nodes[curr];

    if (pred(node.b))
    {
      if (node.leaf())
      {
        const auto ib = node.offset();
        const auto ie = ib + node.length();
        for (auto i = ib; i < ie; ++i)
          if (pred(*(base + i)))
            hit = true;
      }
      else
      {
        si.push(node.right());
        si.push(node.left());
      }
    }
  }

  return hit;
}

template <class BvhT, class DistFuncT, typename DistanceT, class Iter>
inline bool nearest(
  const BvhT &bvh,
  DistFuncT &func,
  DistanceT mind,
  const Iter &base)
{
  const auto &nodes = bvh.nodes();
  if (nodes.empty()) return false;

  bool hit {};
  struct SE { typename BvhT::index_type i; DistanceT d; };
  std::stack<SE> se({ SE { 0, func(nodes.front().b) } });

  while (!se.empty())
  {
    const auto entry = se.top(); se.pop();
    const auto &node = nodes[entry.i];
    if (entry.d >= mind) continue;
    
    if (node.leaf())
    {
      const auto ib = node.offset();
      const auto ie = ib + node.length();

      for (auto i = ib; i < ie; ++i)
      {
        const auto d = func(*(base + i));
        if (mind > d) { mind = d; hit = true; }
      }
    }
    else
    {
      const auto left = node.left();
      const auto right = node.right();
      const auto dl = func(nodes[left].b);
      const auto dr = func(nodes[right].b);

      if (dl > dr)
      {
        se.push({ left, dl });
        se.push({ right, dr });
      }
      else
      {
        se.push({ right, dr });
        se.push({ left, dl });
      }
    }
  }

  return hit;
}

template <class BvhT, class CollideT, class VecT, class Iter>
inline bool intersect(
  const BvhT &bvh,
  CollideT &colli,
  const VecT &dir,
  const Iter &base)
{
  const auto &nodes = bvh.nodes();
  if (nodes.empty()) return false;

  bool hit {};
  std::stack<typename BvhT::index_type> si({ 0 });

  while (!si.empty())
  {
    const auto curr = si.top(); si.pop();
    const auto &node = nodes[curr];

    if (colli(node.b))
    {
      if (node.leaf())
      {
        const auto ib = node.offset();
        const auto ie = ib + node.length();
        for (auto i = ib; i < ie; ++i)
          if (colli(*(base + i)))
            hit = true;
      }
      else
      {
        const auto axis = longest(node.b);

        if (dir[axis] < 0)
        {
          si.push(node.left());
          si.push(node.right());
        }
        else
        {
          si.push(node.right());
          si.push(node.left());
        }
      }
    }
  }

  return hit;
}

////////////////////////////////////////////////////////////////
/// Bvh functor templates
////////////////////////////////////////////////////////////////

/// Bounding box Interfaces:
/// 
/// struct Bound
/// {
///   Aabb operator() (const Value &);
///   ...
/// };
/// 

/// Spatial query Interfaces:
/// 
/// struct Predicate
/// {
///   bool operator() (const Box &);
///   bool operator() (const Value &);
///   ...
/// };
/// 

/// Ray-trace Interfaces:
/// 
/// struct Collide
/// {
///   bool operator() (const Box &);
///   bool operator() (const Value &);
///   ...
/// };
/// 

/// Build Bvh with SAH:
/// 
/// template <class BoundT, class BoxT, class Iter>
/// void build_bvh(
///   Bvh<BoxT> &bvh,
///   const BoundT &bound,
///   Iter begin, Iter end)
/// {
///   SAHSplit<BoundT, BoxT, Iter> split(bound);
///   bvh.build(bound, split, begin, end);
/// }

#endif // !NBVH_HH