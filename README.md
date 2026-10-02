# NBVH
![CMake](https://github.com/IngInx747/nbvh/actions/workflows/.github/workflows/cmake.yml/badge.svg)

A header-only library of N-dimensional Bounding Volume Hierarchy.

## Setup BVH

Suppose your data type is defined as:

```cpp
using Primitive = /* point, triangle, sphere, etc. */;
using Iter = std::vector<Primitive>::iterator;
```

Define the dimension of bounding box:

```cpp
using Box = Aabb<T, N>;
using Bvh = Bvh<Box>;
```

Tell BVH how to create the bounding box per primitive:

```cpp
struct Bound
{
  Box operator() (const Primitive&)
  { /* build the bounding box of the primitive */ }

  /* other necessary attributes */
};
```

Assign a splitting method with your data. There are 3 built-in methods(Middle-point, Equal counts and SAH).

```cpp
Bound bound(/* initializations */);
SAHSplit<Bound, Box, Iter> split(bound);
```

Setup and build the BVH on the given dataset:

```cpp
Bvh bvh ();
std::vector<Primitive> data(/* populated */);
bvh.build(bound, split, data.begin(), data.end());
```

If memory is limited or building time is constrained, use a coarser setting:

```cpp
// stop splitting if #primitives per node is less than the threshold
bvh.build(bound, split, data.begin(), data.begin(), data.end(), 100);
```

Note: The data is reordered after building BVH. If the order matters, consider using indices.

## Range query

Setup your predicate:

```cpp
struct Predicate
{
  bool operator() (const Box&)
  { /* rough query: check if your searching range hit any box(faster) */ }

  bool operator() (const Primitive&)
  { /* fine query: check if your searching range hit any primitive(slower) */ }

  /* you would like to store the results here */
};
```

Query primitives by the predicate:

```cpp
Predicate pred(/* initializations */);

if (query(bvh, pred, data.begin()))
{ /* do something */ }
```

## Nearest primitive

Setup your distance functor:

```cpp
struct Distance
{
  bool operator() (const Box&)
  { /* distance to a box */ }

  bool operator() (const Primitive&)
  { /* distance to a primitive */ }

  /* you would like to store the results here */
};
```

Search for the nearest primitive:

```cpp
Distance func(/* initializations */);

nearest(bvh, func, max_dist, data.begin());

/* minimum distance, nearest primitive, etc. */
```

Bound your search by an estimated maximum distance.

## Ray-trace

Setup your ray collision detector:

```cpp
using Vec3 = VectorN<T, 3>;

struct Collide
{
  bool operator() (const Box&)
  { /* do ray-box collision test */ }

  bool operator() (const Primitive&)
  { /* do ray-primitive collision test */ }

  Vec3 org, dir; // origin and direction
  T dist = +inf; // ray hitting distance

  /* to add more data for your good */
};
```

Trace the ray thru your dataset:

```cpp
Collide collide(/* initializations */);

if (intersect(bvh, collide, dir, data.begin()))
{ /* do something */ }
```

The ray direction is passed to the search for branch pruning.
