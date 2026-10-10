#include "nbvh.hh"

using Vec3 = VectorN<double, 3>;
using Int3 = VectorN<int, 3>;
using Box3 = Aabb<double, 3>;

using Edge = std::array<Vec3, 2>;
using Triangle = std::array<Vec3, 3>;

////////////////////////////////////////////////////////////////
/// Point functors
////////////////////////////////////////////////////////////////

struct PointBound
{
  inline Box3 operator() (const Vec3 &v) const
  { return make_aabb<Box3>(v); }
};

static inline void build_bvh(
  std::vector<Vec3> &vs,
  Bvh<Box3> &bvh)
{
  PointBound bound {};
  SAHSplit<PointBound, Box3, decltype(vs.end())> split(bound);
  bvh.build(bound, split, vs.begin(), vs.end());
}

struct InSphere
{
  InSphere(const std::vector<Vec3> &vs, const Vec3 &p, double r): vs(vs), p(p), r2(r*r) {}

  inline bool operator()(const Box3&) const;
  inline bool operator()(const Vec3&) const;
  
  const std::vector<Vec3> &vs;
  const Vec3 &p;
  const double r2 {};
  mutable int n {};
};

inline bool InSphere::operator()(const Box3 &box) const
{
  const auto d = nearest(box, p) - p;
  return dot(d, d) < r2;
}

inline bool InSphere::operator()(const Vec3 &v) const
{
  const auto d = v - p;
  if (dot(d, d) < r2) { ++n; return true; }
  return false;
}

static int query(
  const std::vector<Vec3> &vs,
  const Bvh<Box3> &bvh,
  const Vec3 &p,
  double r)
{
  InSphere pred(vs, p, r);
  query(bvh, pred, vs.begin());
  return pred.n;
}

static int query(
  const std::vector<Vec3> &vs,
  const Vec3 &p,
  double r)
{
  InSphere pred(vs, p, r);
  for (const auto &v : vs) pred(v);
  return pred.n;
}

////////////////////////////////////////////////////////////////
/// Edge functors
////////////////////////////////////////////////////////////////

struct EdgeBound
{
  EdgeBound(const std::vector<Edge> &es): es(es) {}
  inline Box3 operator() (int) const;
  const std::vector<Edge> &es;
};

inline Box3 EdgeBound::operator()(int edge) const
{
  return make_aabb<Box3>(
    es[edge][0],
    es[edge][1]);
}

static inline void build_bvh(
  const std::vector<Edge> &es,
  std::vector<int> &eids,
  Bvh<Box3> &bvh)
{
  EdgeBound bound(es);
  SAHSplit<EdgeBound, Box3, decltype(eids.end())> split(bound);
  bvh.build(bound, split, eids.begin(), eids.end());
}

inline double hausdorff(const Edge &e, const Vec3 &p)
{
  const auto &a = e[0];
  const auto &b = e[1];
  const auto d0 = p - a;
  const auto d1 = p - b;
  const auto ud = normalize(b - a);
  return dot(d0, ud) < 0 ? norm2(d0)
       : dot(d1,-ud) < 0 ? norm2(d1)
       : norm2(cross(d0, ud));
}

struct EdgeDistance
{
  EdgeDistance(const std::vector<Edge> &es, const Vec3 &p): es(es), p(p) {}

  inline double operator()(const Box3&) const;
  inline double operator()(const int) const;
  
  const std::vector<Edge> &es;
  const Vec3 &p;
  mutable double mind { 1e20 };
  mutable int ec { -1 };
};

inline double EdgeDistance::operator()(const Box3 &box) const
{
  return norm2(nearest(box, p) - p);
}

inline double EdgeDistance::operator()(const int edge) const
{
  const auto d = hausdorff(es[edge], p);
  if (mind > d) { mind = d; ec = edge; }
  return d;
}

static int nearest(
  const std::vector<Edge> &es,
  const std::vector<int> &eids,
  const Bvh<Box3> &bvh,
  const Vec3 &p,
  double &mind)
{
  EdgeDistance distance(es, p);
  nearest(bvh, distance, distance.mind, eids.begin());
  mind = distance.mind;
  return distance.ec;
}

static int nearest(
  const std::vector<Edge> &es,
  const Vec3 &p,
  double &mind)
{
  EdgeDistance distance(es, p);
  for (int i = 0; i < es.size(); ++i)
    distance(i);
  mind = distance.mind;
  return distance.ec;
}

////////////////////////////////////////////////////////////////
/// Triangle functors
////////////////////////////////////////////////////////////////

struct TriangleBound
{
  TriangleBound(const std::vector<Triangle> &fs): fs(fs) {}
  inline Box3 operator() (int) const;
  const std::vector<Triangle> &fs;
};

inline Box3 TriangleBound::operator()(int face) const
{
  return make_aabb<Box3>(
    fs[face][0],
    fs[face][1],
    fs[face][2]);
}

static inline void build_bvh(
  const std::vector<Triangle> &fs,
  std::vector<int> &fids,
  Bvh<Box3> &bvh)
{
  TriangleBound bound(fs);
  SAHSplit<TriangleBound, Box3, decltype(fids.end())> split(bound);
  bvh.build(bound, split, fids.begin(), fids.end());
}

inline bool intersecting(
  const Vec3 &v0,
  const Vec3 &v1,
  const Vec3 &v2,
  const Vec3 &org,
  const Vec3 &dir,
  double &dist,
  bool culling)
{
  constexpr double kEps = std::numeric_limits<double>::epsilon();

  auto v01 = v1 - v0;
  auto v02 = v2 - v0;
  auto pvc = cross(dir, v02); // T
  double det = dot(v01, pvc); // ((P1, V02, V01))

  if (culling && det < 0) return false;
  if (!culling && std::abs(det) < kEps) return false;

  double inv = 1 / det;
  auto tvc = org - v0; // P0 - V0
  double u = dot(tvc, pvc) * inv; // Eq.3
  if (u < 0 || u > 1) return false;

  auto qvc = cross(tvc, v01); // S
  double v = dot(dir, qvc) * inv; // Eq.4
  if (v < 0 || u + v > 1) return false;

  // distance from ray.origin to hit point
  double t = dot(v02, qvc) * inv; // Eq.5

  // update hit distance
  if (t > 0 && dist > t)
  {
    dist = t;
    return true; // ray hit primitive in distance
  }
  else return false; // ray hit primitive out of distance
}

struct TriangleCollide
{
  TriangleCollide(
    const std::vector<Triangle> &fs,
    const Vec3 &org,
    const Vec3 &dir):
    fs(fs),
    org(org),
    dir(dir),
    inv(make_vector<double, 3>(1)/dir) {}

  inline bool operator()(const Box3&) const;
  inline bool operator()(const int) const;
  
  const std::vector<Triangle> &fs;
  const Vec3 &org;
  const Vec3 &dir;
  const Vec3 inv;
  mutable double dist { 1e20 };
  mutable int fc { -1 };
};

inline bool TriangleCollide::operator()(const Box3 &box) const
{
  return intersecting(box, org, inv, dist, bool {});
}

inline bool TriangleCollide::operator()(int face) const
{
  bool hit = intersecting(
    fs[face][0],
    fs[face][1],
    fs[face][2],
    org,
    dir,
    dist,
    false);
  if (hit) fc = face;
  return hit;
}

static int trace(
  const std::vector<Triangle> &fs,
  const std::vector<int> &fids,
  const Bvh<Box3> &bvh,
  const Vec3 &org,
  const Vec3 &dir)
{
  TriangleCollide collide(fs, org, dir);
  trace(bvh, collide, dir, fids.begin());
  return collide.fc;
}

static int trace(
  const std::vector<Triangle> &fs,
  const Vec3 &org,
  const Vec3 &dir)
{
  TriangleCollide collide(fs, org, dir);
  for (int i = 0; i < fs.size(); ++i)
    collide(i);
  return collide.fc;
}

inline double hausdorff(const Triangle &t, const Vec3 &p)
{
  const auto &a = t[0];
  const auto &b = t[1];
  const auto &c = t[2];
  const auto ab = b - a;
  const auto bc = c - b;
  const auto ca = a - c;
  const auto ap = p - a;
  const auto bp = p - b;
  const auto cp = p - c;
  const auto n =  normalize(cross(ab,bc));
  const auto na = normalize(cross(bc, n));
  const auto nb = normalize(cross(ca, n));
  const auto nc = normalize(cross(ab, n));
  if (dot(ab, ap) <= 0 && dot(ca, ap) >= 0) return norm2(ap);
  if (dot(bc, bp) <= 0 && dot(ab, bp) >= 0) return norm2(bp);
  if (dot(ca, cp) <= 0 && dot(bc, cp) >= 0) return norm2(cp);
  if (dot(nc, ap) >= 0 && dot(ab, ap) >= 0 && dot(ab, bp) <= 0) return norm2(cross(ab, ap))/norm2(ab);
  if (dot(na, bp) >= 0 && dot(bc, bp) >= 0 && dot(bc, cp) <= 0) return norm2(cross(bc, bp))/norm2(bc);
  if (dot(nb, cp) >= 0 && dot(ca, cp) >= 0 && dot(ca, ap) <= 0) return norm2(cross(ca, cp))/norm2(ca);
  return fabs(dot(ap, n)); // inside triangle
}

struct TriangleDistance
{
  TriangleDistance(const std::vector<Triangle> &fs, const Vec3 &p): fs(fs), p(p) {}

  inline double operator()(const Box3&) const;
  inline double operator()(const int) const;
  
  const std::vector<Triangle> &fs;
  const Vec3 &p;
  mutable double mind { 1e20 };
  mutable int fc { -1 };
};

inline double TriangleDistance::operator()(const Box3 &box) const
{
  return norm2(nearest(box, p) - p);
}

inline double TriangleDistance::operator()(const int face) const
{
  const auto d = hausdorff(fs[face], p);
  if (mind > d) { mind = d; fc = face; }
  return d;
}

static int nearest(
  const std::vector<Triangle> &fs,
  const std::vector<int> &fids,
  const Bvh<Box3> &bvh,
  const Vec3 &p,
  double &mind)
{
  TriangleDistance distance(fs, p);
  nearest(bvh, distance, distance.mind, fids.begin());
  mind = distance.mind;
  return distance.fc;
}

static int nearest(
  const std::vector<Triangle> &fs,
  const Vec3 &p,
  double &mind)
{
  TriangleDistance distance(fs, p);
  for (int i = 0; i < fs.size(); ++i)
    distance(i);
  mind = distance.mind;
  return distance.fc;
}

////////////////////////////////////////////////////////////////
/// Datasets
////////////////////////////////////////////////////////////////

#include <random>

static double rand(double min, double max)
{
  static std::mt19937_64 _(0);
  return std::uniform_real_distribution<double>(min, max)(_);
}

static void generate_randoms(std::vector<Vec3> &vs, const size_t n)
{
  vs.resize(n);

  for (size_t i = 0; i < n; ++i)
  {
    auto &v = vs[i];
    v[0] = rand(-1, 1);
    v[1] = rand(-1, 1);
    v[2] = rand(-1, 1);
  }
}

static void generate_randoms(std::vector<Edge> &es, const size_t n)
{
  es.resize(n);

  for (size_t i = 0; i < n; ++i)
  {
    auto &e = es[i];

    for (size_t j = 0; j < 2; ++j)
    {
      e[j][0] = rand(-1, 1);
      e[j][1] = rand(-1, 1);
      e[j][2] = rand(-1, 1);
    }
  }
}

static void generate_randoms(std::vector<Triangle> &fs, const size_t n)
{
  fs.resize(n);

  for (size_t i = 0; i < n; ++i)
  {
    auto &f = fs[i];

    for (size_t j = 0; j < 3; ++j)
    {
      f[j][0] = rand(-1, 1);
      f[j][1] = rand(-1, 1);
      f[j][2] = rand(-1, 1);
    }
  }
}

////////////////////////////////////////////////////////////////
/// Tests
////////////////////////////////////////////////////////////////

static int test_query()
{
  std::vector<Vec3> vs {};
  generate_randoms(vs, 100);

  Bvh<Box3> bvh {};
  build_bvh(vs, bvh);
  
  const size_t n = 100;

  for (size_t i = 0; i < n; ++i)
  {
    const double x = rand(-1, 1);
    const double y = rand(-1, 1);
    const double z = rand(-1, 1);
    const Vec3 p { x, y, z };
    const double r = rand(0.01, 0.5);
    const auto n0 = query(vs, p, r);
    const auto n1 = query(vs, bvh, p, r);
    if (n0 != n1)
    {
      printf("test_query %zu: results mismatch %d %d\n", i, n0, n1);
      return 1;
    }
  }

  return 0;
}

static int test_trace()
{
  std::vector<Triangle> fs {};
  generate_randoms(fs, 100);

  std::vector<int> fids(fs.size());
  for (int i = 0; i < fs.size(); ++i) fids[i] = i;

  Bvh<Box3> bvh {};
  build_bvh(fs, fids, bvh);
  
  const Vec3 org { 0, 0, 2 };
  const Vec3 dz { 0, 0, -1 };
  const size_t n = 100;

  for (size_t i = 0; i < n; ++i)
  {
    const double x = rand(-1, 1);
    const double y = rand(-1, 1);
    const auto dir = normalize(dz + Vec3{ x, y, 0 });
    const auto fc0 = trace(fs, org, dir);
    const auto fc1 = trace(fs, fids, bvh, org, dir);
    if (fc0 != fc1)
    {
      printf("test_trace %zu: results mismatch %d %d\n", i, fc0, fc1);
      return 1;
    }
  }
  
  for (size_t i = 0; i < n; ++i)
  {
    const double x = rand(-1, 1);
    const double y = rand(-1, 1);
    const auto eye = org + Vec3{ x, y, 0 };
    const auto fc0 = trace(fs, eye, dz);
    const auto fc1 = trace(fs, fids, bvh, eye, dz);
    if (fc0 != fc1)
    {
      printf("test_trace %zu: results mismatch %d %d\n", i, fc0, fc1);
      return 1;
    }
  }

  return 0;
}

static int test_nearest()
{
  std::vector<Triangle> fs {};
  generate_randoms(fs, 100);

  std::vector<int> fids(fs.size());
  for (int i = 0; i < fs.size(); ++i) fids[i] = i;

  Bvh<Box3> bvh {};
  build_bvh(fs, fids, bvh);
  
  const size_t n = 100;

  for (size_t i = 0; i < n; ++i)
  {
    const double x = rand(-1, 1);
    const double y = rand(-1, 1);
    const double z = rand(-1, 1);
    const Vec3 p { x, y, z };
    double mind0 {};
    double mind1 {};
    const auto fc0 = nearest(fs, p, mind0);
    const auto fc1 = nearest(fs, fids, bvh, p, mind1);
    if (fabs(mind0 - mind1) > 1e-10)
    {
      printf("test_nearest %zu: results mismatch %d %d\n", i, fc0, fc1);
      return 1;
    }
  }

  return 0;
}

int main(int argc, const char **argv)
{
  int err {};

  err = test_query();
  if (err) return err;
  err = test_trace();
  if (err) return err;
  err = test_nearest();
  if (err) return err;

  return 0;
}