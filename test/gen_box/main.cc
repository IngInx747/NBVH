#include <sstream>
#include <fstream>
#include <iomanip>
#include "nbvh.hh"

using Vec3 = VectorN<double, 3>;
using Int2 = VectorN<int, 2>;
using Int3 = VectorN<int, 3>;
using Box3 = Aabb<double, 3>;

struct TriangleBound
{
  TriangleBound(
    const std::vector<Vec3> &vs,
    const std::vector<Int3> &fs):
    vs(vs), fs(fs) {}
  inline Box3 operator() (int) const;
  const std::vector<Vec3> &vs;
  const std::vector<Int3> &fs;
};

inline Box3 TriangleBound::operator()(int fid) const
{
  const auto &f = fs[fid];
  return make_aabb<Box3>(vs[f[0]], vs[f[1]], vs[f[2]]);
}

static int read_obj(
  std::vector<Vec3> &vs,
  std::vector<Int3> &fs,
  const char *filename)
{
  constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

  std::ifstream in(filename, std::ios::in);
  if (!in) return 1;

  for (std::string line; std::getline(in, line); )
  {
    std::stringstream ss(line);
    std::string token {};
    ss >> token;

    if (token.compare("v") == 0) // vertex
    {
      Vec3 v {};
      ss >> v[0]; ss >> v[1]; ss >> v[2];
      vs.push_back(v);
    }
    else if (token.compare("f") == 0) // face
    {
      Int3 vtn[4] = { {-1,-1,-1}, {-1,-1,-1}, {-1,-1,-1}, {-1,-1,-1} };
      std::string ids; // v or v/t or v/t/n

      for (size_t i = 0; ss >> ids; ++i)
      {
        std::stringstream si(ids);
        std::string id {};

        for (size_t j = 0; std::getline(si, id, '/'); ++j)
        {
          std::stringstream sj(id);
          sj >> vtn[i][j];
        }
      }

      fs.push_back(Int3 { vtn[0][0], vtn[1][0], vtn[2][0] } - 1);
      if (vtn[3][0] != -1) fs.push_back(Int3 { vtn[0][0], vtn[2][0], vtn[3][0] } - 1);
    }
    else if (token.compare("vt") == 0) // uv
    {
    }
    else if (token.compare("vn") == 0) // normal
    {
    }
  }

  return 0;
}

static int save_obj(
  const double *vs, const size_t nv,
  const int    *fs, const size_t nf,
  const int    *es, const size_t ne,
  const char *filename,
  const std::streamsize prec = 17)
{
  std::ofstream out(filename, std::ios::out);
  if (!out) return out.bad();

  out << std::defaultfloat << std::setprecision(prec);

  for (size_t i = 0; i < nv; ++i)
    out << "v "
        << vs[i*3 + 0] << " "
        << vs[i*3 + 1] << " "
        << vs[i*3 + 2] << " "
        << "\n";

  for (size_t i = 0; i < nf; ++i)
    out << "f "
        << fs[i*4 + 0] + 1 << " "
        << fs[i*4 + 1] + 1 << " "
        << fs[i*4 + 2] + 1 << " "
        << fs[i*4 + 3] + 1 << " "
        << "\n";

  for (size_t i = 0; i < ne; ++i)
    out << "l "
        << es[i*2 + 0] + 1 << " "
        << es[i*2 + 1] + 1 << " "
        << "\n";

    return 0;
}

static void append(
  std::vector<double> &vs,
  std::vector<int> &fs,
  std::vector<int> &es,
  const Box3 &box)
{
  const auto nv = vs.size()/3;
  const auto nf = fs.size();
  const auto ne = es.size();
  const auto &v0 = box[0];
  const auto &v1 = box[1];

  vs.insert(vs.end(), { v0[0], v0[1], v0[2] });
  vs.insert(vs.end(), { v1[0], v0[1], v0[2] });
  vs.insert(vs.end(), { v0[0], v1[1], v0[2] });
  vs.insert(vs.end(), { v1[0], v1[1], v0[2] });
  vs.insert(vs.end(), { v0[0], v0[1], v1[2] });
  vs.insert(vs.end(), { v1[0], v0[1], v1[2] });
  vs.insert(vs.end(), { v0[0], v1[1], v1[2] });
  vs.insert(vs.end(), { v1[0], v1[1], v1[2] });

  fs.insert(fs.end(), { 0, 2, 3, 1 });
  fs.insert(fs.end(), { 4, 5, 7, 6 });
  fs.insert(fs.end(), { 0, 4, 6, 2 });
  fs.insert(fs.end(), { 1, 3, 7, 5 });
  fs.insert(fs.end(), { 3, 2, 6, 7 });
  fs.insert(fs.end(), { 1, 5, 4, 0 });

  es.insert(es.end(), { 0, 1 });
  es.insert(es.end(), { 2, 3 });
  es.insert(es.end(), { 4, 5 });
  es.insert(es.end(), { 6, 7 });
  es.insert(es.end(), { 0, 2 });
  es.insert(es.end(), { 1, 3 });
  es.insert(es.end(), { 4, 6 });
  es.insert(es.end(), { 5, 7 });
  es.insert(es.end(), { 0, 4 });
  es.insert(es.end(), { 1, 5 });
  es.insert(es.end(), { 2, 6 });
  es.insert(es.end(), { 3, 7 });

  for (size_t i = nf; i < fs.size(); ++i)
    fs[i] += (int)nv;

  for (size_t i = ne; i < es.size(); ++i)
    es[i] += (int)nv;
}

int main(int argc, const char **argv)
{
  int err {};
  if (argc < 2) return 1;

  std::string filename { argv[1] };
  std::string prefix = filename.substr(0, filename.find_last_of("."));
  std::string path = filename.substr(0, filename.find_last_of("/\\"));

  // read mesh data
  std::vector<Vec3> vs {};
  std::vector<Int3> fs {};
  err = read_obj(vs, fs, filename.c_str());
  if (err) return err;

  // build BVH
  Bvh<Box3> bvh {};
  TriangleBound bound(vs, fs);
  std::vector<int> fids {};
  for (int i = 0; i < fs.size(); ++i) fids.push_back(i);
  SAHSplit<TriangleBound, Box3, decltype(fids.end())> split(bound);
  bvh.build(bound, split, fids.begin(), fids.end());

  std::vector<double> vo {};
  std::vector<int> fo {};
  std::vector<int> eo {};
  for (const auto &node : bvh.nodes())
    append(vo, fo, eo, node.b);

  err = save_obj(
    vo.data(), vo.size()/3,
    fo.data(), fo.size()/4,
    eo.data(), eo.size()/2,
    (prefix + ".bvh.obj").c_str());

  return 0;
}