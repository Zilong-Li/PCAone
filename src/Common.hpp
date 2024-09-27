#ifndef PCAone_Common_H
#define PCAone_Common_H

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

#define MAF(a) ((a) > 0.5 ? (1 - a) : (a))

typedef Eigen::MatrixXd MyMatrix;
typedef Eigen::VectorXd MyVector;
typedef Eigen::ArrayXXd MyArrayX;
typedef Eigen::ArrayXd MyArray;
typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long long uint64;
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using Double1D = std::vector<double>;
using String1D = std::vector<std::string>;
using UMapIntInt = std::unordered_map<int, int>;
using UMapIntDouble = std::unordered_map<int, double>;
using UMapIntString = std::unordered_map<int, std::string>;
using Pds = std::pair<double, std::string>;
using UMapIntPds = std::unordered_map<int, Pds>;

inline UMapIntInt vector2map(const Int1D& v) {
  UMapIntInt m;
  int i = 0;
  for (const auto& k : v) m.insert({k, i++});
  return m;
}


#endif  // PCAone_Common_H
