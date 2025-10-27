//
// Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef __RD_POINT_H__
#define __RD_POINT_H__
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <type_traits>
#include <cstddef>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <RDGeneral/Invariant.h>
#include <Numerics/Vector.h>
#include <boost/smart_ptr.hpp>

namespace RDGeom {

class Point {
  // this is the virtual base class, mandating certain functions
 public:
  virtual ~Point() {}

  virtual double operator[](unsigned int i) const = 0;
  virtual double &operator[](unsigned int i) = 0;

  virtual void normalize() = 0;
  virtual double length() const = 0;
  virtual double lengthSq() const = 0;
  virtual unsigned int dimension() const = 0;

  virtual Point *copy() const = 0;
};
#ifndef _MSC_VER
// g++ (at least as of v9.3.0) generates some spurious warnings from here.
// disable them
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// typedef class Point3D Point;
template <typename T = double>
class Point3DT {
  static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Point3DT only supports float or double");

 public:
  T x{0.0};
  T y{0.0};
  T z{0.0};

  Point3DT() = default;
  Point3DT(T xv, T yv, T zv) : x(xv), y(yv), z(zv) {}
  ~Point3DT() = default;
  Point3DT(const Point3DT &) = default;
  Point3DT &operator=(const Point3DT &) = default;

  template <typename U>
  explicit Point3DT(const Point3DT<U> &other)
      : x(static_cast<T>(other.x)),
        y(static_cast<T>(other.y)),
        z(static_cast<T>(other.z)) {}

  inline unsigned int dimension() const { return 3; }

  inline T operator[](unsigned int i) const {
    PRECONDITION(i < 3, "Invalid index on Point3D");
    if (i == 0) {
      return x;
    } else if (i == 1) {
      return y;
    } else {
      return z;
    }
  }

  inline T &operator[](unsigned int i) {
    PRECONDITION(i < 3, "Invalid index on Point3D");
    if (i == 0) {
      return x;
    } else if (i == 1) {
      return y;
    } else {
      return z;
    }
  }

  Point3DT &operator+=(const Point3DT &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Point3DT &operator-=(const Point3DT &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  Point3DT &operator*=(T scale) {
    x *= scale;
    y *= scale;
    z *= scale;
    return *this;
  }

  Point3DT &operator/=(T scale) {
    x /= scale;
    y /= scale;
    z /= scale;
    return *this;
  }

  Point3DT operator-() const {
    Point3DT res(x, y, z);
    res.x *= -1.0;
    res.y *= -1.0;
    res.z *= -1.0;
    return res;
  }

  void normalize() {
    T l = this->length();
    x /= l;
    y /= l;
    z /= l;
  }

  T length() const {
    T res = x * x + y * y + z * z;
    return static_cast<T>(std::sqrt(static_cast<double>(res)));
  }

  T lengthSq() const {
    T res = x * x + y * y + z * z;
    return res;
  }

  T dotProduct(const Point3DT &other) const {
    T res = x * (other.x) + y * (other.y) + z * (other.z);
    return res;
  }

  double angleTo(const Point3DT &other) const {
    double lsq = static_cast<double>(lengthSq()) * static_cast<double>(other.lengthSq());
    double dotProd = static_cast<double>(dotProduct(other));
    dotProd /= std::sqrt(lsq);

    if (dotProd <= -1.0) {
      return M_PI;
    }
    if (dotProd >= 1.0) {
      return 0.0;
    }

    return std::acos(dotProd);
  }

  double signedAngleTo(const Point3DT &other) const {
    double res = this->angleTo(other);
    if ((this->x * other.y - this->y * other.x) < -1e-6) {
      res = 2.0 * M_PI - res;
    }
    return res;
  }

  Point3DT directionVector(const Point3DT &other) const {
    Point3DT res;
    res.x = other.x - x;
    res.y = other.y - y;
    res.z = other.z - z;
    res.normalize();
    return res;
  }

  Point3DT crossProduct(const Point3DT &other) const {
    Point3DT res;
    res.x = y * (other.z) - z * (other.y);
    res.y = -x * (other.z) + z * (other.x);
    res.z = x * (other.y) - y * (other.x);
    return res;
  }

  Point3DT getPerpendicular() const {
    Point3DT res(static_cast<T>(0.0), static_cast<T>(0.0), static_cast<T>(0.0));
    if (x) {
      if (y) {
        res.y = -1 * x;
        res.x = y;
      } else if (z) {
        res.z = -1 * x;
        res.x = z;
      } else {
        res.y = 1;
      }
    } else if (y) {
      if (z) {
        res.z = -1 * y;
        res.y = z;
      } else {
        res.x = 1;
      }
    } else if (z) {
      res.x = 1;
    }
    double l = res.length();
    POSTCONDITION(l > 0.0, "zero perpendicular");
    res /= l;
    return res;
  }
};

template <typename T>
inline Point3DT<T> operator+(const Point3DT<T> &p1, const Point3DT<T> &p2) {
  Point3DT<T> res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  res.z = p1.z + p2.z;
  return res;
}

template <typename T>
inline Point3DT<T> operator-(const Point3DT<T> &p1, const Point3DT<T> &p2) {
  Point3DT<T> res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  res.z = p1.z - p2.z;
  return res;
}

template <typename T>
inline Point3DT<T> operator*(const Point3DT<T> &p1, T v) {
  Point3DT<T> res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  res.z = p1.z * v;
  return res;
}

template <typename T>
inline Point3DT<T> operator/(const Point3DT<T> &p1, T v) {
  Point3DT<T> res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  res.z = p1.z / v;
  return res;
}

template <typename T>
inline double computeDihedralAngle(const Point3DT<T> &pt1,
                                   const Point3DT<T> &pt2,
                                   const Point3DT<T> &pt3,
                                   const Point3DT<T> &pt4) {
  Point3DT<T> begEndVec = pt3 - pt2;
  Point3DT<T> begNbrVec = pt1 - pt2;
  Point3DT<T> crs1 = begNbrVec.crossProduct(begEndVec);

  Point3DT<T> endNbrVec = pt4 - pt3;
  Point3DT<T> crs2 = endNbrVec.crossProduct(begEndVec);

  double ang = crs1.angleTo(crs2);
  return ang;
}

template <typename T>
inline double computeSignedDihedralAngle(const Point3DT<T> &pt1,
                                         const Point3DT<T> &pt2,
                                         const Point3DT<T> &pt3,
                                         const Point3DT<T> &pt4) {
  Point3DT<T> begEndVec = pt3 - pt2;
  Point3DT<T> begNbrVec = pt1 - pt2;
  Point3DT<T> crs1 = begNbrVec.crossProduct(begEndVec);

  Point3DT<T> endNbrVec = pt4 - pt3;
  Point3DT<T> crs2 = endNbrVec.crossProduct(begEndVec);

  double ang = crs1.angleTo(crs2);

  Point3DT<T> crs3 = crs1.crossProduct(crs2);
  double dot = static_cast<double>(crs3.dotProduct(begEndVec));
  if (dot < 0.0) {
    ang *= -1;
  }

  return ang;
}

using Point3D = Point3DT<double>;
using Point3Df = Point3DT<float>;
typedef std::vector<Point3D> POINT3D_VECT;
typedef std::vector<Point3Df> POINT3D_VECT_F;

static_assert(sizeof(Point3Df) == 3 * sizeof(float), "Point3Df must stay tightly packed");
static_assert(alignof(Point3Df) == alignof(float), "Point3Df alignment must match float");

class Point2D : public Point {
 public:
  double x{0.0};
  double y{0.0};

  Point2D() {}
  Point2D(double xv, double yv) : x(xv), y(yv) {}
  ~Point2D() override = default;

  Point2D(const Point2D &other) : Point(other), x(other.x), y(other.y) {}
  //! construct from a Point3D (ignoring the z coordinate)
  /*Point2D(const Point3D &p3d) : Point(p3d), x(p3d.x), y(p3d.y) {}*/

  Point *copy() const override { return new Point2D(*this); }

  inline unsigned int dimension() const override { return 2; }

  inline double operator[](unsigned int i) const override {
    PRECONDITION(i < 2, "Invalid index on Point2D");
    if (i == 0) {
      return x;
    } else {
      return y;
    }
  }

  inline double &operator[](unsigned int i) override {
    PRECONDITION(i < 2, "Invalid index on Point2D");
    if (i == 0) {
      return x;
    } else {
      return y;
    }
  }

  Point2D &operator=(const Point2D &other) {
    x = other.x;
    y = other.y;
    return *this;
  }

  Point2D &operator+=(const Point2D &other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Point2D &operator-=(const Point2D &other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  Point2D &operator*=(double scale) {
    x *= scale;
    y *= scale;
    return *this;
  }

  Point2D &operator/=(double scale) {
    x /= scale;
    y /= scale;
    return *this;
  }

  Point2D operator-() const {
    Point2D res(x, y);
    res.x *= -1.0;
    res.y *= -1.0;
    return res;
  }

  void normalize() override {
    double ln = this->length();
    x /= ln;
    y /= ln;
  }

  void rotate90() {
    double temp = x;
    x = -y;
    y = temp;
  }

  double length() const override {
    // double res = pow(x,2) + pow(y,2);
    double res = x * x + y * y;
    return sqrt(res);
  }

  double lengthSq() const override {
    double res = x * x + y * y;
    return res;
  }

  double dotProduct(const Point2D &other) const {
    double res = x * (other.x) + y * (other.y);
    return res;
  }

  double angleTo(const Point2D &other) const {
    Point2D t1, t2;
    t1 = *this;
    t2 = other;
    t1.normalize();
    t2.normalize();
    double dotProd = t1.dotProduct(t2);
    // watch for roundoff error:
    if (dotProd < -1.0) {
      dotProd = -1.0;
    } else if (dotProd > 1.0) {
      dotProd = 1.0;
    }
    return acos(dotProd);
  }

  double signedAngleTo(const Point2D &other) const {
    double res = this->angleTo(other);
    if ((this->x * other.y - this->y * other.x) < -1e-6) {
      res = 2.0 * M_PI - res;
    }
    return res;
  }

  Point2D directionVector(const Point2D &other) const {
    Point2D res;
    res.x = other.x - x;
    res.y = other.y - y;
    res.normalize();
    return res;
  }
};

class PointND : public Point {
 public:
  typedef boost::shared_ptr<RDNumeric::Vector<double>> VECT_SH_PTR;

  PointND(unsigned int dim) {
    RDNumeric::Vector<double> *nvec = new RDNumeric::Vector<double>(dim, 0.0);
    dp_storage.reset(nvec);
  }

  PointND(const PointND &other) : Point(other) {
    RDNumeric::Vector<double> *nvec =
        new RDNumeric::Vector<double>(*other.getStorage());
    dp_storage.reset(nvec);
  }

  Point *copy() const override { return new PointND(*this); }

#if 0
	template <typename T>
    PointND(const T &vals){
      RDNumeric::Vector<double> *nvec = new RDNumeric::Vector<double>(vals.size(), 0.0);
      dp_storage.reset(nvec);
      unsigned int idx=0;
      typename T::const_iterator it;
      for(it=vals.begin();
          it!=vals.end();
          ++it){
        nvec->setVal(idx,*it);
        ++idx;
      };
    };
#endif

  ~PointND() override = default;

  inline double operator[](unsigned int i) const override {
    return dp_storage.get()->getVal(i);
  }

  inline double &operator[](unsigned int i) override {
    return (*dp_storage.get())[i];
  }

  inline void normalize() override { dp_storage.get()->normalize(); }

  inline double length() const override { return dp_storage.get()->normL2(); }

  inline double lengthSq() const override {
    return dp_storage.get()->normL2Sq();
  }

  unsigned int dimension() const override { return dp_storage.get()->size(); }

  PointND &operator=(const PointND &other) {
    if (this == &other) {
      return *this;
    }

    RDNumeric::Vector<double> *nvec =
        new RDNumeric::Vector<double>(*other.getStorage());
    dp_storage.reset(nvec);
    return *this;
  }

  PointND &operator+=(const PointND &other) {
    (*dp_storage.get()) += (*other.getStorage());
    return *this;
  }

  PointND &operator-=(const PointND &other) {
    (*dp_storage.get()) -= (*other.getStorage());
    return *this;
  }

  PointND &operator*=(double scale) {
    (*dp_storage.get()) *= scale;
    return *this;
  }

  PointND &operator/=(double scale) {
    (*dp_storage.get()) /= scale;
    return *this;
  }

  PointND directionVector(const PointND &other) {
    PRECONDITION(this->dimension() == other.dimension(),
                 "Point dimensions do not match");
    PointND np(other);
    np -= (*this);
    np.normalize();
    return np;
  }

  double dotProduct(const PointND &other) const {
    return dp_storage.get()->dotProduct(*other.getStorage());
  }

  double angleTo(const PointND &other) const {
    double dp = this->dotProduct(other);
    double n1 = this->length();
    double n2 = other.length();
    if ((n1 > 1.e-8) && (n2 > 1.e-8)) {
      dp /= (n1 * n2);
    }
    if (dp < -1.0) {
      dp = -1.0;
    } else if (dp > 1.0) {
      dp = 1.0;
    }
    return acos(dp);
  }

 private:
  VECT_SH_PTR dp_storage;
  inline const RDNumeric::Vector<double> *getStorage() const {
    return dp_storage.get();
  }
};
#ifndef _MSC_VER
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

typedef std::vector<RDGeom::Point *> PointPtrVect;
typedef PointPtrVect::iterator PointPtrVect_I;
typedef PointPtrVect::const_iterator PointPtrVect_CI;

typedef std::vector<RDGeom::Point3D *> Point3DPtrVect;
typedef std::vector<RDGeom::Point2D *> Point2DPtrVect;
typedef Point3DPtrVect::iterator Point3DPtrVect_I;
typedef Point3DPtrVect::const_iterator Point3DPtrVect_CI;
typedef Point2DPtrVect::iterator Point2DPtrVect_I;
typedef Point2DPtrVect::const_iterator Point2DPtrVect_CI;

typedef std::vector<const RDGeom::Point3D *> Point3DConstPtrVect;
typedef Point3DConstPtrVect::iterator Point3DConstPtrVect_I;
typedef Point3DConstPtrVect::const_iterator Point3DConstPtrVect_CI;

typedef std::vector<Point3D> POINT3D_VECT;
typedef std::vector<Point3D>::iterator POINT3D_VECT_I;
typedef std::vector<Point3D>::const_iterator POINT3D_VECT_CI;

typedef std::map<int, Point2D> INT_POINT2D_MAP;
typedef INT_POINT2D_MAP::iterator INT_POINT2D_MAP_I;
typedef INT_POINT2D_MAP::const_iterator INT_POINT2D_MAP_CI;

std::ostream &operator<<(std::ostream &target,
                                                    const RDGeom::Point &pt);

RDGeom::Point2D operator+(const RDGeom::Point2D &p1,
                                                     const RDGeom::Point2D &p2);
RDGeom::Point2D operator-(const RDGeom::Point2D &p1,
                                                     const RDGeom::Point2D &p2);
RDGeom::Point2D operator*(const RDGeom::Point2D &p1,
                                                     double v);
RDGeom::Point2D operator/(const RDGeom::Point2D &p1,
                                                     double v);

RDGeom::PointND operator+(const RDGeom::PointND &p1,
                                                     const RDGeom::PointND &p2);
RDGeom::PointND operator-(const RDGeom::PointND &p1,
                                                     const RDGeom::PointND &p2);
RDGeom::PointND operator*(const RDGeom::PointND &p1,
                                                     double v);
RDGeom::PointND operator/(const RDGeom::PointND &p1,
                                                     double v);
}  // namespace RDGeom

#endif
