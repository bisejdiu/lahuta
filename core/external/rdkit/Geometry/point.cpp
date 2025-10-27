// $Id$
//
// Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "point.h"
//#include <Numerics/Vector.h>

namespace RDGeom {
std::ostream& operator<<(std::ostream& target, const Point& pt) {
  for (unsigned int di = 0; di < pt.dimension(); ++di) {
    target << pt[di] << " ";
  }

  return target;
}

Point2D operator+(const Point2D& p1, const Point2D& p2) {
  Point2D res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  return res;
}

Point2D operator-(const Point2D& p1, const Point2D& p2) {
  Point2D res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  return res;
}

Point2D operator*(const Point2D& p1, double v) {
  Point2D res;
  res.x = p1.x * v;
  res.y = p1.y * v;
  return res;
}

Point2D operator/(const Point2D& p1, double v) {
  Point2D res;
  res.x = p1.x / v;
  res.y = p1.y / v;
  return res;
}

PointND operator+(const PointND& p1, const PointND& p2) {
  unsigned int dim;
  if (p1.dimension() < p2.dimension()) {
    dim = p1.dimension();
  } else {
    dim = p2.dimension();
  }
  PointND res(dim);
  for (unsigned int i = 0; i < dim; ++i) {
    res[i] = p1[i] + p2[i];
  }
  return res;
}
PointND operator-(const PointND& p1, const PointND& p2) {
  unsigned int dim;
  if (p1.dimension() < p2.dimension()) {
    dim = p1.dimension();
  } else {
    dim = p2.dimension();
  }
  PointND res(dim);
  for (unsigned int i = 0; i < dim; ++i) {
    res[i] = p1[i] - p2[i];
  }
  return res;
}

PointND operator*(const PointND& p1, double v) {
  PointND res(p1.dimension());
  for (unsigned int i = 0; i < p1.dimension(); ++i) {
    res[i] = p1[i] * v;
  }
  return res;
}

PointND operator/(const PointND& p1, double v) {
  PointND res(p1.dimension());
  for (unsigned int i = 0; i < p1.dimension(); ++i) {
    res[i] = p1[i] / v;
  }
  return res;
}
}  // namespace RDGeom
