//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class InterfaceMeshCut3DUserObject;

/**
 * Calculate level set values for an interface that is defined by a set of line segments
 */
class MeshCutLevelSetAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MeshCutLevelSetAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /**
   * calculate the signed distance value for a given point.
   * @param p Coordinate of point
   * @return Signed distance
   */
  Real calculateSignedDistance(Point p);

  /**
   * calculate the distance from a point to a line segment.
   * @param x1,x2 Coordinates of line segment end points
   * @param x0 Coordinate of the point
   * @param xp Closest point coordinate on the line segment
   * @return Distance from a point x0 to a line segment defined by x1-x2
   */
  Real pointSegmentDistance(const Point & x0, const Point & x1, const Point & x2, Point & xp);

  /**
   * calculate the distance from a point to triangle.
   * @param x1,x2,x3 Coordinates of triangle vertices
   * @param x0 Coordinate of the point
   * @param xp Closest point coordinate on the triangle
   * @param region The seven regions where the closest point might locate
   * @return distance from a point x0 to a triangle defined by x1-x2-x3
   */

  //        R1
  //         1
  //        *  *
  //   R4  *     * R6
  //     *    R0  *
  //    *           *
  //   2  *  * *  *  3
  // R2       R5       R3

  Real pointTriangleDistance(const Point & x0,
                             const Point & x1,
                             const Point & x2,
                             const Point & x3,
                             Point & xp,
                             unsigned int & region);

  /// Pointer to the InterfaceMeshCut3DUserObject object
  const InterfaceMeshCut3DUserObject * _mesh_cut_uo;

  /// The structural mesh
  MooseMesh & _mesh;

  /// map of pseudo normal
  const std::map<unsigned int, std::array<Point, 7>> * _pseudo_normal;

  /// revert the level set sign
  const bool & _revert_sign;
};
