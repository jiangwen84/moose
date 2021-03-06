//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCutUserObject.h"
#include <array>
/**
 * InterfaceMeshCut3DUserObject: (1) reads in a mesh describing the interface,
 * (2) uses the mesh to do cutting of 3D elements, and
 * (3) grows the mesh based on intreface velocites.
 */

class XFEMMovingInterfaceVelocityBase;

class InterfaceMeshCut3DUserObject : public GeometricCutUserObject
{
public:
  static InputParameters validParams();

  InterfaceMeshCut3DUserObject(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  std::shared_ptr<MeshBase> getCutterMesh() const { return _cutter_mesh; };
  const auto * getPseudoNormal() const { return &_pseudo_normal; };

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes,
                                    Real time) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces,
                                    Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges,
                                     Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces,
                                     Real time) const override;

protected:
  /**
    Check if a line intersects with an element
   */
  virtual bool intersectWithEdge(const Point & p1,
                                 const Point & p2,
                                 const std::vector<Point> & _vertices,
                                 Point & pint) const;

  Point nodeNomal(const unsigned int & node_id);

  /**
    Check if point p is inside the edge p1-p2
   */
  bool isInsideEdge(const Point & p1, const Point & p2, const Point & p) const;

  /**
    Get the relative position of p from p1
   */
  Real getRelativePosition(const Point & p1, const Point & p2, const Point & p) const;

  /**
    Check if point p is inside a plane
   */
  bool isInsideCutPlane(const std::vector<Point> & _vertices, const Point & p) const;

  /// The cutter mesh
  std::shared_ptr<MeshBase> _cutter_mesh;

  /// node to element map of cut mesh
  std::map<dof_id_type, std::vector<dof_id_type>> _node_to_elem_map;

  /// Pointer to XFEMMovingInterfaceVelocityBase object
  const XFEMMovingInterfaceVelocityBase * _interface_velocity;

  /// Map of pseudo normal of element plane, three nodes and three sides of each element
  std::map<unsigned int, std::array<Point, 7>> _pseudo_normal;

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// The points to evaluate at
  std::vector<Point> _points;

  /// Pointer to PointLocatorBase object
  std::unique_ptr<PointLocatorBase> _pl;
};
