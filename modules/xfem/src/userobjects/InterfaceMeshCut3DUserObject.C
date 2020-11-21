//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceMeshCut3DUserObject.h"

#include "XFEMFuncs.h"
#include "MooseError.h"
#include "libmesh/string_to_enum.h"
#include "MooseMesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/plane.h"
#include "Function.h"

#include "libmesh/exodusII_io.h"

registerMooseObject("XFEMApp", InterfaceMeshCut3DUserObject);

InputParameters
InterfaceMeshCut3DUserObject::validParams()
{
  InputParameters params = GeometricCutUserObject::validParams();
  params.addRequiredParam<MeshFileName>(
      "mesh_file",
      "Mesh file for the XFEM geometric cut; currently only the xda type is supported");
  params.addRequiredParam<FunctionName>("velocity", "velocity function for normal direction");
  params.addClassDescription("Creates a UserObject for a mesh cutter in 3D problems");
  return params;
}

// this code does not allow predefined crack growth as a function of time
// all inital cracks are defined at t_start = t_end = 0
InterfaceMeshCut3DUserObject::InterfaceMeshCut3DUserObject(const InputParameters & parameters)
  : GeometricCutUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _velocity(getFunction("velocity"))
{
  // only the xda type is currently supported
  MeshFileName xfem_cut_mesh_file = getParam<MeshFileName>("mesh_file");
  _cut_mesh = libmesh_make_unique<ReplicatedMesh>(_communicator);
  _cut_mesh->read(xfem_cut_mesh_file);

  // test element type; only tri3 elements are allowed
  for (const auto & cut_elem : _cut_mesh->element_ptr_range())
  {
    if (cut_elem->dim() != _cut_elem_dim)
      mooseError("The input cut mesh should have 2D elements only!");
  }

  _exodus_io = std::make_shared<ExodusII_IO>(*_cut_mesh);
}

void
InterfaceMeshCut3DUserObject::initialSetup()
{
}

void
InterfaceMeshCut3DUserObject::initialize()
{
  if (_t_step == 0)
    return;

  std::vector<Point> new_position(_cut_mesh->n_nodes());
  for (const auto & node : _cut_mesh->node_ptr_range())
  {
    Point p = *node;
    p += _dt * findNormalatNode(*node) * _velocity.value(_t, p);
    new_position[node->id()] = p;
  }

  for (const auto & node : _cut_mesh->node_ptr_range())
  {
    _cut_mesh->node_ref(node->id()) = new_position[node->id()];
  }

  std::string name = "interface_cut_mesh_";
  name += std::to_string(_t_step);
  name += ".e";
  //_cut_mesh->write(name);
  // auto exodus_io = std::make_shared<ExodusII_IO>(*_cut_mesh);

  _exodus_io->write(name);
  // exodus_helper->write_timestep(_t_step, _t);
}

Point
InterfaceMeshCut3DUserObject::findNormalatNode(const Node & node)
{
  Point normal(0.0);
  unsigned num_elem_connect_to_node = 0;

  for (const auto & cut_elem : _cut_mesh->element_ptr_range())
  {
    std::vector<Point> vertices{
        cut_elem->node_ref(0), cut_elem->node_ref(1), cut_elem->node_ref(2)};

    for (auto & node_in_elem : cut_elem->node_ref_range())
    {
      if (node_in_elem == node)
      {
        num_elem_connect_to_node++;
        Plane elem_plane(vertices[0], vertices[1], vertices[2]);
        Point point = vertices[0];
        normal += elem_plane.unit_normal(vertices[0]);
        break;
      }
    }
  }

  return normal / num_elem_connect_to_node;
}

bool
InterfaceMeshCut3DUserObject::cutElementByGeometry(const Elem * /*elem*/,
                                                   std::vector<Xfem::CutEdge> & /*cut_edges*/,
                                                   std::vector<Xfem::CutNode> & /*cut_nodes*/,
                                                   Real /*time*/) const
{
  mooseError("invalid method for 3D mesh cutting");
  return false;
}

bool
InterfaceMeshCut3DUserObject::cutElementByGeometry(const Elem * elem,
                                                   std::vector<Xfem::CutFace> & cut_faces,
                                                   Real /*time*/) const
// With the crack defined by a planar mesh, this method cuts a solid element by all elements in
// the planar mesh
// TODO: Time evolving cuts not yet supported in 3D (hence the lack of use of the time variable)
{
  bool elem_cut = false;

  if (elem->dim() != _elem_dim)
    mooseError("The structural mesh to be cut by a surface mesh must be 3D!");

  for (unsigned int i = 0; i < elem->n_sides(); ++i)
  {
    // This returns the lowest-order type of side.
    std::unique_ptr<const Elem> curr_side = elem->side_ptr(i);
    if (curr_side->dim() != 2)
      mooseError("In cutElementByGeometry dimension of side must be 2, but it is ",
                 curr_side->dim());
    unsigned int n_edges = curr_side->n_sides();

    std::vector<unsigned int> cut_edges;
    std::vector<Real> cut_pos;

    for (unsigned int j = 0; j < n_edges; j++)
    {
      // This returns the lowest-order type of side.
      std::unique_ptr<const Elem> curr_edge = curr_side->side_ptr(j);
      if (curr_edge->type() != EDGE2)
        mooseError("In cutElementByGeometry face edge must be EDGE2, but type is: ",
                   libMesh::Utility::enum_to_string(curr_edge->type()),
                   " base element type is: ",
                   libMesh::Utility::enum_to_string(elem->type()));
      const Node * node1 = curr_edge->node_ptr(0);
      const Node * node2 = curr_edge->node_ptr(1);

      for (const auto & cut_elem : _cut_mesh->element_ptr_range())
      {
        std::vector<Point> vertices;

        for (auto & node : cut_elem->node_ref_range())
        {
          Point & this_point = node;
          vertices.push_back(this_point);
        }

        Point intersection;
        if (intersectWithEdge(*node1, *node2, vertices, intersection))
        {
          cut_edges.push_back(j);
          cut_pos.emplace_back(getRelativePosition(*node1, *node2, intersection));
        }
      }
    }

    // if two edges of an element are cut, it is considered as an element being cut
    if (cut_edges.size() == 2)
    {
      elem_cut = true;
      Xfem::CutFace mycut;
      mycut._face_id = i;
      mycut._face_edge.push_back(cut_edges[0]);
      mycut._face_edge.push_back(cut_edges[1]);
      mycut._position.push_back(cut_pos[0]);
      mycut._position.push_back(cut_pos[1]);
      cut_faces.push_back(mycut);
    }
  }
  return elem_cut;
}

bool
InterfaceMeshCut3DUserObject::cutFragmentByGeometry(
    std::vector<std::vector<Point>> & /*frag_edges*/,
    std::vector<Xfem::CutEdge> & /*cut_edges*/,
    Real /*time*/) const
{
  mooseError("invalid method for 3D mesh cutting");
  return false;
}

bool
InterfaceMeshCut3DUserObject::cutFragmentByGeometry(
    std::vector<std::vector<Point>> & /*frag_faces*/,
    std::vector<Xfem::CutFace> & /*cut_faces*/,
    Real /*time*/) const
{
  // TODO: Need this for branching in 3D
  mooseError("cutFragmentByGeometry not yet implemented for 3D mesh cutting");
  return false;
}

bool
InterfaceMeshCut3DUserObject::intersectWithEdge(const Point & p1,
                                                const Point & p2,
                                                const std::vector<Point> & vertices,
                                                Point & pint) const
{
  bool has_intersection = false;

  Plane elem_plane(vertices[0], vertices[1], vertices[2]);
  Point point = vertices[0];
  Point normal = elem_plane.unit_normal(point);

  std::array<Real, 3> plane_point = {{point(0), point(1), point(2)}};
  std::array<Real, 3> planenormal = {{normal(0), normal(1), normal(2)}};
  std::array<Real, 3> edge_point1 = {{p1(0), p1(1), p1(2)}};
  std::array<Real, 3> edge_point2 = {{p2(0), p2(1), p2(2)}};
  std::array<Real, 3> cut_point = {{0.0, 0.0, 0.0}};

  if (Xfem::plane_normal_line_exp_int_3d(
          &plane_point[0], &planenormal[0], &edge_point1[0], &edge_point2[0], &cut_point[0]) == 1)
  {
    Point temp_p(cut_point[0], cut_point[1], cut_point[2]);
    if (isInsideCutPlane(vertices, temp_p) && isInsideEdge(p1, p2, temp_p))
    {
      pint = temp_p;
      has_intersection = true;
    }
  }
  return has_intersection;
}

bool
InterfaceMeshCut3DUserObject::isInsideEdge(const Point & p1,
                                           const Point & p2,
                                           const Point & p) const
{
  Real dotp1 = (p1 - p) * (p2 - p1);
  Real dotp2 = (p2 - p) * (p2 - p1);
  return (dotp1 * dotp2 <= 0.0);
}

Real
InterfaceMeshCut3DUserObject::getRelativePosition(const Point & p1,
                                                  const Point & p2,
                                                  const Point & p) const
{
  Real full_len = (p2 - p1).norm();
  Real len_p1_p = (p - p1).norm();
  return len_p1_p / full_len;
}

bool
InterfaceMeshCut3DUserObject::isInsideCutPlane(const std::vector<Point> & vertices,
                                               const Point & p) const
{
  unsigned int n_node = vertices.size();

  Plane elem_plane(vertices[0], vertices[1], vertices[2]);
  Point normal = elem_plane.unit_normal(vertices[0]);

  bool inside = false;
  unsigned int counter = 0;

  for (unsigned int i = 0; i < n_node; ++i)
  {
    unsigned int iplus1 = (i < n_node - 1 ? i + 1 : 0);
    Point middle2p = p - 0.5 * (vertices[i] + vertices[iplus1]);
    const Point side_tang = vertices[iplus1] - vertices[i];
    Point side_norm = side_tang.cross(normal);
    Xfem::normalizePoint(middle2p);
    Xfem::normalizePoint(side_norm);
    if (middle2p * side_norm <= 0.0)
      counter += 1;
  }
  if (counter == n_node)
    inside = true;
  return inside;
}

const std::vector<Point>
InterfaceMeshCut3DUserObject::getCrackFrontPoints(unsigned int /*num_crack_front_points*/) const
{
  mooseError("getCrackFrontPoints() is not implemented for this object.");
}
