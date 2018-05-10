//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodeValueAtXFEMInterface.h"
#include "MooseVariableFE.h"
#include "XFEM.h"
#include "LineSegmentCutSetUserObject.h"
#include "InterfaceMeshCut3DUserObject.h"

#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

registerMooseObject("XFEMApp", NodeValueAtXFEMInterface);

InputParameters
NodeValueAtXFEMInterface::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addRequiredParam<VariableName>(
      "variable", "The name of the variable that this UserObject operates on");
  params.addParam<UserObjectName>(
      "interface_mesh_cut_userobject",
      "Name of InterfaceMeshCutUserObject that provides cut locations to this UserObject.");
  params.addRequiredParam<VariableName>(
      "level_set_var", "The name of level set variable used to represent the interface");
  params.addParam<bool>("is_3d", false, "3D case");
  params.addClassDescription("Obtain field values and gradients on the interface.");
  return params;
}

NodeValueAtXFEMInterface::NodeValueAtXFEMInterface(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_subproblem.mesh()),
    _var(&_subproblem.getVariable(_tid, parameters.get<VariableName>("variable"))),
    _level_set_var_number(
        _subproblem.getVariable(_tid, parameters.get<VariableName>("level_set_var")).number()),
    _system(_subproblem.getSystem(getParam<VariableName>("level_set_var"))),
    _solution(*_system.current_local_solution.get()),
    _is_3d(getParam<bool>("is_3d"))
{
}

void
NodeValueAtXFEMInterface::initialize()
{
  _pl = _mesh.getPointLocator();
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in NodeValueAtXFEMInterface");

  if (!_is_3d)
  {
    const UserObject * uo =
        &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

    if (dynamic_cast<const LineSegmentCutSetUserObject *>(uo) == nullptr)
      mooseError("UserObject casting to GeometricCutUserObject in XFEMSingleVariableConstraint");

    _geo_cut = dynamic_cast<const LineSegmentCutSetUserObject *>(uo);
    _elem_pairs = _xfem->getXFEMCutElemPairs(_xfem->getGeometricCutID(_geo_cut));
  }
  else
  {
    const UserObject * uo =
        &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

    if (dynamic_cast<const InterfaceMeshCut3DUserObject *>(uo) == nullptr)
      mooseError("UserObject casting to GeometricCutUserObject in XFEMSingleVariableConstraint");

    _geo_cut_3d = dynamic_cast<const InterfaceMeshCut3DUserObject *>(uo);
    _elem_pairs = _xfem->getXFEMCutElemPairs(_xfem->getGeometricCutID(_geo_cut_3d));
  }

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in PointValueAtXFEMInterface");
}

Point
PointValueAtXFEMInterface::getPointCurrentLocation(unsigned int i) const
{
  return _points[i];
};

void
NodeValueAtXFEMInterface::execute()
{
  _values_positive_level_set_side.clear();
  _values_negative_level_set_side.clear();
  _grad_values_positive_level_set_side.clear();
  _grad_values_negative_level_set_side.clear();
  _nodes.clear();

  if (!_is_3d)
  {
    std::vector<Real> cut_data = _geo_cut->getCutData();

    const int line_cut_data_len = 6;
    for (unsigned int i = 0; i < cut_data.size() / line_cut_data_len; ++i)
    {
      _points.push_back(
          Point(cut_data[i * line_cut_data_len + 0], cut_data[i * line_cut_data_len + 1]));
      if (i == cut_data.size() / line_cut_data_len - 1)
        _points.push_back(
            Point(cut_data[i * line_cut_data_len + 2], cut_data[i * line_cut_data_len + 3]));
    }
  }
  else
  {
    std::shared_ptr<MeshBase> cut_mesh = _geo_cut_3d->getCutMesh();

    for (const auto & node : cut_mesh->node_ptr_range())
    {
      Point p = *node;
      _points.push_back(p);
    }
  }

  BoundingBox bbox = _mesh.getInflatedProcessorBoundingBox(0.);

  std::vector<Point> point_vec(1);

  for (MooseIndex(_points) i = 0; i < _points.size(); ++i)
  {
    Point p = _points[i];

    if (bbox.contains_point(p))
    {
      const Elem * elem = getElemContainingPoint(*node, /*positive_level_set = */ true);

      std::cout << "i = " << i << ", p = " << p << std::endl;

      if (elem != nullptr)
      {
        _subproblem.setCurrentSubdomainID(elem, /*_tid */ 0);
        _subproblem.reinitElemPhys(elem, {*node}, 0);

        _values_positive_level_set_side[i] = (dynamic_cast<MooseVariable *>(_var))->sln()[0];
        _grad_values_positive_level_set_side[i] =
            ((dynamic_cast<MooseVariable *>(_var))->gradSln())[0];
      }

      const Elem * elem2 = getElemContainingPoint(*node, false);
      if (elem2 != nullptr)
      {
        _subproblem.setCurrentSubdomainID(elem2, /*_tid */ 0);
        _subproblem.reinitElemPhys(elem2, {*node}, 0);

        _values_negative_level_set_side[i] = (dynamic_cast<MooseVariable *>(_var))->sln()[0];
        _grad_values_negative_level_set_side[i] =
            ((dynamic_cast<MooseVariable *>(_var))->gradSln())[0];
      }
    }
    else // When node is outside of computation domain
    {
      _values_positive_level_set_side[i] = 0;
      _values_negative_level_set_side[i] = 0;
      _grad_values_positive_level_set_side[i] = RealVectorValue(0);
      _grad_values_negative_level_set_side[i] = RealVectorValue(0);
    }
  }
  // Take the value of the x component at only one point
  _grad_x_positive_level_set_side = _grad_values_positive_level_set_side[0](0);
  _grad_x_negative_level_set_side = _grad_values_negative_level_set_side[0](0);

  _current_x = _points[0](0);
}

void
NodeValueAtXFEMInterface::finalize()
{
  _communicator.set_union(_values_positive_level_set_side);
  _communicator.set_union(_grad_values_positive_level_set_side);
  _communicator.set_union(_values_negative_level_set_side);
  _communicator.set_union(_grad_values_negative_level_set_side);
  //_communicator.set_union(_grad_x_positive_level_set_side);
  //_communicator.set_union(_grad_x_negative_level_set_side);
}

const Elem *
NodeValueAtXFEMInterface::getElemContainingPoint(const Node & p, bool positive_level_set)
{
  std::cout << "p = " << p << std::endl;
  const Elem * elem1 = (*_pl)(p);

  if (elem1->processor_id() != processor_id())
    return nullptr;

  const Node * node = elem1->node_ptr(0);

  dof_id_type ls_dof_id = node->dof_number(_system.number(), _level_set_var_number, 0);

  Number ls_node_value = _solution(ls_dof_id);

  bool positive = false;

  if (_xfem->isPointInsidePhysicalDomain(elem1, *node))
  {
    if (ls_node_value > 0.0)
      positive = true;
  }
  else
  {
    if (ls_node_value < 0.0)
      positive = false;
  }

  const Elem * elem2 = nullptr;
  bool found = false;
  for (auto & pair : *_elem_pairs)
  {
    if (pair.first == elem1)
    {
      elem2 = pair.second;
      found = true;
    }
    else if (pair.second == elem1)
    {
      elem2 = pair.first;
      found = true;
    }
  }

  if (!found)
    mooseError("NodeValueAtXFEMInterface: The interface node ",
               p,
               " are not found by element pair locator.");

  if ((positive && positive_level_set) || (!positive && !positive_level_set))
    return elem1;
  else if ((!positive && positive_level_set) || (positive && !positive_level_set))
    return elem2;
  else
    return nullptr;
}
