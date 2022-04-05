//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExtendVelocityLevelSetAux.h"
#include "InterfaceMeshCutUserObjectBase.h"

registerMooseObject("XFEMApp", ExtendVelocityLevelSetAux);

InputParameters
ExtendVelocityLevelSetAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Extends velocity from an interface to a domain.");
  params.addParam<UserObjectName>(
      "qp_point_value_user_object",
      "Name of QpPointValueAtXFEMInterface that gives values at Qp points along an interface.");
  return params;
}

ExtendVelocityLevelSetAux::ExtendVelocityLevelSetAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _qp_value_uo(getUserObjectByName<QpPointValueAtXFEMInterface>(
        getParam<UserObjectName>("qp_point_value_user_object")))
{
  if (!isNodal())
    mooseError("ExtendVelocityLevelSetAux: Aux variable must be nodal variable.");
}

Real
ExtendVelocityLevelSetAux::computeValue()
{
  _values_positive_level_set_side = _qp_value_uo.getValueAtPositiveLevelSet();

  // for (auto const & qp : _values_positive_level_set_side)
  //   std::cout << "value = " << qp.second << std::endl;

  _qp_points = _qp_value_uo.getQpPoint();

  unsigned index = 0;
  Real min_dist = std::numeric_limits<Real>::max();
  for (auto const & qp : _qp_points)
  {
    Real dist = (*_current_node - qp.second).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      index = qp.first;
    }
  }

  return _values_positive_level_set_side[index];
}
