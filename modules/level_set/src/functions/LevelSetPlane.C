//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetPlane.h"
#include "libmesh/utility.h"

registerMooseObject("LevelSetApp", LevelSetPlane);

InputParameters
LevelSetPlane::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Implementation of a level set function to represent a plane.");
  params.addParam<RealVectorValue>("point", RealVectorValue(0, 0, 0), "A point on the plane.");
  params.addParam<RealVectorValue>(
      "normal", RealVectorValue(0, 1, 0), "The normal vector to the plane.");
  return params;
}

LevelSetPlane::LevelSetPlane(const InputParameters & parameters)
  : Function(parameters),
    _point(getParam<RealVectorValue>("point")),
    _normal(getParam<RealVectorValue>("normal"))
{
}

Real
LevelSetPlane::value(Real /*t*/, const Point & p) const
{
  const RealVectorValue unit_normal = _normal / _normal.norm();
  const Real distance_from_orgin = -unit_normal * _point;
  return -(unit_normal * p + distance_from_orgin);
}

RealGradient
LevelSetPlane::gradient(Real /*t*/, const Point & /*p*/) const
{
  const RealVectorValue unit_normal = _normal / _normal.norm();

  RealGradient output;

  for (const auto i : make_range(Moose::dim))
    output(i) = -unit_normal(i);

  return output;
}
