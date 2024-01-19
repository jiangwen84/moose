//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetBubble.h"

registerMooseObject("LevelSetApp", LevelSetBubble);

InputParameters
LevelSetBubble::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Implementation of 'bubble' ranging from 0 to 1.");
  params.addParam<RealVectorValue>(
      "center", RealVectorValue(0.5, 0.5, 0), "The center of the bubble.");
  params.addParam<Real>("radius", 0.15, "The radius of the bubble.");
  return params;
}

LevelSetBubble::LevelSetBubble(const InputParameters & parameters)
  : Function(parameters),
    _center(getParam<RealVectorValue>("center")),
    _radius(getParam<Real>("radius"))
{
}

Real
LevelSetBubble::value(Real /*t*/, const Point & p) const
{
  return (p - _center).norm() - _radius;
}

ADReal
LevelSetBubble::value(const ADReal & /*t*/, const ADPoint & p) const
{
  return (p - _center).norm() - _radius;
}

RealGradient
LevelSetBubble::gradient(Real /*t*/, const Point & p) const
{
  Real norm = (p - _center).norm();
  RealGradient output;

  for (const auto i : make_range(Moose::dim))
    output(i) = (p(i) - _center(i)) / norm;

  return output;
}
