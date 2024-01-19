//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetReinitialization.h"

registerMooseObject("LevelSetApp", LevelSetReinitialization);

InputParameters
LevelSetReinitialization::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The re-initialization equation.");
  params.addRequiredCoupledVar(
      "phi_0", "The level set variable to be reinitialized as signed distance function.");
  params.addRequiredParam<Real>(
      "epsilon", "The epsilon coefficient to be used in the reinitialization calculation.");
  return params;
}

LevelSetReinitialization::LevelSetReinitialization(const InputParameters & parameters)
  : ADKernel(parameters),
    _levelset_0(adCoupledValue("phi_0")),
    _grad_levelset_0(adCoupledGradient("phi_0")),
    _epsilon(parameters.get<Real>("epsilon"))
{
}

ADReal
LevelSetReinitialization::computeQpResidual()
{
  Real ls0 = MetaPhysicL::raw_value(_levelset_0[_qp]);
  Real grad_ls0_norm_sq = MetaPhysicL::raw_value(_grad_levelset_0[_qp]).norm_sq();
  Real sgn = ls0 / std::sqrt(ls0 * ls0 + grad_ls0_norm_sq * _epsilon * _epsilon);

  ADReal phi_norm = _grad_u[_qp].norm();

  return sgn * (phi_norm - 1.0) * _test[_i][_qp];
}
