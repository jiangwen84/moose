//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "ADKernel.h"

/**
 * Implements the re-initialization equation.
 */
class LevelSetReinitialization : public ADKernel
{
public:
  static InputParameters validParams();

  LevelSetReinitialization(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// level set variable at time, \tau = 0.
  const ADVariableValue & _levelset_0;

  /// gradient level set variable at time, \tau = 0.
  const ADVariableGradient & _grad_levelset_0;

  /// regularization
  const Real _epsilon;
};
