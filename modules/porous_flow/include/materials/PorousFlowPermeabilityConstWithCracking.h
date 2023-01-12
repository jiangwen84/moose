//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPermeabilityConst.h"

/**
 * Material designed to provide a constant permeability tensor
 */
class PorousFlowPermeabilityConstWithCracking : public PorousFlowPermeabilityConst
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityConstWithCracking(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Coupled order parameter defining the crack
  const VariableValue & _c;

  /// Gradient of the order parameter defining the crack
  const VectorVariableValue & _grad_c;

  const unsigned int _c_var;

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  /// Mechanical strain material property
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;

  //   /// Number of displacements supplied (1 in 1D, 2 in 2D, 3 in 3D)
  // const unsigned int _ndisp;

  // /// MOOSE variable number of the displacements variables provided
  // std::vector<unsigned int> _disp_var_num;

  /// correction factor for crack surface roughness
  const Real _fc;
};
