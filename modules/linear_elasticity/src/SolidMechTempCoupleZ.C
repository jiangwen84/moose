/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SolidMechTempCoupleZ.h"

template<>
InputParameters validParams<SolidMechTempCoupleZ>()
{
  InputParameters params = validParams<SolidMechTempCouple>();
  return params;
}

SolidMechTempCoupleZ::SolidMechTempCoupleZ(const InputParameters & parameters)
  :SolidMechTempCouple(parameters)
{}

Real
SolidMechTempCoupleZ::computeQpResidual()
{
  if (!_constant_properties)
    recomputeCouplingConstants();

  return -(_c1*(1+2*_c2)*_grad_test[_i][_qp](2)*_thermal_strain[_qp]);
}

Real
SolidMechTempCoupleZ::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (!_constant_properties)
    recomputeCouplingConstants();

  if (jvar == _temp_var)
    return -(_c1*(1+2*_c2)*_grad_test[_i][_qp](2)*_alpha[_qp]*_phi[_j][_qp]);

  return 0.0;
}

