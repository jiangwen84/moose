//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionAux.h"
#include "Function.h"

registerMooseObject("MooseApp", FunctionAux);

defineLegacyParams(FunctionAux);

InputParameters
FunctionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Auxiliary Kernel that creates and updates a field variable by "
                             "sampling a function through space and time.");
  params.addRequiredParam<FunctionName>("function", "The function to use as the value");
  params.addCoupledVar("c",
                       "The name of the temperature variable used in the "
                       "ComputeThermalExpansionEigenstrain.  (Not required for "
                       "simulations without temperature coupling.)");
  return params;
}

FunctionAux::FunctionAux(const InputParameters & parameters)
  : AuxKernel(parameters), _func(getFunction("function")), _c(coupledValue("c"))
{
}

Real
FunctionAux::computeValue()
{
  if (isNodal())
  {
    if ((*_current_node)(0) > 200 && (*_current_node)(0) < 1000 && (*_current_node)(1) > 200 &&
        (*_current_node)(1) < 1000)
    {
      return _c[_qp];
    }
    else
      return 0;
  }
  else
    return _func.value(_t, _q_point[_qp]);
}
