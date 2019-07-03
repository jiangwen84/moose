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
#include "MaterialAuxBase.h"

// Forward declarations
class MaterialRealDenseVectorAux;

template <>
InputParameters validParams<MaterialRealDenseVectorAux>();

/**
 * AuxKernel for outputting a DenseMatrix<Real> material property component to an AuxVariable
 */
class MaterialRealDenseVectorAux : public MaterialAuxBase<DenseVector<Real>>
{
public:
  /**
   * Class constructor
   * @param parameters The input parameters for this AuxKernel
   */
  MaterialRealDenseVectorAux(const InputParameters & parameters);

protected:
  /// Returns the component of the tensor for output
  virtual Real getRealValue() override;

  /// The  index to output
  unsigned int _index;
};
