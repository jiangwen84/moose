//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMMATERIALTENSORMARKERUSEROBJECT_H
#define XFEMMATERIALTENSORMARKERUSEROBJECT_H

#include "XFEMMaterialStateMarkerBase.h"
#include "MaterialTensorCalculator.h"

class XFEMMaterialTensorMarkerUserObject;

template <>
InputParameters validParams<XFEMMaterialTensorMarkerUserObject>();

class XFEMMaterialTensorMarkerUserObject : public XFEMMaterialStateMarkerBase
{
public:
  XFEMMaterialTensorMarkerUserObject(const InputParameters & parameters);
  virtual ~XFEMMaterialTensorMarkerUserObject() {}

protected:
  bool _use_weibull;
  MaterialTensorCalculator _material_tensor_calculator;
  const MaterialProperty<SymmTensor> & _tensor;
  Real _threshold;
  bool _average;
  Real _random_range;

  const MaterialProperty<Real> & _weibull;
  
  virtual bool doesElementCrack(RealVectorValue &direction);
};

#endif // XFEMMATERIALTENSORMARKERUSEROBJECT_H
