//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMEQUALVALUEATINTERFACE_H
#define XFEMEQUALVALUEATINTERFACE_H

#include "ElemElemConstraint.h"

// Forward Declarations
class XFEMEqualValueAtInterface;
class XFEM;

template <>
InputParameters validParams<XFEMEqualValueAtInterface>();

class XFEMEqualValueAtInterface : public ElemElemConstraint
{
public:
  XFEMEqualValueAtInterface(const InputParameters & parameters);
  virtual ~XFEMEqualValueAtInterface();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  // Penalty parameter in penalty's formulation
  Real _alpha;

  /// Value at the interface
  Real _value;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};

#endif /* XFEMEQUALVALUEATINTERFACE_H_ */
