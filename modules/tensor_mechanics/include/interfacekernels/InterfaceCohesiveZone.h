/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef INTERFACECOHESIVEZONE_H
#define INTERFACECOHESIVEZONE_H

#include "InterfaceKernel.h"

// Forward Declarations
class InterfaceCohesiveZone;

template <>
InputParameters validParams<InterfaceCohesiveZone>();

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class InterfaceCohesiveZone : public InterfaceKernel
{
public:
  InterfaceCohesiveZone(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  /// Initial stiffness for the cohesive zone
  Real _stiffness;
  /// Max traction
  Real _max_traction;
  /// Strain energy release rate
  Real _Gc;

  const unsigned int _component;

  const MaterialProperty<Real> & _max_normal_separation;
  const MaterialProperty<Real> & _normal_separation;
};

#endif
