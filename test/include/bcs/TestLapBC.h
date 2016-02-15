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

#ifndef TESTLAPBC_H
#define TESTLAPBC_H

#include "IntegratedBC.h"

//Forward Declarations
class TestLapBC;

template<>
InputParameters validParams<TestLapBC>();

class TestLapBC : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes
   * can be built using the same constructor.
   */
  TestLapBC(const InputParameters & parameters);

  virtual ~TestLapBC() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  VariableSecond & _second_u;
};

#endif //TESTLAPBC_H
