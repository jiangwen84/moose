//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiscreteElementUserObject.h"
#include "NodeValueAtXFEMInterface.h"

class XFEMMovingInterfaceVelocityBase : public DiscreteElementUserObject
{
public:
  static InputParameters validParams();

  XFEMMovingInterfaceVelocityBase(const InputParameters & parameters);
  virtual ~XFEMMovingInterfaceVelocityBase() {}

  virtual void initialize() override;

  /**
   * Compute the interface velocity for a point
   * @param point_id  Point ID
   * @param normal  normal direction at this point
   * @return Real     Interface velocity
   */
  virtual Real computeMovingInterfaceVelocity(dof_id_type point_id,
                                              RealVectorValue normal) const = 0;

  /**
   * Compute total number of points that are used to define an interface
   */
  unsigned int numberPoints() const { return _value_at_interface_uo->numberPoints(); }

protected:
  /// Pointer to NodeValueAtXFEMInterface object
  const NodeValueAtXFEMInterface * _value_at_interface_uo;
};
