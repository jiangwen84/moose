//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VELOCITYOXIDEWEAKMICRO_H
#define XFEMC4VELOCITYOXIDEWEAKMICRO_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMC4VelocityOxideWeakMicro;

template <>
InputParameters validParams<XFEMC4VelocityOxideWeakMicro>();

class XFEMC4VelocityOxideWeakMicro : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityOxideWeakMicro(const InputParameters & parameters);
  virtual ~XFEMC4VelocityOxideWeakMicro() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;

  /** Real getVacancyFlux() const
  {
    return _vacancy_flux;
  }
  */
  //void computeVacancyFlux(unsigned int point_id);

protected:

  /// Diffusivity of oxygen in the Zr alpha phase
  Real _diffusivity_alpha;

  // Initial position of the oxide/metal interface
  Real _x0;

  //Real _vacancy_flux;
};

#endif // XFEMC4VELOCITYOXIDEWEAK_H
