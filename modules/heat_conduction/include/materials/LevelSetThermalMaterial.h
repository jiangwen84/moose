/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "ADMaterial.h"
#include "MeltPoolLevelSetLocation.h"

#define usingLevelSetThermalMateriallMembers usingMaterialMembers;

template <ComputeStage>
class LevelSetThermalMaterial;

declareADValidParams(LevelSetThermalMaterial);

/**
 *
 */
template <ComputeStage compute_stage>
class LevelSetThermalMaterial : public ADMaterial<compute_stage>
{
public:
  LevelSetThermalMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Level set variable
  const ADVariableValue & _ls;

  /// Temperature variable
  const ADVariableValue & _temp;

  /// Enthalpy
  ADMaterialProperty(Real) & _h;

  /// Thermal conductivity
  ADMaterialProperty(Real) & _k;

  /// Specifc heat
  ADMaterialProperty(Real) & _cp;

  /// Gas specific heat
  const Real & _c_g;

  /// Solid specific heat
  const Real & _c_s;

  /// Liquid specific heat
  const Real & _c_l;

  /// Gas thermal conductivity
  const Real & _k_g;

  /// Solid thermal conductivity
  const Real & _k_s;

  /// Liquid thermal conductivity
  const Real & _k_l;

  /// Latent heat
  const Real & _latent_heat;

  /// Solidus temperature
  const Real & _solidus_temperature;

  /// Liquid mass fraction
  const ADMaterialProperty(Real) & _f_l;

  /// Solid mass fraction
  const ADMaterialProperty(Real) & _f_s;

  /// Liquid volume fraction
  const ADMaterialProperty(Real) & _g_l;

  /// Solid volume fraction
  const ADMaterialProperty(Real) & _g_s;

  usingMaterialMembers;
};