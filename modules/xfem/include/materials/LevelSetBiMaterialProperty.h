/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LEVELSETBIMATERIALPROPERTY_H
#define LEVELSETBIMATERIALPROPERTY_H

#include "LevelSetBiMaterialBase.h"

// Forward Declarations
class LevelSetBiMaterialProperty;
class XFEM;

template <>
InputParameters validParams<LevelSetBiMaterialProperty>();

/**
 * Compute material property for bi-materials problem (consisting of two different materials)
 * defined by a level set function
 *
 */
class LevelSetBiMaterialProperty : public LevelSetBiMaterialBase
{
public:
  LevelSetBiMaterialProperty(const InputParameters & parameters);

protected:
  virtual void assignQpPropertiesForLevelSetPositive();
  virtual void assignQpPropertiesForLevelSetNegative();

  /// Material properties for the two separate materials in the bi-material system
  std::vector<const MaterialProperty<Real> *> _bimaterial_material_prop;

  /// Property name
  std::string _prop_name;

  /// Global material property (switch bi-material diffusion coefficient based on level set values)
  MaterialProperty<Real> & _material_prop;
};

#endif // LEVELSETBIMATERIALPROPERTY_H
