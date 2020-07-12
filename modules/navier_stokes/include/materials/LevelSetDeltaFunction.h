#pragma once

#include "ADMaterial.h"
/**
 * This class computes delta function (derivative of the Heaviside function) given by a level set
 */
class LevelSetDeltaFunction : public ADMaterial
{
public:
  static InputParameters validParams();

  LevelSetDeltaFunction(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Gradient of the level set variable
  const ADVectorVariableValue & _grad_c;

  /// Delta function
  ADMaterialProperty<Real> & _delta_function;
};
