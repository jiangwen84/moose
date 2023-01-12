//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPermeabilityConstWithCracking.h"
#include "MathUtils.h"
#include "RankThreeTensor.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityConstWithCracking);

InputParameters
PorousFlowPermeabilityConstWithCracking::validParams()
{
  InputParameters params = PorousFlowPermeabilityConst::validParams();
  params.addClassDescription(
      "This Material calculates the permeability tensor with the sum of constant bulk intrinsic "
      "permeability tensor and anisotropic permeability due to cracking.");
  params.addRequiredCoupledVar("c", "Name of damage variable");
  // params.addRequiredCoupledVar(
  //     "displacements",
  //     "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<Real>(
      "fc", 1, "correction factor related to the roughness of the crack surfaces.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  return params;
}

PorousFlowPermeabilityConstWithCracking::PorousFlowPermeabilityConstWithCracking(
    const InputParameters & parameters)
  : PorousFlowPermeabilityConst(parameters),
    _c(coupledValue("c")),
    _grad_c(coupledGradient("c")),
    _c_var(coupled("c")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    // _ndisp(coupledComponents("displacements")),
    // _disp_var_num(coupledIndices("displacements")),
    _fc(getParam<Real>("fc"))
{
}

void
PorousFlowPermeabilityConstWithCracking::computeQpProperties()
{
  RealVectorValue normal =
      _grad_c[_qp] /
      (_grad_c[_qp] + RealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE)).norm();

  RankTwoTensor iden(RankTwoTensor::initIdentity);
  RankTwoTensor proj;
  proj.vectorOuterProduct(normal, normal);
  proj = iden - proj;

  Real w = _fc *
           std::abs(_current_elem->hmin() * (1.0 + (_mechanical_strain[_qp] * normal) * normal)) *
           MathUtils::heavyside(_c[_qp] - 0.5);

  _permeability_qp[_qp] = _input_permeability + _c[_qp] * _c[_qp] * proj * w * w / 12.0;
  _dpermeability_qp_dvar[_qp].assign(_num_var, RealTensorValue());
  _dpermeability_qp_dgradvar[_qp].resize(LIBMESH_DIM);
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
    _dpermeability_qp_dgradvar[_qp][i].assign(_num_var, RealTensorValue());

  if (_dictator.isPorousFlowVariable(_c_var))
  {
    const unsigned int pvar = _dictator.porousFlowVariableNum(_c_var);
    _dpermeability_qp_dvar[_qp][pvar] = 2.0 * _c[_qp] * proj * w * w / 12.0;

    RankThreeTensor dKdgrad_c =
        -w * w / 12.0 * (iden.mixedProductIkJ(normal) + iden.mixedProductJkI(normal)) /
            sqrt(normal * normal) +
        w * w / 6.0 * (RankTwoTensor::outerProduct(normal, normal)).mixedProductIJk(_grad_c[_qp]) /
            (normal * normal) +
        w / 3.0 * _fc * MathUtils::heavyside(_c[_qp] - 0.5) * _current_elem->hmin() *
            proj.mixedProductIJk(_mechanical_strain[_qp] * normal) / sqrt(normal * normal) -
        w / 3.0 * _fc * MathUtils::heavyside(_c[_qp] - 0.5) * _current_elem->hmin() *
            proj.mixedProductIJk((_mechanical_strain[_qp] * normal) * normal * _grad_c[_qp]) /
            (normal * normal);

    for (unsigned i = 0; i < LIBMESH_DIM; ++i)
      _dpermeability_qp_dgradvar[_qp][i][pvar] = RankTwoTensor(dKdgrad_c(0, 0, i),
                                                               dKdgrad_c(1, 0, i),
                                                               dKdgrad_c(2, 0, i),
                                                               dKdgrad_c(0, 1, i),
                                                               dKdgrad_c(1, 1, i),
                                                               dKdgrad_c(2, 1, i),
                                                               dKdgrad_c(0, 2, i),
                                                               dKdgrad_c(1, 2, i),
                                                               dKdgrad_c(2, 2, i));
  }

  // RankFourTensor dKde = w / 6.0 * _fc * MathUtils::heavyside(_c[_qp] - 0.5) *
  //                       _current_elem->hmin() *
  //                       proj.outerProduct(RankTwoTensor::outerProduct(normal, normal));

  // for (unsigned i = 0; i < _ndisp; ++i)
  //   if (_dictator.isPorousFlowVariable(_disp_var_num[i]))
  //   {
  //     // the i_th displacement is a PorousFlow variable
  //     const unsigned int pvar = _dictator.porousFlowVariableNum(_disp_var_num[i]);
  //   }
}
