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

// MOOSE includes
#include "AddControlAction.h"
#include "FEProblem.h"
#include "Factory.h"
#include "Control.h"

template<>
InputParameters validParams<AddControlAction>()
{
  InputParameters params = validParams<MooseObjectAction>();
  return params;
}

AddControlAction::AddControlAction(InputParameters parameters) :
    MooseObjectAction(parameters)
{
}

void
AddControlAction::act()
{
  _moose_object_pars.addPrivateParam<FEProblem *>("_fe_problem", _problem.get());
  MooseSharedPointer<Control> control = MooseSharedNamespace::static_pointer_cast<Control>(_factory.create(_type, _name, _moose_object_pars));
  _problem->getControlWarehouse().addObject(control);
}
