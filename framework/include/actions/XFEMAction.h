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

#ifndef XFEMACTION_H
#define XFEMACTION_H

#include "Action.h"

class XFEMAction;

template<>
InputParameters validParams<XFEMAction>();


class XFEMAction: public Action
{
public:
  XFEMAction(InputParameters params);

  virtual void act();
protected:
  std::string _xfem_cut_type;
  std::vector<Real> _xfem_cut_data;
  std::vector<Real> _xfem_cut_scale;
  std::vector<Real> _xfem_cut_translate;
  std::string _xfem_qrule;
  std::string _order;
  std::string _family;
  bool _xfem_cut_plane;
};

#endif 
