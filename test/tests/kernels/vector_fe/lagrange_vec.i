# This example reproduces the libmesh vector_fe example 1 results

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 15
  ny = 15
  xmin = -1
  ymin = -1
  elem_type = QUAD9
[]

[Variables]
  [./u]
    family = LAGRANGE_VEC
    order = SECOND
  [../]
[]

[Adaptivity]
  steps = 1
  marker = box
  max_h_level = 2
  initial_steps = 2
  [./Markers]
    [./box]
      bottom_left = '0.3 0.3 0'
      inside = refine
      top_right = '0.6 0.6 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[Kernels]
  [./diff]
    type = VectorDiffusion
    variable = u
  [../]
  [./body_force]
    type = VectorBodyForce
    variable = u
    function_x = 'ffx'
    function_y = 'ffy'
  [../]
[]

[BCs]
  [./bnd]
    type = VectorFunctionDirichletBC
    variable = u
    function_x = 'x_exact_sln'
    function_y = 'y_exact_sln'
    function_z = '0'
    boundary = 'left right top bottom'
  [../]
[]

[Functions]
  [./x_exact_sln]
    type = ParsedFunction
    value = 'cos(.5*pi*x)*sin(.5*pi*y)'
  [../]
  [./y_exact_sln]
    type = ParsedFunction
    value = 'sin(.5*pi*x)*cos(.5*pi*y)'
  [../]
  [./ffx]
    type = ParsedFunction
    value = '.5*pi*pi*cos(.5*pi*x)*sin(.5*pi*y)'
  [../]
  [./ffy]
    type = ParsedFunction
    value = '.5*pi*pi*sin(.5*pi*x)*cos(.5*pi*y)'
  [../]
[]

[Preconditioning]
  [./pre]
    type = SMP
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
