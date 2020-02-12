[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.01
    ymin = 0
    ymax = 0.01
    nx = 64
    ny = 64
    elem_type = QUAD4
  []
[]

[Adaptivity]
  steps = 1
  marker = box
  max_h_level = 1
  initial_steps = 1
  stop_time = 1.0e-10
  [./Markers]
    [./box]
      bottom_left = '0.000 0.004 0'
      inside = refine
      top_right = '0.01 0.006 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[Variables]
  [./ls]
    order = FIRST
  [../]
  [./normal]
    family = LAGRANGE_VEC
  [../]
[]

[AuxVariables]
  [./ls_0]
    order = FIRST
  [../]
[]

[Kernels]
  [./grad_ls]
    type = LevelSetGradientRegulization
    level_set = ls_0
    variable = normal
    normal_regulization = true
  [../]
  [./time]
    type = TimeDerivative
    variable = ls
  [../]
  [./reinit]
    type = LevelSetOlssonReinitialization
    variable = ls
    phi_0 = ls_0
    epsilon = 0.00015 #0.022097
    grad_c = normal
  [../]
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [./arnold]
    type = LevelSetOlssonTerminator
    tol = 0.5
    min_steps = 10
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  num_steps = 10
  nl_abs_tol = 1e-14
  #scheme = crank-nicolson
  line_search = none
  petsc_options_iname = '-pc_type -pc_sub_type'
  petsc_options_value = 'asm      ilu'
  dt = 0.001
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.001
  #   optimal_iterations = 5
  #   growth_factor = 5
  # [../]
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
