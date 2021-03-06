[Mesh]
  file = brick3_mesh.e
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  type = ReferenceResidualProblem
  solution_variables = 'disp_x disp_y disp_z'
  reference_residual_variables = 'saved_x saved_y saved_z'
[]

[Variables]
  [./disp_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./disp_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./disp_z]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./penetration]
    order = SECOND
    family = LAGRANGE
  [../]
  [./saved_x]
    order = SECOND
  [../]
  [./saved_y]
    order = SECOND
  [../]
  [./saved_z]
    order = SECOND
  [../]
  [./diag_saved_x]
    order = SECOND
  [../]
  [./diag_saved_y]
    order = SECOND
  [../]
  [./diag_saved_z]
    order = SECOND
  [../]
  [./inc_slip_x]
    order = SECOND
  [../]
  [./inc_slip_y]
    order = SECOND
  [../]
  [./inc_slip_z]
    order = SECOND
  [../]
  [./accum_slip_x]
    order = SECOND
  [../]
  [./accum_slip_y]
    order = SECOND
  [../]
  [./accum_slip_z]
    order = SECOND
  [../]
[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    save_in_disp_x = saved_x
    save_in_disp_y = saved_y
    save_in_disp_z = saved_z
    diag_save_in_disp_x = diag_saved_x
    diag_save_in_disp_y = diag_saved_y
    diag_save_in_disp_z = diag_saved_z
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yy
    index = 1
  [../]
  [./stress_xy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xy
    index = 3
  [../]
  [./stress_zz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zz
    index = 2
  [../]
  [./inc_slip_x]
    type = PenetrationAux
    variable = inc_slip_x
    execute_on = timestep_begin
    boundary = 4
    paired_boundary = 3
  [../]
  [./inc_slip_y]
    type = PenetrationAux
    variable = inc_slip_y
    execute_on = timestep_begin
    boundary = 4
    paired_boundary = 3
  [../]
  [./accum_slip_x]
    type = PenetrationAux
    variable = accum_slip_x
    execute_on = timestep_end
    boundary = 4
    paired_boundary = 3
  [../]
  [./accum_slip_y]
    type = PenetrationAux
    variable = accum_slip_y
    execute_on = timestep_end
    boundary = 4
    paired_boundary = 3
  [../]
  [./penetration]
    type = PenetrationAux
    variable = penetration
    boundary = 4
    paired_boundary = 3
  [../]
[]

[Postprocessors]
  [./bot_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 1
  [../]
  [./bot_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 1
  [../]
  [./top_react_x]
    type = NodalSum
    variable = saved_x
    boundary = 5
  [../]
  [./top_react_y]
    type = NodalSum
    variable = saved_y
    boundary = 5
  [../]
  [./ref_resid_x]
    type = NodalL2Norm
    execute_on = timestep_end
    variable = saved_x
  [../]
  [./ref_resid_y]
    type = NodalL2Norm
    execute_on = timestep_end
    variable = saved_y
  [../]
  [./sigma_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./sigma_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./disp_x28]
    type = NodalVariableValue
    nodeid = 27
    variable = disp_x
  [../]
  [./disp_x33]
    type = NodalVariableValue
    nodeid = 32
    variable = disp_x
  [../]
  [./disp_y28]
    type = NodalVariableValue
    nodeid = 27
    variable = disp_y
  [../]
  [./disp_y33]
    type = NodalVariableValue
    nodeid = 32
    variable = disp_y
  [../]
  [./_dt]
    type = TimestepSize
  [../]
  [./num_lin_it]
    type = NumLinearIterations
  [../]
  [./num_nonlin_it]
    type = NumNonlinearIterations
  [../]
[]

[BCs]
  [./bot_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./side_x]
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0.0
  [../]
  [./back_z]
    type = DirichletBC
    variable = disp_z
    boundary = 6
    value = 0.0
  [../]
  [./top_press]
    type = Pressure
    variable = disp_y
    boundary = 5
    component = 1
    factor = 109.89
  [../]
[]

[Materials]
  [./bot]
    type = Elastic
    block = 1
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    poissons_ratio = 0.3
    youngs_modulus = 1e6
  [../]
  [./top]
    type = Elastic
    block = 2
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    poissons_ratio = 0.3
    youngs_modulus = 1e6
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  nl_abs_tol = 1e-9
  nl_rel_tol = 1e-8
  l_max_its = 50
  nl_max_its = 100
  dt = 1.0
  end_time = 1.0
  num_steps = 10
  dtmin = 0.1
  l_tol = 1e-5

[]

[VectorPostprocessors]
  [./x_disp]
    type = NodalValueSampler
    variable = disp_x
    boundary = '1 3 4 5'
    sort_by = id
  [../]
  [./y_disp]
    type = NodalValueSampler
    variable = disp_y
    boundary = '1 3 4 5'
    sort_by = id
  [../]
  [./cont_press]
    type = NodalValueSampler
    variable = contact_pressure
    boundary = '3'
    sort_by = id
  [../]
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./exodus]
    type = Exodus
    elemental_as_nodal = true
  [../]
  [./console]
    type = Console
    max_rows = 5
  [../]
  [./chkfile]
    type = CSV
    show = 'bot_react_x bot_react_y disp_x28 disp_y28 disp_x33 disp_y33 ref_resid_x ref_resid_y stress_yy stress_zz top_react_x top_react_y x_disp y_disp cont_press'
    execute_vector_postprocessors_on = timestep_end
  [../]
  [./outfile]
    type = CSV
    delimiter = ' '
    execute_vector_postprocessors_on = none
  [../]
[]

[Contact]
  [./leftright]
    slave = 3
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    master = 4
    order = SECOND
    system = constraint
    normalize_penalty = true
    tangential_tolerance = 1e-3
    penalty = 1e+9
  [../]
[]
