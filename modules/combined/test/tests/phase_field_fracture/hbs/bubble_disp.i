[Mesh]
  type = FileMesh
  file = b101_out.e
[]


#[Adaptivity]
# initial_steps = 1
# max_h_level = 1
# stop_time = 1.0e-10
# initial_marker = err_bnds
#[./Markers]
#   [./err_bnds]
#     type = ErrorFractionMarker
#     coarsen = 0
#     refine = 0.9
#     indicator = ind_bnds
#   [../]
# [../]
# [./Indicators]
#    [./ind_bnds]
#      type = GradientJumpIndicator
#      variable = bnds
#   [../]
# [../]
#[]

[GlobalParams]
  op_num = 9
  var_name_base = gr
  displacements = 'disp_x disp_y'
[]

[UserObjects]
  [./soln]
    type = SolutionUserObject
    mesh = b101_out.e
    timestep = 'LATEST'
    execute_on = 'initial'
  [../]
[]

[MultiApps]
  [damage]
    type = TransientMultiApp
    input_files = 'bubble_c.i'
  []
[]

[Transfers]
  [to_disp_x]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_x'
    variable = 'disp_x'
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_y'
    variable = 'disp_y'
  []
  [from_d]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = from_multiapp
    source_variable = 'd'
    variable = 'd'
  []
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  [d]
  []
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_xy stress_xx strain_xx strain_xy strain_yy strain_zz hydrostatic_stress mid_principal_stress min_principal_stress max_principal_stress'
        decomposition_method = EigenSolution
        save_in = 'force_x force_y'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./bnds]
    [./InitialCondition]
      type = FunctionIC
      function = f_bnds
      variable = bnds
    [../]
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 1.4
      y1 = 1.4
      radius = 0.5
      invalue = 1.0
      outvalue = 0.0
      int_width = 0.05
#    type = SpecifiedSmoothCircleIC
#    invalue = 1.0
#    outvalue = 0
#    int_width = 0.05
#    x_positions = '12.5 27.5'
#    z_positions = '0 0'
#    y_positions = '27.5 12.5'
#    radii = '5 5'
    [../]
  [../]
  [./gr0]
    [./InitialCondition]
      type = FunctionIC
      variable = gr0
      function = f_gr0
    [../]
  [../]
  [./gr1]
    [./InitialCondition]
      type = FunctionIC
      variable = gr1
      function = f_gr1
    [../]
  [../]
  [./gr2]
    [./InitialCondition]
      type = FunctionIC
      variable = gr2
      function = f_gr2
    [../]
  [../]
  [./gr3]
    [./InitialCondition]
      type = FunctionIC
      variable = gr3
      function = f_gr3
    [../]
  [../]
  [./gr4]
    [./InitialCondition]
      type = FunctionIC
      variable = gr4
      function = f_gr4
    [../]
  [../]
  [./gr5]
    [./InitialCondition]
      type = FunctionIC
      variable = gr5
      function = f_gr5
    [../]
  [../]
  [./gr6]
    [./InitialCondition]
      type = FunctionIC
      variable = gr6
      function = f_gr6
    [../]
  [../]
  [./gr7]
    [./InitialCondition]
      type = FunctionIC
      variable = gr7
      function = f_gr7
    [../]
  [../]
  [./gr8]
    [./InitialCondition]
      type = FunctionIC
      variable = gr8
      function = f_gr8
    [../]
  [../]

  # [./euler_angle]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]


[Functions]
  [./tfunc]
    type = ParsedFunction
    value = '0*t'
  [../]
 [./pressure]
   type = PiecewiseLinear
   data_file = 'bubble_pressure_r0.5_gas100_por10_ext0.csv'
   format = columns
 [../]
  # [./pressure]
  #   type = ParsedFunction
  #   value = '100'
  # [../]
  [./f_bnds]
    type = SolutionFunction
    solution = soln
    from_variable = 'bnds'
  [../]
  [./f_gr0]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr0'
  [../]
  [./f_gr1]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr1'
  [../]
  [./f_gr2]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr2'
  [../]
  [./f_gr3]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr3'
  [../]
  [./f_gr4]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr4'
  [../]
  [./f_gr5]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr5'
  [../]
  [./f_gr6]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr6'
  [../]
  [./f_gr7]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr7'
  [../]
  [./f_gr8]
    type = SolutionFunction
    solution = soln
    from_variable = 'gr8'
  [../]

[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = matrix_elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '0.01 1e-3'
  [../]
  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = fracture_pressure
    prop_values = pressure
    factor = 1e-6
  []
  [./gc]
    type = ParsedMaterial
    f_name = gc_prop
    #function = 'if(bnds < 0.75, if(bnds>0.25, 0.5, 2.5), 2.5)'
    function = 'if(bnds < 0.75 & c < 0.5, 0.0012, 0.012)'
    args = 'bnds c'
  [../]

  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l / 3.14159 * 2'
  [../]

  [./stress_void]
    type = ComputeLinearElasticStress
    block = 0
    base_name = void
  [../]
  [./strain_void]
    type = ComputeSmallStrain
    block = 0
    base_name = void
  [../]

  [./strain]
    type = ComputeSmallStrain
    block = 0
    base_name = matrix
  [../]

  [./damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
   decomposition_type = strain_vol_dev
   #decomposition_type = none
    use_snes_vi_solver = true
    base_name = matrix
  [../]
  [./indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'd'
    function = 'd*d'
    derivative_order = 2
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'd'
    function = '((1.0-d)^2+eta)/((1.0-d)^2+d*(1-0.5*d)*(4/3.14159/l*E*gc_prop/sigma^2))'
    material_property_names = 'gc_prop l'
    constant_names       = 'E sigma eta'
    constant_expressions = '385000 130 1e-4'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'd'
    material_property_names = 'gc_prop l'
    function = 'gc_prop/l/3.14159*(2*d-d^2)'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = d
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]
  [./const_stress]
    type = ComputeExtraStressConstant
    block = 0
    base_name = void
    extra_stress_tensor = '-1 -1 -1 0 0 0'
    prefactor = fracture_pressure
  [../]
  [./global_stress]
    type = TwoPhaseStressMaterial
    base_A = matrix
    base_B = void
  [../]
  [./switching]
    type = SwitchingFunctionMaterial
    eta = c
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '385000 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = matrix
  [../]
  [./elasticity_tensor_void]
    type = ComputeElasticityTensor
    C_ijkl = '3.85 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = void
  [../]
[]

[BCs]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
[]

[Postprocessors]
  [./ave_stress_top]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  [../]
  [./disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./react_y_top]
    type = NodalSum
    variable = force_y
    boundary = top
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -snes_type'
  petsc_options_value = 'lu superlu_dist vinewtonrsls'
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 20  ##max nonlinear iterations Previous:50
  start_time=0
  line_search = 'none'
  end_time = 2000
  dtmax = 1
  dtmin = 1e-14
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 10
    iteration_window = 0
    growth_factor = 1.2
    cutback_factor = 0.5
  [../]

  picard_max_its = 100
  picard_rel_tol = 1e-6
  picard_abs_tol = 1e-8
  accept_on_max_picard_iteration = true
[]

[Outputs]
  [exodus]
    type = Exodus
    interval = 1
  []
  csv = true
#gnuplot = true
[]
