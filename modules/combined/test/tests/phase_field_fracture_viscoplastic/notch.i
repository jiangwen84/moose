[Mesh]
  type = FileMesh
  file = notch.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./All]
        add_variables = true
        strain = Finite
        additional_generate_output = stress_yy
        save_in = 'resid_x resid_y resid_z'
      [../]
    [../]
  [../]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = E_el
        mobility = L
        kappa = kappa_op
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
  [./resid_z]
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./peeq]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_work]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
    use_displaced_mesh = true
  [../]
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 2
    c = c
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./peeq]
    type = MaterialRealAux
    variable = peeq
    property = ep_eqv
    execute_on = timestep_end
  [../]
  [./plastic_work]
    type = MaterialRealAux
    variable = plastic_work
    property = plastic_work
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 2
    function = 't'
  [../]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = '1'
    value = 0
  [../]
[]

[UserObjects]
  [./flowstress]
    type = HEVPLinearHardening
    yield_stress = 0.25
    slope = 0
    intvar_prop_name = ep_eqv
  [../]
  [./flowrate]
    type = HEVPFlowRatePowerLawJ2
    reference_flow_rate = 0.0001
    flow_rate_exponent = 10.0
    flow_rate_tol = 1
    strength_prop_name = flowstress
  [../]
  [./ep_eqv]
     type = HEVPEqvPlasticStrain
     intvar_rate_prop_name = ep_eqv_rate
  [../]
  [./ep_eqv_rate]
     type = HEVPEqvPlasticStrainRate
     flow_rate_prop_name = flowrate
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '1 1'
  [../]
  [./pfgc]
    type = GenericFunctionMaterial
    prop_names = 'gc_prop'
    prop_values = '1e-2'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
  [./viscop_damage]
    type = HyperElasticNonSplitPhaseFieldDamage
    #type = HyperElasticPhaseFieldIsoDamage
    resid_abs_tol = 1e-18
    resid_rel_tol = 1e-8
    maxiters = 50
    max_substep_iteration = 1
    flow_rate_user_objects = 'flowrate'
    strength_user_objects = 'flowstress'
    internal_var_user_objects = 'ep_eqv'
    internal_var_rate_user_objects = 'ep_eqv_rate'
    numerical_stiffness = true
    damage_stiffness = 1e-8
    c = c
    F_name = E_el
    use_linear_fracture_energy = true
    use_plastic_work = true
    W0 = 1e-3
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
[]

[Postprocessors]
  [./resid_x]
    type = NodalSum
    variable = resid_x
    boundary = 2
  [../]
  [./resid_y]
    type = NodalSum
    variable = resid_y
    boundary = 2
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = 'ANY_BLOCK_ID'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
#petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
# petsc_options_value = 'asm      31                  preonly       lu           1'
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-8

  l_max_its = 10
  nl_max_its = 15

  dt = 0.0025
  dtmin = 1e-4
  num_steps = 10000
[]

[Outputs]
  exodus = true
  csv = true
[]
