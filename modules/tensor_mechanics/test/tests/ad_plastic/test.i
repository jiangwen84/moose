[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0
    ymin = 0
    zmin = 0
    xmax = 1
    ymax = 1
    zmax = 1
    # second_order = true
  []
  [extra_nodeset]
    type = ExtraNodesetGenerator
    input = gen
    new_boundary = 'fix'
    coord = '0 0 0'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = false
[]

[AuxVariables]
  [hydrostatic_stress]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [hydrostatic_stress]
    type = ADRankTwoScalarAux
    variable = hydrostatic_stress
    rank_two_tensor = stress
    scalar_type = Hydrostatic
  []
[]

[Variables]
  [disp_x]
    order = FIRST
    scaling = 1e-10
  []
  [disp_y]
    order = FIRST
    scaling = 1e-10
  []
  [disp_z]
    order = FIRST
    scaling = 1e-10
  []
[]

[Functions]
  [pull]
    type = PiecewiseLinear
    x = '0 10'
    y = '0 1e-3'
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    generate_output = 'elastic_strain_zz stress_zz'
    use_automatic_differentiation = true
    add_variables = true
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.25
  []
  # [strain]
  #   type = ADComputeIncrementalSmallStrain
  # []

  [elastic_strain]
    type = ADComputeMultipleInelasticStress
    inelastic_models = creep_ten
    internal_solve_full_iteration_history = true
  []

  [creep_ten]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1 #10e-24
    n_exponent = 2
    activation_energy = 0
    base_name = creep_ten
    internal_solve_full_iteration_history = true
    internal_solve_output_on = always
  []
  # [elastic_strain]
  #   type = ADComputeLinearElasticStress
  # []
[]

[BCs]
  [no_disp_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  []

  [no_disp_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  []

  [no_disp_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  []

  [pull_disp_z]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = front
    function = 0.02 #pull
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  solve_type = newton

  petsc_options_iname = -pc_type
  petsc_options_value = lu

  line_search = 'none'
  nl_rel_tol = 1e-10

  nl_forced_its = 2

  num_steps = 1
  dt = 1 #1e-1
[]

[Postprocessors]
  [max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
  []
  [max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
  []
  [max_hydro]
    type = ElementAverageValue
    variable = hydrostatic_stress
  []
  [dt]
    type = TimestepSize
  []
  [num_lin]
    type = NumLinearIterations
    outputs = console
  []
  [num_nonlin]
    type = NumNonlinearIterations
    outputs = console
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
