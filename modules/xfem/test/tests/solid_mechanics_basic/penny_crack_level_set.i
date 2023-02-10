[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 21
  ny = 21
  nz = 21
  xmin = -1.0
  xmax = 1.0
  ymin = -1.0
  ymax = 1.0
  zmin = -1.0
  zmax = 1.0
  elem_type = HEX8
[]

[UserObjects]
  [level_set_cut_uo]
    type = CrackLevelSetCutUserObject
    level_set_phi = phi
    level_set_psi = psi
    heal_always = false
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
  []
[]

[Functions]
  [top_trac_z]
    type = ConstantFunction
    value = 10
  []
  [phi_func]
    type = ParsedFunction
    expression = 'z-0.0'
  []
  [psi_func]
    type = ParsedFunction
    expression = 'sqrt(x*x + y*y) - 0.3-0.1*t'
  []
[]

[AuxKernels]
  [phi_function]
    type = FunctionAux
    variable = phi
    function = phi_func
  []
  [psi_function]
    type = FunctionAux
    variable = psi
    function = psi_func
  []

[]

[AuxVariables]
  [phi]
    order = FIRST
    family = LAGRANGE
  []
  [psi]
    order = FIRST
    family = LAGRANGE
  []
[]

[BCs]
  [top_z]
    type = FunctionNeumannBC
    boundary = front
    variable = disp_z
    function = top_trac_z
  []
  [bottom_x]
    type = DirichletBC
    boundary = back
    variable = disp_x
    value = 0.0
  []
  [bottom_y]
    type = DirichletBC
    boundary = back
    variable = disp_y
    value = 0.0
  []
  [bottom_z]
    type = DirichletBC
    boundary = back
    variable = disp_z
    value = 0.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 207000
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      8'

  line_search = 'none'

  [Predictor]
    type = SimplePredictor
    scale = 1.0
  []

  # controls for linear iterations
  l_max_its = 100
  l_tol = 1e-2

  # controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  # time control
  start_time = 0.0
  dt = 1.0
  end_time = 4.0

  max_xfem_update = 1
[]

[Outputs]
  file_base = penny_crack_out
  execute_on = timestep_end
  exodus = true
  [console]
    type = Console
    output_linear = true
  []
[]
