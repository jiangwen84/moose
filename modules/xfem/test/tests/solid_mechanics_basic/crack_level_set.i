[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 11
  ny = 11
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  elem_type = QUAD4
[]

[UserObjects]
  [level_set_cut_uo]
    type = CrackLevelSetCutUserObject
    level_set_phi = phi
    level_set_psi = psi
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    planar_formulation = plane_strain
    add_variables = true
  []
[]

[Functions]
  [pull]
    type = PiecewiseLinear
    x = '0  50   100'
    y = '0  0.02 0.1'
  []
  [phi_func]
    type = ParsedFunction
    expression = 'y-0.5'
  []
  [psi_func]
    type = ParsedFunction
    expression = 'x-0.3-0.1*t'
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

[BCs]
  [bottomx]
    type = DirichletBC
    boundary = bottom
    variable = disp_x
    value = 0.0
  []
  [bottomy]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  []
  [topx]
    type = DirichletBC
    boundary = top
    variable = disp_x
    value = 0.0
  []
  [topy]
    type = FunctionDirichletBC
    boundary = top
    variable = disp_y
    function = pull
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
    block = 0
  []

  [_elastic_strain]
    type = ComputeFiniteStrainElasticStress
    block = 0
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
  nl_rel_tol = 1e-14
  nl_abs_tol = 1e-9

  # time control
  start_time = 0.0
  dt = 1.0
  end_time = 4.0
  num_steps = 5000

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  [console]
    type = Console
    output_linear = true
  []
[]
