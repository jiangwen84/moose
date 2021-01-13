[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y disp_z'
[]

[XFEM]
  geometric_cut_userobjects = 'cut_mesh'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 6
    ny = 6
    nz = 6
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
    elem_type = HEX8
  []
  [./box1]
    type = BoundingBoxNodeSetGenerator
    new_boundary = box1
    bottom_left = '0 0 0'
    top_right = '0.1 0.1 0.1'
    input = gen
  [../]
  [./box2]
    type = BoundingBoxNodeSetGenerator
    new_boundary = box2
    bottom_left = '0.8 0.8 0'
    top_right = '1.0 1.0 1.0'
    input = box1
  [../]
[]

[UserObjects]
  [./cut_mesh]
    type = InterfaceMeshCut3DUserObject
    mesh_file = sphere.e
    velocity = 0.0
    heal_always = true
  [../]
[]

# [UserObjects]
#   [./level_set_cut_uo]
#     type = LevelSetCutUserObject
#     level_set_var = ls
#   [../]
# []
#
# [AuxVariables]
#   [./ls]
#     order = FIRST
#     family = LAGRANGE
#   [../]
# []
#
# [AuxKernels]
#   [./ls_function]
#     type = FunctionAux
#     variable = ls
#     function = ls_func
#   [../]
# []

[Functions]
  [./ls_func]
    type = ParsedFunction
    value = 'sqrt(x*x+y*y+z*z)-0.4'
  [../]
[]

[Modules/TensorMechanics/Master]
  displacements = 'disp_x disp_y disp_z'
  [./all]
    strain = SMALL
    add_variables = true
    incremental = false
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = MatDiffusion
    variable = u
    diffusivity = 1
  [../]
  [./time_deriv]
    type = TimeDerivative
    variable = u
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 207000
    poissons_ratio = 0.3
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[BCs]
  [./front_u]
    type = DirichletBC
    variable = u
    boundary = front
    value = 0
  [../]
  [./back_u]
    type = DirichletBC
    variable = u
    boundary = back
    value = 1
  [../]
  [box1_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = box1
  []
  [box1_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = box1
  []
  [box1_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = box1
  []
  [box2_x]
    type = FunctionDirichletBC
    variable = disp_x
    function = '0.005*t'
    boundary = box2
  []
  [box2_y]
    type = FunctionDirichletBC
    variable = disp_y
    function = '0.005*t'
    boundary = box2
  []
  [box2_z]
    type = FunctionDirichletBC
    variable = disp_z
    function = '0.005*t'
    boundary = box2
  []
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 20
  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5

  start_time = 0.0
  dt = 1
  end_time = 1

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  csv = true
  perf_graph = true
  [./console]
    type = Console
    output_linear = true
  [../]
[]
