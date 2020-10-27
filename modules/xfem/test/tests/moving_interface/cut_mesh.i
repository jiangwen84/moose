[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  geometric_cut_userobjects = 'cut_mesh'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 1.0
  elem_type = HEX8
[]

[UserObjects]
  [./cut_mesh]
    type = MeshCut3DUserObject
    #mesh_file = mesh_edge_crack.xda
    mesh_file = mesh_in.e
    size_control = 0
    n_step_growth = 0
    function_x = 0
    function_y = 0
    function_z = 0
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

[BCs]
  [./right_u]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      8'

  l_max_its = 20
  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5

  start_time = 0.0
  dt = 1
  end_time = 2

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
