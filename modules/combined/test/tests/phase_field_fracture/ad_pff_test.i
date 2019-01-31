[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 80
  ny = 40
  ymax = 0.5
  elem_type = QUAD4
[]

[MeshModifiers]
  [./noncrack]
    type = BoundingBoxNodeSet
    new_boundary = noncrack
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  #volumetric_locking_correction = true
#order = SECOND
[]

[ADKernels]
  [./stress_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
  [../]
  [./stress_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./c]
  [../]
[]

[ADKernels]
  [./ad_pf]
    type = ADPFFracture
    l_name = l
    gc = gc_prop
    variable = c
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l'
    prop_values = '1e-3 0.03'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
[]

[ADMaterials]
  [./strain]
    type = ADComputeSmallStrain
  [../]
  [./elastic]
    type = ADComputeIsotropicLinearElasticPFFractureStress
    c = c
    kdamage = 1.0e-4
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = 't'
  [../]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = noncrack
    value = 0
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = right
    value = 0
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
  # petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm      31                  preonly       lu           1'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-7
  l_tol = 1e-4
  l_max_its = 50
  nl_max_its = 20

  dtmin = 1e-7
  dt = 5e-5
  num_steps = 200
[]

[Outputs]
  exodus = true
  csv = true
[]
