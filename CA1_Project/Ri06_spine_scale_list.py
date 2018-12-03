# scaling factors
scale_factors = {'no_spine_scale' : 1.0, 'basal_scale' : 2.7, 'med_spine_rad_scale' : 1.2,
		 'med_spine_LM_scale' : 1.1, 'max_spine_rad_scale' : 1.6,
		 'thin_rad_spine_scale' : 2.6, 'thin_LM_spine_scale' : 1.2}

# Primary apical list
prim_apical = [
'somaA',
'dendA5_0',
'dendA5_01',
'dendA5_011',
'dendA5_0111',
'dendA5_01111',
'dendA5_011111',
'dendA5_0111111',
'dendA5_01111111',
'dendA5_011111111',
'dendA5_0111111111',
'dendA5_01111111111',
'dendA5_011111111111',
'dendA5_0111111111111',
'dendA5_01111111111111',
'dendA5_011111111111111',
'dendA5_0111111111111111',
'dendA5_01111111111111111']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
prim_apical = ['_' + item for item in prim_apical]
prim_apical = [item + '_' for item in prim_apical]

# no spine scale list
no_spine_scale = [
'dendA1_0',
'dendA1_01',
'dendA1_010',
'dendA2_0',
'dendA2_01',
'dendA2_010',
'dendA3_0',
'dendA3_00',
'dendA3_01',
'dendA4_0',
'dendA5_0',
'dendA5_01',
'dendA5_011']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
no_spine_scale = ['_' + item for item in no_spine_scale]
no_spine_scale = [item + '_' for item in no_spine_scale]

# basal spine scale list
basal_scale = [
'dendA1_0',
'dendA1_01',
'dendA1_010',
'dendA2_0',
'dendA2_01',
'dendA2_010',
'dendA3_0',
'dendA3_00',
'dendA3_01',
'dendA4_0']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
basal_scale = ['_' + item for item in basal_scale]
basal_scale = [item + '_' for item in basal_scale]

# medium spiny radial scale list
med_spine_rad_scale = [
'dendA5_0111',
'dendA5_01111',
'dendA5_011111',
'dendA5_0111111',
'dendA5_01111111',
'dendA5_011111111',
'dendA5_0111111111',
'dendA5_01111111111',
'dendA5_011111111111']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
med_spine_rad_scale = ['_' + item for item in med_spine_rad_scale]
med_spine_rad_scale = [item + '_' for item in med_spine_rad_scale]

# medium spiny LM scale list
med_spine_LM_scale = [
'dendA5_0111111111111111111',
'dendA5_011111111111111110']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
med_spine_LM_scale = ['_' + item for item in med_spine_LM_scale]
med_spine_LM_scale = [item + '_' for item in med_spine_LM_scale]

# max spiny rad scale list
max_spine_rad_scale = [
'dendA5_0111111111111',
'dendA5_01111111111111',
'dendA5_011111111111111',
'dendA5_0111111111111111',
'dendA5_01111111111111111',
'dendA5_011111111111111111']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
max_spine_rad_scale = ['_' + item for item in max_spine_rad_scale]
max_spine_rad_scale = [item + '_' for item in max_spine_rad_scale]

# thin rad scale list
thin_rad_spine_scale = [
'dendA5_0',
'dendA5_01',
'dendA5_011',
'dendA5_0111',
'dendA5_01111',
'dendA5_011111',
'dendA5_0111111',
'dendA5_01111111',
'dendA5_011111111',
'dendA5_0111111111',
'dendA5_01111111111',
'dendA5_011111111111',
'dendA5_0111111111111',
'dendA5_01111111111111',
'dendA5_011111111111111',
'dendA5_0111111111111111',
'dendA5_01111111111111111',
'dendA5_011111111111111111',
'dendA5_0111111111111111111',
'dendA5_011111111111111110',
'dendA5_01111111111111110',
'dendA5_0111111111111111110',
'dendA5_01111111111111111110',
'dendA5_011111111111111111100',
'dendA5_011111111111111111101',
'dendA5_0111111111111111111010',
'dendA5_0111111111111111111011',
'dendA5_01111111111111111111',
'dendA5_011111111111111111110',
'dendA5_011111111111111111111',
'dendA5_0111111111111111111110',
'dendA5_0111111111111111111111',
'dendA5_0111111111111111100',
'dendA5_0111111111111111101',
'dendA5_01111111111111111010',
'dendA5_01111111111111111011',
'dendA5_011111111111111110110',
'dendA5_011111111111111110111']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
thin_rad_spine_scale = ['_' + item for item in thin_rad_spine_scale]
thin_rad_spine_scale = [item + '_' for item in thin_rad_spine_scale]

# thin LM spine scale list
thin_LM_spine_scale = [
'dendA5_01111111111111110',
'dendA5_0111111111111111110',
'dendA5_01111111111111111110',
'dendA5_011111111111111111100',
'dendA5_011111111111111111101',
'dendA5_0111111111111111111010',
'dendA5_0111111111111111111011',
'dendA5_01111111111111111111',
'dendA5_011111111111111111110',
'dendA5_011111111111111111111',
'dendA5_0111111111111111111110',
'dendA5_0111111111111111111111',
'dendA5_0111111111111111100',
'dendA5_0111111111111111101',
'dendA5_01111111111111111010',
'dendA5_01111111111111111011',
'dendA5_011111111111111110110',
'dendA5_011111111111111110111']

# Add underscores to the beginning and end of each string to aid in proper compartment assignment in main 
# program
thin_LM_spine_scale = ['_' + item for item in thin_LM_spine_scale]
thin_LM_spine_scale = [item + '_' for item in thin_LM_spine_scale]
