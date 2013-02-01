clear all
close all
clc

Use Residual_Stress_3D_MLSEL_v1/residual_stress_3d_mls.m to generate K and Q matrices
Use Residual_Stress_3D_MLSEL_v1/residual_stress_3d.m to generate K and Q matrices
Need to edit inputs (DV GRID / MESH)
Need to vary kappa and BC penalty terms
Output is K and Q
Manually save K and Q
Move K and Q to material folder

FOR No MLS method
CRD of XRD DV DOES NOT MATTER (STEP 3)
HOWEVER, THE ORDERING OF XRD DV WITH RESPECT TO FE MESH MATTERS
WHEN GENERATING A MATRIX IN STEP 1 NOTE THE DV ORDER
MAKE SURE IN STEP4 THE MESH ORDER IS SAME AS DV ORDER
