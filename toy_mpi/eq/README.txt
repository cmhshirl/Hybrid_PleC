All protein diffusion are simulated in ODE, the rest are simulated in SSA
Reactions: CtrA <==> CtrAp are treated as partial equilibrium.
So there is no catalytic reaction type in this case.

when changing grid size, please change # define IJth(y,i,j) y((i-1)*10+j) in files modODE and modSSA

Change K_syn of DivJ to 25.0 instead of 12.5

change ID for reactions:
syn_OrimrDivJ, syn_DupmrDivJ, phos_DivK, deph_DivKp
in modPar.f90 when you change the reaction list in SSA part !!!