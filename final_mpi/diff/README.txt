All protein diffusion are simulated in ODE, the rest are simulated in SSA

!Recompute a0 and propensities every 1 minute to correct the numerical error comes from numerous math operation.

when changing grid size, please change # define IJth(y,i,j) y((i-1)*10+j) in files modODE and modSSA