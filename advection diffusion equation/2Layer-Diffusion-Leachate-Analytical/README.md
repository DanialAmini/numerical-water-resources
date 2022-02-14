# 2Layer-Diffusion-Leachate-Analytical
2 Layer Diffusion Equation of Leachate, Analytical Fourier solution of Li & Cleall 2010

leachate diffusion in two-layer media
<br />landfill on top, layer of impermeable clay of thickness h1 below landfill
<br />normal soil stratum of thickness h2 below clay layer

reference 1: Li, Y. C., & Cleall, P. J. (2010).
<br />Analytical solutions for contaminant diffusion in double-layered porous media.
<br />Journal of Geotechnical and Geoenvironmental Engineering (ASCE)
<br />Volume 136, No. 11, pages 1542-1554.<br />
https://ascelibrary.org/doi/abs/10.1061/(ASCE)GT.1943-5606.0000365

diffusion equation in each layare:
<br />(1/Rd1)* dc_1/dt = Ds1 * d^2c_1/dz^2,    0<z<h1
<br />(1/Rd2)* dc_2/dt = Ds2 * d^2c_2/dz^2,    h1<z<H

boundary condition:
<br />c_1=c0 @ z=0
<br />c_2=0 @ z=H

initial condition
<br />c_1=c_2=0 at t=0

interface boundary condition
<br />n1*Ds1*dc_1/dz=n2*Ds2*dc_2/dz  @ z=h1
<br />c_1=c_2 @ z=h1

H=h1+h2 : overal thickness
<br />Ds1 & Ds2 : diffusion coefficient
<br />n1 & n2: porosities
<br />Rd1 & Rd2: retardation factors
<br /> c1 & c2 : concentration in each layer
  
Implementation of the interface boundary condition (contunuous concentration and concentration mass flux) is done through the method explained in Ref. [2].

Reference 2: Hickson et al. (2011)
<br /> Hickson, R. I., Barry, S. I., Mercer, G. N., & Sidhu, H. S. (2011). Finite difference schemes for multilayer diffusion. Mathematical and Computer Modelling, 54(1-2), 210-220.
<br /> https://www.sciencedirect.com/science/article/pii/S0895717711000938
