# TaylorGreenParticles

<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_555.56.gif" width="900" height="412.5" alt="particles in taylor green vortex"/>
</p>

<p align="center">
Motions of particles in Taylor Green vortex with 
<img src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Ctau_p%3D555"> <br><br>
</p>

**TaylorGreenParticles** is a python script to simulate motions of particles in Taylor Green vortex. The animations can help to understand the phenomenon of preferential concentration of inertial particles in turbulence, see [Preferential concentration](https://en.wikipedia.org/wiki/Preferential_concentration). The fortran code is added by [Guo Chen](None).

<br>

## Solving the fluid phase
The analytical solution for 2D Taylor Green vortex is given by:
 
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbegin%7Balign*%7D%20u%28x%2Cy%2Ct%29%26%3D%20sin%28x%29cos%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20v%28x%2Cy%2Ct%29%26%3D-cos%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20%5Comega%28x%2Cy%2Ct%29%26%3D%20sin%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5Cend%7Balign*%7D">

where ``u``,``v``,``\omega`` are respectively the ``x-velocity``,``y-velocity`` and ``vorticity``.