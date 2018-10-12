# TaylorGreenParticles

<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_555.56.gif" width="700"  alt="particles in taylor green vortex"/>
</p>

<p align="center">
Motions of particles in Taylor Green vortex with 
<img src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Ctau_p%3D555"> <br><br>
</p>

**TaylorGreenParticles** is a single python script to simulate motions of particles in Taylor Green vortex. The animations can help to understand the phenomenon of preferential concentration of inertial particles in turbulence, see [Preferential concentration](https://en.wikipedia.org/wiki/Preferential_concentration). The fortran code is added by [Guo Chen](None).

<br>

## Solve the fluid phase
The analytical solution for 2D Taylor Green vortex is given by:
<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbegin%7Balign*%7D%20u%28x%2Cy%2Ct%29%26%3D%20sin%28x%29cos%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20v%28x%2Cy%2Ct%29%26%3D-cos%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20%5Comega%28x%2Cy%2Ct%29%26%3D%20sin%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5Cend%7Balign*%7D"><br>
</p>

where ``u``,``v``,``\omega`` are respectively the ``x-velocity``,``y-velocity`` and ``vorticity``. In this simulation, we set both ``x`` and ``y`` in range ``[0,2*pi]``. And they are discretized into a mesh with shape of ``64*64`` (change it as you wish).

## Solve the particle phase
The particle dispersed in fluid experiences a hydrodynamic force, and this force will drive the particle's motion following the Newton's Second Law. The procedure to caculate the force lists as follows.

- First, caculate the particle Reynolds number ![re_p](http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20Re_p):

<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20Re_p%3D%7C%5Cmathbf%7BU%7D_f-%5Cmathbf%7BV%7D_p%7Cd_p/%5Cnu"><br>
</p>

- Next, caculate the drag coefficient ![drag_coff](http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbeta) :

<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbeta%3D1+0.15Re_p%5E%7B0.687%7D"><br>
</p>

- If we define \tau_p as ![tau_p](http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Ctau_%7Bp%7D%3D%5Crho_%7Bp%7Dd_%7Bp%7D%5E2/18%5Cnu), we can get:

<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cfrac%7Bd%5Cmathbf%7BV%7D_p%7D%7Bdt%7D%3D%5Cfrac%7B%5Cbeta%7D%7B%5Ctau_p%7D%28%5Cmathbf%7BU%7D_f-%5Cmathbf%7BV%7D_p%29"><br>
</p>





