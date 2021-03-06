# TaylorGreenParticles

<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_555.56.gif" width="700"  alt="particles in taylor green vortex"/>
</p>

<p align="center">
Motions of particles in Taylor Green vortex with 
<img src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Ctau_p%3D555"> <br><br>
</p>

**TaylorGreenParticles** is a single python script to simulate motions of particles in Taylor Green vortex. The animations can help to understand the phenomenon of preferential concentration of inertial particles in turbulence, see [Preferential concentration](https://en.wikipedia.org/wiki/Preferential_concentration). The fortran code is added by [Guo Chen](https://github.com/chenguo960627).

<br>

## Solve the fluid phase
The analytical solution for 2D Taylor Green vortex is given by:
<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbegin%7Balign*%7D%20u%28x%2Cy%2Ct%29%26%3D%20sin%28x%29cos%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20v%28x%2Cy%2Ct%29%26%3D-cos%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5C%5C%20%5Comega%28x%2Cy%2Ct%29%26%3D%20sin%28x%29sin%28y%29e%5E%7B-2.0%5Cnu%20t%7D%20%5Cend%7Balign*%7D"><br>
</p>

where ``u``,``v``,``\omega`` are respectively the ``x-velocity``,``y-velocity`` and ``vorticity``. In this simulation, we set both ``x`` and ``y`` in range ``[0,2*pi]``. And they are discretized into a **PERIODIC** mesh with shape of ``64*64`` (change it as you wish).

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

- If we define the particle relaxation time as ![tau_p](http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Ctau_%7Bp%7D%3D%5Crho_%7Bp%7Dd_%7Bp%7D%5E2/18%5Cnu), we can get:

<p align="center">
<img align="center" height="30" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cfrac%7Bd%5Cmathbf%7BV%7D_p%7D%7Bdt%7D%3D%5Cfrac%7B%5Cbeta%7D%7B%5Ctau_p%7D%28%5Cmathbf%7BU%7D_f-%5Cmathbf%7BV%7D_p%29"><br>
</p>

- Finally, the particle translation motion (location) is governed by :
<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cfrac%7Bd%5Cmathbf%7BX%7D_p%7D%7Bdt%7D%3D%5Cmathbf%7BV%7D_p"><br>
</p>

- In practice, the particle velocity and location is updated by:

<p align="center">
<img align="center" src = "http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cdpi%7B120%7D%20%5Cfn_cm%20%5Cbegin%7Balign*%7D%20%5Cmathbf%7BV%7D_p%5E%7Bnew%7D%26%3D%5Cmathbf%7BV%7D_p%5E%7Bold%7D+%5CDelta%20t%20%5Cfrac%7B%5Cbeta%7D%7B%5Ctau_p%7D%28%5Cmathbf%7BU%7D_f%5E%7Bold%7D-%5Cmathbf%7BV%7D_p%5E%7Bold%7D%29%20%5C%5C%20%5Cmathbf%7BX%7D_p%5E%7Bnew%7D%26%3D%5Cmathbf%7BX%7D_p%5E%7Bold%7D+0.5%5CDelta%20t%28%5Cmathbf%7BV%7D_p%5E%7Bnew%7D+%5Cmathbf%7BV%7D_p%5E%7Bold%7D%29%20%5Cend%7Balign*%7D"><br>
</p>

# Results
The particle relaxation time (or particle inertia) can be changed by adjusting either ``dp`` or ``rho_p``. Here we give five typical results.

- If ``tau_p = 5.56``, it means particle has small inertia, so it will exactly follow the motion of fluid.
<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_5.56.gif" width="700"  alt="particles in taylor green vortex"/><br>
</p>

- If we increase ``tau_p`` to ``555.56``, the particle will have similar inertia with the fluid. There will be a strong competition between the particle and fluid vortex. It turns out that the particle's power is much smaller than the vortex core, so they choose to stay in the edge of the vortexs. 
<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_555.56.gif" width="700"  alt="particles in taylor green vortex"/><br>
</p>

- Further increase ``tau_p`` to ``55555.56``, the particles now can cross the outer part of the vortex, but they can still not penetrate into the vortex core. 
<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_55555.56.gif" width="700"  alt="particles in taylor green vortex"/><br>
</p>

- Finally, ``tau_p`` is increased to ``277777.78``, now they have much higher inertia than the vortex, so they just freely move in the whole domain. Particles win!
<p align="center">
    <img src="https://github.com/hietwll/TaylorGreenParticles/raw/master/gif/tau_p_277777.78.gif" width="700"  alt="particles in taylor green vortex"/><br>
</p>

# Contact

Send me an email at zjuwangzhuo@zju.edu.cn (Zhuo Wang) or 11827039@zju.edu.cn (Guo Chen) if you have any suggestion.



