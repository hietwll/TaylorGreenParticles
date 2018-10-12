#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Simulate pa
# @Date    : 2018-10-12 11:25:18
# @Author  : Zhuo Wang (zjuwangzhuo@zju.edu.cn)

import numpy as np
import os
import time
import shutil
from math import pi
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams["font.size"] = 8.0
plt.rcParams['axes.linewidth'] = "0.2"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"



class TaylorGreen():
    def __init__(self,arg1,arg2,arg3):
        self.ng = 64 #number of grids
        self.npar = 1000 #number of particles
        self.vis = 1.0/1000.0 # vicosity = 1.0/Re
        self.dp = arg1 # particle diameter
        self.rho = arg2 # particle density
        self.t = 0.0 # time
        self.h = 2.0*pi/(self.ng-1) # grid width
        self.taup = self.rho*self.dp**2/18.0/self.vis #particle relax time
        self.cfl = arg3 #cfl number
        self.dt = self.h/1.0*self.cfl # initial time step
        self.u = np.zeros((self.ng,self.ng)) # fluid velocity
        self.v = np.zeros((self.ng,self.ng))
        self.omega = np.zeros((self.ng,self.ng))
        self.up = np.random.random(self.npar) # particle velocity
        self.vp = np.random.random(self.npar)
        self.rep = np.zeros(self.npar)
        self.upf = np.zeros(self.npar)
        self.vpf = np.zeros(self.npar)
        self.up_old = np.zeros(self.npar)
        self.vp_old = np.zeros(self.npar)
        self.xp = np.random.random(self.npar)*2.0*pi # particle location
        self.yp = np.random.random(self.npar)*2.0*pi
        self.xmesh,self.ymesh = np.meshgrid(np.linspace(0,2.0*pi,self.ng),np.linspace(0,2.0*pi,self.ng))
        self.u = np.sin(self.xmesh)*np.cos(self.ymesh)*np.exp(-2.0*self.t*self.vis)
        self.v = -np.cos(self.xmesh)*np.sin(self.ymesh)*np.exp(-2.0*self.t*self.vis)
        self.omega = (self.deriv(self.v,1)-self.deriv(self.u,0))*0.5
        # self.omega = np.sin(self.xmesh)*np.sin(self.ymesh)*np.exp(-2.0*self.t*self.vis) # this is analytical solution

        # draw figure
        self.fig , axe = plt.subplots(nrows=1,ncols=2,figsize=(6.0,2.75),dpi=300)
        self.axes1,self.axes2 = axe.flat
        plt.subplots_adjust(left=0.075,hspace=0.1,wspace=0.5,bottom=0.075,right=0.85, top=0.9)
 
    def ShowInit(self,axes,var,ti):
        axes.set_aspect(1.0)
        axes.set_xlim((0.0,2.0*pi))
        axes.set_ylim((0.0,2.0*pi))
        la1 = [0.0,pi,2.0*pi]
        la2 = [0.5*pi,1.5*pi]
        la3 = np.linspace(-1.0,1.2,12,endpoint=True)
        tickla = [r'0',r'$\pi$',r'2$\pi$']
        tickla2 = list(map(lambda x: str(int(x*10)/10.0),la3))
        axes.set_xticks(la1,minor=False)
        axes.set_xticks(la2,minor=True)
        axes.set_xticklabels(tickla,minor=False)
        axes.set_yticks(la1,minor=False)
        axes.set_yticks(la2,minor=True)
        axes.set_yticklabels(tickla,minor=False)
        axes.tick_params(axis='both',direction='in',which='both',right=True,top=True,width=0.2,length=2.0)
        cf = axes.contourf(self.xmesh,self.ymesh,var,np.linspace(-1.0,1.0,50,endpoint=True),cmap='rainbow',extend='both')
        sc = axes.scatter(self.xp,self.yp,color='k',s=5,marker='.',edgecolors='none')
        p = axes.get_position().get_points().flatten()
        ax_cbar = self.fig.add_axes([p[2]+0.02,p[1],0.016,p[3]-p[1]])
        cbar = plt.colorbar(cf,cax=ax_cbar,orientation='vertical')
        cbar.ax.set_xlabel(ti)
        cbar.set_ticks(la3)
        cbar.ax.set_yticklabels(tickla2)
        axes.set_xlabel(r'X')
        axes.set_ylabel(r'Y')
        return cf,sc

    def deriv(self,vin,ax):
        #ax 0 = y, 1 = x
        vout = (np.roll(vin,-1,axis=ax)-np.roll(vin,1,axis=ax))/self.h/2.0
        return vout #一定要用返回形式才能传出

    def UpdateFluid(self): # update the fluid velocity
        self.u = np.sin(self.xmesh)*np.cos(self.ymesh)*np.exp(-2.0*self.t*self.vis)
        self.v = -np.cos(self.xmesh)*np.sin(self.ymesh)*np.exp(-2.0*self.t*self.vis)
        self.omega = (self.deriv(self.v,1)-self.deriv(self.u,0))*0.5

    def UpdatePar(self): # update partilce velocity and postion
        #fluid velocity at partilce position
        self.upf = np.sin(self.xp)*np.cos(self.yp)*np.exp(-2.0*self.t*self.vis)
        self.vpf = -np.cos(self.xp)*np.sin(self.yp)*np.exp(-2.0*self.t*self.vis)

        # particle reynolds number
        self.rep = ((self.upf-self.up)**2+(self.vpf-self.vp)**2)**0.5*self.dp/self.vis

        # drag coffecient(reuse)
        self.rep = 1.0 + 0.15*self.rep**0.687

        #old velocity
        self.up_old = self.up
        self.vp_old = self.vp

        #update velocity
        self.up = self.up_old+self.dt*18.0/self.rho/self.dp**2.0*self.rep*(self.upf-self.up_old) # note we assure fluid density is 1.0
        self.vp = self.vp_old+self.dt*18.0/self.rho/self.dp**2.0*self.rep*(self.vpf-self.vp_old)

        #update location
        self.xp = self.xp + (self.up_old+self.up)*0.5*self.dt
        self.yp = self.yp + (self.vp_old+self.vp)*0.5*self.dt

        # check boundary
        self.xp[self.xp>2.0*pi] = self.xp[self.xp>2.0*pi] - 2.0*pi
        self.xp[self.xp<0] = self.xp[self.xp<0] + 2.0*pi

        self.yp[self.yp>2.0*pi] = self.yp[self.yp>2.0*pi] - 2.0*pi
        self.yp[self.yp<0] = self.yp[self.yp<0] + 2.0*pi

    def UpdateDt(self,i,j): # Update timestep based on CFL number
        max_u = max(np.max(np.abs(self.up)),np.max(np.abs(self.vp)))
        assert max_u < 1000 ,'Divegence reached!!! Max velocity : {:.3f}'.format(max_u)
        self.dt = self.h/max_u*self.cfl
        print('max u = {:.3f}  dt = {:.3f}  i = {:d}  j = {:d}'.format(max_u,self.dt,i,j))

    def update(self,i,cf1,cf2,sc1,sc2,nsave,substep):
        for j in range(substep):
            self.UpdateDt(i,j)
            self.t = self.t + self.dt
            self.UpdatePar()
        self.UpdateFluid()
        labels = r'$\tau_p$ : {:.3f}   $\Delta t$  :  {:.3e}   t : {:.3f}    i : {:d}  '.format(self.taup,self.dt,self.t,i)
        title = self.fig.suptitle(labels)
        cf1.set_array(self.u)
        cf2.set_array(self.omega)
        sc1.set_offsets(np.c_[self.xp,self.yp])
        sc2.set_offsets(np.c_[self.xp,self.yp])
        # save fig
        if i%nsave == 0:
            dirs = 'tau_p_{:.2f}'.format(self.taup)
            plt.savefig(dirs+'/'+str(i).zfill(4)+'.png',dpi=150)
        return title


    def ShowAnimation(self,nt,nsave,substep):
        cf1,sc1=self.ShowInit(self.axes1,self.u,r'$u$')
        cf2,sc2=self.ShowInit(self.axes2,self.omega,r'$\omega$')
        # figure out put dir
        dirs = 'tau_p_{:.2f}'.format(self.taup)
        if os.path.isdir(dirs):
            shutil.rmtree(dirs)
        os.makedirs(dirs)
        ani = animation.FuncAnimation(fig=self.fig,func=self.update,fargs=(cf1,cf2,sc1,sc2,nsave,substep, ),frames=nt+1,repeat=False,interval=1,save_count=1,blit=0)
        plt.show()

        # here we do not use the default save, because we want to skip some frames
        # ani.save('test.gif')
        self.CreateGif(dirs)

    def CreateGif(self,dirs):
        ow = os.walk(dirs)
        _,_,files = next(ow)
        tms = sorted(list(map(lambda x: int(x.split('.')[0]),files)))
        with imageio.get_writer('tau_p_{:.2f}'.format(self.taup)+'.gif', fps=4) as writer:
            for i in tms:
                image = imageio.imread(dirs+"\\"+str(i).zfill(4)+'.png')
                writer.append_data(image)
        writer.close()

tg = TaylorGreen(0.1,10.0,0.01) # dp ,rho, cfl
tg.ShowAnimation(200,4,50) # total frame, freq to save frame, sub step for integration









