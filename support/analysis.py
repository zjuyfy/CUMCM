from simpson import *
from calcdiff34 import *
import scipy.integrate as intode
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time
import scipy.optimize as opt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('svg')
import pandas as pd


from optimize34 import *

def equt(t,x,nu,nuM,f,L,w,k,km):
    xv=x[0]
    xf=x[1]
    tv=x[2]
    tf=x[3]
    vv=x[4]
    vf=x[5]
    wv=x[6]
    wf=x[7]
    dxv=vv
    dxf=vf
    dtv=wv
    dtf=wf
    dwf=(-rm*tf-muM*wf+nuM*wv+km*tv+L*np.cos(w*t))/(If+Ia)
    Iv=Iv0+mv*xv**2
    dwv=(-nuM*wv-km*tv+mv*g*xv*np.sin(tv+tf))/Iv-dwf
    m=mf+ma+mv*np.sin(tv+tf)**2
    dvf=(k*(xv-l+mv*g/k)*np.cos(tv+tf)+nu*vv*np.cos(tv+tf)-mu*vf+f*np.cos(w*t)-r*xf+np.sin(tv+tf)*(mv*(xv*dwv+2*vv*wv+xv*dwf+2*vv*wf)-mv*g*np.sin(tv+tf)))/m
    dvv=(-k*(xv-l)-nu*vv-mv*g*np.cos(tv+tf)-mv*(-xv*wv**2+dvf*np.cos(tv+tf)-xv*wf**2-2*xv*wf*wv))/mv
    return np.array([dxv,dxf,dtv,dtf,dvv,dvf,dwv,dwf])

def calcPow(nu,nuM,f,L,w,k,km):
    print(x0)
    sol=intode.solve_ivp(lambda t,x:equt(t,x,nu,nuM,f,L,w,k,km),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
    xres=sol.y
    tres=sol.t
    p=calcP(nu,nuM,xres)
    index=xres[1]/np.cos(xres[3])+np.tan(xres[3])
    return Dint(p,titv)/tmax,max(index)

nu0,num0=res.x[0],res.x[1]

ff=range(0,2500,40)
mm=[0 for i in ff]
pp=mm[:]
for i in range(len(ff)):
    temp=calcPow(nu0,num0,ff[i],L,w,k,km)
    mm[i]=temp[1]
    pp[i]=temp[0]

plt.plot(ff,mm)
plt.ylim(0.69,1.2)
plt.xlim(-100,2600)
plt.plot(range(-1000,9000,1000),[1 for i in range(-1000,9000,1000)],color='r')
plt.title('激励力-浮子振幅指数')
plt.xlabel(r'激励力($N$)')
plt.ylabel('浮子振幅指数(m)')
xxx=2289.24
plt.vlines([xxx],[0],[1],linestyles='dashed',colors=['grey'])
plt.text(xxx+50,0.7,r'x='+str(xxx))

ll=range(0,7500,100)
mm=[0 for i in ll]
pp=mm[:]
for i in range(len(ll)):
    temp=calcPow(nu0,num0,f,ll[i],w,k,km)
    mm[i]=temp[1]
    pp[i]=temp[0]
    
plt.plot(ll,mm)
plt.ylim(0.88,1.02)
plt.xlim(-100,7500)
plt.plot(range(-1000,9000,1000),[1 for i in range(-1000,9000,1000)],color='r')
plt.title('激励力矩-浮子振幅指数')
plt.xlabel(r'激励力矩($N\cdot m$)')
plt.ylabel('浮子振幅指数(m)')
xxx=6932.31
plt.vlines([xxx],[0],[1],linestyles='dashed',colors=['grey'])
plt.text(xxx+50,0.883,r'x='+str(xxx))