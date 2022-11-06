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

solvemethod='Radau'
problem=2
plt.rcParams['font.sans-serif'] = ['STSong']
plt.rcParams['axes.unicode_minus'] = False

k=80000
r=1025*9.8*np.pi
mv=2433
mf=4866
if(problem==1):
    ma=1335.535
    mu=656.3616
    w=1.4005
    f=6250
elif(problem==2):
    ma=1165.992
    mu=167.8395
    w=2.2143
    f=4890
prd=2*np.pi/w
tmax=prd*40
titv=0.2
nu1=10000
n1=0.5

start=time.time()
try:
    x0=np.array([xres1[i][-1] for i in range(4)])
    x1=np.array([xres2[i][-1] for i in range(4)])
except:
    x0=np.array([0,0,0,0])
    x1=np.array([0,0,0,0])
def equ1(t,x,nu):
    xv=x[0]
    xf=x[1]
    dxv=x[2]
    dxf=x[3]
    fxv=dxv
    fxf=dxf
    fdxv=(-k*xv+k*xf-nu*(dxv-dxf))/mv
    m=mf+ma
    fdxf=(-k*xf+k*xv+nu*(dxv-dxf)-mu*dxf+f*np.cos(w*t)-r*xf)/m
    return np.array([fxv,fxf,fdxv,fdxf])
def fun1(t,x):
    return equ1(t,x,nu1)
def equ2(t,x,nu,n):
    xv=x[0]
    xf=x[1]
    dxv=x[2]
    dxf=x[3]
    fxv=dxv
    fxf=dxf
    fdxv=(-k*xv+k*xf-nu*(dxv-dxf)*(np.abs(dxv-dxf))**n)/mv
    m=mf+ma
    fdxf=(-k*xf+k*xv+nu*(dxv-dxf)*(np.abs(dxv-dxf))**n-mu*dxf+f*np.cos(w*t)-r*xf)/m
    return np.array([fxv,fxf,fdxv,fdxf])
def fun2(t,x):
    return equ2(t,x,nu1,n1)
print(x0)
print(x1)
te=np.arange(start=0,stop=tmax,step=titv)
sol1 = intode.solve_ivp(fun1,(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
sol2 = intode.solve_ivp(fun2,(te[0],te[-1]),x1,method=solvemethod,t_eval=te)
end=time.time()
print(end-start)

xres1=sol1.y
xres2=sol2.y
tres=sol1.t

data=pd.DataFrame(np.transpose(xres1))
writer=pd.ExcelWriter('xres1-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()
data=pd.DataFrame(np.transpose(xres2))
writer=pd.ExcelWriter('xres2-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()
data=pd.DataFrame(np.transpose(tres))
writer=pd.ExcelWriter('t-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()

x1vres=xres1[0]
x1fres=xres1[1]
v1vres=xres1[2]
v1fres=xres1[3]
x2vres=xres2[0]
x2fres=xres2[1]
v2vres=xres2[2]
v2fres=xres2[3]

def draw1():
    plt.figure(figsize=(15,5))
    plt.suptitle(r'$\nu\in[0,100000],n=0$时最大功率下浮子与振子运动时间图')
    plt.subplots_adjust(wspace=0.25,hspace=0.45)
    plt.subplot(2,2,1)
    plt.plot(tres,x1vres)
    plt.title('振子位置')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$x_{v}$(m)')
    plt.subplot(2,2,2)
    plt.plot(tres,v1vres)
    plt.title('振子速度')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$\dot x_{v}$(m/s)')
    plt.subplot(2,2,3)
    plt.plot(tres,x1fres)
    plt.title('浮子位置')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$x_{f}$(m)')
    plt.subplot(2,2,4)
    plt.plot(tres,v1fres)
    plt.title('浮子速度')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$\dot x_{f}$(m/s)')    
    
def draw2():
    plt.figure(figsize=(15,5))
    plt.suptitle(r'$\nu\in[0,100000],n\in[0,1]$时最大功率下浮子与振子稳定运动时间图')
    plt.subplots_adjust(wspace=0.25,hspace=0.45)
    plt.subplot(2,2,1)
    plt.plot(tres,x2vres)
    plt.title('振子位置')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$x_{v}$(m)')
    plt.subplot(2,2,2)
    plt.plot(tres,v2vres)
    plt.title('振子速度')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$\dot x_{v}$(m/s)')
    plt.subplot(2,2,3)
    plt.plot(tres,x2fres)
    plt.title('浮子位置')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$x_{f}$(m)')
    plt.subplot(2,2,4)
    plt.plot(tres,v2fres)
    plt.title('浮子速度')
    plt.xlabel(r't(s)')
    plt.ylabel(r'$\dot x_{f}$(m/s)')

draw1()
plt.show()