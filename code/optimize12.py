from simpson import *
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

def errequ2(nu,n):
    sol=intode.solve_ivp(lambda t,x:equ2(t,x,nu,n),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
    xres=sol.y
    dxv1res=xres[2]
    dxf1res=xres[3]
    p1=nu*(dxv1res-dxf1res)**2*np.abs(dxv1res-dxf1res)**n
    p1average=Dint(p1,titv)/tmax
    print(nu,'&&',n,'->',p1average)
    return -p1average
def calcP(nu,n,xres):
    dxv1res=xres[2]
    dxf1res=xres[3]
    p1=nu*(dxv1res-dxf1res)**2*np.abs(dxv1res-dxf1res)**n
    return Dint(p1,titv)/tmax

aa=np.linspace(start=0,stop=100000,num=200)
cc=np.linspace(start=0,stop=1,num=100)
X,Y=np.meshgrid(aa,cc)
dd=[[-errequ2(i,j) for j in cc] for i in aa]

Z=np.array(dd)
Z=np.transpose(Z)
fig=plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,linewidth=0, antialiased=False)
plt.title('阻尼系数-幂次-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'幂次$n$')
ax.set_zlabel(r'功率$p(N\cdot m/s)$')

cs=plt.contourf(X,Y,Z,50,extend='max')
plt.title('阻尼系数-幂次-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'幂次$n$')
plt.colorbar(cs,label=r'功率$p(N\cdot m/s)$')
plt.contour(X,Y,Z>228.2,1,colors=['None','Black'])

nun0=[40000,0.8]
ptstore=np.array(nun0)
def drawpt(x):
    global ptstore
    plt.scatter(x[0],x[1],color='b')
    plt.plot([x[0],ptstore[0]],[x[1],ptstore[1]],color='b')
    ptstore=x
cs=plt.contourf(X,Y,Z,50,extend='max')

plt.title('阻尼系数-幂次-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'幂次$n$')
plt.colorbar(cs,label=r'功率$p(N\cdot m/s)$')

res2=opt.minimize(lambda x:errequ2(x[0],x[1]),nun0,method='Nelder-Mead',bounds=[(0,100000),(0,1)],tol=0.0001,callback=drawpt)

print(res2)
print(-errequ2(res2.x[0],res2.x[1]))