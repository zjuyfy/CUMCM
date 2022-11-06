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

#参数
solvemethod='Radau'
problem=4
k=80000
km=250000
r=1025*9.8*np.pi
rm=8890.7
mv=2433
mf=4866
if(problem==3):
    ma=1028.876
    Ia=7001.914
    mu=683.4558
    muM=654.3383
    w=1.7152
    f=3640
    L=1690
elif(problem==4):
    ma=1091.099
    Ia=7142.493
    mu=528.5018
    muM=1655.909
    w=1.9806
    f=1760
    L=2140
elif(problem==5):
    ma=1091.099
    Ia=7142.493
    mu=528.5018
    muM=1655.909
    w=1.9806
    f=1760
    L=0
prd=2*np.pi/w
tmax=prd*40
titv=0.2
nu1=10000
num1=1000
g=9.8
l=0.5

#计算转动惯量
R=1
Hcone=0.8
Hclnd=3
M=mf
LL=np.sqrt(R**2+Hcone**2)
temp=np.pi*R*LL+2*np.pi*R*Hclnd
mcone=np.pi*R*LL*M/temp
mclnd=2*np.pi*R*Hclnd*M/temp
If=mclnd*R*R/2+mclnd*Hclnd**2/3+mcone*R*R/4+mcone*Hcone*Hcone/6
rr,h=0.5,0.5
Iv0=mv*rr**2/4+mv*h**2/12

#微分方程
def equ(t,x,nu,nuM):
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

#求解
te=np.arange(start=0,stop=tmax,step=titv)
try:
    x0=np.array([xres[i][-1] for i in range(8)])
except:
    x0=np.array([l-mv*g/k,0,0,0,0,0,0,0])
print(x0)

#计算功率
def calcP(nu,nuM,res):
    vv=res[4]
    wv=res[6]
    
    pv=nu*vv**2
    pw=nuM*wv**2
    p=pv+pw
    return p

#误差函数
def errequ(nu,num):
    sol=intode.solve_ivp(lambda t,x:equ(t,x,nu,num),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
    xres=sol.y
    p=calcP(nu,num,xres)
    paverage=Dint(p,titv)/tmax
    print(nu,'&&',num,'->',paverage)
    return -paverage

#做初步图像
aa=np.linspace(start=0,stop=100000,num=200)
cc=np.linspace(start=0,stop=100000,num=100)
X,Y=np.meshgrid(aa,cc)
dd=[[-errequ(i,j) for j in cc] for i in aa]
Z=np.array(dd)
Z=np.transpose(Z)
fig=plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,linewidth=0, antialiased=False)
plt.title('阻尼系数-转动阻尼系数-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'转动阻尼系数$\mu_m$')
ax.set_zlabel(r'功率$p(N\cdot m/s)$')
cs=plt.contourf(X,Y,Z,50,extend='max')
plt.title('阻尼系数-转动阻尼系数-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'转动阻尼系数$\nu_m(N\cdot m)$')
plt.colorbar(cs,label=r'功率$p(N\cdot m/s)$')
plt.contour(X,Y,Z>731.17,1,colors=['None','Black'])

#优化
nunum0=[10000,50000]
ptstore=np.array(nunum0)
def drawpt(x):
    global ptstore
    plt.scatter(x[0],x[1],color='b')
    plt.plot([x[0],ptstore[0]],[x[1],ptstore[1]],color='b')
    ptstore=x
cs=plt.contourf(X,Y,Z,50,extend='max')
plt.title('阻尼系数-转动阻尼系数-功率')
plt.xlabel(r'阻尼系数$\nu(N\cdot s/m)$')
plt.ylabel(r'转动阻尼系数$\nu_m(N\cdot m)$')
plt.colorbar(cs,label=r'功率$p(N\cdot m/s)$')
res=opt.minimize(lambda x:errequ(x[0],x[1]),nunum0,method='Nelder-Mead',bounds=[(0,100000),(0,100000)],tol=0.0001,callback=drawpt)
print(res)
print(-errequ(res.x[0],res.x[1]))