# %%
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


# %%
solvemethod='Radau'
problem=2

plt.rcParams['font.sans-serif'] = ['STSong']
plt.rcParams['axes.unicode_minus'] = False

# %%
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


# %%

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


# %%
xres1=sol1.y
xres2=sol2.y
tres=sol1.t


# %%
data=pd.DataFrame(np.transpose(xres1))
writer=pd.ExcelWriter('xres1-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()

# %%
data=pd.DataFrame(np.transpose(xres2))
writer=pd.ExcelWriter('xres2-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()

# %%
data=pd.DataFrame(np.transpose(tres))
writer=pd.ExcelWriter('t-temp.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()


# %%
x1vres=xres1[0]
x1fres=xres1[1]
v1vres=xres1[2]
v1fres=xres1[3]

# %%
x2vres=xres2[0]
x2fres=xres2[1]
v2vres=xres2[2]
v2fres=xres2[3]

# %%
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



# %%
plt.figure(figsize=(15,5))
plt.suptitle(r'$\nu\in[0,100000],n\in[0,1]$时最大功率下浮子与振子运动时间图')
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



# %%
def Dint(arr,h):
    #复合辛普森积分
    n=len(arr)
    sum=0
    if (n%2==0):
        sum+=(arr[-1]+arr[-2])*h/2
        n-=1
    tsum=arr[0]+arr[n-1]
    for i in range(1,n-1):
        if(i%2):
            tsum+=4*arr[i]
        else:
            tsum+=2*arr[i]
    sum+=tsum*h/3
    return sum

# %%

def draw1():
    plt.subplot(2,1,1)

    xres=sol1.y
    xv1res=xres[0]
    dxv1res=xres[2]
    xf1res=xres[1]
    dxf1res=xres[3]
    tres=sol1.t

    p1=nu1*(dxv1res-dxf1res)**2
    p1total=np.sum(p1)*titv

    print('1平均功率为',p1total/tmax)


    plt.plot(tres,xv1res)
    plt.plot(tres,xf1res)

    # plt.plot(tres,xv1res-xf1res)

    # plt.plot(tres,p1)

    plt.subplot(2,1,2)

    xres=sol2.y
    xv2res=xres[0]
    dxv2res=xres[2]
    xf2res=xres[1]
    dxf2res=xres[3]
    tres=sol2.t

    p2=nu1*(dxv2res-dxf2res)**2
    p2total=np.sum(p2)*titv

    print('2平均为',p2total/tmax)

    # plt.plot(tres,xv1res-xv2res)
    # plt.plot(tres,xf1res-xf2res)
    # plt.plot(tres,xv2res-xf2res)
    # plt.plot(tres,p1)

draw1()

# %%

def errequ1(value):
    value=np.abs(value)
    sol=intode.solve_ivp(lambda t,x:equ1(t,x,value),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
    xres=sol.y
    dxv1res=xres[2]
    dxf1res=xres[3]
    p1=value*(dxv1res-dxf1res)**2
    p1average=Dint(p1,titv)/tmax
    print(value,'->',p1average)
    return -p1average


# %%
def errequ2(nu,n):
    sol=intode.solve_ivp(lambda t,x:equ2(t,x,nu,n),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
    xres=sol.y
    dxv1res=xres[2]
    dxf1res=xres[3]
    p1=nu*(dxv1res-dxf1res)**2*np.abs(dxv1res-dxf1res)**n
    p1average=Dint(p1,titv)/tmax
    print(nu,'&&',n,'->',p1average)
    return -p1average
    

# %%



aa=np.arange(start=0,stop=100000,step=1000)
bb=[-errequ1(i) for i in aa]



# %%
aa=np.linspace(start=0,stop=100000,num=200)
cc=np.linspace(start=0,stop=1,num=100)
X,Y=np.meshgrid(aa,cc)
dd=[[-errequ2(i,j) for j in cc] for i in aa]


# %%
Z=np.array(dd)
Z=np.transpose(Z)

# %%


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

# %%

plt.rcParams['font.sans-serif'] = ['STSong']
plt.rcParams['axes.unicode_minus'] = False
plt.plot(aa,bb)
plt.xlabel(r' 阻尼系数$\nu(N/m)$')
plt.ylabel(r' 平均功率$P(N/m^2)$')
plt.title('各个阻尼系数下的平均功率')


# %%
res=opt.minimize(errequ1,[10000],method='Nelder-Mead',bounds=[(0,100000)])

print(res)
print(-errequ1(res.x[0]))

# %%
nu1=res.x[0]

# %%
errequ2(49875,0),errequ1(49875.06)

# %%
res2=opt.minimize(lambda x:errequ2(x[0],x[1]),[40000,0.2],method='Nelder-Mead',bounds=[(0,100000),(0,1)],tol=0.0001)

print(res2)
print(-errequ2(res2.x[0],res2.x[1]))

# %%
res2=opt.minimize(lambda x:errequ2(x[0],x[1]),[40000,0.6],method='Nelder-Mead',bounds=[(0,100000),(0,1)],tol=0.0001)

print(res2)
print(-errequ2(res2.x[0],res2.x[1]))

# %%
res2=opt.differential_evolution(lambda x:errequ2(np.abs(x[0]),x[1]),x0=[40000,0.2],bounds=[(0,100000),(0,1)],tol=0.001)

# %%
nu1=res2.x[0]
n1=res2.x[1]


