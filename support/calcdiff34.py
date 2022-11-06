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
sol=intode.solve_ivp(lambda t,x:equ(t,x,nu1,num1),(te[0],te[-1]),x0,method=solvemethod,t_eval=te)
xres=sol.y
tres=sol.t

#输出结果
data=pd.DataFrame(np.transpose(xres))
writer=pd.ExcelWriter('3-xres.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()
data=pd.DataFrame(np.transpose(tres))
writer=pd.ExcelWriter('3-t.xlsx')
data.to_excel(writer,'page_1',float_format='%.5f')
writer.save()
writer.close()
xvres=xres[0]
xfres=xres[1]
tvres=xres[2]
tfres=xres[3]
vvres=xres[4]
vfres=xres[5]
wvres=xres[6]
wfres=xres[7]
plt.rcParams['font.sans-serif'] = ['STSong']
plt.rcParams['axes.unicode_minus'] = False
fig = plt.figure(figsize=(15, 8))
plt.suptitle(r'$\nu\in[0,100000],\nu_m\in[0,100000]$中最优解时浮子与振子稳定运动时间图')
plt.subplots_adjust(wspace=0.25,hspace=0.25)
ax1 = fig.add_subplot(2,2,1)
ax1.set_ylim(-0.4,0.4)
ha11,=ax1.plot(tres,xvres,color='b')
ax1.set_ylabel('位移x(m)')
ax1.set_title("振子位移-时间图")  
ax1.set_xlabel('时间t(s)')
bx1 = ax1.twinx()
bx1.set_ylim(-0.01,0.01)
hb11,=bx1.plot(tres,wvres,color='g')
bx1.set_ylabel(r'角位移$\theta(rad)$')
ax1.legend([ha11,hb11,],['振子位移','振子角位移'])
ax2 = fig.add_subplot(2,2,3)
ax2.set_ylim(-2,2)
ha2,=ax2.plot(tres,xfres,color='b')
ax2.set_ylabel('位移x(m)')
ax2.set_title("浮子位移-时间图")  
ax2.set_xlabel('时间t(s)')
bx2 = ax2.twinx()
bx2.set_ylim(-0.1,0.1)
hb2,=bx2.plot(tres,wfres,color='g')
bx2.set_ylabel(r'角位移$\theta(rad)$')
ax2.legend([ha2,hb2],['浮子位移','浮子角位移'])
ax3 = fig.add_subplot(2,2,2)
ax3.set_ylim(-0.5,0.5)
ha3,=ax3.plot(tres,vvres,color='b')
ax3.set_ylabel('速度v(m/s)')
ax3.set_title("振子速度-时间图")  
ax3.set_xlabel('时间t(s)')
bx3 = ax3.twinx()
bx3.set_ylim(-0.005,0.005)
hb3,=bx3.plot(tres,wvres,color='g')
bx3.set_ylabel(r'角速度$\omega(rad/s)$')
ax3.legend([ha3,hb3],['振子速度','振子角速度'])
ax4 = fig.add_subplot(2,2,4)
ax4.set_ylim(-2,2)
ha4,=ax4.plot(tres,vfres,color='b')
ax4.set_ylabel('速度v(m/s)')
ax4.set_title("浮子速度-时间图")  
ax4.set_xlabel('时间t(s)')
bx4 = ax4.twinx()
bx4.set_ylim(-0.1,0.1)
hb4,=bx4.plot(tres,wfres,color='g')
bx4.set_ylabel(r'角速度$\omega(rad/s)$')
ax4.legend([ha4,hb4],['浮子速度','浮子角速度'])

plt.show()