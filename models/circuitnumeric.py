#module
#numerical solution and misc for asymmetric circuit

from numpy import *

tauu=1.
dt=tauu/10.
ttot=60 #seconds
itot=int(ttot*1000./dt)

def ccc(si):
    S1=P['J1']*P[si]['c1']**P['etaS1']+P['S01']
    S2=P['J2']*P[si]['c2']**P['etaS2']+P['S02']
    u1_=[]
    u2_=[]
    a1_=[]
    a2_=[]
    u1=random.uniform()
    u2=random.uniform()
    a1=0.001*random.uniform()
    a2=0.001*random.uniform()
    for i in range(itot):
        x=(S1+P['alpha1']*u1-P['beta2']*u2-P['gamma1']*a1)
        x=((x>0)*x)**P['etag1']
        u1+=dt*(-u1+x)
        x=(S2+P['alpha2']*u2-P['beta1']*u1-P['gamma2']*a2)
        x=((x>0)*x)**P['etag2']
        u2+=dt*(-u2+x)
        u1_.append(u1)
        u2_.append(u2)
        a1+=(-a1+u1)*dt/(1000.*P['taua1'])
        a2+=(-a2+u2)*dt/(1000.*P['taua2'])
        a1_.append(a1)
        a2_.append(a2)
    u1_=array(u1_)
    u2_=array(u2_)
    a1_=array(a1_)
    a2_=array(a2_)
    urec=zeros((2,itot))
    urec[0,:]=u1_
    urec[1,:]=u2_
    arec=zeros((2,itot))
    arec[0,:]=a1_
    arec[1,:]=a2_
    # tdA,tdB=calcTNpool(urec,dt,2)
    return urec,arec

def calcTNpool(urec,dt=dt,Npools=2): #activity of any pool over the sum > three
    xA=urec[0,:]/sum(urec,0)
    n=len(xA)
    xB=urec[-1,:]/sum(urec,0)
    if Npools>2:
        xO=urec[1:-1,:]/sum(urec,0)
    else:
        xO=zeros((1,n))
    s_=-100*ones(n) #percept state-time vector
    for i in range(n):
        if xA[i]>.5:
            s_[i]=1
        elif xB[i]>.5:
            s_[i]=-1
        elif max(xO[:,i])>.5:
            s_[i]=0
        else:
            s_[i]=-100
    T={}
    T[1]=[]
    T[0]=[]
    T[-1]=[]
    T[-100]=[] #not dom or mixed
    ind=0
    s=s_[ind] #state
    sbuff_=zeros(3)
    sbuff_[0]=s
    rcount=0 #number of reversions
    scount=0 #number of switches
    gate=0
    first_dom=0
    if s!=0 and s!=-100 and gate==0:
        first_dom=s
        gate=1
    for i in range(1,n):
            if s!=0 and s!=-100 and gate==0:
                first_dom=s
                gate=1
            if s_[i]!=s:#if state changes then record length of last epoch
                T[s].append((i-ind)*dt)
                s=s_[i]
                ind=i
                scount+=1
                sbuff_=roll(sbuff_,1)
                sbuff_[0]=s
                if scount>1 and sbuff_[1]==0:
                    if sbuff_[0]==sbuff_[2]:
                        rcount+=1
    # rf=1.*rcount/(scount-1)
    # T[s].append((i-ind)*dt)
    if len(T[1][2:])>0:
        tdA=mean(T[1][2:])/1000.
    else:
        tdA=0
    if len(T[-1][2:])>0:
        tdB=mean(T[-1][2:])/1000.
    else:
        tdB=0
    return tdA,tdB




# def calcuinf(S1S2):
#     si=S1S2[2]
#     urec=ccc(S1S2)[-1]
#     du=urec[0,:]-urec[1,:]
#     x=where(diff(sign(du))!=0)[0]
#     # print(len(x))
#     if len(x)<2: #wta
#         uinf1=urec[0,-1]
#         uinf2=urec[1,-1]
#     else:
#         a0=urec[0,x[-1]-1000]
#         a1=urec[0,x[-2]-1000]
#         uinf1=max((a0,a1))
#         a0=urec[1,x[-1]-1000]
#         a1=urec[1,x[-2]-1000]
#         uinf2=max((a0,a1))
#     # print(uinf1)
#     # print(uinf2)
#     P[si]['uinf1']=uinf1
#     P[si]['uinf2']=uinf2
#     return
#
# def checkgap(S1S2,cmax):
#     gap=0
#     si=S1S2[2]
#     urec=ccc(S1S2)[-1]
#     du=urec[0,:]-urec[1,:]
#     theta=0.1*max(du[int(1000/dt):])
#     du=du*(absolute(du)>theta)
#     s0=sign(du[0]) #state
#     count=0
#     for i in range(1,len(du)):
#         s1=sign(du[i])
#         if s1!=s0:
#             if s1==0:
#                 count+=1
#             else:
#                 s0=s1 #update state
#                 count=0
#         if count>cmax:
#             gap=1
#             break
#     return gap
