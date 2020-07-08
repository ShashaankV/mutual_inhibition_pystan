#module
#mostly analytic solutions using wilson approximation that a(T)=uinf
#
import numpy as np

def thetaids():
    thetaids=['etag1','etag2','J1','J2','etaS1','etaS2','S01','S02','alpha1','alpha2','beta1','beta2','gamma1','gamma2','taua1','taua2']
    return thetaids


def barrier1(si):
    #calc u1 escape distance
    S1=DF["J1"]*DF[si]["c1"]**DF["etaS1"]+DF["S01"]
    return S1-DF["beta2"]*DF[si]["uinf2"]

def barrier2(si):
    #calc u2 escape distance
    S2=DF["J2"]*DF[si]["c2"]**DF["etaS2"]+DF["S02"]
    return S2-DF["beta1"]*DF[si]["uinf1"]

def rivalry_error():
    err=0
    for si in range(DF['ctot']):
        err+=0.5*((DF[si]['T1hat']-DF[si]['T1mu'])**2)/(DF[si]['T1se']**2)
        err+=0.5*((DF[si]['T2hat']-DF[si]['T2mu'])**2)/(DF[si]['T2se']**2)
    return err


# def calcuquad(S,alpha,gamma):
#     a=1
#     b=-(alpha-gamma)
#     c=-S
#     up=(-b+(b**2-4*a*c)**0.5)/(2*a)
#     un=(-b-(b**2-4*a*c)**0.5)/(2*a)
#     return up*(1-np.iscomplex(up)),un*(1-np.iscomplex(un))
#
#
# def calcuinf_quad(si):
#     for i in range(1,3):
#         i=str(i)
#         S=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         upos,uneg=calcuquad(S,DF['alpha'+i],DF['gamma'+i])
#         DF[si]['uinf'+i]=max([upos,uneg,0])
#     return

# def uinferr(u,S,a,g,etag):
#     x=(S+(a-g)*u)**etag
#     err=(u-(x*(x>0)))**2
#     return err
#
# uinfjac=grad(uinferr,0)
# uinfhess=grad(grad(uinferr,0),0)

def calcuinf_newton(si):
    u1=1e-8
    u2=1e-8
    a=1e-1#newton learning rate
    S1=DF['J1']*DF[si]['c1']**DF['etaS1']+DF['S01']
    S2=DF['J2']*DF[si]['c2']**DF['etaS2']+DF['S02']
    for i in range(100):
        x=(S1+(DF['alpha1']-DF['gamma1'])*u1)
        if x>0:
            f1=u1-((x>0)*(x))**DF['etag1']
            f1prime=1-DF['etag1']*((x>0)*x)**(DF['etag1']-1)*(DF['alpha1']-DF['gamma1'])
            u1+=-a*f1/f1prime
        x=(S2+(DF['alpha2']-DF['gamma2'])*u2)
        if x>0:
            f2=u2-((x>0)*(x))**DF['etag2']
            f2prime=1-DF['etag2']*((x>0)*x)**(DF['etag2']-1)*(DF['alpha2']-DF['gamma2'])
            u2+=-a*f2/f2prime
    DF[si]['uinf1']=u1
    DF[si]['uinf2']=u2
    return


#state:= whether the pool can ever be on (in the absence of the other pool),
#however, the deterministic system could be stuck in one parity if gamma is too small

def checkwta(si):
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    ind=0
    i=i_[ind]
    j=j_[ind]
    Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
    wta1=int((Si-DF['beta'+j]*DF[si]['uinf'+j])<=0)
    ind=1
    i=i_[ind]
    j=j_[ind]
    Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
    wta2=int((Si-DF['beta'+j]*DF[si]['uinf'+j])<=0)
    return wta1,wta2
#
def calcdudtheta(si):
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    k_=["1","2"]
    ind=-1
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        for k in k_:
            S=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
            x=S+DF['alpha'+i]*DF[si]['uinf'+i]-DF['gamma'+i]*DF[si]['uinf'+i]
            if x<=0:
                x=1e-12
                gain=0
                dgain=0
            else:
                gain=x**DF['etag'+i]
                dgain=DF['etag'+i]*x**(DF['etag'+i]-1)
            #do not overwrite x
            DF[si]['dudalpha_'+i+k]=(i==k)*(dgain*DF[si]['uinf'+i])/(1-dgain*(DF['alpha'+i]-DF['gamma'+i]))
            DF[si]['dudgamma_'+i+k]=-DF[si]['dudalpha_'+i+k]
            DF[si]['dudJ_'+i+k]=(i==k)*dgain*DF[si]['c'+i]**DF['etaS'+i]/(1-dgain*(DF['alpha'+i]-DF['gamma'+i]))
            DF[si]['dudS0_'+i+k]=(i==k)*dgain/(1-dgain*(DF['alpha'+i]-DF['gamma'+i]))
            DF[si]['dudetaS_'+i+k]=(i==k)*dgain*DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]*np.log(DF[si]['c'+i])/(1-dgain*(DF['alpha'+i]-DF['gamma'+i]))
            DF[si]['dudetag_'+i+k]=(i==k)*np.log(x)*gain/(1-dgain*(DF['alpha'+i]-DF['gamma'+i]))
    return




####barriers...give barrier w.t. to barrier index (indb)
#
def bgrad_alphak(si,k,indb):
    #k is the pool-parameter index for partial
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    grad_alphak=-DF['beta'+j]*DF[si]['dudalpha_'+j+k]
    return grad_alphak
#
def bgrad_gammak(si,k,indb):
    #k is the pool-parameter index for partial
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    grad_gammak=-DF['beta'+j]*DF[si]['dudgamma_'+j+k]
    return grad_gammak
#
def bgrad_betak(si,k,indb):
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    grad_betak=-(j==k)*DF[si]['uinf'+j]
    return grad_betak
#
def bgrad_tauak(si,k,indb):
    return 0
#
def bgrad_S0k(si,k,indb):
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    return (i==k)-DF['beta'+j]*DF[si]['dudS0_'+j+k]

def bgrad_Jk(si,k,indb):
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    return (i==k)*DF[si]['c'+i]**DF['etaS'+i]-DF['beta'+j]*DF[si]['dudJ_'+j+k]

def bgrad_etaSk(si,k,indb):
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    return (i==k)*DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]*np.log(DF[si]['c'+i])-DF['beta'+j]*DF[si]['dudetaS_'+j+k]

def bgrad_etagk(si,k,indb):
    #k is the pool-parameter index for partial
    k=str(k)
    i_=["1","2"]
    j_=["2","1"]
    i=i_[indb]
    j=j_[indb]
    grad_etagk=-DF['beta'+j]*DF[si]['dudetag_'+j+k]
    return grad_etagk

####rivalry


def calcT(si):
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        x=(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j])
        DF[si]['T'+i+'hat']=-DF['taua'+j]*np.log(x)
    return

def rivgrad_alphak(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_alphak=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b=DF[si]['dudalpha_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudalpha_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_alphak+=(a*(b+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_alphak


def rivgrad_gammak(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_gammak=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b0=(j==k)*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]**2*DF[si]['uinf'+j])
        b1=DF[si]['dudgamma_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudgamma_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_gammak+=(a*(b0+b1+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_gammak

def rivgrad_betak(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_betak=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=(i==k)*DF['taua'+j]*DF[si]['uinf'+i]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        rivgrad_betak+=a*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_betak

def rivgrad_tauak(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_tauak=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF[si]['dTdtaua_'+i+k]=(j==k)*-np.log((Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]))
        rivgrad_tauak+=a*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_tauak


def rivgrad_Jk(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_Jk=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b0=-(DF[si]['c'+k]**DF['etaS'+k]*(j==k))/(DF['gamma'+j]*DF[si]['uinf'+j])
        b1=DF[si]['dudJ_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudJ_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_Jk+=(a*(b0+b1+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_Jk


def rivgrad_S0k(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_S0k=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b0=-(j==k)/(DF['gamma'+j]*DF[si]['uinf'+j])
        b1=DF[si]['dudS0_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudS0_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_S0k+=(a*(b0+b1+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_S0k

def rivgrad_etaSk(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_etaSk=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b0=-DF['J'+k]*DF[si]['c'+k]**DF['etaS'+k]*np.log(DF[si]['c'+k])*(j==k)/(DF['gamma'+j]*DF[si]['uinf'+j])
        b1=DF[si]['dudetaS_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudetaS_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_etaSk+=(a*(b0+b1+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_etaSk

def rivgrad_etagk(si,k):
    k=str(k)
    #i:= self, j:=other, k:=parameter index
    i_=["1","2"]
    j_=["2","1"]
    ind=-1
    rivgrad_etagk=0
    for ind in range(2):
        i=i_[ind]
        j=j_[ind]
        Sj=DF['J'+j]*DF[si]['c'+j]**DF['etaS'+j]+DF['S0'+j]
        a=DF['taua'+j]*DF['gamma'+j]*DF[si]['uinf'+j]/(Sj-DF['beta'+i]*DF[si]['uinf'+i])
        b=DF[si]['dudetag_'+j+k]*(Sj-DF['beta'+i]*DF[si]['uinf'+i])/(DF['gamma'+j]*DF[si]['uinf'+j]**2)
        c=DF[si]['dudetag_'+i+k]*DF['beta'+i]/(DF['gamma'+j]*DF[si]['uinf'+j])
        rivgrad_etagk+=(a*(b+c))*(1.0/DF[si]['T'+i+'se']**2*(DF[si]['T'+i+'hat']-DF[si]['T'+i+'mu']))
    return rivgrad_etagk


# ####wta..NEED TO FIX RETURN
# def wtagrad_alpha1(si):
#     #k is the pool-parameter index for partial
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdalphak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdalphak+=x*-DF['beta'+j]*DF[si]['dudalpha_'+j+k]
#     return wtagrad_alphak
#
# def wtagrad_alpha2(si):
#     #k is the pool-parameter index for partial
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdalphak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdalphak+=x*-DF['beta'+j]*DF[si]['dudalpha_'+j+k]
#     return wtagrad_alphak
#
# def wtagrad_gamma1(si):
#     #k is the pool-parameter index for partial
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdgammak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdgammak+=x*-DF['beta'+j]*DF[si]['dudgamma_'+j+k]
#     return wtagrad_gammak

# def wtagrad_gamma2(si):
#     #k is the pool-parameter index for partial
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdgammak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdgammak+=x*-DF['beta'+j]*DF[si]['dudgamma_'+j+k]
#     return wtagrad_gammak
#
# def wtagrad_beta1(si):
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdbetak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdbetak+=x*-(j==k)*DF[si]['uinf'+j]
#     return wtagrad_betak
#
# def wtagrad_beta2(si):
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdbetak=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdbetak+=x*-(j==k)*DF[si]['uinf'+j]
#     return wtagrad_betak
#
# def wtagrad_taua1(si):
#     return 0
#
# def wtagrad_taua2(si):
#     return 0
#
# def wtagrad_J1(si):
#     #k is the pool-parameter index for partial
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdJk=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdJk+=x*((i==k)*DF[si]['c'+i]**DF['etaS'+i]-DF['beta'+j]*DF[si]['dudJ_'+j+k])
#     return wtagrad_Jk
#
# def wtagrad_J2(si):
#     #k is the pool-parameter index for partial
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdJk=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdJk+=x*((i==k)*DF[si]['c'+i]**DF['etaS'+i]-DF['beta'+j]*DF[si]['dudJ_'+j+k])
#     return wtagrad_Jk
#
#
# def wtagrad_S01(si):
#     #k is the pool-parameter index for partial
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdS0k=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdS0k+=x*((i==k)-DF['beta'+j]*DF[si]['dudS0_'+j+k])
#     return wtagrad_S0k
#
# def wtagrad_S02(si):
#     #k is the pool-parameter index for partial
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdS0k=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=-int(x<0)
#         derrdS0k+=x*((i==k)-DF['beta'+j]*DF[si]['dudS0_'+j+k])
#     return wtagrad_S0k

# def calcderrdetaS1(si):
#     #k is the pool-parameter index for partial
#     k="1"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdetaSk=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=x*(x<0)
#         derrdetaSk+=2*x*((i==k)*DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]*np.log(DF[si]['c'+i])-DF['beta'+j]*DF[si]['dudJ_'+j+k])
#     return derrdetaSk
#
# def calcderrdetaS2(si):
#     #k is the pool-parameter index for partial
#     k="2"
#     i_=["1","2"]
#     j_=["2","1"]
#     ind=-1
#     derrdetaSk=0 #add derivatives for pool1 and pool2 wta'ness
#     for ind in range(2):
#         i=i_[ind]
#         j=j_[ind]
#         Si=DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]+DF['S0'+i]
#         x=Si-DF['beta'+j]*DF[si]['uinf'+j]
#         x=x*(x<0)
#         derrdetaSk+=2*x*((i==k)*DF['J'+i]*DF[si]['c'+i]**DF['etaS'+i]*np.log(DF[si]['c'+i])-DF['beta'+j]*DF[si]['dudJ_'+j+k])
#     return derrdetaSk


# def calcuinf_human_newtons(S1S2):
#     si=S1S2[-1]
#     S1=float(S1S2[0])
#     S2=float(S1S2[1])
#     u1inf=0
#     u1inf0=u1inf
#     S=DF['JS1']*S1+DF['Soff1']
#     x1=S+(DF['alpha1']-DF['gamma1'])*u1inf
#     u2inf=0
#     u2inf0=u2inf
#     S=DF['JS2']*S2+DF['Soff2']
#     x2=S+(DF['alpha2']-DF['gamma2'])*u2inf
#     for i in range(10):
#         while x1<0:
#             u1inf=u1inf0+(u1inf-u1inf0)/2.
#             S=DF['JS1']*S1+DF['Soff1']
#             x1=S+(DF['alpha1']-DF['gamma1'])*u1inf
#         f1=u1inf-x1**.5
#         fprime1=1-.5*x1**(-.5)*(DF['alpha1']-DF['gamma1'])
#         while x2<0:
#             u2inf=u2inf0+(u2inf-u2inf0)/2.
#             S=DF['JS2']*S2+DF['Soff2']
#             x2=S+(DF['alpha2']-DF['gamma2'])*u2inf
#         f2=u2inf-x2**.5
#         fprime2=1-.5*x2**(-.5)*(DF['alpha2']-DF['gamma2'])
#         u1inf0=u1inf
#         u1inf+=-f1/fprime1
#         S=DF['JS1']*S1+DF['Soff1']
#         x1=S+(DF['alpha1']-DF['gamma1'])*u1inf
#         # f1=u1inf-x**.5
#         # fprime1=1-.5*x**(-.5)*(DF['alpha1']-DF['gamma1'])
#
#         u2inf0=u2inf
#         u2inf+=-f2/fprime2
#         S=DF['JS2']*S2+DF['Soff2']
#         x2=S+(DF['alpha2']-DF['gamma2'])*u2inf
#         # f2=u2inf-x**.5
#         # fprime2=1-.5*x**(-.5)*(DF['alpha2']-DF['gamma2'])
#     uinf1=absolute(u1inf)*(u1inf>=0)
#     uinf2=absolute(u2inf)*(u2inf>=0)
#     DF[si]['uinf1']=uinf1
#     DF[si]['uinf2']=uinf2
#
#     #check stability, call it rhythmif greater than one than unstable
#     # Jac1=array([[DF['alpha1']*fprime1-1,-DF['gamma1']*fprime1],[1,-1]])
#     # Jac1[1,:]*=1/DF['taua1']
#     gprime1 = .5*x1**(-.5)
#     DF[si]['rhythm1'] = (DF['alpha1'] * gprime1) #if greater than one than unstable
#
#     # Jac2=array([[DF['alpha2']*fprime2-1,-DF['gamma1']*fprime1],[1,-1]])
#     # Jac1[1,:]*=1/DF['taua1']
#     gprime2 = .5*x2**(-.5)
#     DF[si]['rhythm2'] = (DF['alpha2'] * gprime2)
#
#     return
