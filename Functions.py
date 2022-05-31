from numpy import pi, cos, sin, exp, log, sqrt, arctan
mt = 172.76
vev = 246.22
mhsm = 125.10
mw = 80.38
mz = 91.19
sm = 540.00
mh1, mhpm, tb = tuple(map(float,input("Double value for pseudoscalar, charged Higgs, tan(β):").split()))
mh2 = (sm**4 - mh1**4 - 2*mhpm**4)**(1/4)
print("Mass of scalar at tree level: ", int(mh2), "GeV")
lambda3init = 1/vev**2*(2*mhpm**2 - mh2**2)


steps = 20000

scanrange = .3
bias = 2
lambda3list = [lambda3init + bias*i/(steps) for i in range(steps+1)]

delt = 3*mt**4/(pi**2*vev**2)

def MixNeutral(lam3):
    l5ncpv = -mh1**2/vev**2
    l4 = -(2*mhpm**2 - mh1**2)/vev**2
    mh2 = sqrt(-vev**2*(lam3 + l4 + l5ncpv))
    A = 1/(64*pi**2*vev**4)*(mh2**4*(-3/2 + 2*log(mh2/vev)) + mh1**4*(-3/2 + 2*log(mh1/vev)) + 2*mhpm**4*(-3/2 + 2*log(mhpm/vev)) + \
     6*mw**4*(-5/6 + 2*log(mw/vev)) + 3*mz**4*(-5/6+2*log(mz/vev)) - 12*mt**4*(-1 + 2*log(mt/vev)))

    B = 1/(64*pi**2*vev**4)*(mh1**4 + mh2**4 + 2*mhpm**4 + 6*mw**4 + 3*mz**4 - 12*mt**4)

    lval = vev*exp(A/(2*B) + 1/4)
    lgw = lval

    mh2prime = -vev**2*(lam3 + l4 + l5ncpv) + 1/4*delt*(1 - 2*tb**2 - 4*log(mt/lgw))
    mH20 = -vev**2*(lam3 + l4 + l5ncpv)

    mixang = 1/2*arctan(delt*tb/(mh2prime - mhsm**2))
    mh10r = mhsm**2
    pertm = delt*tb

    mh1sq = cos(mixang)**2*mh10r + sin(mixang)**2*mh2prime - sin(2*mixang)*pertm/2 
    mh2sq = sin(mixang)**2*mh10r + cos(mixang)**2*mh2prime + sin(2*mixang)*pertm/2

    return sqrt(mh1sq), sqrt(mh2sq), mixang, sqrt(mH20) 




def mixcouple(List1, List2, Angle):
    Mixh1 = List1 
    Mixh2 = List2
    Couple1 = []
    Couple2 = [] 
    for j in range(min(len(Mixh1), len(Mixh2))):
        ang = Angle[j]
        q1, q2 = abs(cos(ang)**2), abs(sin(ang)**2)
        if Mixh2[j]>Mixh1[j]:
            Couple2.append(q2)
            Couple1.append(q1)
        else:
            Couple2.append(q1)
            Couple1.append(q2)
    return Couple1, Couple2


def masshiar(List1, List2):
    Mixh1 = List1 
    Mixh2 = List2
    Couple1 = []
    Couple2 = [] 
    for j in range(min(len(Mixh1), len(Mixh2))):
        a,b = Mixh2[j],Mixh1[j]
        Couple2.append(max(a,b))
        Couple1.append(min(a,b))
    return Couple1, Couple2

def call(arg):

    Testi = [MixNeutral(k) for k in lambda3list]
    N1p = []
    N2p = []
    N3p = []
    N4p = []
    Lnep = [N1p, N2p, N3p, N4p]
    for i in range(len(Testi)):
        q=0
        for j in Lnep:
            j.append(Testi[i][q])
            q+=1
    if arg == "mass":
        N1o, N2o = masshiar(N1p, N2p)
        return N1o, N2o, N4p
    if arg == "hVV":
        ang1p, ang2p = mixcouple(N1p, N2p, N3p)
        return ang1p, ang2p 
