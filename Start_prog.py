from sympy import *
from time import time 
init_printing()
N = 6
L = zeros(N,1)
#Lp = zeros(N-1,1)
for i in range(1,N):
    L[i-1,0] = Symbol('lambda_' + '{}'.format(i), real = True)
L = Matrix(L)
l5 = Symbol('lambda_5', real = True)
l1, l2, l3, l4 = L[0], L[1], L[2], L[3]
l1,l2,l3,l4,l5

D1, D2,  = symbols('Delta_top Delta_mu')

varname = zeros(4,2)
vev = Symbol('nu', real = True)
angphase = Symbol('theta', real = True)
angvev = Symbol('beta', real = True)
angcp = angphase
phasefac = exp(I*angphase)

phases= Matrix([[cos(angvev), sin(angvev)*phasefac, conjugate(sin(angvev)*phasefac)]])
for k in [1,2]:
    w=0
    for q in ['phi__+','varphi__0', 'a__0', 'phi__-']:
        varname[w,k-1] = Symbol(q + '_{}'.format(k), real = True)
        w+=1





H1 = 1/sqrt(2)*Matrix([[sqrt(2)*varname[0,0]],
            [vev*phases[0] + varname[1,0] + I*varname[2,0]]])
H2 = 1/sqrt(2)*Matrix([[sqrt(2)*varname[0,1]],
            [vev*phases[1] + varname[1,1] + I*varname[2,1]]])

H1dag = 1/sqrt(2)*Matrix([[sqrt(2)*varname[3,0], vev*phases[0] + varname[1,0] - I*varname[2,0]]])
H2dag = 1/sqrt(2)*Matrix([[sqrt(2)*varname[3,1], vev*phases[2] + varname[1,1] - I*varname[2,1]]])
x1 = (H1dag*H1)[0]
x2 = (H2dag*H2)[0]
x3 = (H1dag*H2)[0]
x4 = (H2dag*H1)[0]
V = Symbol('V')

x1,x2,x3,x4

z1, z2, z3, z4 = symbols('z_1 z_2 z_3 z_4', real = True)

w1, w2, w3, w4 = symbols('w_1 w_2 w_3 w_4', real = True)

Z1, Z2, Z3, Z4, Z5, Z6, Z7 = symbols('Z_1 Z_2 Z_3 Z_4 Z_5 Z_6 Z_7')

U = Matrix([[cos(angvev), sin(angvev)*exp(-I*angcp)],
           [-sin(angvev)*exp(I*angcp), cos(angvev)]])
A = U.det()
Ma1 = trigsimp(conjugate(Matrix((U.inv()*Matrix([z3,z4]))).T))

Ma2 = trigsimp(U.inv()*Matrix([z1,z2]))

#y1 = Ma1[0]*Ma2[0]
#y2 = Ma1[1]*Ma2[1]
#y3 = Ma1[0]*Ma2[1]
#y4 = Ma1[1]*Ma2[0]

y1 = w1*cos(angvev)**2 + w2*sin(angvev)**2 - w4*sin(angvev)*cos(angvev)*exp(I*angcp) - w3*sin(angvev)*cos(angvev)*exp(-I*angcp)
y2 = w1*sin(angvev)**2 + cos(angvev)*sin(angvev)*(w4*exp(I*angcp) + w3*exp(-I*angcp)) + w2*cos(angvev)**2

y3 = sin(angvev)*cos(angvev)*(w1 - w2)*exp(I*angcp) + cos(angvev)**2*w3 - sin(angvev)**2*w4*exp(2*I*angcp)

y4 = exp(-I*angphase)*sin(angvev)*cos(angvev)*(w1 - w2) + cos(angvev)**2*w4 - sin(angvev)**2*w3*exp(-2*I*angcp)


portHBT = 1/2*l1*y1**2 + 1/2*l2*y2**2 + l3*y1*y2 + l4*y4*y3 + 1/2*(conjugate(l5)*y4**2 + l5*y3**2)


Hp, Hm, Gp, Gm, g0, a0, p1, p2 = symbols('H__+ H__- G__+ G__- g__0 a__0 pi_1__0 pi_2__0')

varnameHBT = Matrix([[Hp, Hm, Gp, Gm, g0, a0, p1, p2]])

HB1 = 1/sqrt(2)*Matrix([[sqrt(2)*Gp],
            [vev + p1 + I*g0]])
HB2 = 1/sqrt(2)*Matrix([[sqrt(2)*Hp],
            [p2 + I*a0]])

HB1dag = 1/sqrt(2)*Matrix([[sqrt(2)*Gm, vev + p1 - I*g0]])
HB2dag = 1/sqrt(2)*Matrix([[sqrt(2)*Hm, p2 - I*a0]])
HB1, HB2, HB1dag, HB2dag
q1 = (HB1dag*HB1)[0]
q2 = (HB2dag*HB2)[0]
q3 = (HB1dag*HB2)[0]
q4 = (HB2dag*HB1)[0]
#q1, q2, q3, q4

def f(x,y):
    a=x
    for i in y:
        a=limit(a,i,0)
    return a

mH2 = Symbol('mu_H__2')
mHc2 = Symbol('mu_Hc__2')
mZ2 = Symbol('mu_Z__2')
mA2 = Symbol('mu_A__2')
mW2 = Symbol('mu_W__2')
mt2 = Symbol('mu_t__2')
RG = Symbol('Lambda_WG__2')
hI = Symbol('h_I')
mH = Symbol('mu__2', real = True)
RG = Symbol('Lambda_WG__2')
lpp = Symbol('kappa_1')
ksi = Symbol('xi_1')
lpp2 = Symbol('kappa_2')
ksi2 = Symbol('xi_2')

mH = Symbol('mu__2', real = True)
RG = Symbol('Lambda_WG__2')
lp = Symbol('kappa')

potHB = 1/2*Z1*q1**2 + 1/2*Z2*q2**2 + Z3*q1*q2 + Z4*q4*q3 + 1/2*(conjugate(Z5)*q4**2 + Z5*q3**2) + q1*(Z6*q3 + conjugate(Z6)*q4) + q2*(Z7*q3 + conjugate(Z7)*q4)

a4,a5,a6 = symbols('alpha_23 alpha_24 alpha_34')
angv = [a4,a5,a6]



y1p = q1*cos(angvev)**2 + q2*sin(angvev)**2 - q4*sin(angvev)*cos(angvev)*exp(I*angcp) - q3*sin(angvev)*cos(angvev)*exp(-I*angcp)

y2p = q1*sin(angvev)**2 + q2*cos(angvev)**2 + q4*sin(angvev)*cos(angvev)*exp(I*angcp) + q3*sin(angvev)*cos(angvev)*exp(-I*angcp)
Z2 = [y1p,y2p]

evenfield = int(input("Input which Phi is even:"))

inpfield = Z2[evenfield-1]
def t(x):
    li = expand(lpp*(q1 + q2))
    li2 = expand(lpp2*inpfield)

    
    newitems = [q for q in varnameHBT[:]]
    
    newitems.remove(x)
            
    li = li.subs([(k,0) for k in newitems])
    li2 = li2.subs([(k,0) for k in newitems])

    o = 0
    masses = [mH, mt2]
    cc = [li,li2]
    poteff = [mH**2*(ksi + log(mH/RG)), mt2**2*(ksi2 + log(mt2/RG))]
    for i in range(2):
        o+=diff(poteff[i], masses[i])*diff(cc[i],x)
    return o.subs(x,0)


def g(x,y):
    li = lpp*(q1 + q2)
    li2 = lpp2*inpfield
    
    newitems = [q for q in varnameHBT[:]]
    for i in [x,y]:
        if i in newitems:
            newitems.remove(i)
        else:
            break
    li,li2 = li.subs([(j,0) for j in newitems]), li2.subs([(j,0) for j in newitems])
    o = 0
    masses = [mH, mt2]
    cc = [li,li2]
    poteff = [mH**2*(ksi + log(mH/RG)), mt2**2*(ksi2 + log(mt2/RG))]
    for i in range(2):
        o+=diff(poteff[i],masses[i])*diff(cc[i],x,y) + diff(cc[i],y)*diff(cc[i],x)*diff(poteff[i],masses[i],2)
    return o.subs([(x,0), (y,0)])


mHp, mG, mg, ma1, mh1, mh2 = symbols('mu_H mu_G mu_g mu_a1 mu_h1 mu_h2')


def Treetadpole():
    tad1, tad2, tad3,tad4 = f(diff(potHB, p1),varnameHBT), f(diff(potHB, p2),varnameHBT), f(diff(potHB, a0),varnameHBT), f(diff(potHB, g0),varnameHBT)
    return tad1, tad2, tad3, tad4

def OneLoopTadPole():
    tad1p, tad2p, tad3p,tad4p = f(diff(potHB, p1),varnameHBT) + t(p1), f(diff(potHB, p2),varnameHBT) + t(p2), f(diff(potHB, a0),varnameHBT) + t(a0), f(diff(potHB, g0),varnameHBT) + t(g0)
    return tad1p, tad2p, tad3p, tad4p


def TreeMassmatrix(B = True):
    subvarq = [p1,p2,a0]
    Ad = zeros(len(subvarq), len(subvarq))
    Ap = zeros(len(subvarq), len(subvarq))
    q=0
    if B == True:
        pote = potHB
    if B == False:
        pote = potHB.subs([(Z1, 0), (Z6,0)])
    for i in subvarq:
        w=0
        for j in subvarq:
            Ap[q,w] = Derivative(V,i,j)
            Ad[q,w] = f(diff(pote,i,j), varnameHBT)
            w+=1
        q+=1
    
    
    subvarp = [Hp]
    subvarm = [Hm]
    A1 = zeros(len(subvarp), len(subvarp))
    A1p = zeros(len(subvarp), len(subvarp))
    q=0
    for i in subvarp:
        w=0
        for j in subvarm:
            A1p[q,w] = Derivative(V,i,j)
            A1[q,w] = f(diff(pote,i,j), varnameHBT)
            w+=1
        q+=1
    return  Ad, Ap,A1

def LoopMassMatrix():
    subvarq = [p1, p2,a0]
    C = zeros(len(subvarq), len(subvarq))
    Cp = zeros(len(subvarq), len(subvarq))
    q=0
    for i in subvarq:
        w=0
        for j in subvarq:
            Cp[q,w] = Derivative(V,i,j)
            C[q,w] = trigsimp(g(i,j)) 
            w+=1
        q+=1
    return C


def ChargedLoopMatrix():
    subvarp = [Hp]
    subvarm = [Hm]
    C = zeros(len(subvarp), len(subvarp))
    Cp = zeros(len(subvarp), len(subvarp))
    q=0
    for i in subvarp:
        w=0
        for j in subvarm:
            Cp[q,w] = Derivative(V,i,j)
            C[q,w] = trigsimp(f(g(i,j), varnameHBT)) 
            w+=1
        q+=1
    return C


def ZtoLambda(n):
    Zlist = [Z1, Z2, Z3, Z4, Z5, Z6, Z7]
#    n = int(input("Insert k value:"))
    items = [(w1,w1), (w2,w2), (w1,w2), (w3,w4), (w3,w3), (w1,w3), (w2,w3)]
    a,b = items[n-1][0], items[n-1][1]
    return diff(portHBT, a, b)

    
def ShortJacobiSol(i,j,w):
    c23, c24, c34, s23, s24, s34 = symbols('c_23 c_24 c_34 s_23 s_24 s_34')
    cosym = [c23, c24, c34]
    sisym = [s23, s24, s34]

    romat = eye(4) 

    romat[i-1,i-1] = cosym[w]
    romat[j-1,j-1] = cosym[w]
    romat[i-1,j-1] = -sisym[w]
    romat[j-1,i-1] = sisym[w]
    
    return romat 

#Rotmatrix = GivenMatrix(a6,3,4)*GivenMatrix(a5,2,4)*GivenMatrix(a4,2,3)


st, ct, sb, cb, c2b = symbols('s_Theta c_Theta s_beta c_beta c_A')
substitutes = [(sin(angcp-angphase), st), (cos(angcp-angphase), ct), (sin(angvev), sb), (cos(angvev), cb), (cos(2*angvev), c2b)]