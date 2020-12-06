#!/usr/bin/env python3.5

from SUSYalgebra import *
#from sympy import solve

print("""
Fields
******
""")

Phi = SuperField("chiral")
Phibar = SuperField("antichiral")

phi = Phi.lowest()
phiC = Phibar.lowest()
phi.print()
phiC.print()
chi = [ -half * Phi.component[a][0] for a in range(1,3) ]
chiC = [ half * Phibar.component[0][b] for b in range(1,3) ]
chi[0].print()
chiC[0].print()
chi[1].print()
chiC[1].print()
F = Phi.Fterm()
Fc = Phibar.Fbarterm()
F.print()
Fc.print()
print()



print("Coupling:")
G = SuperField("chiralcoupling")
Gbar = SuperField("antichiralcoupling")
g = half * sqrt(2) * (G + Gbar).lowest()
g.sort()
g.print()
theta = -I * half * sqrt(2) * (G - Gbar).lowest()
theta.sort()
theta.print()
xi = [ -half * G.component[a][0] for a in range(1,3) ]
xiC = [ half * Gbar.component[0][b] for b in range(1,3) ]
xi[0].print()
xiC[0].print()
xi[1].print()
xiC[1].print()
f = G.Fterm()
f.print()
fC = Gbar.Fbarterm()
fC.print()
print()


print("""
Lagrangian
**********
""")

dL = half * sqrt(2) * (G + Gbar) * Phibar * Phi
dL.sort()

L = dL.Dterm()

print("""
L = g d_mu phi^* d^mu phi
    - 1/2 (d_mu d^mu g) phi^* phi
    + i/2 d_mu theta ( d^mu phi^* phi - phi^* d^mu phi )
    + g F^* F
    + 1/sqrt(2) ( f phi F^* + f^* phi^* F )
    + ...
""")
printbegin("Check : ")
X = (
    sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g * phiC.d(mu[0]) * phi.d(mu[1]),
        0, 0, 2
    )
    - half * g.box() * phiC * phi
    - half * I * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.d(mu[0]) * phiC * phi.d(mu[1]),
        0, 0, 2
    )
    + half * I * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.d(mu[0]) * phiC.d(mu[1]) * phi,
        0, 0, 2
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * chiC[b[0]] * chi[a[0]].d(mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * chiC[b[0]].d(mu[0]) * chi[a[0]],
        1, 1, 1
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            theta.d(mu[0]) * chiC[b[0]] * chi[a[0]],
        1, 1, 1
    )
    + g * Fc * F
    + half * sqrt(2) * f * phi * Fc
    + half * sqrt(2) * fC * phiC * F
    - I * sqrt(2) * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * phiC * chi[a[0]].d(mu[0]),
        1, 1, 1
    )
    - I * sqrt(2) * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * phi * chiC[b[0]].d(mu[0]),
        1, 1, 1
    )
    - sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * chiC[b[1]] * F,
        0, 2, 0
    )
    - sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * chi[a[1]] * Fc,
        2, 0, 0
    )
    # boundary terms
    - quarter * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            (g * phiC.d(mu[0]) * phi).d(mu[1]),
        0, 0, 2
    )
    - quarter * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            (g * phiC * phi.d(mu[0])).d(mu[1]),
        0, 0, 2
    )
    + quarter * sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        (g.d(mu[0]) * phiC * phi).d(mu[1]),
        0, 0, 2
    )
    + I * half * sqrt(2) * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (xiC[b[0]] * phiC * chi[a[0]]).d(mu[0]),
        1, 1, 1
    )
    + I * half * sqrt(2) * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (xi[a[0]] * phi * chiC[b[0]]).d(mu[0]),
        1, 1, 1
    )
    - L
)
X.sort()
printend([ X.iszero() ])



print("""
Field equations
***************
""")

eomPhi = -quarter * (half * sqrt(2) * (G + Gbar) * Phi).covD2()
eomPhibar = -quarter * (half * sqrt(2) * (G + Gbar) * Phibar).covDbar2()

eomF = eomPhi.lowest().bosonicpart()
eomFc = eomPhibar.lowest().bosonicpart()

eomphi = eomPhi.Fbarterm().bosonicpart()
eomphiC = eomPhibar.Fterm().bosonicpart()


eomG = -quarter * (half * sqrt(2) * Phi * Phibar).covD2()
eomGbar = -quarter * (half * sqrt(2) * Phi * Phibar).covDbar2()

eomf = eomG.lowest().bosonicpart()
eomfC = eomGbar.lowest().bosonicpart()

eomg1 = eomG.Fbarterm().bosonicpart()
eomg2 = eomGbar.Fterm().bosonicpart()

eomg = half * sqrt(2) * (eomg1 + eomg2)
eomg.sort()
eomtheta = -half * sqrt(2) * I * (eomg1 - eomg2)
eomtheta.sort()


print("""
Field equations in components:
    sqrt(2) g F + f phi = 0
    sqrt(2) g F^* + f^* phi^* = 0
    g d_mu d^mu phi + d_mu (g + i theta) d^mu phi
        + 1/2 d_mu d^mu (g + i theta) phi - sqrt(2)/2 f^* F = 0
    ...
""")

printbegin("Check : ")
X = [
    eomF - (g * F + half * sqrt(2) * f * phi),
    eomFc - (g * Fc + half * sqrt(2) * fC * phiC),
    eomphi - (
        -1 * g * phi.box()
        - sum(
            lambda a, b, mu: eta[mu[0]][mu[1]],
            lambda a, b, mu: (g + I * theta).d(mu[0]) * phi.d(mu[1]),
            0, 0, 2
        )
        - half * (g + I * theta).box() * phi
        + half * sqrt(2) * fC * F
    ),
    eomphiC - (
        -1 * g * phiC.box()
        - sum(
            lambda a, b, mu: eta[mu[0]][mu[1]],
            lambda a, b, mu: (g - I * theta).d(mu[0]) * phiC.d(mu[1]),
            0, 0, 2
        )
        - half * (g - I * theta).box() * phiC
        + half * sqrt(2) * f * Fc
    ),
    eomf - half * sqrt(2) * phiC * F,
    eomfC - half * sqrt(2) * phi * Fc,
    eomg - (
        - half * (phiC * phi.box() + phiC.box() * phi)
        + Fc * F
    ),
    eomtheta - (
        half * I * (phiC * phi.box() - phiC.box() * phi)
    )
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])




print("""
Supercurrent
************
""")


printbegin("Computing the conserved J_ab...")
J = [[
    - I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar.dd(a,b) * Phi
    + I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi.dd(a,b)
    + Rational(1,48) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi.covD(a)
    - Rational(1,48) * sqrt(2) * G.covD(a) * Phibar.covDbar(b) * Phi
    + Rational(1,48) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi.covD(a)
    + I * Rational(1,24) * sqrt(2) * G.dd(a,b) * Phibar * Phi
    - I * Rational(1,24) * sqrt(2) * Gbar.dd(a,b) * Phibar * Phi
    for b in range(2)] for a in range(2)]
for Jrow in J:
    for j in Jrow:
        j.sort()
printend()


printbegin("Computing D^a J_ab...")
DJ = [ sum(lambda a, bb, mu: epsilon[a[0]][a[1]],
           lambda a, bb, mu: J[a[0]][b].covD(a[1]),
           2, 0, 0)
       for b in range(2) ]
for dj in DJ:
    dj.sort()
printend()
printbegin("Computing Dbar^b J_ab...")
DbarJ = [ sum(lambda aa, b, mu: epsilon[b[0]][b[1]],
              lambda aa, b, mu: J[a][b[0]].covDbar(b[1]),
              0, 2, 0)
          for a in range(2) ]
for dj in DbarJ:
    dj.sort()
printend()

printbegin("D^a J_ab = 0 (up to field equations) : ")
X = [
    DJ[b]
    + Rational(1,12) * Phibar * eomPhi.covDbar(b)
    - Rational(1,6) * Phibar.covDbar(b) * eomPhi
    - quarter * Gbar.covDbar(b) * eomG
    for b in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])

printbegin("Dbar^b J_ab = 0 (up to field equations) : ")
X = [
    DbarJ[a]
    + Rational(1,12) * Phi * eomPhibar.covD(a)
    - Rational(1,6) * Phi.covD(a) * eomPhibar
    - quarter * G.covD(a) * eomGbar
    for a in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])


printbegin("Computing J^mu...")
Jmu = [sum(lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
           lambda a, b, rho: J[a[0]][b[0]],
           1, 1, 0)
       for mu in range(4)]
for j in Jmu:
    j.sort()
printend()


print("""
J^mu = i/6 g ( d^mu phi^* phi - phi^* d^mu phi )
       + 1/6 d_mu theta phi^* phi
       + ...
""")
printbegin("Check: ")
X = [
    - Rational(1,6) * theta.d(mu) * phiC * phi
    - Rational(1,6) * I * g * phiC.d(mu) * phi
    + Rational(1,6) * I * g * phiC * phi.d(mu)
    - sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: Jmu[rho[0]].lowest().bosonicpart(),
        0, 0, 1
    )
    for mu in range(4)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])

printbegin("Computing d_mu J^mu...")
dJ = sum(lambda a, b, mu: 1,
         lambda a, b, mu: Jmu[mu[0]].d(mu[0]),
         0, 0, 1)
dJ.sort()
printend()

printbegin("d_mu J^mu = 0 (up to field equations) : ")
X = (
    dJ
    + Rational(1,12) * I * Phibar.covDbar2() * eomPhi
    + Rational(1,24) * I * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Phibar.covDbar(b[0]) * eomPhi.covDbar(b[1]),
        0, 2, 0
    )
    - Rational(1,24) * I * Phibar * eomPhi.covDbar2()
    - Rational(1,12) * I * Phi.covD2() * eomPhibar
    - Rational(1,24) * I * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: Phi.covD(a[0]) * eomPhibar.covD(a[1]),
        2, 0, 0
    )
    + Rational(1,24) * I * Phi * eomPhibar.covD2()
    + Rational(1,8) * I * Gbar.covDbar2() * eomG
    + Rational(1,8) * I * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Gbar.covDbar(b[0]) * eomG.covDbar(b[1]),
        0, 2, 0
    )
    - Rational(1,8) * I * G.covD2() * eomGbar
    - Rational(1,8) * I * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: G.covD(a[0]) * eomGbar.covD(a[1]),
        2, 0, 0
    )
)
X.sort()
printend([X.iszero()])


print("""
Stress tensor
*************
""")

Tasymmetric = [
    [
        sum(lambda a, b, rho:
                sigmabar[mu][b[0]][a[0]],
            lambda a, b, rho:
                Jmu[nu].component[a[0] + 1][b[0] + 1],
            1, 1, 0)
        for nu in range(4)
    ]
    for mu in range(4)
]

T = [[half * Tasymmetric[mu][nu] + half * Tasymmetric[nu][mu]
      for nu in range(4)]
     for mu in range(4)]
for mu in range(4):
    for nu in range(4):
        T[mu][nu].sort()

print("""
T^munu =
    ...
""")
printbegin("Check: ")
X = [
    -Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g * phiC.d(rho[0]) * phi.d(rho[1]),
        0, 0, 2
    )
    + Rational(2,3) * g * phiC.d(mu) * phi.d(nu)
    + Rational(2,3) * g * phiC.d(nu) * phi.d(mu)
    - Rational(1,3) * g * phiC * phi.d(mu).d(nu)
    - Rational(1,3) * g * phiC.d(mu).d(nu) * phi
    - Rational(1,6) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * phiC * phi.d(rho[1]),
        0, 0, 2
    )
    + Rational(1,6) * g.d(mu) * phiC * phi.d(nu)
    + Rational(1,6) * g.d(nu) * phiC * phi.d(mu)
    - Rational(1,6) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * phiC.d(rho[1]) * phi,
        0, 0, 2
    )
    + Rational(1,6) * g.d(mu) * phiC.d(nu) * phi
    + Rational(1,6) * g.d(nu) * phiC.d(mu) * phi
    + Rational(1,6) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC * phi.d(rho[1]),
        0, 0, 2
    )
    - half * I * theta.d(mu) * phiC * phi.d(nu)
    - half * I * theta.d(nu) * phiC * phi.d(mu)
    - Rational(1,6) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC.d(rho[1]) * phi,
        0, 0, 2
    )
    + half * I * theta.d(mu) * phiC.d(nu) * phi
    + half * I * theta.d(nu) * phiC.d(mu) * phi
    - Rational(1,3) * eta[mu][nu] * g * Fc * F
    - Rational(1,6) * sqrt(2) * eta[mu][nu] * f * phi * Fc
    - Rational(1,6) * sqrt(2) * eta[mu][nu] * fC * phiC * F
    - Rational(1,3) * g.d(mu).d(nu) * phiC * phi
    - sum(
        lambda a, b, rho: eta[mu][rho[0]] * eta[nu][rho[1]],
        lambda a, b, rho: T[rho[0]][rho[1]].bosonicpart(),
        0, 0, 2
    )
    for mu in range(4)
    for nu in range(4)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])



