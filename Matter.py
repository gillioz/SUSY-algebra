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




print("""
Lagrangian
**********
""")

dL = Phibar * Phi
dL.sort()

L = dL.Dterm()

print("""
L = d_mu phi^* d^mu phi
    - i chi^* sigma^mu d_mu chi + h.c.
    + F^* F
    + ...
""")
printbegin("Check : ")
X = (
    sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            phiC.d(mu[0]) * phi.d(mu[1]),
        0, 0, 2
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            chiC[b[0]] * chi[a[0]].d(mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            chiC[b[0]].d(mu[0]) * chi[a[0]],
        1, 1, 1
    )
    + Fc * F
    # boundary terms
    - quarter * (phiC * phi).box()
    - L
)
X.sort()
printend([ X.iszero() ])



print("""
Field equations
***************
""")

eomPhi = -quarter * Phi.covD2()
eomPhibar = -quarter * Phibar.covDbar2()

eomF = eomPhi.lowest()
eomFc = eomPhibar.lowest()

eomphi = eomPhi.Fbarterm()
eomphiC = eomPhibar.Fterm()


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
    eomF - F,
    eomFc - Fc,
    eomphi + phi.box(),
    eomphiC + phiC.box()
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
    - I * Rational(1,12) * Phibar.dd(a,b) * Phi
    + I * Rational(1,12) * Phibar * Phi.dd(a,b)
    + Rational(1,24) * Phibar.covDbar(b) * Phi.covD(a)
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
       + 1/6 chi^* sigma^mu chi
""")
printbegin("Check: ")
X = [
    - Rational(1,6) * I * (phiC.d(mu) * phi - phiC * phi.d(mu))
    + Rational(1,6) * sum(
        lambda a, b, rho:
            eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho:
            chiC[b[0]] * chi[a[0]],
        1, 1, 1
    )
    - sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: Jmu[rho[0]].lowest(),
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
        lambda a, b, rho: phiC.d(rho[0]) * phi.d(rho[1]),
        0, 0, 2
    )
    + Rational(2,3) * phiC.d(mu) * phi.d(nu)
    + Rational(2,3) * phiC.d(nu) * phi.d(mu)
    - Rational(1,3) * phiC * phi.d(mu).d(nu)
    - Rational(1,3) * phiC.d(mu).d(nu) * phi
    + Rational(1,3) * I * eta[mu][nu] * sum(
        lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * chi[a[0]].d(rho[0]),
        1, 1, 1
    )
    - Rational(1,3) * I * eta[mu][nu] * sum(
        lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]].d(rho[0]) * chi[a[0]],
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * chi[a[0]].d(nu),
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * chi[a[0]].d(mu),
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]].d(nu) * chi[a[0]],
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]].d(mu) * chi[a[0]],
        1, 1, 1
    )
    - Rational(1,3) * eta[mu][nu] * Fc * F
    - sum(
        lambda a, b, rho: eta[mu][rho[0]] * eta[nu][rho[1]],
        lambda a, b, rho: T[rho[0]][rho[1]],
        0, 0, 2
    )
    for mu in range(4)
    for nu in range(4)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])



