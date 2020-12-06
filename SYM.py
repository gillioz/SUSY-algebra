#!/usr/bin/env python3.5

from SUSYalgebra import *

print("""
Fields
******
""")

print("Vector superfield:")

V = SuperField("WZvector")

W = [ -quarter * V.covD(a).covDbar2() for a in range(2) ]
Wbar = [ -quarter * V.covDbar(b).covD2() for b in range(2) ]

e = symbols("e")
exp2V = 2 * e * V + 2 * e**2 * V * V
exp2V.sort()

print("lambda:")
lam = [-2 * W[a].lowest() for a in range(2)]
lamC = [-2 * Wbar[b].lowest() for b in range(2)]
lam[0].print()
lam[1].print()
lamC[0].print()
lamC[1].print()
print()


print("A^mu:")
AA = [ [V.component[i][j]
       for j in range(1, 3)]
      for i in range(1, 3)]
A = [
    sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: AA[a[0]][b[0]],
        1, 1, 0)
    for mu in range(4)]
for x in A:
    x.sort()
A[0].print()
A[1].print()
A[2].print()
A[3].print()
print()

FF = [[sum(lambda a, b, rho: eta[nu][rho[0]],
          lambda a, b, rho: A[rho[0]].d(mu),
          0, 0 , 1)
      -sum(lambda a, b, rho: eta[mu][rho[0]],
           lambda a, b, rho: A[rho[0]].d(nu),
           0, 0 , 1)
      for nu in range(4)]
      for mu in range(4)]
for mu in range(4):
    for nu in range(4):
        FF[mu][nu].sort()
FFdual = [[ 0 * FF[mu][nu] for nu in range(4)] for mu in range(4)]
for mu in range(4):
    for nu in range(4):
        X = half * sum(
            lambda a, b, rho: eta[mu][rho[0]] * eta[nu][rho[1]] * LeviCivita(rho),
            lambda a, b, rho: FF[rho[2]][rho[3]],
            0, 0, 4
        )
        if X != 0:
            X.sort()
            FFdual[mu][nu] = X

print("D:")
D = quarter * sum(
    lambda a, b, mu: epsilon[a[0]][a[1]],
    lambda a, b, mu: V.covD(a[0]).covDbar2().covD(a[1]),
    2, 0, 0
).lowest()
D.sort()
D.print()
print()


print("Chiral superfield:")

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


# covariant derivatives
def Dphi(mu):
    return phi.d(mu) - half * e * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: phi * A[rho[0]],
        0, 0, 1
    )
def DphiC(mu):
    return phiC.d(mu) + half * e * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: phiC * A[rho[0]],
        0, 0, 1
    )
def Dchi(a, mu):
    return chi[a].d(mu) - half * e * I * sum(
        lambda aa, bb, rho: eta[mu][rho[0]],
        lambda aa, bb, rho: chi[a] * A[rho[0]],
        0, 0, 1
    )
def DchiC(b, mu):
    return chiC[b].d(mu) + half * e * I * sum(
        lambda aa, bb, rho: eta[mu][rho[0]],
        lambda aa, bb, rho: chiC[b] * A[rho[0]],
        0, 0, 1
    )
def D2phi(mu, nu):
    return Dphi(nu).d(mu) - half * e * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: Dphi(nu) * A[rho[0]],
        0, 0, 1
    )
def D2phiC(mu,nu):
    return DphiC(nu).d(mu) + half * e * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: DphiC(nu) * A[rho[0]],
        0, 0, 1
    )

print("""
Lagrangian
**********
""")

printbegin("Computation of the Lagrangian (1/4)...")
dL1 = 2 * sum(
    lambda a, b, mu: epsilon[a[1]][a[0]],
    lambda a, b, mu: W[a[0]] * W[a[1]],
    2, 0, 0
)
L1 = dL1.Fterm()
L1.sort()
printend()

printbegin("Computation of the Lagrangian (2/4)...")
dL2 = 2 * sum(
    lambda a, b, mu: epsilon[b[0]][b[1]],
    lambda a, b, mu: Wbar[b[0]] * Wbar[b[1]],
    0, 2, 0
)
L2 = dL2.Fbarterm()
L2.sort()
printend()


printbegin("Computation of the Lagrangian (3/4)...")
dL3 = Phibar * Phi + exp2V * Phibar * Phi
L3 = dL3.Dterm()
L3.sort()
printend()

printbegin("Computation of the Lagrangian (4/4)...")
L = L1 + L2 + L3
L.sort()
printend()


print("""
L = -1/2 F_munu F^munu
    - i lambda^* sigma^mu d_mu lambda + h.c.
    + D^2
    + D_mu phi^* D^mu phi
    - i chi^* sigma^mu D_mu chi + h.c.
    + e ( phi^* chi lambda + h.c.)
    + e/2 D phi^* phi
    + F^* F
""")
printbegin("Check : ")
X = (
    -half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
        0, 0, 4
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            lamC[b[0]] * lam[a[0]].d(mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            lamC[b[0]].d(mu[0]) * lam[a[0]],
        1, 1, 1
    )
    + D * D
    + sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: DphiC(mu[0]) * Dphi(mu[1]),
        0, 0, 2
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            chiC[b[0]] * Dchi(a[0], mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            DchiC(b[0], mu[0]) * chi[a[0]],
        1, 1, 1
    )
    + e * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            phiC * chi[a[0]] * lam[a[1]],
        2, 0, 0
    )
    + e * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            phi * chiC[b[0]] * lamC[b[1]],
        0, 2, 0
    )
    + half * e * D * phiC * phi
    + Fc * F
    # boundary terms
    -quarter * (phiC * phi).box()
    - L
)
X.sort()
printend([ X.iszero() ])



print("""
Field equations
***************
""")

printbegin("Computing the field equation for W...")
eomV = (
    -1 * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: W[a[1]].covD(a[0]),
        2, 0, 0
    )
    -1 * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Wbar[b[1]].covDbar(b[0]),
        0, 2, 0
    )
    + half * e * (Phibar * Phi + exp2V * Phibar * Phi)
)
eomV.sort()
eomD = eomV.lowest()
eomA = [
    half * sum(
        lambda aa, bb, rho: sigmabar[mu][bb[0]][aa[0]],
        lambda aa, bb, rho: eomV.component[1 + aa[0]][1 + bb[0]],
        1, 1, 0)
for mu in range(4)]
for x in eomA:
    x.sort()
printend()


printbegin("Computing the field equation for Phi...")
eomPhi = -quarter * (Phi + exp2V * Phi).covD2()
eomPhibar = -quarter * (Phibar + exp2V * Phibar).covDbar2()

eomF = eomPhi.lowest()
eomFc = eomPhibar.lowest()

eomphi = eomPhi.Fbarterm()
eomphiC = eomPhibar.Fterm()
printend()


print("""
Field equations in components:
    ...
""")
printbegin("Check: ")
X = (
    [eomD - (2 * D + half * e * phiC * phi)]
    + [eomA[mu]
        + 2 * sum(
            lambda aa, bb, rho: eta[mu][rho[0]] * eta[rho[1]][rho[2]],
            lambda aa, bb, rho: FF[rho[0]][rho[1]].d(rho[2]),
            0, 0, 3
        )
        - e * half * I * sum(
            lambda aa, bb, rho: eta[mu][rho[0]],
            lambda aa, bb, rho: phiC * Dphi(rho[0]) - DphiC(rho[0]) * phi,
            0, 0, 1
        )
        + e * sum(
            lambda a, b, rho:
                sigmabar[mu][b[0]][a[0]],
            lambda a, b, rho:
                chiC[b[0]] * chi[a[0]],
            1, 1, 0
        )
       for mu in range(4)]
    + [
        eomF - F,
        eomFc - Fc,
        eomphi - (
            half * e * D * phi
            - sum(
                lambda a, b, rho:
                    eta[rho[0]][rho[1]],
                lambda a, b, rho:
                    D2phi(rho[0],rho[1]),
                0, 0, 2
            )
            + e * sum(
                lambda a, b, mu:
                    epsilon[a[1]][a[0]],
                lambda a, b, mu:
                    chi[a[0]] * lam[a[1]],
                2, 0, 0
            )
        ),
        eomphiC - (
            half * e * D * phiC
            - sum(
                lambda a, b, rho:
                    eta[rho[0]][rho[1]],
                lambda a, b, rho:
                    D2phiC(rho[0],rho[1]),
                0, 0, 2
            )
            + e * sum(
                lambda a, b, mu:
                    epsilon[b[0]][b[1]],
                lambda a, b, mu:
                    chiC[b[0]] * lamC[b[1]],
                0, 2, 0
            )
        )
    ]
)
for x in X:
    x.sort()
printend([x.iszero() for x in X])




print("""
Supercurrent
************
""")

printbegin("Computing the conserved J_ab...")
J = [[
    2 * W[a] * Wbar[b]
    - I * Rational(1,12) * Phibar.dd(a, b) * Phi
    - I * Rational(1,12) * Phibar.dd(a, b) * Phi * exp2V
    + I * Rational(1,12) * Phibar * Phi.dd(a,b)
    + I * Rational(1,12) * Phibar * Phi.dd(a,b) * exp2V
    + Rational(1,24) * Phibar.covDbar(b) * Phi.covD(a)
    + Rational(1,24) * Phibar.covDbar(b) * Phi.covD(a) * exp2V
    - e * Rational(1,12) * Phibar * Phi.covD(a) * V.covDbar(b)
    - e * Rational(1,12) * Phibar * Phi.covD(a) * V.covDbar(b) * exp2V
    + e * Rational(1,12) * Phibar.covDbar(b) * Phi * V.covD(a)
    + e * Rational(1,12) * Phibar.covDbar(b) * Phi * V.covD(a) * exp2V
    + e * Rational(1,12) * Phibar * Phi * V.covDbar(b).covD(a)
    + e * Rational(1,12) * Phibar * Phi * V.covDbar(b).covD(a) * exp2V
    - e * Rational(1,12) * Phibar * Phi * V.covD(a).covDbar(b)
    - e * Rational(1,12) * Phibar * Phi * V.covD(a).covDbar(b) * exp2V
    + e**2 * Rational(1,6) * Phibar * Phi * V.covDbar(b) * V.covD(a)
    + e**2 * Rational(1,6) * Phibar * Phi * V.covDbar(b) * V.covD(a) * exp2V
    for b in range(2)]
    for a in range(2)]
for Jrow in J:
    for j in Jrow:
        j.sort()
printend()

printbegin("Computing D^a J_ab...")
DJ = [ sum(lambda a, bb, mu: epsilon[a[0]][a[1]],
           lambda a, bb, mu: J[a[0]][b].covD(a[1]),
           2, 0, 0)
       for b in range(2) ]
printend()
printbegin("Computing Dbar^b J_ab...")
DbarJ = [ sum(lambda aa, b, mu: epsilon[b[0]][b[1]],
              lambda aa, b, mu: J[a][b[0]].covDbar(b[1]),
              0, 2, 0)
          for a in range(2) ]
printend()

printbegin("D^a J_ab = 0 (up to field equations) : ")
X = [
    DJ[b]
    + Wbar[b] * eomV
    + Rational(1,12) * Phibar * eomPhi.covDbar(b)
    - Rational(1,6) * Phibar.covDbar(b) * eomPhi
    - e * half * V.covDbar(b) * Phibar * eomPhi
    for b in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])
printbegin("Dbar^b J_ab = 0 (up to field equations) : ")
X = [
    DbarJ[a]
    + W[a] * eomV
    + Rational(1,12) * Phi * eomPhibar.covD(a)
    - Rational(1,6) * Phi.covD(a) * eomPhibar
    - e * half * V.covD(a) * Phi * eomPhibar
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
J^mu = -i/6 ( phi^* D^mu phi - D^mu phi^* phi )
       -1/2 lambda^* (sigmabar^mu) lambda
       +1/6 chi^* (sigmabar^mu) chi
""")
printbegin("Check: ")
X = [
    Rational(1,6) * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: phiC * Dphi(rho[0]),
        0, 0, 1
    )
    - Rational(1,6) * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: DphiC(rho[0]) * phi,
        0, 0, 1
    )
    - half * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]] * lam[a[0]],
        1, 1, 0
    )
    + Rational(1,6) * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * chi[a[0]],
        1, 1, 0
    )
    - Jmu[mu].lowest()
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
    + half * I * sum(
        lambda a, b, rho: epsilon[a[1]][a[0]],
        lambda a, b, rho: W[a[0]] * eomV.covD(a[1]),
        2, 0, 0
    )
    - half * I * sum(
        lambda a, b, rho: epsilon[b[0]][b[1]],
        lambda a, b, rho: Wbar[b[0]] * eomV.covDbar(b[1]),
        0, 2, 0
    )
    + Rational(1,12) * I * Phibar.covDbar2() * eomPhi
    + Rational(1,24) * I * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Phibar.covDbar(b[0]) * eomPhi.covDbar(b[1]),
        0, 2, 0
    )
    - Rational(1,24) * I * Phibar * eomPhi.covDbar2()
    + quarter * I * e * V.covDbar2() * Phibar * eomPhi
    + quarter * I * e * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: V.covDbar(b[0]) * Phibar * eomPhi.covDbar(b[1]),
        0, 2, 0
    )
    + quarter * I * e * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: V.covDbar(b[0]) * Phibar.covDbar(b[1]) * eomPhi,
        0, 2, 0
    )
    - Rational(1,12) * I * Phi.covD2() * eomPhibar
    - Rational(1,24) * I * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: Phi.covD(a[0]) * eomPhibar.covD(a[1]),
        2, 0, 0
    )
    + Rational(1,24) * I * Phi * eomPhibar.covD2()
    - quarter * I * e * V.covD2() * Phi * eomPhibar
    - quarter * I * e * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: V.covD(a[0]) * Phi * eomPhibar.covD(a[1]),
        2, 0, 0
    )
    - quarter * I * e * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: V.covD(a[0]) * Phi.covD(a[1]) * eomPhibar,
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
T^munu = ...
""")
printbegin("Check: ")
X = [
    half * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho: FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
        0, 0, 4
    )
    - 2 * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: FF[mu][rho[0]] * FF[nu][rho[1]],
        0, 0, 2
    )
    + eta[mu][nu] * D * D
    - half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]] * lam[a[0]].d(nu),
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]] * lam[a[0]].d(mu),
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]].d(nu) * lam[a[0]],
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]].d(mu) * lam[a[0]],
        1, 1, 1
    )
    -Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: DphiC(rho[0]) * Dphi(rho[1]),
        0, 0, 2
    )
    + Rational(2,3) * DphiC(mu) * Dphi(nu)
    + Rational(2,3) * DphiC(nu) * Dphi(mu)
    - Rational(1,6) * phiC * D2phi(mu, nu)
    - Rational(1,6) * phiC * D2phi(nu, mu)
    - Rational(1,6) * D2phiC(mu, nu) * phi
    - Rational(1,6) * D2phiC(nu, mu) * phi
    + Rational(1,3) * I * eta[mu][nu] * sum(
        lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * Dchi(a[0], rho[0]),
        1, 1, 1
    )
    - Rational(1,3) * I * eta[mu][nu] * sum(
        lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: DchiC(b[0], rho[0]) * chi[a[0]],
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * Dchi(a[0], nu),
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: chiC[b[0]] * Dchi(a[0], mu),
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: DchiC(b[0], nu) * chi[a[0]],
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
        lambda a, b, rho: DchiC(b[0], mu) * chi[a[0]],
        1, 1, 1
    )
    - Rational(1,3) * eta[mu][nu] * Fc * F
    + Rational(1,3) * e * eta[mu][nu] * D * phiC * phi
    + Rational(1,6) * e * eta[mu][nu] * sum(
        lambda a, b, rho: epsilon[a[1]][a[0]],
        lambda a, b, rho: phiC * lam[a[0]] * chi[a[1]],
        2, 0, 0
    )
    + Rational(1,6) * e * eta[mu][nu] * sum(
        lambda a, b, rho: epsilon[b[0]][b[1]],
        lambda a, b, rho: phi * lamC[b[0]] * chiC[b[1]],
        0, 2, 0
    )
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


