#!/usr/bin/env python3.5

from SUSYalgebra import *

print("""
Fields
******
""")

V = SuperField("vector")

W = [ -quarter * V.covD(a).covDbar2() for a in range(2) ]
Wbar = [ -quarter * V.covDbar(b).covD2() for b in range(2) ]


print("lambda:")
lam = [-2 * W[a].lowest() for a in range(2)]
lamC = [-2 * Wbar[b].lowest() for b in range(2)]
lam[0].print()
lam[1].print()
lamC[0].print()
lamC[1].print()
print()


print("A^mu:")
A = [ [V.component[i][j]
       for j in range(1, 3)]
      for i in range(1, 3)]
Amu = [
    sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: A[a[0]][b[0]],
        1, 1, 0)
    for mu in range(4)]
for x in Amu:
    x.sort()
Amu[0].print()
Amu[1].print()
Amu[2].print()
Amu[3].print()
print()

FF = [[sum(lambda a, b, rho: eta[nu][rho[0]],
          lambda a, b, rho: Amu[rho[0]].d(mu),
          0, 0 , 1)
      -sum(lambda a, b, rho: eta[mu][rho[0]],
           lambda a, b, rho: Amu[rho[0]].d(nu),
           0, 0 , 1)
      for nu in range(4)]
      for mu in range(4)]

print("D:")
D = quarter * sum(
    lambda a, b, mu: epsilon[a[0]][a[1]],
    lambda a, b, mu: V.covD(a[0]).covDbar2().covD(a[1]),
    2, 0, 0
).lowest()
D.sort()
D.print()
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

printbegin("Computation of the Lagrangian (1/3)...")
dL1 = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[a[1]][a[0]],
    lambda a, b, mu: G * W[a[0]] * W[a[1]],
    2, 0, 0
)
L1 = dL1.Fterm()
L1.sort()
printend()

printbegin("Computation of the Lagrangian (2/3)...")
dL2 = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[b[0]][b[1]],
    lambda a, b, mu: Gbar * Wbar[b[0]] * Wbar[b[1]],
    0, 2, 0
)
L2 = dL2.Fbarterm()
L2.sort()
printend()

printbegin("Computation of the Lagrangian (3/3)...")
L = L1 + L2
L.sort()
printend()


print("""
L = - g/4 F_{mu nu} F^{mu nu}
    - theta/8 eps^{mu nu rho tau} F_{mu nu} F_{rho tau}
    + i/2 (g + i theta) d_ab lambda^*_b lambda_a
    - i/2 (g - i theta) lambda^*_b d_ab lambda_a
    + ...
    + g/2 D^2
""")
printbegin("Check : ")
X = (
    -half * g * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
        0, 0, 4
    )
    - quarter * theta * sum(
        lambda a, b, mu: LeviCivita(mu),
        lambda a, b, mu: FF[mu[0]][mu[1]] * FF[mu[2]][mu[3]],
        0, 0, 4
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * lamC[b[0]] * lam[a[0]].d(mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * lamC[b[0]].d(mu[0]) * lam[a[0]],
        1, 1, 1
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            theta.d(mu[0]) * lamC[b[0]] * lam[a[0]],
        1, 1, 1
    )
    + half * sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]] * epsilon[a[0]][a[1]],
        lambda a, b, mu:
            xiC[b[0]] * lamC[b[1]] * FF[mu[0]][mu[1]],
        2, 2, 2
    )
    + half * sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]] * epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xi[a[0]] * lam[a[1]] * FF[mu[0]][mu[1]],
        2, 2, 2
    )
    + half * sqrt(2) * f * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: lam[a[0]] * lam[a[1]],
        2, 0, 0
    )
    + half * sqrt(2) * fC * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: lamC[b[0]] * lamC[b[1]],
        0, 2, 0
    )
    + g * D * D
    + sqrt(2) * D * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: xi[a[0]] * lam[a[1]],
        2, 0, 0
    )
    + sqrt(2) * D * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: xiC[b[0]] * lamC[b[1]],
        0, 2, 0
    )
    # boundary terms
    - sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (theta * lamC[b[0]] * lam[a[0]]).d(mu[0]),
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

printbegin("D^a W_a = Dbar_b W^b : ")
DW = sum(
    lambda a, b, mu: -quarter * epsilon[a[1]][a[0]],
    lambda a, b, mu: W[a[1]].covD(a[0]),
    2, 0, 0
)
DW.sort()
DbarWbar = sum(
    lambda a, b, mu: -quarter * epsilon[b[0]][b[1]],
    lambda a, b, mu: Wbar[b[1]].covDbar(b[0]),
    0, 2, 0
)
DbarWbar.sort()
DbarWbar.check()
X = DW - DbarWbar
X.sort()
printend([ X.iszero() ])


printbegin("Computing the field equation for W...")
eomV = (
    - sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: (G * W[a[1]]).covD(a[0]),
        2, 0, 0
    )
    - sqrt(2) * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: (Gbar * Wbar[b[1]]).covDbar(b[0]),
        0, 2, 0
    )
)
eomV.sort()
eomD = eomV.lowest().bosonicpart()
eomA = [
    half * sum(
        lambda aa, bb, rho: sigmabar[mu][bb[0]][aa[0]],
        lambda aa, bb, rho: eomV.component[1 + aa[0]][1 + bb[0]].bosonicpart(),
        1, 1, 0)
for mu in range(4)]
for x in eomA:
    x.sort()
printend()


printbegin("Computing the field equation for g...")
eomG = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[b[0]][b[1]],
    lambda a, b, mu: Wbar[b[0]] * Wbar[b[1]],
    0, 2, 0
)
eomG.sort()
eomGbar = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[a[1]][a[0]],
    lambda a, b, mu: W[a[0]] * W[a[1]],
    2, 0, 0
)
eomGbar.sort()

eomf = eomG.lowest().bosonicpart()
eomfC = eomGbar.lowest().bosonicpart()

eomg1 = eomG.Fbarterm().bosonicpart()
eomg2 = eomGbar.Fterm().bosonicpart()

eomg = half * sqrt(2) * (eomg1 + eomg2)
eomg.sort()
eomtheta = -half * sqrt(2) * I * (eomg1 - eomg2)
eomtheta.sort()

printend()


print("""
Field equations in components (bosonic part only):
    g D = 0
    d_nu ( g F^munu)
        + 1/2 epsilon^munurhotau d_nu theta F_rhotau = 0
    -1/4 F_munu F^munu + 1/2 D^2 = 0
    1/8 epsilon^munurhotau F_munu F_rhotau = 0
""")
printbegin("Check: ")
X = (
    [eomD - 2 * g * D]
    + [eomA[nu]
       + 2 * sum(
            lambda aa, bb, mu: eta[nu][mu[0]] * eta[mu[1]][mu[2]],
            lambda aa, bb, mu: (g * FF[mu[0]][mu[1]]).d(mu[2]),
            0, 0, 3)
       + sum(
            lambda aa, bb, mu: LeviCivita([nu] + mu),
            lambda aa, bb, mu: theta.d(mu[0]) * FF[mu[1]][mu[2]],
            0, 0, 3)
       for nu in range(4)]
    + [
        eomf,
        eomfC,
        eomg + half * sum(
            lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
            lambda a, b, mu: FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
        0, 0, 4)
        - D * D,
        eomtheta + quarter * sum(
            lambda a, b, mu: LeviCivita(mu),
            lambda a, b, mu: FF[mu[0]][mu[1]] * FF[mu[2]][mu[3]],
        0, 0, 4)
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
J = [[sqrt(2) * (G + Gbar) * W[a] * Wbar[b]
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
    - quarter * Gbar.covDbar(b) * eomG
    for b in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])
printbegin("Dbar^b J_ab = 0 (up to field equations) : ")
X = [
    DbarJ[a]
    + W[a] * eomV
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
J^mu = lambda^* (sigmabar^mu) lambda
""")
printbegin("Check: ")
X = [
    - half * g * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: lamC[b[0]] * lam[a[0]],
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
    + Rational(1,8) * I * Gbar.covDbar2() * eomG
    - Rational(1,8) * I * G.covD2() * eomGbar
    - Rational(1,8) * I * sum(
        lambda a, b, rho: epsilon[a[1]][a[0]],
        lambda a, b, rho: G.covD(a[0]) * eomGbar.covD(a[1]),
        2, 0, 0
    )
    + Rational(1,8) * I * sum(
        lambda a, b, rho: epsilon[b[0]][b[1]],
        lambda a, b, rho: Gbar.covDbar(b[0]) * eomG.covDbar(b[1]),
        0, 2, 0
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
        lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
        0, 0, 4
    )
    - 2 * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g * FF[mu][rho[0]] * FF[nu][rho[1]],
        0, 0, 2
    )
    + eta[mu][nu] * g * D * D
#    - I * sum(
#        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: g * lamC[b[0]] * lam[a[0]].d(nu),
#        1, 1, 1
#    )
#    - I * sum(
#        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: g * lamC[b[0]] * lam[a[0]].d(mu),
#        1, 1, 1
#    )
#    + I * sum(
#        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: g * lamC[b[0]].d(nu) * lam[a[0]],
#        1, 1, 1
#    )
#    + I * sum(
#        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: g * lamC[b[0]].d(mu) * lam[a[0]],
#        1, 1, 1
#    )
#    + sum(
#        lambda a, b, rho: eta[mu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: theta.d(nu) * lamC[b[0]] * lam[a[0]],
#        1, 1, 1
#    )
#    + sum(
#        lambda a, b, rho: eta[nu][rho[0]] * sigmabar[rho[0]][b[0]][a[0]],
#        lambda a, b, rho: theta.d(mu) * lamC[b[0]] * lam[a[0]],
#        1, 1, 1
#    )
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



print("""
Static configuration:
T^00 = g/4 F_ij F_ij + 8 g D^2
""")
printbegin("Check: ")
X = (
    half * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
        0, 0, 4
    )
    + g * D * D
    - T[0][0].bosonicpart()
)
X.sort()
printend([X.staticpart().iszero()])
