#!/usr/bin/env python3.5

from SUSYalgebra import *
#from sympy import solve


print("""
Fields
******
""")


print("Vector superfield:")

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


b0, c1, c2 = symbols("b0 c1 c2")


print("""
Lagrangian
**********
""")


printbegin("Computation of the Lagrangian (gauge part)...")
dLgauge1 = b0 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[a[1]][a[0]],
    lambda a, b, mu: G * W[a[0]] * W[a[1]],
    2, 0, 0
)
Lgauge1 = dLgauge1.Fterm()
Lgauge1.sort()
dLgauge2 = b0 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[b[0]][b[1]],
    lambda a, b, mu: Gbar * Wbar[b[0]] * Wbar[b[1]],
    0, 2, 0
)
Lgauge2 = dLgauge2.Fbarterm()
Lgauge2.sort()
printend()

printbegin("Computation of the Lagrangian (matter part)...")
dLmatter = b0 * quarter * sqrt(2) * (G + Gbar) * Phibar * Phi
dLmatter.sort()
Lmatter = dLmatter.Dterm()
printend()

printbegin("Computation of the Lagrangian (vacuum part)...")
dLvac1 = c1 * sum(
    lambda a, b, mu: eta[mu[0]][mu[1]],
    lambda a, b, mu: Gbar.d(mu[0]) * G.d(mu[1]),
    0, 0, 2
)
dLvac1.sort()
Lvac1 = dLvac1.Dterm()
dLvac2 = c2 * Rational(1,32) * sum(
    lambda a, b, mu:
        epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
    lambda a, b, mu:
        G.covD(a[0]) * G.covD(a[1]) * Gbar.covDbar(b[0]) * Gbar.covDbar(b[1]),
    2, 2, 0
)
dLvac2.sort()
Lvac2 = dLvac2.Dterm()
printend()


printbegin("Computation of the Lagrangian (combination)...")
L = Lgauge1 + Lgauge2 + Lmatter + Lvac1 + Lvac2
L.sort()
printend()


print("""
L = ...
""")
printbegin("Check : ")
H = (
    half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g.d(mu[0]) * g.d(mu[1]),
        0, 0, 2
    )
    + half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.d(mu[0]) * theta.d(mu[1]),
        0, 0, 2
    )
    + fC * f
)
X = (
    -b0 * quarter * g * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
        0, 0, 4
    )
    - b0 * quarter * theta * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: FF[mu[0]][mu[2]] * FFdual[mu[1]][mu[3]],
        0, 0, 4
    )
    + b0 * half * g * D * D
    + b0 * half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g * phiC.d(mu[0]) * phi.d(mu[1]),
        0, 0, 2
    )
    - b0 * quarter * g.box() * phiC * phi
    - b0 * quarter * I * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.d(mu[0]) * phiC * phi.d(mu[1]),
        0, 0, 2
    )
    + b0 * quarter * I * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.d(mu[0]) * phiC.d(mu[1]) * phi,
        0, 0, 2
    )
    + b0 * half * g * Fc * F
    + b0 * quarter * sqrt(2) * f * phi * Fc
    + b0 * quarter * sqrt(2) * fC * phiC * F
    + c1 * half * g.box() * g.box()
    + c1 * half * theta.box() * theta.box()
    + c1 * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: fC.d(mu[0]) * f.d(mu[1]),
        0, 0, 2
    )
    + c2 * half * H * H
    - c2 * half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]) * g.d(mu[1]) * theta.d(mu[2]) * theta.d(mu[3]),
        0, 0, 4
    )
    + c2 * half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]) * g.d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
        0, 0, 4
    )
    # boundary terms
    - b0 * Rational(1,8) * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            (g * phiC.d(mu[0]) * phi).d(mu[1]),
        0, 0, 2
    )
    - b0 * Rational(1,8) * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            (g * phiC * phi.d(mu[0])).d(mu[1]),
        0, 0, 2
    )
    + b0 * Rational(1,8) * sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        (g.d(mu[0]) * phiC * phi).d(mu[1]),
        0, 0, 2
    )
    - c1 * half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: (g.d(mu[0]) * g.box()).d(mu[1]),
        0, 0, 2
    )
    + c1 * quarter * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: (g.d(mu[0]) * g.d(mu[1]).d(mu[2])).d(mu[3]),
        0, 0, 4
    )
    - c1 * half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: (theta.d(mu[0]) * theta.box()).d(mu[1]),
        0, 0, 2
    )
    + c1 * quarter * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: (theta.d(mu[0]) * theta.d(mu[1]).d(mu[2])).d(mu[3]),
        0, 0, 4
    )
    - L.bosonicpart()
)
X.sort()
printend([ X.iszero() ])


print("""
Field equations
***************
""")

printbegin("Computing the field equation for W...")
eomV = (
    - b0 * half * sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: (G * W[a[1]]).covD(a[0]),
        2, 0, 0
    )
    - b0 * half * sqrt(2) * sum(
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



printbegin("Computing the field equation for Phi...")
eomPhi = -quarter * (b0 * quarter * sqrt(2) * (G + Gbar) * Phi).covD2()
eomPhibar = -quarter * (b0 * quarter * sqrt(2) * (G + Gbar) * Phibar).covDbar2()

eomF = eomPhi.lowest().bosonicpart()
eomFc = eomPhibar.lowest().bosonicpart()

eomphi = eomPhi.Fbarterm().bosonicpart()
eomphiC = eomPhibar.Fterm().bosonicpart()
printend()


printbegin("Computing the field equation for the coupling...")
eomG = (
    b0 * sqrt(2) * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Wbar[b[0]] * Wbar[b[1]],
        0, 2, 0
    )
    - quarter * (b0 * quarter * sqrt(2) * Phi * Phibar).covD2()
    + quarter * c1 * G.box().covD2()
    + quarter * c2 * Rational(1,16) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
        lambda a, b, mu:
            (G.covD(a[0]) * G.covD(a[1]) * Gbar.covDbar(b[1])).covDbar(b[0]),
        2, 2, 0).covD2()
)
eomG.sort()

eomGbar = (
    b0 * sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: W[a[0]] * W[a[1]],
        2, 0, 0
    )
    - quarter * (b0 * quarter * sqrt(2) * Phi * Phibar).covDbar2()
    + quarter * c1 * Gbar.box().covDbar2()
    + quarter * c2 * Rational(1,16) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
        lambda a, b, mu:
            (G.covD(a[1]) * Gbar.covDbar(b[0]) * Gbar.covDbar(b[1])).covD(a[0]),
        2, 2, 0).covDbar2()
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
    ...
""")
printbegin("Check: ")
X = (
    [
        eomD - b0 * g * D
    ] + [
        eomA[nu] - (
            - b0 * sum(
                lambda aa, bb, mu: eta[nu][mu[0]] * eta[mu[1]][mu[2]],
                lambda aa, bb, mu: (g * FF[mu[0]][mu[1]]).d(mu[2]),
                0, 0, 3)
            - b0 * sum(
                lambda aa, bb, mu: eta[nu][mu[0]] * eta[mu[1]][mu[2]],
                lambda aa, bb, mu: theta.d(mu[1]) * FFdual[mu[0]][mu[2]],
                0, 0, 3)
        )
        for nu in range(4)
    ] + [
        eomF - (b0 * half * g * F + b0 * quarter * sqrt(2) * f * phi),
        eomFc - (b0 * half * g * Fc + b0 * quarter * sqrt(2) * fC * phiC),
        eomphi - (
            - b0 * half * g * phi.box()
            - b0 * half * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]],
                lambda a, b, mu: (g + I * theta).d(mu[0]) * phi.d(mu[1]),
                0, 0, 2
            )
            - b0 * quarter * (g + I * theta).box() * phi
            + b0 * quarter * sqrt(2) * fC * F
        ),
        eomphiC - (
            - b0 * half * g * phiC.box()
            - b0 * half * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]],
                lambda a, b, mu: (g - I * theta).d(mu[0]) * phiC.d(mu[1]),
                0, 0, 2
            )
            - b0 * quarter * (g - I * theta).box() * phiC
            + b0 * quarter * sqrt(2) * f * Fc
        ),
        eomf - (
            b0 * quarter * sqrt(2) * phiC * F
            - c1 * f.box()
            + c2 * f * H
        ),
        eomfC - (
            b0 * quarter * sqrt(2) * phi * Fc
            - c1 * fC.box()
            + c2 * fC * H
        ),
        eomg - (
            - b0 * quarter * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu: FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
            0, 0, 4)
            + b0 * half * D * D
            - b0 * quarter * (phiC * phi.box() + phiC.box() * phi)
            + b0 * half * Fc * F
            + c1 * g.box().box()
            - c2 * half * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]],
                lambda a, b, mu:
                    g.box() * g.d(mu[0]) * g.d(mu[1]),
                0, 0, 2
            )
            + c2 * half * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]],
                lambda a, b, mu:
                    g.box() * theta.d(mu[0]) * theta.d(mu[1]),
                0, 0, 2
            )
            - c2 * g.box() * fC * f
            - c2 * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu:
                    g.d(mu[0]).d(mu[2]) * g.d(mu[1]) * g.d(mu[3]),
                0, 0, 4
            )
            - c2 * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu:
                    g.d(mu[0]).d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
                0, 0, 4
            )
            - c2 * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]],
                lambda a, b, mu:
                    theta.box() * g.d(mu[0]) * theta.d(mu[1]),
                0, 0, 2
            )
            - c2 * sum(
                lambda a, b, mu:
                    eta[mu[0]][mu[1]],
                lambda a, b, mu:
                    g.d(mu[0]) * (fC * f).d(mu[1]),
                0, 0, 2
            )
        ),
        eomtheta - (
            - b0 * quarter * sum(
                lambda a, b, mu: eta[mu[0]][mu[2]] * eta[mu[1]][mu[3]],
                lambda a, b, mu: FF[mu[0]][mu[1]] * FFdual[mu[2]][mu[3]],
            0, 0, 4)
            + b0 * quarter * I * (phiC * phi.box() - phiC.box() * phi)
            + c1 * theta.box().box()
            - c2 * half * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]],
                lambda a, b, mu:
                theta.box() * theta.d(mu[0]) * theta.d(mu[1]),
                0, 0, 2
            )
            + c2 * half * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]],
                lambda a, b, mu:
                theta.box() * g.d(mu[0]) * g.d(mu[1]),
                0, 0, 2
            )
            - c2 * theta.box() * fC * f
            - c2 * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu:
                theta.d(mu[0]).d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
                0, 0, 4
            )
            - c2 * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu:
                theta.d(mu[0]).d(mu[2]) * g.d(mu[1]) * g.d(mu[3]),
                0, 0, 4
            )
            - c2 * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]],
                lambda a, b, mu:
                g.box() * g.d(mu[0]) * theta.d(mu[1]),
                0, 0, 2
            )
            - c2 * sum(
                lambda a, b, mu:
                eta[mu[0]][mu[1]],
                lambda a, b, mu:
                theta.d(mu[0]) * (fC * f).d(mu[1]),
                0, 0, 2
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


printbegin("Computing J_ab...")
J = [[
    # gauge part
    b0 * half * sqrt(2) * (G + Gbar) * W[a] * Wbar[b]
    # matter part
    - b0 * I * Rational(1,48) * sqrt(2) * (G + Gbar) * Phibar.dd(a,b) * Phi
    + b0 * I * Rational(1,48) * sqrt(2) * (G + Gbar) * Phibar * Phi.dd(a,b)
    + b0 * Rational(1,96) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi.covD(a)
    - b0 * Rational(1,96) * sqrt(2) * G.covD(a) * Phibar.covDbar(b) * Phi
    + b0 * Rational(1,96) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi.covD(a)
    + b0 * I * Rational(1,48) * sqrt(2) * G.dd(a,b) * Phibar * Phi
    - b0 * I * Rational(1,48) * sqrt(2) * Gbar.dd(a,b) * Phibar * Phi
    # coupling part 1
    - c1 * Rational(1,24) * Gbar.covDbar(b).box() * G.covD(a)
    + c1 * Rational(1,24) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]).dd(a,b) * G.covD(aa[1]),
        2, 2, 0)
    - c1 * Rational(1,6) * I * Gbar.box() * G.dd(a,b)
    + c1 * Rational(1,12) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,b).dd(aa[0],bb[0]) * G.dd(aa[1],bb[1]),
        2, 2, 0)
    + c1 * Rational(1,64) * I * Gbar.covDbar2().dd(a,b) * G.covD2()
    + c1 * Rational(1,384) * Gbar.covDbar2().covD(a) * G.covD2().covDbar(b)
    - c1 * Rational(1,12) * sum(lambda aa, bb, mu:
            eta[mu[0]][mu[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(b).d(mu[0]) * G.covD(a).d(mu[1]),
        0, 0, 2)
    - c1 * Rational(1,12) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]) * G.covD(a).dd(aa[1],b),
        2, 2, 0)
    + c1 * Rational(1,12) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(b).dd(a,bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]),
        2, 2, 0)
    + c1 * Rational(1,6) * I * Gbar.dd(a,b) * G.box()
    - c1 * Rational(1,12) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(aa[0],bb[0]) * G.dd(a,b).dd(aa[1],bb[1]),
        2, 2, 0)
    - c1 * Rational(1,64) * I * Gbar.covDbar2() * G.covD2().dd(a,b)
    - c1 * Rational(1,24) * Gbar.covDbar(b) * G.covD(a).box()
    - c1 * Rational(1,24) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]).dd(a,b),
        2, 2, 0)
    # coupling part 2
    + c2 * Rational(1,64) * I * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a, bb[0]).covDbar(bb[1]) * Gbar.covDbar(b)
            * G.covD(aa[0]) * G.covD(aa[1]),
        2, 2, 0)
    - c2 * Rational(1,128) * I * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]],
        lambda aa, bb, mu:
            Gbar.dd(a, b) * Gbar.covDbar2()
            * G.covD(aa[0]) * G.covD(aa[1]),
        2, 0, 0)
    + c2 * Rational(1,256) * Gbar.covDbar2() * Gbar.covDbar(b) * G.covD2() * G.covD(a)
    - c2 * Rational(1,16) * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,b) * Gbar.covDbar(bb[0])
            * G.dd(aa[0],bb[1]) * G.covD(aa[1]),
        2, 2, 0)
    + c2 * Rational(1,16) * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(aa[0],bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(a,b) * G.covD(aa[1]),
        2, 2, 0)
    - c2 * Rational(1,16) * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(aa[0],b) * G.covD(aa[1]),
        2, 2, 0)
    - c2 * Rational(1,128) * I * sum(
        lambda aa, bb, mu:
            epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(a, b) * G.covD2(),
        0, 2, 0)
    - c2 * Rational(1,64) * I * sum(
        lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(aa[0], b).covD(aa[1]) * G.covD(a),
        2, 2, 0)
    for b in range(2)] for a in range(2)]
printend()

printbegin("Simplification of J_ab...")
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
    + Wbar[b] * eomV
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
    + W[a] * eomV
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
J^mu = i/12 g ( phi^* d^mu phi - d^mu phi^* phi )
       - 1/12 d_mu theta phi^* phi
       - 1/3 d_nu ( d^mu g d^nu theta - d^nu g d^mu theta )
       + i/2 ( d^mu f^* f - f^* d^mu f )
       + fermions
""")
printbegin("Check: ")
X = [
    b0 * Rational(1,12) * I * g * phiC * phi.d(mu)
    - b0 * Rational(1,12) * I * g * phiC.d(mu) * phi
    - b0 * Rational(1,12) * theta.d(mu) * phiC * phi
    - c1 * Rational(1,3) * theta.box() * g.d(mu)
    - c1 * Rational(1,3) * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * g.d(mu).d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(1,3) * theta.d(mu) * g.box()
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(mu).d(rho[0]) * g.d(rho[1]),
        0, 0, 2
    )
    + c1 * half * I * fC.d(mu) * f
    - c1 * half * I * fC * f.d(mu)
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
    - Rational(1,12) * I * Phi.covD2() * eomPhibar
    - Rational(1,24) * I * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: Phi.covD(a[0]) * eomPhibar.covD(a[1]),
        2, 0, 0
    )
    + Rational(1,24) * I * Phi * eomPhibar.covD2()
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

printbegin("Computing T^munu...")
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
printend()



print("""
T^munu =
    ...
""")
printbegin("Check: ")
X = [
    b0 * quarter * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
        0, 0, 4
    )
    - b0 * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g * FF[mu][rho[0]] * FF[nu][rho[1]],
        0, 0, 2
    )
    + b0 * half * eta[mu][nu] * g * D * D
    - b0 * Rational(1,6) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g * phiC.d(rho[0]) * phi.d(rho[1]),
        0, 0, 2
    )
    + b0 * Rational(1,3) * g * phiC.d(mu) * phi.d(nu)
    + b0 * Rational(1,3) * g * phiC.d(nu) * phi.d(mu)
    - b0 * Rational(1,6) * g * phiC * phi.d(mu).d(nu)
    - b0 * Rational(1,6) * g * phiC.d(mu).d(nu) * phi
    - b0 * Rational(1,6) * g.d(mu).d(nu) * phiC * phi
    - b0 * Rational(1,12) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * phiC * phi.d(rho[1]),
        0, 0, 2
    )
    + b0 * Rational(1,12) * g.d(mu) * phiC * phi.d(nu)
    + b0 * Rational(1,12) * g.d(nu) * phiC * phi.d(mu)
    - b0 * Rational(1,12) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * phiC.d(rho[1]) * phi,
        0, 0, 2
    )
    + b0 * Rational(1,12) * g.d(mu) * phiC.d(nu) * phi
    + b0 * Rational(1,12) * g.d(nu) * phiC.d(mu) * phi
    + b0 * Rational(1,12) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC * phi.d(rho[1]),
        0, 0, 2
    )
    - b0 * quarter * I * theta.d(mu) * phiC * phi.d(nu)
    - b0 * quarter * I * theta.d(nu) * phiC * phi.d(mu)
    - b0 * Rational(1,12) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC.d(rho[1]) * phi,
        0, 0, 2
    )
    + b0 * quarter * I * theta.d(mu) * phiC.d(nu) * phi
    + b0 * quarter * I * theta.d(nu) * phiC.d(mu) * phi
    - b0 * Rational(1,6) * eta[mu][nu] * g * Fc * F
    - b0 * Rational(1,12) * sqrt(2) * eta[mu][nu] * f * phi * Fc
    - b0 * Rational(1,12) * sqrt(2) * eta[mu][nu] * fC * phiC * F
    - c1 * half * eta[mu][nu] * g.box() * g.box()
    + c1 * 2 * g.d(mu).d(nu) * g.box()
    + c1 * Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]).d(rho[2]) * g.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    - c1 * Rational(4,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(mu).d(rho[0]) * g.d(nu).d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.box().d(rho[1]),
        0, 0, 2
    )
    - c1 * g.d(mu) * g.box().d(nu)
    - c1 * g.d(nu) * g.box().d(mu)
    + c1 * Rational(2,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(mu).d(nu).d(rho[1]),
        0, 0, 2
    )
    - c1 * half * eta[mu][nu] * theta.box() * theta.box()
    + c1 * 2 * theta.d(mu).d(nu) * theta.box()
    + c1 * Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            theta.d(rho[0]).d(rho[2]) * theta.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    - c1 * Rational(4,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(mu).d(rho[0]) * theta.d(nu).d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.box().d(rho[1]),
        0, 0, 2
    )
    - c1 * theta.d(mu) * theta.box().d(nu)
    - c1 * theta.d(nu) * theta.box().d(mu)
    + c1 * Rational(2,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(mu).d(nu).d(rho[1]),
        0, 0, 2
    )
    - c1 * Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            fC.d(rho[0]) * f.d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(2,3) * fC.d(mu) * f.d(nu)
    + c1 * Rational(2,3) * fC.d(nu) * f.d(mu)
    + c1 * Rational(1,3) * eta[mu][nu] * fC * f.box()
    - c1 * Rational(1,3) * fC * f.d(mu).d(nu)
    + c1 * Rational(1,3) * eta[mu][nu] * fC.box() * f
    - c1 * Rational(1,3) * fC.d(mu).d(nu) * f
    - c2 * Rational(1,8) * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
        0, 0, 4
    )
    + c2 * half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(nu) * g.d(rho[0]) * g.d(rho[1]),
        0, 0, 2
    )
    - c2 * Rational(1,8) * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    + c2 * half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        theta.d(mu) * theta.d(nu) * theta.d(rho[0]) * theta.d(rho[1]),
        0, 0, 2
    )
    + c2 * quarter * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    - c2 * half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(nu) * theta.d(rho[0]) * theta.d(rho[1]),
        0, 0, 2
    )
    - c2 * half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * theta.d(mu) * theta.d(nu),
        0, 0, 2
    )
    - c2 * half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
        0, 0, 4
    )
    + c2 * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(rho[0]) * theta.d(nu) * theta.d(rho[1]),
        0, 0, 2
    )
    + c2 * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(nu) * g.d(rho[0]) * theta.d(mu) * theta.d(rho[1]),
        0, 0, 2
    )
    - c2 * half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + c2 * g.d(mu) * g.d(nu) * fC * f
    - c2 * half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + c2 * theta.d(mu) * theta.d(nu) * fC * f
    - c2 * half * eta[mu][nu] * fC * fC * f * f
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



printbegin("T^mu_mu = 0 (upon field equations): ")
X = (
    sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: T[mu[0]][mu[1]].bosonicpart(),
        0, 0, 2
    )
    - (
        2 * D * eomD
        - Rational(2,3) * Fc * eomF
        - Rational(2,3) * F * eomFc
        + Rational(1,3) * phiC * eomphi
        + Rational(1,3) * phi * eomphiC
        - fC * eomf
        - f * eomfC
    )
)
X.sort()
printend([X.iszero()])


printbegin("d_nu T^munu = 0 (upon field equations): ")
X = [
    sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: T[rho[0]][rho[1]].d(rho[1]).bosonicpart(),
        0, 0, 2
    )
    - (
        D * eomD.d(mu)
        - sum(
            lambda a, b, rho: 1,
            lambda a, b, rho: FF[mu][rho[0]] * eomA[rho[0]],
            0, 0, 1
        )
        + Rational(1,3) * Fc * eomF.d(mu)
        + Rational(1,3) * F * eomFc.d(mu)
        - Rational(2,3) * Fc.d(mu) * eomF
        - Rational(2,3) * F.d(mu) * eomFc
        - Rational(2,3) * phiC.d(mu) * eomphi
        - Rational(2,3) * phi.d(mu) * eomphiC
        + Rational(1,3) * phiC * eomphi.d(mu)
        + Rational(1,3) * phi * eomphiC.d(mu)
        - f.d(mu) * eomfC
        - fC.d(mu) * eomf
        - g.d(mu) * eomg
        - theta.d(mu) * eomtheta
    )
    for mu in range(4)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])



print("""
Static configuration:
T^00 =
    ...
""")
printbegin("Check: ")
X = (
    b0 * quarter * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
        0, 0, 4
    )
    + b0 * half * g * D * D
    + b0 * Rational(1,6) * sum(
        lambda a, b, rho: -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho: g * phiC.d(rho[0]) * phi.d(rho[1]),
        0, 0, 2
    )
    + b0 * Rational(1,12) * sum(
        lambda a, b, rho: -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * (phiC * phi).d(rho[1]),
        0, 0, 2
    )
    - b0 * Rational(1,12) * I * sum(
        lambda a, b, rho: -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC * phi.d(rho[1]),
        0, 0, 2
    )
    + b0 * Rational(1,12) * I * sum(
        lambda a, b, rho: -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC.d(rho[1]) * phi,
        0, 0, 2
    )
    - b0 * Rational(1,6) * g * Fc * F
    - b0 * Rational(1,12) * sqrt(2) * f * phi * Fc
    - b0 * Rational(1,12) * sqrt(2) * fC * phiC * F
    - c1 * half * g.box() * g.box()
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]).d(rho[2]) * g.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.box().d(rho[1]),
        0, 0, 2
    )
    - c1 * half * theta.box() * theta.box()
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            theta.d(rho[0]).d(rho[2]) * theta.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.box().d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            fC.d(rho[0]) * f.d(rho[1]),
        0, 0, 2
    )
    - c1 * Rational(1,3) * (-1) * fC * f.box()
    - c1 * Rational(1,3) * (-1) * fC.box() * f
    - c2 * Rational(1,8) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
        0, 0, 4
    )
    - c2 * Rational(1,8) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    + c2 * quarter * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    - c2 * half * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
        0, 0, 4
    )
    + c2 * half * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + c2 * half * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
        0, 0, 2
    )
    - c2 * half * fC * fC * f * f
    - T[0][0].bosonicpart()
)
X.sort()
printend([X.staticpart().iszero()])


# in this form, all dependence on matter fields has been removed
printbegin("Check: ")
B = [
    b0 * Rational(1, 12) * g * (phiC * phi).d(mu)
    - b0 * Rational(1, 6) * g.d(mu) * phiC * phi
    - b0 * quarter * I * theta * (phiC * phi.d(mu) - phiC.d(mu) * phi)
    + c1 * (-1) * g * g.d(mu).box()
    - c1 * (-1) * g.d(mu) * g.box()
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]] ,
        lambda a, b, rho:
            g.d(rho[0]) * g.d(mu).d(rho[1]),
        0, 0, 2
    )
    + c1 * (-1) * theta * theta.d(mu).box()
    - c1 * (-1) * theta.d(mu) * theta.box()
    + c1 * Rational(1,3) * sum(
        lambda a, b, rho:
        -1 * eta[rho[0]][rho[1]] ,
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(mu).d(rho[1]),
        0, 0, 2
    )
    + c1 * Rational(2,3) * (fC * f).d(mu)
    - c2 * half * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            ((g * g.d(mu) - theta * theta.d(mu))
            * (g.d(rho[0]) * g.d(rho[1]) - theta.d(rho[0]) * theta.d(rho[1]))),
        0, 0, 2
    )
    - c2 * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g * g.d(rho[0]) * theta.d(mu) * theta.d(rho[1]),
        0, 0, 2
    )
    - c2 * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(mu) * g.d(rho[0]) * theta * theta.d(rho[1]),
        0, 0, 2
    )
    + c2 * (g * g.d(mu) + theta * theta.d(mu)) * fC * f
    for mu in range(4)
]
X = (
    c1 * half * g.box() * g.box()
    + c1 * half * theta.box() * theta.box()
    - c1 * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            fC.d(rho[0]) * f.d(rho[1]),
        0, 0, 2
    )
    ######
    + c2 * Rational(3,8) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
        0, 0, 4
    )
    + c2 * Rational(3,8) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    - c2 * Rational(3,4) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    + c2 * Rational(3,2) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
        0, 0, 4
    )
    - c2 * Rational(3,2) * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(rho[1]) * fC * f,
        0, 0, 2
    )
    - c2 * Rational(3,2) * sum(
        lambda a, b, rho:
            -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + c2 * Rational(3,2) * fC * fC * f * f
    ######
    + sum(
        lambda a, b, rho: -1 * eta[rho[0]][rho[1]],
        lambda a, b, rho: B[rho[1]].d(rho[0]),
        0, 0, 2
    )
    ######
    + D * eomD
    - g * eomg
    - theta * eomtheta
    - (fC * eomf + f * eomfC)
    + Rational(1,3) * (phiC * eomphi + phi * eomphiC)
    + Rational(1,3) * (Fc * eomF + F * eomFc)
    - T[0][0].bosonicpart()
)
X.sort()
printend([X.staticpart().iszero()])





# general formulation of the energy
#
# we would like:
# 1/2 < alpha < 1 (e.g. alpha = 1)
# 1/4 < gamma < 1/2
# 2 alpha + 2 gamma - 1 < 0 (this one is impossible given the above two)
#
# alpha, gamma = symbols("alpha gamma")
# printbegin("Check: ")
# X = (
#     (1 - alpha) * b0 * quarter * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
#         0, 0, 4
#     )
#     + (1 - alpha) * b0 * half * g * D * D
#     ######
#     + (alpha - gamma) * b0 * g * Fc * F
#     ######
#     + (2 * alpha - 1) * c1 * half * g.box() * g.box()
#     + (2 * alpha - 1) * c1 * half * theta.box() * theta.box()
#     + (1 - 2 * gamma) * c1 * sum(
#         lambda a, b, rho:
#             -1 * eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             fC.d(rho[0]) * f.d(rho[1]),
#         0, 0, 2
#     )
#     ######
#     + (4 * alpha - 1) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
#         0, 0, 4
#     )
#     + (4 * alpha - 1) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (4 * alpha - 1) * c2 * quarter * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     + (4 * alpha - 1) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (2 * alpha + 2 * gamma - 1) * c2 * half * sum(
#         lambda a, b, rho:
#             -1 * eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     - (2 * alpha + 2 * gamma - 1) * c2 * half * sum(
#         lambda a, b, rho:
#             -1 * eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     + (4 * gamma - 1) * c2 * half * fC * fC * f * f
#     ######
#     - b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g * phiC * phi.d(rho[1])).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g * phiC.d(rho[1]) * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - (1 - 3 * alpha) * b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g.d(rho[1]) * phiC * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + alpha * b0 * quarter * I * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (theta * phiC * phi.d(rho[1])).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * b0 * quarter * I * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (theta * phiC.d(rho[1]) * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g * g.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g.d(rho[1]) * g.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (g.d(rho[2]) * g.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta * theta.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta.d(rho[1]) * theta.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (theta.d(rho[2]) * theta.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + (1 - 3 * gamma) * c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (fC * f).d(rho[1]).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * g.d(mu[2]) * g.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[2]) * theta.d(mu[3]) * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * theta.d(mu[2]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[1]) * g.d(mu[2]) * theta * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[2]) * g.d(mu[3]) * theta * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#         eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     ######
#     + alpha * D * eomD
#     - alpha * g * eomg
#     - alpha * theta * eomtheta
#     - gamma * (fC * eomf + f * eomfC)
#     - (1 - 3 * alpha) * Rational(1,6) * (phiC * eomphi + phi * eomphiC)
#     - (1 + 3 * alpha - 6 * gamma) * Rational(1,6) * (Fc * eomF + F * eomFc)
#     - T[0][0].bosonicpart()
# )
# X.sort()
# printend([X.staticpart().iszero()])



#alpha, beta, gamma = symbols("alpha beta gamma")
# printbegin("Check: ")
# X = (
#     (1 - alpha) * b0 * quarter * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
#         0, 0, 4
#     )
#     + (1 + alpha) * b0 * half * g * D * D
#     ######
#     + (alpha - beta) * b0 * half * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: g * phiC.d(rho[0]) * phi.d(rho[1]),
#         0, 0, 2
#     )
#     - (alpha - beta) * b0 * quarter * g.box() * phiC * phi
#     + (alpha + beta - 2 * gamma) * b0 * half * g * Fc * F
#     ######
#     - (1 - 2 * alpha) * c1 * half * g.box() * g.box()
#     - (1 - 2 * beta) * c1 * half * theta.box() * theta.box()
#     - (1 - 2 * gamma) * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             fC.d(rho[0]) * f.d(rho[1]),
#         0, 0, 2
#     )
#     ######
#     - (1 - 4 * alpha) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 4 * beta) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     + (1 - 2 * alpha - 2 * beta) * c2 * quarter * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 2 * alpha - 2 * beta) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 2 * alpha - 2 * gamma) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     - (1 - 2 * beta - 2 * gamma) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     - (1 - 4 * gamma) * c2 * half * fC * fC * f * f
#     ######
#     - (1 + 3 * alpha - 3 * beta) * b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g * phiC * phi.d(rho[1])).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - (1 + 3 * alpha - 3 * beta) * b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g * phiC.d(rho[1]) * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - (1 - 3 * alpha) * b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (g.d(rho[1]) * phiC * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + beta * b0 * quarter * I * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (theta * phiC * phi.d(rho[1])).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - beta * b0 * quarter * I * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: (theta * phiC.d(rho[1]) * phi).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g * g.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g.d(rho[1]) * g.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (g.d(rho[2]) * g.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + beta * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta * theta.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - beta * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta.d(rho[1]) * theta.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (theta.d(rho[2]) * theta.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + (1 - 3 * gamma) * c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (fC * f).d(rho[1]).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * g.d(mu[2]) * g.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[2]) * theta.d(mu[3]) * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * theta.d(mu[2]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     - beta * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - beta * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[1]) * g.d(mu[2]) * theta * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + beta * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[2]) * g.d(mu[3]) * theta * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - beta * c2 * sum(
#         lambda a, b, mu:
#         eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     ######
#     - alpha * g * eomg
#     - beta * theta * eomtheta
#     - gamma * (fC * eomf + f * eomfC)
#     - (1 - 3 * beta) * Rational(1,6) * (phiC * eomphi + phi * eomphiC)
#     - (1 + 3 * beta - 6 * gamma) * Rational(1,6) * (Fc * eomF + F * eomFc)
#     - T[0][0].bosonicpart()
# )
# X.sort()
# printend([X.staticpart().iszero()])


# Most general use of field equations in energy:
# #
# alpha, beta, gamma, delta = symbols("alpha beta gamma delta")
# printbegin("Check: ")
# X = (
#     (1 - alpha) * b0 * quarter * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho: g * FF[rho[0]][rho[2]] * FF[rho[1]][rho[3]],
#         0, 0, 4
#     )
#     + (1 + alpha) * b0 * half * g * D * D
#     ######
#     - b0 * Rational(1,6) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: g * phiC.d(rho[0]) * phi.d(rho[1]),
#         0, 0, 2
#     )
#     - (alpha + 2 * delta) * b0 * quarter * g * (phiC * phi.box() + phiC.box() * phi)
#     - (1 + 6 * delta) * b0 * Rational(1,12) * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: g.d(rho[0]) * (phiC * phi).d(rho[1]),
#         0, 0, 2
#     )
#     - 6 * delta * b0 * Rational(1,12) * g.box() * phiC * phi
#     + (1 - 6 * delta) * b0 * Rational(1,12) * I * sum(
#         lambda a, b, rho: eta[rho[0]][rho[1]],
#         lambda a, b, rho: theta.d(rho[0]) * (phiC * phi.d(rho[1]) - phiC.d(rho[1]) * phi),
#         0, 0, 2
#     )
#     + beta * b0 * quarter * I * theta * (phiC * phi.box() - phiC.box() * phi)
#     + (1 + 3 * alpha - 6 * gamma - 6 * delta) * b0 * Rational(1,3) * half * g * Fc * F
#     ######
#     - (1 - 2 * alpha) * c1 * half * g.box() * g.box()
#     - (1 - 2 * beta) * c1 * half * theta.box() * theta.box()
#     - (1 - 2 * gamma) * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             fC.d(rho[0]) * f.d(rho[1]),
#         0, 0, 2
#     )
#     ######
#     - (1 - 4 * alpha) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 4 * beta) * c2 * Rational(1,8) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     + (1 - 2 * alpha - 2 * beta) * c2 * quarter * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 2 * alpha - 2 * beta) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
#         0, 0, 4
#     )
#     - (1 - 2 * alpha - 2 * gamma) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             g.d(rho[0]) * g.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     - (1 - 2 * beta - 2 * gamma) * c2 * half * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
#         0, 0, 2
#     )
#     - (1 - 4 * gamma) * c2 * half * fC * fC * f * f
#     ######
#     + alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g * g.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (g.d(rho[1]) * g.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (g.d(rho[2]) * g.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + beta * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta * theta.d(rho[1]).box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - beta * c1 * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (theta.d(rho[1]) * theta.box()).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     + c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
#         lambda a, b, rho:
#             (theta.d(rho[2]) * theta.d(rho[1]).d(rho[3])).d(rho[0]), # boundary term
#         0, 0, 4
#     )
#     + (1 - 3 * gamma) * c1 * Rational(1,3) * sum(
#         lambda a, b, rho:
#             eta[rho[0]][rho[1]],
#         lambda a, b, rho:
#             (fC * f).d(rho[1]).d(rho[0]), # boundary term
#         0, 0, 2
#     )
#     - alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * g.d(mu[2]) * g.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[2]) * theta.d(mu[3]) * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + alpha * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * theta.d(mu[2]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - alpha * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (g * g.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     - beta * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - beta * c2 * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[1]) * g.d(mu[2]) * theta * theta.d(mu[3])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     + beta * c2 * half * sum(
#         lambda a, b, mu:
#             eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
#         lambda a, b, mu:
#             (g.d(mu[2]) * g.d(mu[3]) * theta * theta.d(mu[1])).d(mu[0]), # boundary term
#         0, 0, 4
#     )
#     - beta * c2 * sum(
#         lambda a, b, mu:
#         eta[mu[0]][mu[1]],
#         lambda a, b, mu:
#             (theta * theta.d(mu[1]) * fC * f).d(mu[0]), # boundary term
#         0, 0, 2
#     )
#     ######
#     - alpha * g * eomg
#     - beta * theta * eomtheta
#     - gamma * (fC * eomf + f * eomfC)
#     - delta * (phiC * eomphi + phi * eomphiC)
#     - (1 - 3 * gamma - 3 * delta) * Rational(1,3) * (Fc * eomF + F * eomFc)
#     - T[0][0].bosonicpart()
# )
# X.sort()
# printend([X.staticpart().iszero()])

