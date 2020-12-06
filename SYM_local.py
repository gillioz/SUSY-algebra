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

printbegin("Computation of the Lagrangian (gauge part)...")
dLgauge1 = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[a[1]][a[0]],
    lambda a, b, mu: G * W[a[0]] * W[a[1]],
    2, 0, 0
)
Lgauge1 = dLgauge1.Fterm()
dLgauge2 = 2 * sqrt(2) * sum(
    lambda a, b, mu: epsilon[b[0]][b[1]],
    lambda a, b, mu: Gbar * Wbar[b[0]] * Wbar[b[1]],
    0, 2, 0
)
Lgauge2 = dLgauge2.Fbarterm()
printend()


printbegin("Computation of the Lagrangian (matter part)...")
dLmatter = half * sqrt(2) * (G + Gbar) \
           * (Phibar * Phi + exp2V * Phibar * Phi)
dLmatter.sort()
Lmatter = dLmatter.Dterm()
printend()


printbegin("Computation of the Lagrangian (combination)...")
L = Lgauge1 + Lgauge2 + Lmatter
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
    + ...
""")
printbegin("Check : ")
X = (
    # # terms in g
    -half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: g * FF[mu[0]][mu[2]] * FF[mu[1]][mu[3]],
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
    + g * D * D
    + sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: g * DphiC(mu[0]) * Dphi(mu[1]),
        0, 0, 2
    )
    - half * g.box() * phiC * phi
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * chiC[b[0]] * Dchi(a[0], mu[0]),
        1, 1, 1
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            g * DchiC(b[0], mu[0]) * chi[a[0]],
        1, 1, 1
    )
    + e * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            g * phiC * chi[a[0]] * lam[a[1]],
        2, 0, 0
    )
    + e * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            g * phi * chiC[b[0]] * lamC[b[1]],
        0, 2, 0
    )
    + half * e * g * D * phiC * phi
    + g * Fc * F
    # terms in theta
    -half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: theta * FF[mu[0]][mu[2]] * FFdual[mu[1]][mu[3]],
        0, 0, 4
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            theta.d(mu[0]) * lamC[b[0]] * lam[a[0]],
        1, 1, 1
    )
    - half * I * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: theta.d(mu[0]) * (phiC * Dphi(mu[1]) - DphiC(mu[1]) * phi),
        0, 0, 2
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            theta.d(mu[0]) * chiC[b[0]] * chi[a[0]],
        1, 1, 1
    )
    # terms in xi
    - half * sqrt(2) * I * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xiC[b[0]] * lamC[b[1]] * FF[mu[0]][mu[1]],
        2, 2, 2
    )
    + half * sqrt(2) * I * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]] * sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xi[a[0]] * lam[a[1]] * FF[mu[0]][mu[1]],
        2, 2, 2
    )
    + sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * lam[a[1]] * D,
        2, 0, 0
    )
    + sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * lamC[b[1]] * D,
        0, 2, 0
    )
    + half * sqrt(2) * e * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * lam[a[1]] * phiC * phi,
        2, 0, 0
    )
    + half * sqrt(2) * e * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * lamC[b[1]] * phiC * phi,
        0, 2, 0
    )
    - sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * DchiC(b[0], mu[0]) * phi,
        1, 1, 1
    )
    - sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * Dchi(a[0], mu[0]) * phiC,
        1, 1, 1
    )
    - sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * chi[a[1]] * Fc,
        2, 0, 0
    )
    - sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * chiC[b[1]] * F,
        0, 2, 0
    )
    # terms in f
    + half * sqrt(2) * f * phi * Fc
    + half * sqrt(2) * fC * phiC * F
    + half * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            f * lam[a[0]] * lam[a[1]],
        2, 0, 0
    )
    + half * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            fC * lamC[b[0]] * lamC[b[1]],
        0, 2, 0
    )
    # # boundary terms
    - quarter * (g * phiC * phi).box()
    + half * sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        (g.d(mu[1]) * phiC * phi).d(mu[0]),
        0, 0, 2
    )
    - sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (theta * lamC[b[0]] * lam[a[0]]).d(mu[0]),
        1, 1, 1
    )
    + half * sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (xi[a[0]] * chiC[b[0]] * phi).d(mu[0]),
        1, 1, 1
    )
    + half * sqrt(2) * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            (xiC[b[0]] * chi[a[0]] * phiC).d(mu[0]),
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

printbegin("Computing the field equation for W...")
eomV = (
    -sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: (G * W[a[1]]).covD(a[0]),
        2, 0, 0
    )
    -sqrt(2) * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: (Gbar * Wbar[b[1]]).covDbar(b[0]),
        0, 2, 0
    )
    + quarter * sqrt(2) * e * (G + Gbar) * (Phibar * Phi + exp2V * Phibar * Phi)
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
eomPhi = -quarter * ( half * sqrt(2) * (G + Gbar)
    * (Phi + exp2V * Phi) ).covD2()
eomPhibar = -quarter * ( half * sqrt(2) * (G + Gbar)
    * (Phibar + exp2V * Phibar) ).covDbar2()

eomF = eomPhi.lowest()
eomFc = eomPhibar.lowest()

eomphi = eomPhi.Fbarterm()
eomphiC = eomPhibar.Fterm()
printend()



printbegin("Computing the field equation for the coupling...")
eomG = (
    2 * sqrt(2) * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: Wbar[b[0]] * Wbar[b[1]],
        0, 2, 0
    )
    - quarter * (half * sqrt(2) * Phi * Phibar
                 + half * sqrt(2) * Phi * Phibar * exp2V).covD2()
)
eomG.sort()
eomGbar = (
    2 * sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: W[a[0]] * W[a[1]],
        2, 0, 0
    )
    - quarter * (half * sqrt(2) * Phi * Phibar
                 + half * sqrt(2) * Phi * Phibar * exp2V).covDbar2()
)
eomGbar.sort()

eomf = eomG.lowest()
eomfC = eomGbar.lowest()

eomg1 = eomG.Fbarterm()
eomg2 = eomGbar.Fterm()

eomg = half * sqrt(2) * (eomg1 + eomg2)
eomg.sort()
eomtheta = -half * sqrt(2) * I * (eomg1 - eomg2)
eomtheta.sort()

printend()


print("""
Field equations in components:
    ...
""")
printbegin("Check: ")
X = (
    [eomD - (
        2 * g * D
        + half * e * g * phiC * phi
        + sqrt(2) * sum(
            lambda a, b, mu:
                epsilon[a[1]][a[0]],
            lambda a, b, mu:
                xi[a[0]] * lam[a[1]],
            2, 0, 0
        )
        + sqrt(2) * sum(
            lambda a, b, mu:
                epsilon[b[0]][b[1]],
            lambda a, b, mu:
                xiC[b[0]] * lamC[b[1]],
            0, 2, 0
        )
    )]
    + [eomA[mu] - (
        - 2 * sum(
            lambda aa, bb, rho: eta[mu][rho[0]] * eta[rho[1]][rho[2]],
            lambda aa, bb, rho: (g * FF[rho[0]][rho[1]]).d(rho[2]),
            0, 0, 3
        )
        - 2 * sum(
            lambda aa, bb, rho: eta[mu][rho[0]] * eta[rho[1]][rho[2]],
            lambda aa, bb, rho: theta.d(rho[2]) * FFdual[rho[0]][rho[1]],
            0, 0, 3
        )
        + e * half * I * sum(
            lambda aa, bb, rho: eta[mu][rho[0]],
            lambda aa, bb, rho: g * (phiC * Dphi(rho[0]) - DphiC(rho[0]) * phi),
            0, 0, 1
        )
        - half * e * sum(
            lambda aa, bb, rho: eta[mu][rho[0]],
            lambda aa, bb, rho: theta.d(rho[0]) * phiC * phi,
            0, 0, 1
        )
        - e * sum(
            lambda a, b, rho:
                sigmabar[mu][b[0]][a[0]],
            lambda a, b, rho:
                g * chiC[b[0]] * chi[a[0]],
            1, 1, 0
        )
        + half * sqrt(2) * e * sum(
            lambda a, b, rho:
                sigmabar[mu][b[0]][a[0]],
            lambda a, b, rho:
                xi[a[0]] * chiC[b[0]] * phi,
            1, 1, 0
        )
        - half * sqrt(2) * e * sum(
            lambda a, b, rho:
                sigmabar[mu][b[0]][a[0]],
            lambda a, b, rho:
                xiC[b[0]] * chi[a[0]] * phiC,
            1, 1, 0
        )
        + half * sqrt(2) * I * sum(
            lambda a, b, rho:
                epsilon[b[0]][b[1]] * (sigmabar[mu][b[0]][a[0]] * sigmabar[rho[0]][b[1]][a[1]]
                 - sigmabar[rho[0]][b[0]][a[0]] * sigmabar[mu][b[1]][a[1]]),
            lambda a, b, rho:
                (xi[a[0]] * lam[a[1]]).d(rho[0]),
            2, 2, 1
        )
        - half * sqrt(2) * I * sum(
            lambda a, b, rho:
                epsilon[a[1]][a[0]] * (sigmabar[mu][b[0]][a[0]] * sigmabar[rho[0]][b[1]][a[1]]
                 - sigmabar[rho[0]][b[0]][a[0]] * sigmabar[mu][b[1]][a[1]]),
            lambda a, b, rho:
                (xiC[b[0]] * lamC[b[1]]).d(rho[0]),
            2, 2, 1
        )
    )
    for mu in range(4)]
    + [
        eomF - (
            g * F
            + half * sqrt(2) * f * phi
            - sqrt(2) * sum(
                lambda a, b, mu:
                    epsilon[a[1]][a[0]],
                lambda a, b, mu:
                    xi[a[0]] * chi[a[1]],
                2, 0, 0
            )
        ),
        eomFc - (
            g * Fc
            + half * sqrt(2) * fC * phiC
            - sqrt(2) * sum(
                lambda a, b, mu:
                    epsilon[b[0]][b[1]],
                lambda a, b, mu:
                    xiC[b[0]] * chiC[b[1]],
                0, 2, 0
            )
        ),
        eomphi - (
            half * e * g * D * phi
            - sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: g * D2phi(rho[0], rho[1]),
                0, 0, 2
            )
            - sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: g.d(rho[0]) * Dphi(rho[1]),
                0, 0, 2
            )
            - half * g.box() * phi
            + e * sum(
                lambda a, b, mu: epsilon[a[1]][a[0]],
                lambda a, b, mu: g * chi[a[0]] * lam[a[1]],
                2, 0, 0
            )
            - I * sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: theta.d(rho[0]) * Dphi(rho[1]),
                0, 0, 2
            )
            - half * I * theta.box() * phi
            + half * sqrt(2) * e * sum(
                lambda a, b, rho: epsilon[a[1]][a[0]],
                lambda a, b, rho: xi[a[0]] * lam[a[1]] * phi,
                2, 0, 0
            )
            + half * sqrt(2) * e * sum(
                lambda a, b, rho: epsilon[b[0]][b[1]],
                lambda a, b, rho: xiC[b[0]] * lamC[b[1]] * phi,
                0, 2, 0
            )
            - sqrt(2) * I * sum(
                lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
                lambda a, b, rho: xiC[b[0]] * Dchi(a[0], rho[0]),
                1, 1, 1
            )
            + half * sqrt(2) * fC * F
        ),
        eomphiC - (
            half * e * g * D * phiC
            - sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: g * D2phiC(rho[0], rho[1]),
                0, 0, 2
            )
            - sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: g.d(rho[0]) * DphiC(rho[1]),
                0, 0, 2
            )
            - half * g.box() * phiC
            + e * sum(
                lambda a, b, mu: epsilon[b[0]][b[1]],
                lambda a, b, mu: g * chiC[b[0]] * lamC[b[1]],
                0, 2, 0
            )
            + I * sum(
                lambda a, b, rho: eta[rho[0]][rho[1]],
                lambda a, b, rho: theta.d(rho[0]) * DphiC(rho[1]),
                0, 0, 2
            )
            + half * I * theta.box() * phiC
            + half * sqrt(2) * e * sum(
                lambda a, b, rho: epsilon[a[1]][a[0]],
                lambda a, b, rho: xi[a[0]] * lam[a[1]] * phiC,
                2, 0, 0
            )
            + half * sqrt(2) * e * sum(
                lambda a, b, rho: epsilon[b[0]][b[1]],
                lambda a, b, rho: xiC[b[0]] * lamC[b[1]] * phiC,
                0, 2, 0
            )
            - sqrt(2) * I * sum(
                lambda a, b, rho: sigmabar[rho[0]][b[0]][a[0]],
                lambda a, b, rho: xi[a[0]] * DchiC(b[0], rho[0]),
                1, 1, 1
            )
            + half * sqrt(2) * f * Fc
        ),
        eomf - (
            half * sqrt(2) * phiC * F
            + half * sqrt(2) * sum(
                lambda a, b, mu: epsilon[b[0]][b[1]],
                lambda a, b, mu: lamC[b[0]] * lamC[b[1]],
                0, 2, 0
            )
        ),
        eomfC - (
            half * sqrt(2) * phi * Fc
            + half * sqrt(2) * sum(
                lambda a, b, mu: epsilon[a[1]][a[0]],
                lambda a, b, mu: lam[a[0]] * lam[a[1]],
                2, 0, 0
            )
        ),
        eomg - (
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
            - half * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]],
                lambda a, b, mu: D2phiC(mu[0], mu[1]) * phi,
                0, 0, 2
            )
            - half * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]],
                lambda a, b, mu: phiC * D2phi(mu[0], mu[1]),
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
        ),
        eomtheta - (
            -half * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
                lambda a, b, mu: FF[mu[0]][mu[2]] * FFdual[mu[1]][mu[3]],
                0, 0, 4
            )
            - sum(
                lambda a, b, mu: sigmabar[mu[0]][b[0]][a[0]],
                lambda a, b, mu: (lamC[b[0]] * lam[a[0]]).d(mu[0]),
                1, 1, 1
            )
            + half * I * sum(
                lambda a, b, mu: eta[mu[0]][mu[1]],
                lambda a, b, mu: phiC * D2phi(mu[0],mu[1]) - D2phiC(mu[0], mu[1]) * phi,
                0, 0, 2
            )
            - sum(
                lambda a, b, mu: sigmabar[mu[0]][b[0]][a[0]],
                lambda a, b, mu: (chiC[b[0]] * chi[a[0]]).d(mu[0]),
                1, 1, 1
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

# printbegin("Computing the conserved J_ab...")
# J = [[
#     sqrt(2) * (G + Gbar) * W[a] * Wbar[b]
#     - I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar.dd(a,b) * Phi
#     - I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar.dd(a,b) * Phi * exp2V
#     + I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi.dd(a,b)
#     + I * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi.dd(a,b) * exp2V
#     + Rational(1,48) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi.covD(a)
#     + Rational(1,48) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi.covD(a) * exp2V
#     - e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi.covD(a) * V.covDbar(b)
#     - e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi.covD(a) * V.covDbar(b) * exp2V
#     + e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi * V.covD(a)
#     + e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar.covDbar(b) * Phi * V.covD(a) * exp2V
#     + e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covDbar(b).covD(a)
#     + e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covDbar(b).covD(a) * exp2V
#     - e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covD(a).covDbar(b)
#     - e * Rational(1,24) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covD(a).covDbar(b) * exp2V
#     + e**2 * Rational(1,12) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covDbar(b) * V.covD(a)
#     + e**2 * Rational(1,12) * sqrt(2) * (G + Gbar) * Phibar * Phi * V.covDbar(b) * V.covD(a) * exp2V
#     - Rational(1,48) * sqrt(2) * G.covD(a) * Phibar.covDbar(b) * Phi
#     - Rational(1,48) * sqrt(2) * G.covD(a) * Phibar.covDbar(b) * Phi * exp2V
#     - e * Rational(1,24) * sqrt(2) * G.covD(a) * Phibar * Phi * V.covDbar(b)
#     - e * Rational(1,24) * sqrt(2) * G.covD(a) * Phibar * Phi * V.covDbar(b) * exp2V
#     + Rational(1,48) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi.covD(a)
#     + Rational(1,48) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi.covD(a) * exp2V
#     + e * Rational(1,24) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi * V.covD(a)
#     + e * Rational(1,24) * sqrt(2) * Gbar.covDbar(b) * Phibar * Phi * V.covD(a) * exp2V
#     + I * Rational(1,24) * sqrt(2) * G.dd(a,b) * Phibar * Phi
#     + I * Rational(1,24) * sqrt(2) * G.dd(a,b) * Phibar * Phi * exp2V
#     - I * Rational(1,24) * sqrt(2) * Gbar.dd(a,b) * Phibar * Phi
#     - I * Rational(1,24) * sqrt(2) * Gbar.dd(a,b) * Phibar * Phi * exp2V
#     for b in range(2)]
#     for a in range(2)]
# for Jrow in J:
#     for j in Jrow:
#         j.sort()
# printend()


printbegin("Computing the conserved J_ab...")
PhibarPhi = Phibar * Phi + Phibar * Phi * exp2V
PhibarPhi.sort()
Phibarexp2V = Phibar + Phibar * exp2V
Phibarexp2V.sort()
Phiexp2V = Phi + Phi * exp2V
Phiexp2V.sort()
exp2Vinv = -2 * e * V + 2 * e**2 * V * V
exp2Vinv.sort()
J = [[
    sqrt(2) * (G + Gbar) * W[a] * Wbar[b]
    + Rational(1,48) * sqrt(2) * (G + Gbar) * PhibarPhi.covDbar(b).covD(a)
    - Rational(1,48) * sqrt(2) * (G + Gbar) * PhibarPhi.covD(a).covDbar(b)
    + Rational(1,16) * sqrt(2) * (G + Gbar) * Phibarexp2V.covDbar(b) * Phiexp2V.covD(a)
    + Rational(1,16) * sqrt(2) * (G + Gbar) * Phibarexp2V.covDbar(b) * Phiexp2V.covD(a) * exp2Vinv
    - Rational(1,48) * sqrt(2) * G.covD(a) * PhibarPhi.covDbar(b)
    + Rational(1,48) * sqrt(2) * Gbar.covDbar(b) * PhibarPhi.covD(a)
    - Rational(1,48) * sqrt(2) * G.covD(a).covDbar(b) * Phibar * Phi
    - Rational(1,48) * sqrt(2) * G.covD(a).covDbar(b) * Phibar * Phi * exp2V
    + Rational(1,48) * sqrt(2) * Gbar.covDbar(b).covD(a) * Phibar * Phi
    + Rational(1,48) * sqrt(2) * Gbar.covDbar(b).covD(a) * Phibar * Phi * exp2V
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
    # + Rational(1,12) * Phibar * eomPhi.covDbar(b)
    # - Rational(1,6) * Phibar.covDbar(b) * eomPhi
    # - e * half * V.covDbar(b) * Phibar * eomPhi
    + Rational(1,12) * (Phibar * eomPhi).covDbar(b)
    - quarter * Phibarexp2V.covDbar(b) * eomPhi
    - quarter * Phibarexp2V.covDbar(b) * eomPhi * exp2Vinv
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
    # + Rational(1,12) * Phi * eomPhibar.covD(a)
    # - Rational(1,6) * Phi.covD(a) * eomPhibar
    # - e * half * V.covD(a) * Phi * eomPhibar
    + Rational(1,12) * (Phi * eomPhibar).covD(a)
    - quarter * Phiexp2V.covD(a) * eomPhibar
    - quarter * Phiexp2V.covD(a) * eomPhibar * exp2Vinv
    - quarter * G.covD(a) * eomGbar
    for a in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])



X = (
    sqrt(2) * sum(
        lambda a, b, mu: epsilon[a[1]][a[0]],
        lambda a, b, mu: (G * W[a[1]] * Wbar[0]).covD(a[0]),
        2, 0, 0
    )
    + sqrt(2) * Wbar[0] * sum(
        lambda a, b, mu: epsilon[b[0]][b[1]],
        lambda a, b, mu: (Gbar * Wbar[b[1]]).covDbar(b[0]),
        0, 2, 0
    )
    + quarter * Gbar.covDbar(0) * (
        2 * sqrt(2) * sum(
            lambda a, b, mu: epsilon[b[0]][b[1]],
            lambda a, b, mu: Wbar[b[0]] * Wbar[b[1]],
            0, 2, 0
        )
    )
    - quarter * sqrt(2) * e * Wbar[0] * (G + Gbar) * (Phibar * Phi + exp2V * Phibar * Phi)
    - Rational(1,32) * sqrt(2) * Gbar.covDbar(0) * (
        Phi * Phibar + Phi * Phibar * exp2V).covD2()
    + Wbar[0] * eomV
    - quarter * Gbar.covDbar(0) * eomG
)
X.sort()
print(X.iszero())


printbegin("Computing J^mu...")
Jmu = [sum(lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
           lambda a, b, rho: J[a[0]][b[0]],
           1, 1, 0)
       for mu in range(4)]
for j in Jmu:
    j.sort()
printend()

print("""
J^mu = ...
""")
printbegin("Check: ")
X = [
    Rational(1,6) * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: g * phiC * Dphi(rho[0]),
        0, 0, 1
    )
    - Rational(1,6) * I * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: g * DphiC(rho[0]) * phi,
        0, 0, 1
    )
    - half * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: g * lamC[b[0]] * lam[a[0]],
        1, 1, 0
    )
    + Rational(1,6) * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: g * chiC[b[0]] * chi[a[0]],
        1, 1, 0
    )
    - Rational(1,6) * sum(
        lambda a, b, rho: eta[mu][rho[0]],
        lambda a, b, rho: theta.d(rho[0]) * phiC * phi,
        0, 0, 1
    )
    + Rational(1,12) * sqrt(2) * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: xiC[b[0]] * chi[a[0]] * phiC,
        1, 1, 0
    )
    - Rational(1,12) * sqrt(2) * sum(
        lambda a, b, rho: sigmabar[mu][b[0]][a[0]],
        lambda a, b, rho: xi[a[0]] * chiC[b[0]] * phi,
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
T^munu =
    ...
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
    + Rational(1,3) * e * eta[mu][nu] * g * D * phiC * phi
    - Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g * DphiC(rho[0]) * Dphi(rho[1]),
        0, 0, 2
    )
    + Rational(2,3) * g * DphiC(mu) * Dphi(nu)
    + Rational(2,3) * g * DphiC(nu) * Dphi(mu)
    - Rational(1,6) * g * phiC * D2phi(mu, nu)
    - Rational(1,6) * g * phiC * D2phi(nu, mu)
    - Rational(1,6) * g * D2phiC(mu, nu) * phi
    - Rational(1,6) * g * D2phiC(nu, mu) * phi
    - Rational(1,3) * g.d(mu).d(nu) * phiC * phi
    - Rational(1,6) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * phiC * Dphi(rho[1]),
        0, 0, 2
    )
    + Rational(1,6) * g.d(mu) * phiC * Dphi(nu)
    + Rational(1,6) * g.d(nu) * phiC * Dphi(mu)
    - Rational(1,6) * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: g.d(rho[0]) * DphiC(rho[1]) * phi,
        0, 0, 2
    )
    + Rational(1,6) * g.d(mu) * DphiC(nu) * phi
    + Rational(1,6) * g.d(nu) * DphiC(mu) * phi
    + Rational(1,6) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * phiC * Dphi(rho[1]),
        0, 0, 2
    )
    - half * I * theta.d(mu) * phiC * Dphi(nu)
    - half * I * theta.d(nu) * phiC * Dphi(mu)
    - Rational(1,6) * I * eta[mu][nu] * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * DphiC(rho[1]) * phi,
        0, 0, 2
    )
    + half * I * theta.d(mu) * DphiC(nu) * phi
    + half * I * theta.d(nu) * DphiC(mu) * phi
    - Rational(1,3) * eta[mu][nu] * g * Fc * F
    - Rational(1,6) * sqrt(2) * eta[mu][nu] * f * phi * Fc
    - Rational(1,6) * sqrt(2) * eta[mu][nu] * fC * phiC * F
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
