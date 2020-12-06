#!/usr/bin/env python3.5

from SUSYalgebra import *
#from sympy import solve


print("""
Fields
******
""")

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

dL = Rational(1,32) * sum(
    lambda a, b, mu:
        epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
    lambda a, b, mu:
        G.covD(a[0]) * G.covD(a[1]) * Gbar.covDbar(b[0]) * Gbar.covDbar(b[1]),
    2, 2, 0
)
dL.sort()

L = dL.Dterm()


print("""
L = 1/2 ( 1/2 d_mu g d^mu g + 1/2 d_mu theta d^mu theta + f^* f )^2
    - 1/2 d_mu g d^mu g d_nu theta d^nu theta
    + 1/2 d_mu g d_nu g d^mu theta d^nu theta
    + fermions
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
K = [[
    half * g.d(mu) * g.d(nu)
    + half * theta.d(mu) * theta.d(nu)
    + eta[mu][nu] * fC * f
    for nu in range(4)]
    for mu in range(4)
]
R = [[
    I * (g.d(mu) * theta.d(nu) - theta.d(mu) * g.d(nu))
    for nu in range(4)]
    for mu in range(4)
]
X = (
    half * H * H
    - half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]) * g.d(mu[1]) * theta.d(mu[2]) * theta.d(mu[3]),
        0, 0, 4
    )
    + half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]) * g.d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
        0, 0, 4
    )
    - half * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]].d(mu[0]) * H,
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[0]) * xi[a[0]] * H,
        1, 1, 1
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]].d(mu[1]) * K[mu[0]][mu[2]],
        1, 1, 3
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[1]) * xi[a[0]] * K[mu[0]][mu[2]],
        1, 1, 3
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]].d(mu[1]) * R[mu[0]][mu[2]],
        1, 1, 3
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[1]) * xi[a[0]] * R[mu[0]][mu[2]],
        1, 1, 3
    )
    - half * I * sum(
        lambda a, b, mu:
            LeviCivita(mu[0:4]) * eta[mu[3]][mu[4]] * sigmabar[mu[4]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]].d(mu[0]) * g.d(mu[1]) * theta.d(mu[2]),
        1, 1, 5
    )
    + half * I * sum(
        lambda a, b, mu:
            LeviCivita(mu[0:4]) * eta[mu[3]][mu[4]] * sigmabar[mu[4]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[0]) * xi[a[0]] * g.d(mu[1]) * theta.d(mu[2]),
        1, 1, 5
    )
    + half * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]] * g.box() * theta.d(mu[0]),
        1, 1, 1
    )
    - half * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]] * g.d(mu[0]) * theta.box(),
        1, 1, 1
    )
    - quarter * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]],
        lambda a, b, mu:
            xi[a[0]] * xi[a[1]] * (g - I * theta).box() * fC,
        2, 0, 0
    )
    - quarter * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]] * (g + I * theta).box() * f,
        0, 2, 0
    )
    + quarter * sqrt(2) * I * sum(
        lambda a, b, mu:
            (LeviCivita(mu[0:4]) * eta[mu[2]][mu[4]] * eta[mu[3]][mu[5]]
            * sigmabar[mu[0]][b[0]][a[0]] * epsilon[b[0]][b[1]] * sigmabar[mu[1]][b[1]][a[1]]),
        lambda a, b, mu:
            xi[a[0]] * xi[a[1]].d(mu[4]) * g.d(mu[5]) * fC,
        2, 2, 6
    )
    + quarter * sqrt(2) * sum(
        lambda a, b, mu:
            (LeviCivita(mu[0:4]) * eta[mu[2]][mu[4]] * eta[mu[3]][mu[5]]
            * sigmabar[mu[0]][b[0]][a[0]] * epsilon[b[0]][b[1]] * sigmabar[mu[1]][b[1]][a[1]]),
        lambda a, b, mu:
            xi[a[0]] * xi[a[1]].d(mu[4]) * theta.d(mu[5]) * fC,
        2, 2, 6
    )
    + quarter * sqrt(2) * I * sum(
        lambda a, b, mu:
            (LeviCivita(mu[0:4]) * eta[mu[2]][mu[4]] * eta[mu[3]][mu[5]]
            * sigmabar[mu[0]][b[0]][a[0]] * epsilon[a[0]][a[1]] * sigmabar[mu[1]][b[1]][a[1]]),
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].d(mu[4]) * g.d(mu[5]) * f,
        2, 2, 6
    )
    - quarter * sqrt(2) * sum(
        lambda a, b, mu:
            (LeviCivita(mu[0:4]) * eta[mu[2]][mu[4]] * eta[mu[3]][mu[5]]
            * sigmabar[mu[0]][b[0]][a[0]] * epsilon[a[0]][a[1]] * sigmabar[mu[1]][b[1]][a[1]]),
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].d(mu[4]) * theta.d(mu[5]) * f,
        2, 2, 6
    )
    + quarter * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * eta[mu[0]][mu[1]],
        lambda a, b, mu:
            xi[a[0]] * xi[a[1]] * (g - I * theta).d(mu[0]) * fC.d(mu[1]),
        2, 0, 2
    )
    + quarter * sqrt(2) * sum(
        lambda a, b, mu:
            epsilon[b[0]][b[1]] * eta[mu[0]][mu[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]] * (g + I * theta).d(mu[0]) * f.d(mu[1]),
        0, 2, 2
    )
    - half * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]] * fC.d(mu[0]) * f,
        1, 1, 1
    )
    + half * I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]],
        lambda a, b, mu:
            xiC[b[0]] * xi[a[0]] * fC * f.d(mu[0]),
        1, 1, 1
    )
    + quarter * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]] * xi[a[0]] * xi[a[1]].box(),
        2, 2, 0
    )
    - quarter * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]] * eta[mu[0]][mu[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]] * xi[a[0]].d(mu[0]) * xi[a[1]].d(mu[1]),
        2, 2, 2
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]] * xi[a[0]].d(mu[0]) * xi[a[1]].d(mu[1]),
        2, 2, 2
    )
    + sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]] * eta[mu[0]][mu[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].d(mu[0]) * xi[a[0]] * xi[a[1]].d(mu[1]),
        2, 2, 2
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].d(mu[0]) * xi[a[0]] * xi[a[1]].d(mu[1]),
        2, 2, 2
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].d(mu[1]) * xi[a[0]] * xi[a[1]].d(mu[0]),
        2, 2, 2
    )
    + quarter * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
        lambda a, b, mu:
            xiC[b[0]] * xiC[b[1]].box() * xi[a[0]] * xi[a[1]],
        2, 2, 0
    )
    - quarter * sum(
        lambda a, b, mu:
            epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]] * eta[mu[0]][mu[1]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[0]) * xiC[b[1]].d(mu[1]) * xi[a[0]] * xi[a[1]],
        2, 2, 2
    )
    + sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * sigmabar[mu[1]][b[1]][a[1]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[0]) * xiC[b[1]].d(mu[1]) * xi[a[0]] * xi[a[1]],
        2, 2, 2
    )
    - L
)
X.sort()
printend([ X.iszero() ])




print("""
Field equations
***************
""")

eomG = quarter * Rational(1,16) * sum(
    lambda a, b, mu:
        epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
    lambda a, b, mu:
        (G.covD(a[0]) * G.covD(a[1]) * Gbar.covDbar(b[1])).covDbar(b[0]),
    2, 2, 0).covD2()
eomGbar = quarter * Rational(1,16) * sum(
    lambda a, b, mu:
        epsilon[a[1]][a[0]] * epsilon[b[0]][b[1]],
    lambda a, b, mu:
        (G.covD(a[1]) * Gbar.covDbar(b[0]) * Gbar.covDbar(b[1])).covD(a[0]),
    2, 2, 0).covDbar2()

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
    f (d_mu g d^mu g + d_mu theta d^mu theta + f^* f) = 0
    f^* (d_mu g d^mu g + d_mu theta d^mu theta + f^* f) = 0
    ...
""")
printbegin("Check : ")
X = [
    eomf - f * H,
    eomfC - fC * H,
    eomg
    + half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g.box() * g.d(mu[0]) * g.d(mu[1]),
        0, 0, 2
    )
    - half * sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g.box() * theta.d(mu[0]) * theta.d(mu[1]),
        0, 0, 2
    )
    + g.box() * fC * f
    + sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]).d(mu[2]) * g.d(mu[1]) * g.d(mu[3]),
        0, 0, 4
    )
    + sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
            g.d(mu[0]).d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
        0, 0, 4
    )
    + sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            theta.box() * g.d(mu[0]) * theta.d(mu[1]),
        0, 0, 2
    )
    + sum(
        lambda a, b, mu:
            eta[mu[0]][mu[1]],
        lambda a, b, mu:
            g.d(mu[0]) * (fC * f).d(mu[1]),
        0, 0, 2
    ),
    eomtheta
    + half * sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        theta.box() * theta.d(mu[0]) * theta.d(mu[1]),
        0, 0, 2
    )
    - half * sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        theta.box() * g.d(mu[0]) * g.d(mu[1]),
        0, 0, 2
    )
    + theta.box() * fC * f
    + sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
        theta.d(mu[0]).d(mu[2]) * theta.d(mu[1]) * theta.d(mu[3]),
        0, 0, 4
    )
    + sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu:
        theta.d(mu[0]).d(mu[2]) * g.d(mu[1]) * g.d(mu[3]),
        0, 0, 4
    )
    + sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        g.box() * g.d(mu[0]) * theta.d(mu[1]),
        0, 0, 2
    )
    + sum(
        lambda a, b, mu:
        eta[mu[0]][mu[1]],
        lambda a, b, mu:
        theta.d(mu[0]) * (fC * f).d(mu[1]),
        0, 0, 2
    )
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])


print("""
Supercurrent
************
""")

# Building blocks for G:
#   0 derivative:
#       G
#   1 derivative:
#       G.covD(a)
#   2 derivatives
#       G.covD2()
#       G.dd(a,b)
#   3 derivatives
#       G.covD2().covDbar(b)
#       G.covD(a).dd(a,b)
#   4 derivatives
#       G.box()
#       G.dd(a,b).dd(a,b)
#   5 derivatives
#       G.box().covD(a)
#       G.dd(a,b).dd(a,b).covD(a)
#   6 derivatives
#       G.box().dd(a,b)
#
# Building blocks for Gbar:
#   0 derivative:
#       Gbar
#   1 derivative:
#       Gbar.covDbar(b)
#   2 derivatives
#       Gbar.covDbar2()
#       Gbar.dd(a,b)
#   3 derivatives
#       Gbar.covDbar2().covD(a)
#       Gbar.covDbar(b).dd(a,b)
#   4 derivatives
#       Gbar.box()
#       Gbar.dd(a,b).dd(a,b)
#   5 derivatives
#       Gbar.box().covDbar(b)
#       Gbar.dd(a,b).dd(a,b).covDbar(b)
#   6 derivatives
#       Gbar.box().dd(a,b)

# the terms that have been commented out have been found to
# printbegin("Computing the most general J_ab...")
# Jbasis = [[[
#     # number of derivatives acting resp. on ((Gbar, Gbar), (G, G))
#     # ((6,0), (0,0))
# #    Gbar.dd(a,b).box() * Gbar * G * G,
#     # ((4,2), (0,0))
# #    Gbar.box() * Gbar.dd(a, b) * G * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a,b).dd(aa[0],bb[0]) * Gbar.dd(aa[1],bb[1]) * G * G,
# #        2, 2, 0),
#     # ((5,0), (1,0))
# #    Gbar.covDbar(b).box() * Gbar * G.covD(a) * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]).dd(a,b) * Gbar
# #            * G.covD(aa[1]) * G,
# #        2, 2, 0),
#     # ((4,1), (1,0))
# #    Gbar.box() * Gbar.covDbar(b) * G.covD(a) * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[0]).dd(a,b) * Gbar.covDbar(bb[1])
# #            * G.covD(aa[1]) * G,
# #        2, 2, 0),
#     # ((3,2), (1,0))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a,b).covDbar(bb[0]) * Gbar.dd(aa[0],bb[1])
# #            * G.covD(aa[1]) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[0]).covDbar(bb[1]) * Gbar.dd(a,b)
# #            * G.covD(aa[1]) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[0]).covDbar(b) * Gbar.dd(aa[1],bb[1])
# #            * G.covD(a) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[0]).covDbar(bb[1]) * Gbar.dd(aa[1],b)
# #            * G.covD(a) * G,
# #        2, 2, 0),
#     # ((4,0), (2,0))
# #    Gbar.box() * Gbar * G.dd(a,b) * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a,b).dd(aa[0],bb[0]) * Gbar
# #            * G.dd(aa[1],bb[1]) * G,
# #        2, 2, 0),
# #    Gbar.covDbar2().dd(a,b) * Gbar * G.covD2() * G,
#     # ((4,0), (1,1))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar2().dd(a, b) * Gbar
# #            * G.covD(aa[0]) * G.covD(aa[1]),
# #        2, 0, 0),
#     # ((3,1), (2,0))
# #    sum(lambda aa, bb, mu:
# #            epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a, b).covDbar(bb[0]) * Gbar.covDbar(bb[1])
# #            * G.covD2() * G,
# #        0, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a, bb[0]).covDbar(bb[1]) * Gbar.covDbar(b)
# #            * G.covD2() * G,
# #        0, 2, 0),
#     # ((3,1), (1,1))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a, b).covDbar(bb[0]) * Gbar.covDbar(bb[1])
# #            * G.covD(aa[0]) * G.covD(aa[1]),
# #        2, 2, 0),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(a, bb[0]).covDbar(bb[1]) * Gbar.covDbar(b)
#             * G.covD(aa[0]) * G.covD(aa[1]),
#         2, 2, 0),
#     # ((2,2), (2,0))
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.d(mu[0]) * Gbar.d(mu[1])
# #            * G.dd(a, b) * G,
# #        0, 0, 2),
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.d(mu[0]) * Gbar.dd(a,b)
# #            * G.d(mu[1]) * G,
# #        0, 0, 2),
# #    Gbar.dd(a, b) * Gbar.covDbar2() * G.covD2() * G,
#     # ((2,2), (1,1))
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(a, b) * Gbar.covDbar2()
#             * G.covD(aa[0]) * G.covD(aa[1]),
#         2, 0, 0),
#     # ((3,0), (2,1))
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(b).d(mu[0]) * Gbar
# #            * G.d(mu[1]) * G.covD(a),
# #        0, 0, 2),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]).dd(a,b) * Gbar
# #            * G.dd(aa[0],bb[1]) * G.covD(aa[1]),
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]) * Gbar
# #            * G.dd(a,b) * G.covD(aa[1]),
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]).dd(a,bb[1]) * Gbar
# #            * G.dd(aa[0],b) * G.covD(aa[1]),
# #        2, 2, 0),
#     # ((3,0), (3,0))
# #    Gbar.covDbar2().covD(a) * Gbar * G.covD2().covDbar(b) * G,
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(b).d(mu[0]) * Gbar
# #            * G.covD(a).d(mu[1]) * G,
# #        0, 0, 2),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]) * Gbar
# #            * G.covD(a).dd(aa[1],b) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(b).dd(a,bb[0]) * Gbar
# #            * G.covD(aa[0]).dd(aa[1],bb[1]) * G,
# #        2, 2, 0),
#     # ((2,1), (2,1))
#     Gbar.covDbar2() * Gbar.covDbar(b) * G.covD2() * G.covD(a),
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.d(mu[0]) * Gbar.covDbar(b)
# #            * G.d(mu[1]) * G.covD(a),
# #        0, 0, 2),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(a,b) * Gbar.covDbar(bb[0])
#             * G.dd(aa[0],bb[1]) * G.covD(aa[1]),
#         2, 2, 0),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(aa[0],bb[0]) * Gbar.covDbar(bb[1])
#             * G.dd(a,b) * G.covD(aa[1]),
#         2, 2, 0),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(a,bb[0]) * Gbar.covDbar(bb[1])
#             * G.dd(aa[0],b) * G.covD(aa[1]),
#         2, 2, 0),
#     # ((2,1), (3,0))
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.d(mu[0]) * Gbar.covDbar(b)
# #            * G.d(mu[1]).covD(a) * G,
# #        0, 0, 2),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a,b) * Gbar.covDbar(bb[0])
# #            * G.covD(aa[1]).dd(aa[0],bb[1]) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[1]) * Gbar.covDbar(bb[0])
# #            * G.covD(aa[1]).dd(a,b) * G,
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a,bb[1]) * Gbar.covDbar(bb[0])
# #            * G.covD(aa[1]).dd(aa[0],b) * G,
# #        2, 2, 0),
#     # ((1,1), (2,2))
#     sum(lambda aa, bb, mu:
#             epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
#             * G.dd(a, b) * G.covD2(),
#         0, 2, 0),
#     # ((2,0), (2,2))
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(a, b) * Gbar
# #            * G.d(mu[0]) * G.d(mu[1]),
# #        0, 0, 2),
# #    sum(lambda aa, bb, mu:
# #            eta[mu[0]][mu[1]],
# #        lambda aa, bb, mu:
# #            Gbar.d(mu[0]) * Gbar
# #            * G.d(mu[1]) * G.dd(a,b),
# #        0, 0, 2),
# #    Gbar.covDbar2() * Gbar * G.dd(a, b) * G.covD2(),
#     # ((1,1), (3,1))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
# #            * G.dd(a, b).covD(aa[0]) * G.covD(aa[1]),
# #        2, 2, 0),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
#             * G.dd(aa[0], b).covD(aa[1]) * G.covD(a),
#         2, 2, 0),
#     # ((2,0), (3,1))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar2() * Gbar
# #            * G.dd(a, b).covD(aa[0]) * G.covD(aa[1]),
# #        2, 0, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar2() * Gbar
# #            * G.dd(aa[0], b).covD(aa[1]) * G.covD(a),
# #        2, 0, 0),
#     # ((1,1), (4,0))
# #    sum(lambda aa, bb, mu:
# #            epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
# #            * G.covD2().dd(a, b) * G,
# #        0, 2, 0),
#     # ((2,0), (4,0))
# #    Gbar.dd(a,b) * Gbar * G.box() * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.dd(aa[0],bb[0]) * Gbar
# #            * G.dd(a,b).dd(aa[1],bb[1]) * G,
# #        2, 2, 0),
# #    Gbar.covDbar2() * Gbar * G.covD2().dd(a,b) * G,
#     # ((1,0), (3,2))
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar
# #            * G.dd(a,b).covD(aa[0]) * G.dd(aa[1],bb[1]),
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar
# #            * G.dd(aa[0],bb[1]).covD(aa[1]) * G.dd(a,b),
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(b) * Gbar
# #            * G.dd(aa[0],bb[0]).covD(a) * G.dd(aa[1],bb[1]),
# #        2, 2, 0),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(b) * Gbar
# #            * G.dd(aa[0],bb[0]).covD(aa[1]) * G.dd(a,bb[1]),
# #        2, 2, 0),
#     # ((1,0), (4,1))
# #    Gbar.covDbar(b) * Gbar * G.box() * G.covD(a),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar
# #            * G.dd(aa[0],bb[1]).dd(a,b) * G.covD(aa[1]),
# #        2, 2, 0),
#     # ((1,0), (5,0))
# #    Gbar.covDbar(b) * Gbar * G.covD(a).box() * G,
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar.covDbar(bb[0]) * Gbar
# #            * G.covD(aa[0]).dd(aa[1],bb[1]).dd(a,b) * G,
# #        2, 2, 0),
#     # ((0,0), (4,2))
# #    Gbar * Gbar * G.box() * G.dd(a, b),
# #    sum(lambda aa, bb, mu:
# #            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
# #        lambda aa, bb, mu:
# #            Gbar * Gbar * G.dd(a,b).dd(aa[0],bb[0]) * G.dd(aa[1],bb[1]),
# #        2, 2, 0),
#     # ((0,0), (6,0))
# #    Gbar * Gbar * G.dd(a,b).box() * G
# ]
# for b in range(2)]
# for a in range(2)]
#
# #c = list(symbols("c6000, c4200a, c4200b, c5010a, c5010b, "
# #                 "c4110a, c4110b, c3210a, c3210b, c3210c, c3210d, "
# #                 "c4020a, c4020b, c4020c, c4011, c3120a, c3120b, "
# #                 "c3111a, c3111b, c2220a, c2220b, c2220c, c2211, "
# #                 "c3021a, c3021b, c3021c, c3021d, "
# #                 "c3030a, c3030b, c3030c, c3030d, "
# #                 "c2121a, c2121b, c2121c, c2121d, c2121e, "
# #                 "c2130a, c2130b, c2130c, c2130d, "
# #                 "c1122, c2022a, c2022b, c2022c, c1131a, c1131b, "
# #                 "c2031a, c2031b, c1140, c2040a, c2040b, c2040c, "
# #                 "c1032a, c1032b, c1032c, c1032d, c1041a, c1041b, "
# #                 "c1050a, c1050b, c0042a, c0042b, c0060"))
# c = list(symbols("c3111, c2211, "
#                  "c2121a, c2121b, c2121c, c2121d, "
#                  "c1122, c1131"))
#
# J = [[c[0] * Jbasis[a][b][0]
#       for b in range(2)]
#      for a in range(2)]
# for a in range(2):
#     for b in range(2):
#         for i in range(1,len(c)):
#             J[a][b] = J[a][b] + c[i] * Jbasis[a][b][i]
# printend()

printbegin("Computing J_ab...")
J = [[
    Rational(1,64) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a, bb[0]).covDbar(bb[1]) * Gbar.covDbar(b)
            * G.covD(aa[0]) * G.covD(aa[1]),
        2, 2, 0)
    - Rational(1,128) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]],
        lambda aa, bb, mu:
            Gbar.dd(a, b) * Gbar.covDbar2()
            * G.covD(aa[0]) * G.covD(aa[1]),
        2, 0, 0)
    + Rational(1,256) * Gbar.covDbar2() * Gbar.covDbar(b) * G.covD2() * G.covD(a)
    - Rational(1,16) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,b) * Gbar.covDbar(bb[0])
            * G.dd(aa[0],bb[1]) * G.covD(aa[1]),
        2, 2, 0)
    + Rational(1,16) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(aa[0],bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(a,b) * G.covD(aa[1]),
        2, 2, 0)
    - Rational(1,16) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(aa[0],b) * G.covD(aa[1]),
        2, 2, 0)
    - Rational(1,128) * I * sum(lambda aa, bb, mu:
            epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(a, b) * G.covD2(),
        0, 2, 0)
    - Rational(1,64) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * Gbar.covDbar(bb[1])
            * G.dd(aa[0], b).covD(aa[1]) * G.covD(a),
        2, 2, 0)
    for b in range(2)]
    for a in range(2)]
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

# printbegin("Computing a basis of terms vanishing "
#            "by the field equations (1/2)...")
# zerobasis1 = [ [Gbar * eomG.covDbar(b),
#                 Gbar.covDbar(b) * eomG]
#                for b in range(2) ]
# for basis in zerobasis1:
#     for element in basis:
#         element.sort()
# printend()
# printbegin("Computing a basis of terms vanishing "
#            "by the field equations (2/2)...")
# zerobasis2 = [ [G * eomGbar.covD(a),
#                 G.covD(a) * eomGbar]
#                for a in range(2) ]
# for basis in zerobasis2:
#     for element in basis:
#         element.sort()
# printend()
#
# eq = []
# printbegin("Finding the linear combination(s) of terms \n")
# printbegin("    vanishing by the field equations (1/2)...")
# for b in range(2):
#     eq += DJ[b].linearcombinationcoeffs(zerobasis1[b], c)
# printend()
# printbegin("Finding the linear combination(s) of terms \n")
# printbegin("    vanishing by the field equations (2/2)...")
# for a in range(2):
#     eq += DbarJ[a].linearcombinationcoeffs(zerobasis2[a], c)
# printend()
#
# printbegin("Solution: ")
# sol = solve(eq, c)
# printend(sol)
#
# printbegin("Number of free coefficients: ")
# printend(len(c) - len(sol))
# if len(c) - len(sol) != 1:
#     raise Error("FreeChiral.py",
#                 "Number of free coefficients not equal to one", "")
#
# printbegin("Explicit choice of coefficients: ")
# #eq0 = eq + [c[18] - 1]
# eq0 = eq + [c[0] + I * Rational(1,32)]
# sol0 = solve(eq0, c)
# c0 = [sol0[i] for i in c]
# printend(c0)


printbegin("D^a J_ab = 0 (up to field equations) : ")
X = [
    DJ[b]
    - quarter * Gbar.covDbar(b) * eomG
    for b in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])

printbegin("Dbar^b J_ab = 0 (up to field equations) : ")
X = [
    DbarJ[a]
    - quarter * G.covD(a) * eomGbar
    for a in range(2)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])


#
# printbegin("Computing a basis of terms vanishing "
#            "by the field equations (1/2)...")
# zerobasis1 = [ [Gbar * eomG.covDbar(b),
#                 Gbar.covDbar(b) * eomG]
#                for b in range(2) ]
# for basis in zerobasis1:
#     for element in basis:
#         element.sort()
# printend()
# printbegin("Computing a basis of terms vanishing "
#            "by the field equations (2/2)...")
# zerobasis2 = [ [G * eomGbar.covD(a),
#                 G.covD(a) * eomGbar]
#                for a in range(2) ]
# for basis in zerobasis2:
#     for element in basis:
#         element.sort()
# printend()
#
# eq = []
# printbegin("Finding the linear combination(s) of terms \n")
# printbegin("    vanishing by the field equations (1/2)...")
# for b in range(2):
#     eq += DJ[b].linearcombinationcoeffs(zerobasis1[b], c)
# printend()
# printbegin("Finding the linear combination(s) of terms \n")
# printbegin("    vanishing by the field equations (2/2)...")
# for a in range(2):
#     eq += DbarJ[a].linearcombinationcoeffs(zerobasis2[a], c)
# printend()
#
# printbegin("Solution: ")
# sol = solve(eq, c)
# printend(sol)
#
# printbegin("Number of free coefficients: ")
# printend(len(c) - len(sol))
# if len(c) - len(sol) != 1:
#     raise Error("FreeChiral.py",
#                 "Number of free coefficients not equal to one", "")
#
# printbegin("Explicit choice of coefficients: ")
# eq0 = eq + [c[1] + Rational(1,8)]
# sol0 = solve(eq0, c)
# c0 = [sol0[i] for i in c]
# printend(c0)
#
# printbegin("Computing the conserved J_ab...")
# J0 = [[c0[0] * Jbasis[a][b][0]
#       for b in range(2)]
#      for a in range(2)]
# for a in range(2):
#     for b in range(2):
#         for i in range(1,16):
#             J0[a][b] = J0[a][b] + c0[i] * Jbasis[a][b][i]
# printend()
#
# printbegin("Simplification of J_ab...")
# for Jrow in J0:
#     for j in Jrow:
#         j.sort()
# printend()
#
# printbegin("Re-computing D^a J_ab...")
# DJ0 = [ sum(lambda a, bb, mu: epsilon[a[0]][a[1]],
#             lambda a, bb, mu: J0[a[0]][b].covD(a[1]),
#             2, 0, 0)
#        for b in range(2) ]
# for dj in DJ0:
#     dj.sort()
# printend()
# printbegin("Re-computing Dbar^b J_ab...")
# DbarJ0 = [ sum(lambda aa, b, mu: epsilon[b[0]][b[1]],
#                lambda aa, b, mu: J0[a][b[0]].covDbar(b[1]),
#                0, 2, 0)
#           for a in range(2) ]
# for dj in DbarJ0:
#     dj.sort()
# printend()
#
# printbegin("D^a J_ab = 0 (up to field equations) : ")
# printend([ DJ0[b].islinearcombination(zerobasis1[b])
#            for b in range(2) ])
#
# printbegin("Dbar^b J_ab = 0 (up to field equations) : ")
# printend([ DbarJ0[a].islinearcombination(zerobasis2[a])
#            for a in range(2) ])

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
    Jmu[mu].lowest().bosonicpart()
    for mu in range(4)
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])


printbegin("Computing d_mu J^mu...")
dJ = sum(lambda a, b, mu: 1,
         lambda a, b, mu: Jmu[mu[0]].d(mu[0]),
         0, 0, 1)
printend()


printbegin("d_mu J^mu = 0 (up to field equations) : ")
X = (
    dJ
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
T^munu = ...
""")
printbegin("Check: ")
X = [
    - Rational(1,8) * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * g.d(rho[2]) * g.d(rho[3]),
        0, 0, 4
    )
    + half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(nu) * g.d(rho[0]) * g.d(rho[1]),
        0, 0, 2
    )
    + quarter * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    - half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(nu) * theta.d(rho[0]) * theta.d(rho[1]),
        0, 0, 2
    )
    - half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * theta.d(mu) * theta.d(nu),
        0, 0, 2
    )
    - half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[2]) * theta.d(rho[1]) * theta.d(rho[3]),
        0, 0, 4
    )
    + sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(mu) * g.d(rho[0]) * theta.d(nu) * theta.d(rho[1]),
        0, 0, 2
    )
    + sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(nu) * g.d(rho[0]) * theta.d(mu) * theta.d(rho[1]),
        0, 0, 2
    )
    - Rational(1,8) * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
        theta.d(rho[0]) * theta.d(rho[1]) * theta.d(rho[2]) * theta.d(rho[3]),
        0, 0, 4
    )
    + half * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        theta.d(mu) * theta.d(nu) * theta.d(rho[0]) * theta.d(rho[1]),
        0, 0, 2
    )
    - half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        g.d(rho[0]) * g.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + g.d(mu) * g.d(nu) * fC * f
    - half * eta[mu][nu] * sum(
        lambda a, b, rho:
        eta[rho[0]][rho[1]],
        lambda a, b, rho:
        theta.d(rho[0]) * theta.d(rho[1]) * fC * f,
        0, 0, 2
    )
    + theta.d(mu) * theta.d(nu) * fC * f
    - half * eta[mu][nu] * fC * fC * f * f
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