#!/usr/bin/env python3.5

from SUSYalgebra import *
from sympy import solve


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

dL = sum(
    lambda a, b, mu: eta[mu[0]][mu[1]],
    lambda a, b, mu: Gbar.d(mu[0]) * G.d(mu[1]),
    0, 0, 2
)
dL.sort()

L = dL.Dterm()

print("""
L = 1/2 (d_mu d^mu g)^2
    + 1/2 (d_mu d^mu theta)^2
    + d_mu f^* d^mu f
    + fermions
""")
printbegin("Check : ")
X = (
    half * g.box() * g.box()
    + half * theta.box() * theta.box()
    + sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: fC.d(mu[0]) * f.d(mu[1]),
        0, 0, 2
    )
    - I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[1]) * xi[a[0]].d(mu[0]).d(mu[2]),
        1, 1, 3
    )
    + I * sum(
        lambda a, b, mu:
            sigmabar[mu[0]][b[0]][a[0]] * eta[mu[1]][mu[2]],
        lambda a, b, mu:
            xiC[b[0]].d(mu[1]).d(mu[0]) * xi[a[0]].d(mu[2]),
        1, 1, 3
    )
    # boundary terms
    - half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: (g.d(mu[0]) * g.box()).d(mu[1]),
        0, 0, 2
    )
    + quarter * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: (g.d(mu[0]) * g.d(mu[1]).d(mu[2])).d(mu[3]),
        0, 0, 4
    )
    - half * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]],
        lambda a, b, mu: (theta.d(mu[0]) * theta.box()).d(mu[1]),
        0, 0, 2
    )
    + quarter * sum(
        lambda a, b, mu: eta[mu[0]][mu[1]] * eta[mu[2]][mu[3]],
        lambda a, b, mu: (theta.d(mu[0]) * theta.d(mu[1]).d(mu[2])).d(mu[3]),
        0, 0, 4
    )
    - L
)
X.sort()
printend([ X.iszero() ])


print("""
Field equations
***************
""")

eomG = quarter * G.box().covD2()
eomGbar = quarter * Gbar.box().covDbar2()

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
    d_mu d^mu f = 0
    d_mu d^mu f^* = 0
    d_mu d^mu d_nu d^nu g = 0
    d_mu d^mu d_nu d^nu theta = 0
""")

printbegin("Check : ")
X = [
    eomf + f.box(),
    eomfC + fC.box(),
    eomg - g.box().box(),
    eomtheta - theta.box().box()
]
for x in X:
    x.sort()
printend([x.iszero() for x in X])


print("""
Supercurrent
************
""")


# printbegin("Computing the most general J_ab...")
# Jbasis = [[[
#     # (6,0) : number of cov. derivatives acting resp. on Gbar and G
#     Gbar.dd(a,b).box() * G,
#     # (5,1)
#     Gbar.covDbar(b).box() * G.covD(a),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(bb[0]).dd(aa[0],bb[1]).dd(a,b) * G.covD(aa[1]),
#         2, 2, 0),
#     # (4,2)
#     Gbar.box() * G.dd(a,b),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(a,b).dd(aa[0],bb[0]) * G.dd(aa[1],bb[1]),
#         2, 2, 0),
#     Gbar.covDbar2().dd(a,b) * G.covD2(),
#     # (3,3)
#     Gbar.covDbar2().covD(a) * G.covD2().covDbar(b),
#     sum(lambda aa, bb, mu:
#             eta[mu[0]][mu[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(b).d(mu[0]) * G.covD(a).d(mu[1]),
#         0, 0, 2),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(bb[0]).dd(aa[0],bb[1]) * G.covD(a).dd(aa[1],b),
#         2, 2, 0),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(b).dd(a,bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]),
#         2, 2, 0),
#     # (2,4)
#     Gbar.dd(a,b) * G.box(),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.dd(aa[0],bb[0]) * G.dd(a,b).dd(aa[1],bb[1]),
#         2, 2, 0),
#     Gbar.covDbar2() * G.covD2().dd(a,b),
#     # (1,5)
#     Gbar.covDbar(b) * G.covD(a).box(),
#     sum(lambda aa, bb, mu:
#             epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
#         lambda aa, bb, mu:
#             Gbar.covDbar(bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]).dd(a,b),
#         2, 2, 0),
#     # (0,6)
#     Gbar * G.dd(a,b).box()
#     ]
#     for b in range(2)]
#     for a in range(2)]
#
# c = list(symbols("c6, c5a, c5b, c4a, c4b, c4c, c3a, c3b, c3c, c3d, c2a, c2b, c2c, c1a, c1b, c0"))
# J = [[c[0] * Jbasis[a][b][0]
#       for b in range(2)]
#      for a in range(2)]
# for a in range(2):
#     for b in range(2):
#         for i in range(1,16):
#             J[a][b] = J[a][b] + c[i] * Jbasis[a][b][i]
# printend()


printbegin("Computing J_ab...")
J = [[
    - Rational(1,24) * Gbar.covDbar(b).box() * G.covD(a)
    + Rational(1,24) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]).dd(a,b) * G.covD(aa[1]),
        2, 2, 0)
    - Rational(1,6) * I * Gbar.box() * G.dd(a,b)
    + Rational(1,12) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(a,b).dd(aa[0],bb[0]) * G.dd(aa[1],bb[1]),
        2, 2, 0)
    + Rational(1,64) * I * Gbar.covDbar2().dd(a,b) * G.covD2()
    + Rational(1,384) * Gbar.covDbar2().covD(a) * G.covD2().covDbar(b)
    - Rational(1,12) * sum(lambda aa, bb, mu:
            eta[mu[0]][mu[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(b).d(mu[0]) * G.covD(a).d(mu[1]),
        0, 0, 2)
    - Rational(1,12) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]).dd(aa[0],bb[1]) * G.covD(a).dd(aa[1],b),
        2, 2, 0)
    + Rational(1,12) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(b).dd(a,bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]),
        2, 2, 0)
    + Rational(1,6) * I * Gbar.dd(a,b) * G.box()
    - Rational(1,12) * I * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.dd(aa[0],bb[0]) * G.dd(a,b).dd(aa[1],bb[1]),
        2, 2, 0)
    - Rational(1,64) * I * Gbar.covDbar2() * G.covD2().dd(a,b)
    - Rational(1,24) * Gbar.covDbar(b) * G.covD(a).box()
    - Rational(1,24) * sum(lambda aa, bb, mu:
            epsilon[aa[0]][aa[1]] * epsilon[bb[0]][bb[1]],
        lambda aa, bb, mu:
            Gbar.covDbar(bb[0]) * G.covD(aa[0]).dd(aa[1],bb[1]).dd(a,b),
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
J^mu = 2/3 d_nu ( d^mu g d^nu theta - d^nu g d^mu theta )
       - i/2 ( d^mu f^* f - f^* d^mu f )
       + ...
""")
printbegin("Check: ")
X = [
    - Rational(1,3) * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(rho[0]) * g.d(mu).d(rho[1]),
        0, 0, 2
    )
    + Rational(1,3) * theta.d(mu) * g.box()
    - Rational(1,3) * theta.box() * g.d(mu)
    + Rational(1,3) * sum(
        lambda a, b, rho: eta[rho[0]][rho[1]],
        lambda a, b, rho: theta.d(mu).d(rho[0]) * g.d(rho[1]),
        0, 0, 2
    )
    + half * I * fC.d(mu) * f
    - half * I * fC * f.d(mu)
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
    - half * eta[mu][nu] * g.box() * g.box()
    + 2 * g.d(mu).d(nu) * g.box()
    + Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            g.d(rho[0]).d(rho[2]) * g.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    - Rational(4,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(mu).d(rho[0]) * g.d(nu).d(rho[1]),
        0, 0, 2
    )
    + Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.box().d(rho[1]),
        0, 0, 2
    )
    - g.d(mu) * g.box().d(nu)
    - g.d(nu) * g.box().d(mu)
    + Rational(2,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            g.d(rho[0]) * g.d(mu).d(nu).d(rho[1]),
        0, 0, 2
    )
    - half * eta[mu][nu] * theta.box() * theta.box()
    + 2 * theta.d(mu).d(nu) * theta.box()
    + Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]] * eta[rho[2]][rho[3]],
        lambda a, b, rho:
            theta.d(rho[0]).d(rho[2]) * theta.d(rho[1]).d(rho[3]),
        0, 0, 4
    )
    - Rational(4,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(mu).d(rho[0]) * theta.d(nu).d(rho[1]),
        0, 0, 2
    )
    + Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.box().d(rho[1]),
        0, 0, 2
    )
    - theta.d(mu) * theta.box().d(nu)
    - theta.d(nu) * theta.box().d(mu)
    + Rational(2,3) * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            theta.d(rho[0]) * theta.d(mu).d(nu).d(rho[1]),
        0, 0, 2
    )
    - Rational(1,3) * eta[mu][nu] * sum(
        lambda a, b, rho:
            eta[rho[0]][rho[1]],
        lambda a, b, rho:
            fC.d(rho[0]) * f.d(rho[1]),
        0, 0, 2
    )
    + Rational(2,3) * fC.d(mu) * f.d(nu)
    + Rational(2,3) * fC.d(nu) * f.d(mu)
    + Rational(1,3) * eta[mu][nu] * fC * f.box()
    - Rational(1,3) * fC * f.d(mu).d(nu)
    + Rational(1,3) * eta[mu][nu] * fC.box() * f
    - Rational(1,3) * fC.d(mu).d(nu) * f
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


