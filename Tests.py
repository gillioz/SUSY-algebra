#!/usr/bin/env python3.5

from SUSYalgebra import *

print("""
------------
Spin algebra
------------
""")

print("epsilon^ab epsilon_bc = delta^a_c :",
      [
            sum(
                  lambda aa, b, mu:
                        epsilon[a1][aa[0]] * epsilon[a2][aa[0]],
                  lambda aa, b, mu: 1,
                  1, 0, 0
            ) == delta[a1][a2]
            for a1 in range(0, 2)
            for a2 in range(0, 2)
      ])
print()

print("(sigma^mu)_ab (sigmabar^nu)^ba = 2 eta^{mu nu} :")
print([
      sum(
            lambda a, b, rho:
                  sigma[mu][a[0]][b[0]]
                  * sigmabar[nu][b[0]][a[0]],
            lambda a, b, rho: 1,
            1, 1, 0
      ) == 2 * eta[mu][nu]
      for mu in range(0, 4)
      for nu in range(0, 4)
])
print()

print("eta_{mu nu} (sigma^mu)_ab (sigma^nu)_cd"
      " = 2 epsilon_ac epsilon_bd :")
print([
      sum(
            lambda aa, bb, mu:
                  eta[mu[0]][mu[1]]
                  * sigma[mu[0]][a][b]
                  * sigma[mu[1]][c][d],
            lambda aa, bb, mu: 1,
            0, 0, 2
      ) == 2 * epsilon[c][a] * epsilon[d][b]
      for a in range(0, 2)
      for b in range(0, 2)
      for c in range(0, 2)
      for d in range(0, 2)
])
print()

print("""
-----------
SuperFields
-----------
""")

# definition of a field
phi = Field("\\phi", dim=1)
phi.check()

# check the validity of its derivatives
phi.d(0).check()
phi.box().check()

# definition of a spinor field
psi = Field("\\psi", [1], [], dim = 1 + half)
psi.check()

# definition of a superfield
S = SuperField("")
S.check()

# check the validity of its derivatives
S.d(1).check()
S.box().check()
S.covD(0).check()
S.covD(1).check()
S.covDbar(0).check()
S.covDbar(1).check()

#print("Components of an unconstrained superfield:")
#S.printall()
#print()

print("Chiral and anti-chiral superfields:")
Phi = SuperField("chiral")
Phi.check()
print("  Dbar_b Phi = 0 :",
      [Phi.covDbar(0).iszero(),
       Phi.covDbar(1).iszero()])
Phi.printlowest()
Phi.covD(0).printlowest()
Phi.covD(1).printlowest()
(-quarter * Phi.covD2()).printlowest()
Phibar = SuperField("antichiral")
Phibar.check()
print("  D_a Phibar = 0 :",
      [Phibar.covD(0).iszero(),
       Phibar.covD(1).iszero()])
Phibar.printlowest()
Phibar.covDbar(0).printlowest()
Phibar.covDbar(1).printlowest()
(-quarter * Phibar.covDbar2()).printlowest()
print()

print("Vector superfield:")
V = SuperField("vector")
V.check()
W = [ -quarter * V.covD(a).covDbar2() for a in range(0, 2) ]
W[0].check()
W[1].check()
print("  Dbar_b W_a = 0 :",
      [ W[a].covDbar(b).iszero() for a in range(0, 2)
        for b in range(0, 2)])
Wbar = [ -quarter * V.covDbar(b).covD2() for b in range(0, 2) ]
Wbar[0].check()
Wbar[1].check()
print("  D_a Wbar_b = 0 :",
      [ Wbar[b].covD(a).iszero() for a in range(0, 2)
        for b in range(0, 2)])
W[0].printlowest()
W[1].printlowest()
Wbar[0].printlowest()
Wbar[1].printlowest()
DW = sum(
    lambda a, b, mu: -quarter * epsilon[a[1]][a[0]],
    lambda a, b, mu: W[a[1]].covD(a[0]),
    2, 0, 0
)
DW.sort()
DW.check()
DbarWbar = sum(
    lambda a, b, mu: -quarter * epsilon[b[0]][b[1]],
    lambda a, b, mu: Wbar[b[1]].covDbar(b[0]),
    0, 2, 0
)
DbarWbar.sort()
DbarWbar.check()
X = DW - DbarWbar
X.sort()
print("  D^a W_a = Dbar_b W^b :", [ X.iszero() ])
DW.printlowest()
#V.printall()
#W[0].printall()
#Wbar[1].printall()
print()

print("Commutation relations for superfields:")
S1 = SuperField(1)
S2 = SuperField(2)
X = S1 * S2 - S2 * S1
X.sort()
print("  S1 * S2 - S2 * S1 = 0 :",
      [X.iszero()])
X = S1 * S2.covDbar(0) - S2.covDbar(0) * S1
X.sort()
print("  S1 * Dbar_1 S2 - Dbar_1 S2 * S1 = 0 :",
      [X.iszero()])
X = S1.covD(1) * S2.d(0) - S2.d(0) * S1.covD(1)
X.sort()
print("  D_2 S1 * d_0 S2 - d_0 S2 * D_2 S1 = 0 :",
      [X.iszero()])
X = S1.covD(0) * S2.covDbar(1) + S2.covDbar(1) * S1.covD(0)
X.sort()
print("  D_1 S1 * Dbar_2 S2 + Dbar_2 S2 * D_1 S1 = 0 :",
      [X.iszero()])
print()

print("Associativity of the product:")
S3 = SuperField(3)
X = (S1 * S2) * S3 - S1 * (S2 * S3)
X.sort()
print("  (S1 * S2) * S3 - S1 * (S2 * S3) = 0 :",
      [X.iszero()])
print()

print("""
---------------------
Covariant derivatives
---------------------
""")

X = [ S.covD(a1).covD(a2) + S.covD(a2).covD(a1)
      for a1 in range(0, 2)
      for a2 in range(0, 2)]
for x in X:
    x.sort()
print("{D_a1, D_a2} = 0 :", [ x.iszero() for x in X ])
X = [ S.covDbar(b1).covDbar(b2) + S.covDbar(b2).covDbar(b1)
      for b1 in range(0, 2)
      for b2 in range(0, 2)]
for x in X:
    x.sort()
print("{Dbar_b1, Dbar_b2} = 0 :", [ x.iszero() for x in X ])
X = [ S.covD(a).covDbar(b) + S.covDbar(b).covD(a)
      + 2 * I * S.dd(a, b)
      for a in range(0, 2)
      for b in range(0, 2)]
for x in X:
    x.sort()
print("{D_a, Dbar_b} = -2 i d_ab :", [ x.iszero() for x in X ])
print()

X = [ S.covD(a1).covD(a2).covD(a3)
      for a1 in range(0, 2)
      for a2 in range(0, 2)
      for a3 in range(0, 2)]
for x in X:
    x.sort()
print("D_a D_b D_c = 0 :", [ x.iszero() for x in X ])
X = [ S.covDbar(b1).covDbar(b2).covDbar(b3)
      for b1 in range(0, 2)
      for b2 in range(0, 2)
      for b3 in range(0, 2)]
for x in X:
    x.sort()
print("Dbar_a Dbar_b Dbar_c = 0 :", [ x.iszero() for x in X ])
print()

X = [ S.covD(a2).covD(a1) - half * epsilon[a2][a1] * S.covD2()
      for a1 in range(0, 2)
      for a2 in range(0, 2)]
for x in X:
      x.sort()
print("D_a1 D_a2 = 1/2 epsilon_{a1 a2} D^2 :", [ x.iszero() for x in X ])
X = [ S.covDbar(b2).covDbar(b1) + half * epsilon[b2][b1] * S.covDbar2()
      for b1 in range(0, 2)
      for b2 in range(0, 2)]
for x in X:
      x.sort()
print("Dbar_b1 Dbar_b2 = 1/2 epsilon_{b1 b2} Dbar^2 :", [ x.iszero() for x in X ])
print()

X = [ S.covDbar2().covD(a) - S.covD(a).covDbar2()
      + sum(
            lambda aa, b, mu: 4 * I * epsilon[b[0]][b[1]],
            lambda aa, b, mu: S.covDbar(b[1]).dd(a, b[0]),
            0, 2, 0
      )
      for a in range(0, 2)]
for x in X:
      x.sort()
print("[ D_a, Dbar^2 ] = -4 i d_ab Dbar^b :", [ x.iszero() for x in X ])
X = [ S.covD2().covDbar(b) - S.covDbar(b).covD2()
      + sum(
            lambda a, bb, mu: -4 * I * epsilon[a[0]][a[1]],
            lambda a, bb, mu: S.covD(a[1]).dd(a[0], b),
            2, 0, 0
      )
      for b in range(0, 2)]
for x in X:
      x.sort()
print("[ Dbar_b, D^2 ] = 4 i d_ab D^a :", [ x.iszero() for x in X ])
print()

X = S.covD2().covDbar2().covD2() + 16 * S.covD2().box()
X.sort()
print("D^2 Dbar^2 D^2 = - 16 box D^2 :", [ X.iszero() ])
X = S.covDbar2().covD2().covDbar2() + 16 * S.covDbar2().box()
X.sort()
print("Dbar^2 D^2 Dbar^2 = - 16 box Dbar^2 :", [ X.iszero() ])
print()


