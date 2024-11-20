# BLS12-381 parameters (https://pkg.go.dev/github.com/consensys/gnark-crypto/ecc/bls12-381)
k = 12
q = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513

h = (q^12 - 1) // r

from sage.modules.free_module_integer import IntegerLattice
M = matrix([[r, -q, -q^2, -q^3],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
L = IntegerLattice(M.transpose()) # transpose because Sage uses columns for lattices
v = L.shortest_vector()

# the shortest vector is (1, 0, -1, -15132376222941642752), which means that lambda = 1 - q^2 + x * q^3
assert(v[0] == 1)
assert(v[1] == 0)
assert(v[2] == -1)
assert(v[3] == -15132376222941642752)
assert(len(v) == 4)

x = -15132376222941642752
lambda_value = 1 - q^2 + x * q^3
m = lambda_value // r



# Tonelli-Shanks parameters

import random

quadratic_non_residue = 1
while kronecker(quadratic_non_residue, (q ** 12 - 1)) != -1:
    quadratic_non_residue = random.sample(range(99999999), 1)


def find_s_t(n):
    for s in range(1, 50):
        if n % (2**s) == 0:
            t = n / 2**s
            assert t.is_integer()
            if not ((t-1)/2).is_integer():
                continue
            return s, t


print("h: ", h)
print("lambda: ", lambda_value)
print("m: ", m)
print("quadratic_non_residue: ", quadratic_non_residue)
print("s, t: ", find_s_t(q ** 12 - 1))