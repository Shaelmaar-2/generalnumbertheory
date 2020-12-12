from fastmod import fastm
from random import randint, getrandbits


# run n rabin millers on a
# if all fail return true (probably prime)
def decent_rabin_miller(a, n):
    q = a - 1
    k = 0

    while q % 2 == 0:
        k += 1
        q /= 2
    # a-1 = 2^k * q

    for _ in range(n):
        x = randint(2, a - 2)
        if fastm(x, q, a) != 1 and all(fastm(x, q*2**i, a) != a-1 for i in range(1, k)):
            return False
    return True


def newprime(bits, rounds):
    p = getrandbits(bits)
    p |= (1 << bits - 1) | 1
    while not decent_rabin_miller(p, rounds):
        p = getrandbits(bits)
        p |= (1 << bits - 1) | 1
    return p
