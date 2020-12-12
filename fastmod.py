def fastm(a, k, m):
    b = 1
    while k >= 1:
        if (k % 2) == 1:
            b = (b * a) % m
        a = (a**2) % m
        k = k//2
    return b
