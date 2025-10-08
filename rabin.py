from typing import List

def rabin_karp_double(text: str, pattern: str) -> List[int]:
    if len(pattern) == 0:
        return list(range(len(text)+1))
    n, m = len(text), len(pattern)
    base1, mod1 = 257, 2**61-1  # big prime-like modulus using Mersenne trick
    base2, mod2 = 263, 10**9+7

    def modmul(a,b,mod):
        return (a*b) % mod

    # precompute base^m
    power1 = pow(base1, m, mod1)
    power2 = pow(base2, m, mod2)

    ph1 = ph2 = 0
    th1 = th2 = 0
    for i in range(m):
        ph1 = (ph1*base1 + ord(pattern[i])) % mod1
        ph2 = (ph2*base2 + ord(pattern[i])) % mod2
        th1 = (th1*base1 + ord(text[i])) % mod1
        th2 = (th2*base2 + ord(text[i])) % mod2

    res = []
    if (ph1,ph2) == (th1,th2) and text[:m] == pattern:
        res.append(0)
    for i in range(m, n):
        # remove leading char, add trailing char
        th1 = (th1*base1 - ord(text[i-m])*power1 + ord(text[i])) % mod1
        th2 = (th2*base2 - ord(text[i-m])*power2 + ord(text[i])) % mod2
        if (th1,th2) == (ph1,ph2):
            start = i-m+1
            if text[start:start+m] == pattern:
                res.append(start)
    return res

# Example
if __name__ == "__main__":
    t = "ababcabcababc"
    p = "abc"
    print(rabin_karp_double(t,p))  # [2,5,9]
