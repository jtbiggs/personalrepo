def fact(n):
    if n<2:
        return 1
    else: return n*fact(n-1)

def comb(n,r):
    return fact(n)/(fact(r)*fact(n-r))

a=9+6+12
print((comb(48,1)*comb(4,2))/comb(52,3)+(comb(48,2)*comb(4,1))/comb(52,3)+(comb(48,0)*comb(4,3))/comb(52,3))