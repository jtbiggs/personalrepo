#Simple algorithm for checking primality
def isprime(x: int):
    if x>0:
        xold=x-1
        while xold>1:
            if x%xold==0:
                return False
            else: xold=xold-1
        return True
    else:
        xold=x+1
        while xold<-1:
            if x%xold==0:
                return False
            else: xold=xold+1
        return True

#Determines parity
def isodd(x: int):
    if x%2!=0:
        return 'odd'
    else: return 'even'

#Returns max number in a tuple
def largestnumber(x,y,z):
   return max(x,y,z)

#Returns max number in a tuple
def maxproduct(x,y,z):
    S=[x,y,z]
    max1=max(tuple(S))
    S.remove(max1)
    max2=max(tuple(S))
    return max1*max2

# Main
if __name__=="__main__":
    x=int(input('Enter a number:\t'))
    y=int(input('Enter a number:\t'))
    z=int(input('Enter a number:\t'))

    print('--------------------------')
    print('Num\tParity\tPrime')
    print('--------------------------')
    print('{}\t{}\t{}'.format(x,isodd(x),isprime(x)))
    print('{}\t{}\t{}'.format(y,isodd(y),isprime(y)))
    print('{}\t{}\t{}'.format(z,isodd(z),isprime(z)))
    print('--------------------------')
    print('The largest number is {}'.format(largestnumber(x,y,z)))
    print('The largest possible product is {}'.format(maxproduct(x,y,z)))