def pynum(n: int):
    if n==1:
        return 1
    else: return n**2+pynum(n-1)

if __name__=="__main__":
    num=int(input("How many levels will your pyramid have? "))
    print("For {} levels, you will need {} blocks.".format(num,pynum(num)))