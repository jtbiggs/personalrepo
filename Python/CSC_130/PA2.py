def namefunc():
    name=input('Please enter your name: ')
    return name

def tscore():
    t1=int(input('Enter your score: '))
    return t1
    
def tavg(name,t1,t2):
    return 'Hi {}, your average score is: '.format(name,) +str((t1+t2)/2)

if __name__=="__main__":
    name=namefunc()
    t1=tscore()
    t2=tscore()
    avg=print(tavg(name,t1,t2))

