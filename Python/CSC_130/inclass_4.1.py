#Grade Calculator

grade=100
letter_grade=''
if grade>89.5:
    letter_grade='A'
elif grade>79.5:
    letter_grade='B'
elif grade>69.5:
    letter_grade='C'
elif grade>59.5:
    letter_grade='D'
else: letter_grade='F'
#print(letter_grade)

if (letter_grade=='A'):
    print('Congrats!')

#Repetition
counter=0
sum=0

while counter<10:
    sum+=2
    counter+=1
#print(sum)
i=0
while(i<9):
    i+=1
    j=0
    while j<9:
        j+=1
        print(str(i)+'*'+str(j)+'='+str(i*j))

bottles =99

while bottles>0:
    if bottles>1:
        print(str(bottles)+' bottles of beer on the wall.')
        print(str(bottles)+' bottles of beer.')
        print('Take one down pass it around.')
        bottles -=1
        print(str(bottles)+' bottles of beer on the wall!')

    else:
        print(str(bottles)+' bottle of beer on the wall!')
        print(str(bottles)+' bottle of beer.')
        print('Take one down pass it around.')
        print('0 bottles of beer on the wall.')
        
