#Justin Biggs - Lab 1 Submission
#6/11/2024


x <- c(1, 3, 2, 5)
x

x = c(1, 6, 2)  
x

y = c(1, 4, 3)

?c

length(x)
length(y)

x+y

ls()

rm(x, y)

ls()
rm(list=ls())

?matrix
x = matrix(data=c(1, 2, 3, 4), nrow=2, ncol=2)
x

x = matrix(c(1, 2, 3, 4), 2, 2)
matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE)

sqrt(x)
x^2

x = rnorm(50)  # x ~ N(0, 1^2)
y = x + rnorm(50, mean=50, sd=.1)

cor(x,y)

