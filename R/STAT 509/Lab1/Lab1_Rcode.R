### Chapter 2 Lab: Introduction to R

## 2.3.1 Basic Commands --------------------------------------------------------

# ******************************************************************************
# Example. Joint together the numbers 1, 3, 2, and 5, and to save thme as 
#          a vector named x.
# ******************************************************************************

x <- c(1, 3, 2, 5)
x

# We can also save things "=" rather than "<-"
x = c(1, 6, 2)  
x
y = c(1, 4, 3)

# Type "?" to open a new help file window
?c


# ******************************************************************************
# Example. Add two sets of numbers together.
# ******************************************************************************

# Check their length using the "length()"
length(x)
length(y)

# Add two sets of numbers together. However, x and y should be the same length.
x+y


# ******************************************************************************
# Example. Delete both x and y
# ******************************************************************************

# Look at a list of the objects, show as data and functions.
ls()

# Delete both x and y
rm(x, y)

# Remove all objects at once
ls()
rm(list=ls())


# ******************************************************************************
# Example. Create a simple matrix 
# ******************************************************************************

?matrix
x = matrix(data=c(1, 2, 3, 4), nrow=2, ncol=2)
x

x = matrix(c(1, 2, 3, 4), 2, 2)
matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE)


# ******************************************************************************
# Example. Take the square root of x and square of x
# ******************************************************************************

sqrt(x)
x^2


# ******************************************************************************
# Example. Compute the correlation between x and y
# ******************************************************************************

# Randomly generate 50 observations in a normal distribution.
x = rnorm(50)  # x ~ N(0, 1^2)
y = x + rnorm(50, mean=50, sd=.1)  # y = x + e, where e ~ N(50, 0.1^2)

# compute the correlation between x and y
cor(x, y)


# ******************************************************************************
# Example. Reproduce the exact same set of random numbers
# ******************************************************************************

set.seed(1303)
rnorm(50)


# ******************************************************************************
# Example. Find the mean, the variance and the standard deviation for y
# ******************************************************************************

set.seed(3)
y = rnorm(100)

# Compute the mean of y
mean(y)

# Compute the variance of y
var(y)

# Compute the standard deviation of y
sqrt(var(y))
sd(y)


## 2.3.2 Graphics --------------------------------------------------------------

# ******************************************************************************
# Example. Draw a scatterplot using x and y
# ******************************************************************************

x = rnorm(100)  # x ~ N(0, 1)
y = rnorm(100)  # y ~ N(0, 1)

# Draw a simple scatterplot using x and y
plot(x, y)

# Draw a scatterplot, and specify the x-axis and the y-axis label.
# Also, specify the main title for the plot.
plot(x, y, xlab="this is the x-axis",
           ylab="this is the y-axis",
           main="Plot of X vs Y")


# ******************************************************************************
# Example. Export the plot to a PDF file
# ******************************************************************************

# Create a pdf file named Figure.pdf
pdf("Figure.pdf")  # You can use jpeg() to create a jpeg file.
plot(x, y, col="green")  # Specify the dot color as green.
dev.off()


# ******************************************************************************
# Example. Create a sequence of numbers
# ******************************************************************************

x = seq(1, 10)
x
x = 1:10  # This is a short cut for seq(1, 10)
x

# Create a sequence with a length of 50 between -3.14 and 3.14
x = seq(-pi, pi, length=50)


# ******************************************************************************
# Example. Create a contour plot
# ******************************************************************************

y = x
f = outer(x, y, function(x,y)cos(y)/(1+x^2))
contour(x, y, f)
contour(x, y, f, nlevels=45, add=T)
fa = (f-t(f))/2
contour(x, y, fa, nlevels=15)

# The image() works the same way as contour(), expect that it produces 
# a color-coded plot whose colors depend on the z value.
image(x, y, fa)  # This is known as a heatmap


# ******************************************************************************
# Example. Produce a three-dimensional plot. The arguments theta and phi 
#          control the angles at which the plot is viewed.
# ******************************************************************************

persp(x, y, fa)
persp(x, y, fa, theta=30)
persp(x, y, fa, theta=30, phi=20)
persp(x, y, fa, theta=30, phi=70)
persp(x, y, fa, theta=30, phi=40)


## 2.3.3 Indexing Data ---------------------------------------------------------

# ******************************************************************************
# Example. Examine part of a set of data.
# ******************************************************************************

A = matrix(1:16, 4, 4)
A

# We want to check the element corresponding to 
# the 2nd row and the 3rd column.
A[2,3]

# Check the element corresponding to the combination of 
# the 1st, the 3rd rows, and the 2nd, the 4th columns.
A[c(1,3),c(2,4)]

# More options
A[1:3, 2:4]
A[1:2, ]
A[, 1:2]
A[1, ]

# Check all elements except the 1st row and the 3rd row.
A[-c(1,3), ]

# More example.
A[-c(1,3), -c(1,3,4)]

# Check the dimension of the matrix A.
dim(A)


## 2.3.4 Loading Data ----------------------------------------------------------

# ******************************************************************************
# Example. Load the data in your computer
# ******************************************************************************

# Let R search your data under your data file
setwd("C:/Users/gdsdl/Desktop/UCFClassFiles/STA 4241/Class note/2_3")

# Read the table using read.table()
Auto = read.table("Auto.data")  # The data will be stored as a data frame.

# Use the fix() function to view your data in a spreadsheet like window
fix(Auto)

# The options: header and na.strings.
# header: Use the first line of the file as the variable names.
# na.string: Treat the "?" as a missing element of the data matrix.
Auto = read.table("Auto.data", header=T, na.strings="?")
fix(Auto)

# ******************************************************************************
# Example. Load a csv file (comma separated value) file 
# ******************************************************************************

Auto = read.csv("Auto.csv", header=T, na.strings="?")
fix(Auto)
dim(Auto)
Auto[1:4, ]

# Remove the rows (the observations) with any missing variables
Auto = na.omit(Auto)
dim(Auto)
fix(Auto)

# Check the variable names in our data
names(Auto)


## 2.3.5 Additional Graphical and Numerical Summaries --------------------------

# ******************************************************************************
# Example. Plot the cylinders and the mpg variables from the Auto data
# ******************************************************************************

plot(cylinders, mpg)  # Obiviously, this cannot work
plot(Auto$cylinders, Auto$mpg)  # This is a regular scatterplot

# We can use attach() to make the variable in the data frame (Auto) avaliable 
# by name.
attach(Auto)
plot(cylinders, mpg)

# We want to treat cylinder as a categorical variable.
# (In R, its called factor.)
cylinders = as.factor(cylinders)

# Now, since the x-axis is categorical, the plot is changed into a boxplot.
plot(cylinders, mpg)
plot(cylinders, mpg, col="red")  # Change the color.

# "varwidth=T": The box width will be proportional to the  square root of the 
#               number of observations in the data. 
plot(cylinders, mpg, col="red", varwidth=T) 

# "horizontal=T": Creates horizontal boxes
plot(cylinders, mpg, col="red", varwidth=T, horizontal=T)

# Create text labels for x-axis and y-axis
plot(cylinders, mpg, col="red", varwidth=T, xlab="cylinders", ylab="MPG")


# ******************************************************************************
# Example. Plot a histogram using the mpg variable
# ******************************************************************************

hist(mpg)
hist(mpg, col=2)  # Change the color into red.
hist(mpg, col=2, breaks=15)


# ******************************************************************************
# Example. Use pair() to create a scatterplot matrix
# ******************************************************************************

pairs(Auto[, -9])  # We can only use non-character variables.

# We can also produce scatterplots for just a subset of the variables
pairs(~ mpg + displacement + horsepower + weight + acceleration, Auto)


# ******************************************************************************
# Example. identify() provides a useful interactive method for identifying
#          the value for a particular variable for points on a plot.
# ******************************************************************************

plot(horsepower, mpg)

# Three inputs for the identify function:
#   x: The x-axis. 
#   y: The y-axis. 
#   labels: the variable whose values we would like to see printed 
#           for each point.
identify(x = horsepower, y = mpg, labels = name)


# ******************************************************************************
# Example. Use summary() function to produce a numerical summary for each 
#          variable in the Auto data
# ******************************************************************************

summary(Auto)

# Can also produce a summary of just a single variable.
summary(mpg)
