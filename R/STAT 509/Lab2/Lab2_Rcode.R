### Chapter 3 Lab: Linear Regression

## 3.6.1 Libraries -------------------------------------------------------------

library(MASS)  # contains Boston data set.

# Need to download the ISLR library. on a Windows system, select
# Packages -> Install package(s) -> Select any location you want 
# (US location is preferred) -> ISLR
library(ISLR)


## 3.6.2 Simple Linear Regression. ---------------------------------------------

# ******************************************************************************
# Example. Need to load the Boston data set in the MASS library,
#          - Number of observations: 506.
#          - Response: medv
#          - Predictor: 13 predictors such as rm, lstat,...,etc.
# ******************************************************************************

fix(Boston)  # Display the Boston data.
names(Boston)  # Display the variable for the Boston data.


# ******************************************************************************
# Example. Perform a linear regression of
#          medv = beta_0 + beta_1(lstat) + epsilon.
# ******************************************************************************

lm.fit = lm(medv ~ lstat)  # Cause error (why?)
lm.fit = lm(medv ~ lstat, data = Boston)  


attach(Boston)
lm.fit = lm(medv ~ lstat)  # Now it does not cause any problem (Why?)


# ******************************************************************************
# Example.  Some very basic information about the linear regression.
# ******************************************************************************

lm.fit
summary(lm.fit)  # Summary, including F-statistic, R-squared,...,etc.
coef(lm.fit)  # Display the estimates only.
confint(lm.fit)  # Display the confidence intervals for the estimates.

# We can use the names() function in order to find out what other pieces of
# information are stored in lm.fit.
names(lm.fit)

# Similar to the coef(). But is much safer to use the function.
lm.fit$coefficients 


# ******************************************************************************
# Example. The predict() function can be used to produce confidence interval and 
#          prediction interval for the prediction of medv for a given value of 
#          lstat.
# ******************************************************************************

# The 95% PI is wider than the 95% CI.
predict(lm.fit, data.frame(lstat = (c(5, 10, 15))), interval = "confidence")
predict(lm.fit, data.frame(lstat = (c(5, 10, 15))), interval = "prediction")

# The 99% CI.
predict(lm.fit, data.frame(lstat = (c(5, 10, 15))), interval = "confidence", 
        level = 0.99)


# ******************************************************************************
# Example. Plot medv and lstat along with the least sqAuares regression line
# ******************************************************************************

plot(lstat, medv)
abline(lm.fit)
abline(lm.fit, lwd = 3)  # Line width = 3
abline(lm.fit, lwd = 3, col = "red")  # line color = "red"

# Similarly, we can design the color of our dots.
plot(lstat, medv, col = "red")

# Different plotting symbols
plot(lstat,medv, pch = 20)
plot(lstat,medv, pch = "+")
plot(1:20, 1:20, pch = 1:20)

# Using par() divides the plotting region into a 2X2 grid of panels
par(mfrow = c(2, 2))
plot(lm.fit)


par(mfrow = c(1, 2))  # 1X2 grid of panels

# The residuals() will return the residuals.
plot(predict(lm.fit), residuals(lm.fit))

# The rstudent() will return the studentized residuals (scaled residuals)
plot(predict(lm.fit), rstudent(lm.fit))


par(mfrow = c(1, 1))  # 1X1 grid of panels

# Using hatvalues() to display the leverage statistics.
plot(hatvalues(lm.fit))

# Using which.max() to determine the largest element of a vector 
which.max(hatvalues(lm.fit))  # The 375th element is the largest
max(hatvalues(lm.fit))


## 3.6.3 Multiple Linear Regression. -------------------------------------------

# ******************************************************************************
# Example. Perform a linear regression of
#          medv = beta_0 + beta_1(lstat) + beta_2(age) + epsilon.
# ******************************************************************************

lm.fit = lm(medv ~ lstat+age, data = Boston)
summary(lm.fit)


# ******************************************************************************
# Example. Perform a linear regression of
#          medv = beta_0 + beta_1(lstat) + beta_2(age) +...
#                 + \beta_13(crim) + epsilon.
# ******************************************************************************

lm.fit = lm(medv ~ ., data = Boston)  # Simply use "." to represent the rest.
summary(lm.fit)
summary(lm.fit)$r.sq  # Only display the R-squares for the model


# Need to download the car and all other libraries. 
install.packages("car") 
install.packages("carData")
install.packages("hms")
install.packages("tibble")
install.packages("openxlsx")
install.packages("readxl")


# If this is your first time to use R, you will have to install
# more libraries than I displayed above.
# 
# Unfortunately, you will have to manually install these packages before you 
# can use the car library. Please check your red error message to install those 
# libraries. 
library(car)

# The vif() function is in the car library. 
# Most of the VIF is less than 10. Hence there is no collinearity problem for 
# our model.
vif(lm.fit)


# ******************************************************************************
# Example. Perform a linear regression of for the Boston using all the variables
#          but the age.
# ******************************************************************************

coefficients(lm.fit)

lm.fit1 = lm(medv ~ . -age, data=Boston) 
coefficients(lm.fit1)

lm.fit1 = update(lm.fit, ~ . -age)  # Similarly, you can use update().
coefficients(lm.fit1)


## 3.6.4 Interaction Terms. ----------------------------------------------------

# ******************************************************************************
# Example. Perform a linear regression of
#          medv = beta_0 + beta_1(lstat) + beta_2(age) 
#                 + beta_3(lstat : age) + epsilon.
# ******************************************************************************

lm.fit2 = lm(medv ~ lstat*age, data = Boston)
summary(lm.fit2)


## 3.6.5 Non-linear Transformations of the Predictors --------------------------

# ******************************************************************************
# Example. Perform a linear regression of
#          medv = beta_0 + beta_1(lstat) + beta_2(lstat)^2 + epsilon.
# ******************************************************************************

# To do so, you will need to use I() function inside of the lm() function.
lm.fit2=lm(medv ~ lstat + I(lstat^2), data = Boston)
summary(lm.fit2)  # R-sq = 0.6047

# Compare it with the original model
lm.fit = lm(medv ~ lstat, data = Boston)
summary(lm.fit)  # R-sq = 0.5441


# Use anova function to further quantify the extent to which the quadratic fit
# is superior to the linear fit.
anova(lm.fit, lm.fit2)
# Since the p-value is very small, adding the quaratic fit can largely increase
# the performance of the model.


par(mfrow=c(2,2))
plot(lm.fit2)  # Still a little discernible pattern in the residuals.

# Perform the linear regression with polynomial term up to fifth order. i.e.,
# Y = beta_0 + beta_1(X) + beta_2(X)^2 +...+ beta_5(X)^5 + epsilon.
lm.fit5 = lm(medv ~ poly(lstat, 5), data = Boston)
summary(lm.fit5)  # R-sq = 0.6817

# Perform a log transformation.
summary(lm(medv ~ log(rm), data = Boston))  # R-sq = 0.4358


## 3.6.6 Qualitative Predictors ------------------------------------------------

# ******************************************************************************
# Example. Need to load the Carseats data set.
#          - Number of observations: 400.
#          - Response: Sales
#          - Predictor: 10 predictors such as CompPrice, Income,...,etc.
# 
#          The Carseat contains some qualitative predictors.
#          - Shelveloc: Bad/Medium/Good
#          - Urban: Yes/No
#          - US: Yes/No
# ******************************************************************************

Carseats = read.csv(file = "Lab2_Carseats.csv", header = T)
fix(Carseats)
summary(Carseats)
names(Carseats)

# Perform a simple linear regression includes 2 interaction terms.
lm.fit=lm(Sales ~ . + Income:Advertising + Price:Age, data = Carseats)

# Now, check the ShelveLocGood/ShelveLocMedium, UrbanYes, and USYes. 
# What is your interpretation?
summary(lm.fit)


# ******************************************************************************
# Example. Use the contrasts() function to check the dummy variables
# ******************************************************************************

attach(Carseats)

# The as.factor() function treat the character variable as a factor variable.
contrasts(as.factor(ShelveLoc))
contrasts(ShelveLoc)  # This will cause an error.
?Contrasts()


## 3.6.7 Writing Functions -----------------------------------------------------

# ******************************************************************************
# Example. Try to create a function to load ISLR and MASS libraries and return
#          a message.
# ******************************************************************************

LoadLibraries()  # There is no such function by default.

# Lets create the function.
LoadLibraries = function(){
    library(ISLR)
    library(MASS)
    print("The libraries have been loaded.")
}


LoadLibraries
LoadLibraries()  # Now, the libraries are successfully loaded.


