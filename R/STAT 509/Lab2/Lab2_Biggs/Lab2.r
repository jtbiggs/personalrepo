#Justin Biggs - Lab 2 Submission

#set working directory
setwd("/home/justinbiggs/Documents/Code/personalrepo/R/STAT 509/Lab2")

#import table data
SEATS=read.csv("Lab2_Carseats.csv",header=T, na.strings="?")

lm.fit=lm(Sales ~ CompPrice+Income+ShelveLoc+ShelveLoc*CompPrice, data=SEATS)

plot(lm.fit)
