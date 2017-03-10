#Comments are always preceded by a # sign, this tells the R interpreter to
# ignore anything following it on the same line

#To execute a command press ctrl + enter
#New objects you create will appear in the upper right box under the 
# "Environment" tab.  You can click on them to see more information

install.packages(c("tidyverse", "memisc")) #This line installs two libraries (packages) that you will need

library(tidyverse)# this loads several very useful packages so that you can easily call their functions

#The following line loads the spss data set.  
#You may need to point R to the right directory if 
getwd() #if this does not show the same directory your data is in
#then you will need to set the correct directory using 
#setwd("C:/THe directory and parent directories of your data")
data <- memisc::as.data.set( memisc::spss.system.file("College coping data (complete).sav"))

#This line converts the "data.set" object to a much more useful data.frame object
df <- as.data.frame(data)

#Now you have a data frame with rows of observations and columns of variables...
