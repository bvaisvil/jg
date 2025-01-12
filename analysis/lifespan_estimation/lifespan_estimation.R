###############################################
###                                         ###
###           Lifespan estimation           ###
###   from CpG density with BLAST Promoters ###
###                                         ###
###############################################
# Mayne, B., Berry, O., Davies, C. et al. A genomic predictor of lifespan in vertebrates. Sci Rep 9, 17866 (2019). https://doi.org/10.1038/s41598-019-54447-w
### Load the required libraries
library(data.table)
library(Biostrings)
library(dplyr)
### Read in the coefficients for the promoters into the linear model
datPromoters <- read.csv("Estimate_Lifepan_AllSpecies.csv", row.names = 1)
args = commandArgs(trailingOnly=TRUE)
#####################
### PGLS Function ###
#####################

pgls <- function(x, class){
  
  ### x is the raw prediction value (log)
  ### class can either be Aves, Fish, Mammalia, Reptilia
  
  ### -4.38996 + 2.57328x + ax + b (PGLS adjustment formula)
  
  if(class == "Aves"){
    
    lifespan = exp(-4.38996 + (2.57328 * x) + (-0.90323 * x) + 2.14857)
    
  } else if(class == "Fish"){
    
    lifespan = exp(-4.38996 + (2.57328 * x) + (2.14632 * x) + -6.58228)
    
  } else if(class == "Mammalia"){
    
    lifespan = exp(-4.38996 + (2.57328 * x) + (-0.92888 * x) + 2.33508)
    
  } else if(class == "Reptilia"){
    
    lifespan = exp(-4.38996 + (2.57328 * x) + (-0.48958 * x) + 1.17281)
    
  }
  return(lifespan)
}

estimate_Lifespan <- function(dat){
  # sum of the product of coefficient weights multiplied by the respective CpG densities and coefficient intercept
  dat$Product <- dat$Coefficient * dat$CG_Density;
  return(sum(dat$Product) + 3.16778948197575)
}


### Read in the BLAST Output and determine the CpG density
datAnimal <- read.table(args[1])
datAnimal$Length <- "Length"
datAnimal$CGs <- "CGs"

for(i in 1:nrow(datAnimal)){
  
  datAnimal$Length[i] <- nchar(as.character(datAnimal[i,10]))
  datAnimal$CGs[i] <- dinucleotideFrequency(DNAString(datAnimal[i,10]))["CG"]
}

datAnimal$CG_Density <- as.numeric(datAnimal$CGs) / as.numeric(datAnimal$Length)
row.names(datAnimal) <- datAnimal$V1
datAnimal$Promoters <- datAnimal$V1
datAnimalSubset <- subset(datAnimal, datAnimal$Promoters %in% datPromoters$Promoters[-1]) # subset data to those in the LE promoter list.
x <- left_join(datAnimalSubset, datPromoters, by="Promoters") # join promoters so we can compute using their coeffecients
write.table(as.data.frame(x), "out.tab", sep="\t")
print(as.numeric(exp(estimate_Lifespan(as.data.frame(x)))))  
print(pgls(x = as.numeric(estimate_Lifespan(as.data.frame(x))), class = args[2])) # class needs to be changed to fit input data.

