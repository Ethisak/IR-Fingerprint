library(readxl)
library(gdata)
library(rcdk)
library(proxy)
library(corrplot)
library(stringr)

d <- read_xlsx("Database.xlsx")
dt <- read.csv("FP-Theoretical.csv")
dt <- dt[2:102]

All <- 1:nrow(d)


cluster_split <- function(x,y,z){
  
  xx <- c()
  for(i in 1:length(y)){
    xx <- c(xx, which(dt[,y[i]] == 1))
  }
  
  if("TRUE" %in% duplicated(xx)){
    xx <- xx[-which(duplicated(xx) == "TRUE")]
  }else{}
  
  xxx <- c(x,xx)
  positive <- sort(xxx[which(duplicated(xxx) == "TRUE")])
  negative <- sort(x[!x %in% positive])
  
  assign(paste("Remaining-", z, sep = ""),negative,.GlobalEnv)
  assign(paste(z),positive,.GlobalEnv)
}


cluster_split2 <- function(x,y,z){
  
  xx <- 1:101
  xx <- xx[!xx %in% c(1,y)]
  yy <- c()
  
  for(i in x){
    bits <- sum(dt[i,xx])
    if (bits == 0) {
      yy <- c(yy,i)
    }else{}
  }
 
  positive <- yy
  negative <- x[!x %in% yy]
  
  assign(paste("Remaining-", z, sep = ""),negative,.GlobalEnv)
  assign(paste(z),positive,.GlobalEnv)
}

#Acids
cluster_split(x = All, y = c(7,47,81,75,76,77), z = "Organic-Acids")

cluster_split(x = `Organic-Acids`, y = c(7,47,81), z = "Carboxylic-Acids")

mv("Remaining-Carboxylic-Acids","Sulfones-Sulfonates")



#Acids derivatives

cluster_split(x = `Remaining-Organic-Acids`, y = c(44,82,43,42,9,10,48,49,62,63,86,87,64,90), z = "Acid-Derivatives")
rm(`Remaining-Organic-Acids`)

cluster_split(x = `Acid-Derivatives`, y = c(44,82), z = "Esters")

cluster_split(x = `Remaining-Esters`, y = 43, z = "Lactones")
rm(`Remaining-Esters`)

cluster_split(x = `Remaining-Lactones`, y = 42, z = "Acyl-Chloride")
rm(`Remaining-Lactones`)

cluster_split(x = `Remaining-Acyl-Chloride`, y = c(9, 10, 48, 49, 62, 63, 86, 87), z = "Amides")
rm(`Remaining-Acyl-Chloride`)

cluster_split(x = `Remaining-Amides`, y = 64, z = "Lactams")
rm(`Remaining-Amides`)

cluster_split(x = `Remaining-Lactams`, y = 90, z = "Thioamides")
rm(`Remaining-Lactams`,`Remaining-Thioamides`)


#Nitriles and isocyanate

cluster_split(x = `Remaining-Acid-Derivatives`, y = c(29,30,34,36), z = "Nitriles-Isocyanate")
rm(`Remaining-Acid-Derivatives`)

cluster_split(x = `Nitriles-Isocyanate`, y = c(29,30), z = "Nitriles")

cluster_split(x = `Remaining-Nitriles`, y = c(34), z = "Isocyanate")
rm(`Remaining-Nitriles`)

cluster_split(x = `Remaining-Isocyanate`, y = c(36), z = "Isothiocyanate")
rm(`Remaining-Isocyanate`,`Remaining-Isothiocyanate`)


#Aldehydes

cluster_split(x = `Remaining-Nitriles-Isocyanate`, y = c(27,45), z = "Aldehydes")
rm(`Remaining-Nitriles-Isocyanate`)


#Ketones and sulfoxides

cluster_split(x = `Remaining-Aldehydes`, y = c(46,78,79), z = "Ketones-Sulfoxides")
rm(`Remaining-Aldehydes`)

cluster_split(x = `Ketones-Sulfoxides`, y = c(46), z = "Ketones")

cluster_split(x = `Remaining-Ketones`, y = c(78,79), z = "Sulfoxides")
rm(`Remaining-Ketones`,`Remaining-Sulfoxides`)


#Alcohol and thiole

cluster_split(x = `Remaining-Ketones-Sulfoxides`, y = c(4,5,28,83), z = "Alcohols-Thiols")
rm(`Remaining-Ketones-Sulfoxides`)

cluster_split(x = `Alcohols-Thiols`, y = c(4,5,83), z = "Alcohols")

cluster_split(x = `Remaining-Alcohols`, y = 28, z = "Thiols")
rm(`Remaining-Alcohols`,`Remaining-Thiols`)


#Amines

cluster_split(x = `Remaining-Alcohols-Thiols`, y = c(11,12,60,61,88,89,65,6,66,94), z = "Amines-Oximes")
rm(`Remaining-Alcohols-Thiols`)

cluster_split(x = `Amines-Oximes`, y = c(11,12,13,14,60,61,88,89), z = "Amines")

cluster_split(x = `Remaining-Amines`, y = 65, z = "Imines")

cluster_split(x = `Remaining-Imines`, y = c(6,66,94), z = "Oximes")
rm(`Remaining-Amines`,`Remaining-Imines`,`Remaining-Oximes`)


#Other acid substituents

cluster_split(x = `Remaining-Amines-Oximes`, y = c(57,58,93), z = "Acid-Substituents")
rm(`Remaining-Amines-Oximes`)

cluster_split(x = `Acid-Substituents`, y = c(57,58), z = "Nitro-Compounds")

cluster_split(x = `Remaining-Nitro-Compounds`, y = 93, z = "Aliphatic-Phosphates")
rm(`Remaining-Aliphatic-Phosphates`,`Remaining-Nitro-Compounds`)


#Halocarbons

cluster_split(x = `Remaining-Acid-Substituents`, y = c(91,98,100,101), z = "Halocarbons")
rm(`Remaining-Acid-Substituents`)

cluster_split(x = Halocarbons, y = 98, z = "Chlorinated-Hydrocarbons")

cluster_split(x = `Remaining-Chlorinated-Hydrocarbons`, y = 100, z = "Brominated-Hydrocarbons")

cluster_split(x = `Remaining-Brominated-Hydrocarbons`, y = 101, z = "Iodated-Hydrocarbons")

cluster_split(x = `Remaining-Iodated-Hydrocarbons`, y = 91, z = "Fluorinated-Hydrocarbons")
rm(`Remaining-Iodated-Hydrocarbons`,`Remaining-Brominated-Hydrocarbons`,`Remaining-Chlorinated-Hydrocarbons`,`Remaining-Fluorinated-Hydrocarbons`)


#Hydrocarbons

cluster_split2(x = `Remaining-Halocarbons`, y = c(15,16,17,18,19,20,21,22,23,24,25,26,39,40,51,52,53,55,67,68,69,70,71,72), z = "Hydrocarbons")


cluster_split(x = `Hydrocarbons`, y = c(20,21,22,23,24), z = "Aromatic-Heterocycles")

cluster_split(x = `Remaining-Aromatic-Heterocycles`, y = c(19,67), z = "Aromatic-Rings")

cluster_split(x = `Remaining-Aromatic-Rings`, y = c(18,53,69), z = "Cyclic-Hydrocarbons")

cluster_split(x = `Remaining-Cyclic-Hydrocarbons`, y = c(16,17,39,40,52,72), z = "Unsaturated-Hydrocarbons")

cluster_split(x = `Remaining-Unsaturated-Hydrocarbons`, y = c(25,26,70,71), z = "Saturated-Hydrocarbons")
rm(`Remaining-Halocarbons`,`Remaining-Aromatic-Heterocycles`,`Remaining-Aromatic-Rings`,
   `Remaining-Cyclic-Hydrocarbons`,`Remaining-Saturated-Hydrocarbons`,`Remaining-Unsaturated-Hydrocarbons`)

mv("Remaining-Hydrocarbons", "Remaining-Compounds")
