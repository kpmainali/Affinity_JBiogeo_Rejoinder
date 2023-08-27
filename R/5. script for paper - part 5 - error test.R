library(BiasedUrn)
library(CooccurrenceAffinity)


rm(list = ls())

# the following three functions should be obtained from the Science Advances paper SI
# (under "Auxiliary Supplementary Materials and Other Supporting Files").
# dowloaded from https://www.science.org/doi/10.1126/sciadv.abj9204
# the functions have also been uploaded to this GitHub repo.

# AlphInt <- dget("R/ScienceAdvances_functions/AlphInt function.R")
# AlphMLE <- dget("R/ScienceAdvances_functions/AlphMLE.R")
# NewAlph <- dget("R/ScienceAdvances_functions/NewAlph.R")

# Alternatively, you can directly access the functions from the GitHub
devtools::source_url("https://raw.githubusercontent.com/kpmainali/Affinity_JBiogeo_Rejoinder/main/R/ScienceAdvances_functions/AlphInt%20function.R")
devtools::source_url("https://raw.githubusercontent.com/kpmainali/Affinity_JBiogeo_Rejoinder/main/R/ScienceAdvances_functions/AlphMLE.R")
devtools::source_url("https://raw.githubusercontent.com/kpmainali/Affinity_JBiogeo_Rejoinder/main/R/ScienceAdvances_functions/NewAlph.R")

ls()

# USSG indicated correctly in their GitHub that following computation fails.
NewAlph(4, mA=7, mB=10, N=20)$AlphMLE
# Error in pFNCHypergeo(t, mA, N - mA, mB, exp(alp[i])) :
#   Inconsistency. mean = 16, lower limit = 0, upper limit = 7

# USSG used this function from supplementary material of our Science Advances paper.
# The source of this error was the BiasedUrn package that we used for various computations of Extended Hypergeometric
# (also called Fisher noncentral hypergeometric) distribution.
# We did not encounter this problem while analyzing data for our Science Advances paper (Mainali et al., 2022).
# This problem from the BiasedUrn package was later resolved in our package revision on July 12, 2022.
# When USSG published their commentary in October 2022, the package had resolved this problem for several months.
# see below:
ML.Alpha(4, c(7, 10, 20))


# USSG provided R code that iterates through various combination of input values
# to check the stability of our function available as supplement to Science Advances paper.
# We performed ten runs of their code (with some corrections) simulating input data for 100K times in each run.
# We found that the indicated error happened 5-13 times in 100,000 input values.
# In the run that we saved, this problem happened 7 times.

# As we show below, our package function does NOT fail for those input values.

# erdf <- read.csv("data/error_parameters_locked.csv"); erdf

# Alternatively, directly load it from GitHub
url <- "https://raw.githubusercontent.com/kpmainali/Affinity_JBiogeo_Rejoinder/main/data/error_parameters_locked.csv"
erdf <- read.csv(url)


nrow(erdf)
cranfunction <- sciadvfunction <- c()

# alpha estimate from CRAN function
for(i in 1:nrow(erdf)) {
  X <- erdf$cooccurrence[i]
  mA <- erdf$entityAfreq[i]
  mB <- erdf$entityBfreq[i]
  N <- erdf$totalN[i]

  cranfunction <- c(cranfunction, ML.Alpha(X, c(mA, mB, N))$est)
  if(i == nrow(erdf)) print(cranfunction)
}

# alpha estimate from Science Advances function
for(i in 1:nrow(erdf)) {
  X <- erdf$cooccurrence[i]
  mA <- erdf$entityAfreq[i]
  mB <- erdf$entityBfreq[i]
  N <- erdf$totalN[i]

  sciadvfunction <- c(sciadvfunction, NewAlph(X, mA, mB, N)$AlphMLE)
  if(i == nrow(erdf)) print(sciadvfunction)
}

# compare the alpha outputs from CRAN vs Science Advances functions
df <- cbind(cranfunction, sciadvfunction); df

