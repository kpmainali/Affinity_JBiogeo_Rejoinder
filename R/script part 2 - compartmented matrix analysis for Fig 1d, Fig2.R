# load the libraries and some functions at the beginning of script part 1 first


# An artificial matrix of 22 species occupying various number of 50 sites
# -----------------------------------------------------------------------

# compartmented matrix used in USSG's Fig 1B, C

par(mfrow=c(3,4))

mat_file <- 'data/compartmented.csv'
mat <- t(read.csv(mat_file,header=T,row.names=1)) # 50 x 22 -- Commentary says 22 species, 50 sites
head(mat)


head(mat) #compartmented.csv supplied by authors
mA = apply(mat,2,sum); mA        ## species site-counts range from 18:39
Sp.Count = apply(mat,1,sum); Sp.Count        ## sites   species-counts range from 6:20

## Commentary article states that over all pairs of species,
## counts of sites occupied by at least  one species ranges from 27 to 49

# compute alpha and other quantities for the compartmented.csv table
mat.anal = array(0, c(22,22,9), dimnames=list(paste0("Spec1.",1:22),
                                              paste0("Spec2.",1:22),c("Sites-Occ","X","mA","mB","Aff.MLE",
                                                                      "Aff.Lo","Aff.Hi","pval2sid","p_gt")))
head(mat.anal)
head(mat)
i <- 1
j <- 2
for(i in 1:21) for(j in (i+1):22) {
  X = sum(mat[,i]*mat[,j])
  tmp = ML.Alpha(X,c(mA[i],mA[j],50))
  mat.anal[i,j,] = c(sum(mat[,i]+mat[,j]>0), X, mA[i], mA[j], tmp$est,
                     tmp$CI.Blaker, tmp$pval, 1-phyper(X-1,mA[i],50-mA[i],mA[j])) }
head(mat.anal)
aux = c(mat.anal[,,1]); aux
range(aux[aux>0])
# [1] 27 49               ## OK, as stated in USSG Commentary; most between 36,45
hist(aux[aux>0])

sum(mat.anal[,,3]>0 & (mat.anal[,,2] == pmin(mat.anal[,,3],mat.anal[,,4]) |
                         +          mat.anal[,,2] == pmax(0, mat.anal[,,3]+mat.anal[,,4]-50)))
# [1] 0                   ## OK, no extreme X's this time

# make a single dataframe from a list of several square ones above
choose(22,2)
Exmp1bc = array(0, c(231,9), dimnames=list(NULL,dimnames(mat.anal)[[3]]))
ctr=0
for(i in 1:21) {
  Exmp1bc[ctr+(1:(22-i)),] = mat.anal[i,(i+1):22,]
  ctr=ctr+22-i }

head(Exmp1bc)
dim(Exmp1bc)
class(Exmp1bc)
Exmp2bc <- data.frame(Exmp1bc)

hist(Exmp2bc$Aff.MLE)
table(cut(Exmp2bc$Aff.MLE, breaks = c(-10,-1,1.5,10)))
table(cut(Exmp2bc$Aff.MLE, breaks = c(-10,-1,1.5,10)))/length(Exmp2bc$Aff.MLE)
# 87% of the pairs have affinity between -1 and 1.5

plot(Exmp2bc$p_gt, Exmp2bc$Aff.MLE, xlab="p_gt", ylab="Affinity")

plot(Exmp2bc$Sites.Occ, Exmp2bc$Aff.MLE)   ## Affinity versus occupied sites
### does show definite quadratic bend


# Let's generate four sets of X's as medians of Extended Hypergeometric with same mA[i],mA[j],N
# as in USSG's Compartmented Matrix but with alpha values successively taken as log(2)=0.69,
# log(2.75)=1.01, log(3.5)=1.25, log(4.25)=1.45. The picture, with points for medians
# from different alphas plotted in different colors, shows that Affinity does track nearly
# along a simple curve for X's from each alpha. For smaller alpha, Affinity and p_gt show
# similar ability to discriminate between X values in the tail of the hypergeomnetric.
# However, for larger alpha, Affinity shows progressively greater ability to discriminate
# between values far out in the tail of the hypergeometric distribution that are compressed together by p_gt.


## Now for 4 different alpha values, generate X's as median using qFNCHypergeo
## using same mA{i],mA[j],N as in Exmp1bc
Xarr = array(0, c(231,4,3), dimnames=list(NULL, paste0("exp.alph",c(2,2.75,3.5,4.25)),
                                          c("Xmed","p_gt","alphMLE")))
head(Xarr)
for(i in 1:231) {
  mA=Exmp1bc[i,3]; mB=Exmp1bc[i,4]
  Xarr[i,,1] = c(qFNCHypergeo(0.5,mA,50-mA,mB,2),
                 qFNCHypergeo(0.5,mA,50-mA,mB,2.75),
                 qFNCHypergeo(0.5,mA,50-mA,mB,3.5),
                 qFNCHypergeo(0.5,mA,50-mA,mB,4.25))
  Xarr[i,,2] = 1-phyper(Xarr[i,,1]-1,mA,50-mA,mB)
  for(j in 1:4) Xarr[i,j,3] = ML.Alpha(Xarr[i,j,1],c(mA,mB,50))$est }


# ============= Fig 1d =============
colrs=c("black","blue","red","orange")
plot(1,2, xlab="Veech_p", ylab="Affinity", type="n", xlim=c(0,0.4), ylim=c(0.3,2),main=
       paste0("Plots of Affinity vs p_gt for X's from Different Alphas","\n",
              "(X = median for Alpha, Prevalence as in Compartmented Matrix)"))
for(j in 1:4) points(Xarr[,j,2], Xarr[,j,3], pch=18,col=colrs[j])
legend(0.2,2, legend=paste0("alpha=log(",c("2","2.75","3.5","4.25"),")"),
       pch=rep(18,4), col=colrs, cex=1.1)



#------------------------------------------------------------------
## scatterplotting  Affinity vs null-standardized cooccurrence in the original compartmented data with moderate affinity

null.mean = Exmp1bc[,3]*Exmp1bc[,4]/50
null.sdev = sqrt(null.mean*(50-Exmp1bc[,3])*(50-Exmp1bc[,4])/(50*49))
plot( (Exmp1bc[,2]-null.mean)/null.sdev,Exmp1bc[,5], xlab="null-standardized Co-occ",ylab="Affinity",
      main="Affinity vs Null-Stdizd X \nin Artificial Compartmented-Matrix Data")


# ============= Fig 2a =============
## Affinity vs null-standardized cooccurrence
## in the settting of the four-color picture with X's generated as medians
## according to Extended Hypergeometric with four different alpha's

plot(0,0, xlab="Null-standardized Co-Occ", ylab="Affinity",
     type="n", xlim=c(0.5,3), ylim=c(0.5,2), main=paste0(
       "Affinity vs Null-Stdizd X in Compartmented-Matrix Data","\n",
       "Generated as Median Etended-Hypergeometric with 4 alpha's"))
for(j in 1:4)
  points((Xarr[,j,1]-null.mean)/null.sdev, Xarr[,j,3], pch=18,col=colrs[j])
legend(0.5,2, legend=paste0("alpha=log(",c("2","2.75","3.5","4.25"),")"),
       pch=rep(18,4), col=colrs, cex=1)


# "Number of Sites" in USSG's Fig.1c was mA[i]+mA[j]-X, given in array Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,,1]
# and # "Affinity" is given in Xarr[,,3]. We plot these next:

plot(c(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,,1]),c(Xarr[,,3]), xlab="Number of sites",
     ylab="Affinity", type="n", ylim=c(0.25, 3) )
for(i in 1:4) points(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,i,1],Xarr[,i,3],
                     pch=18, col=colrs[i])

## For each of these 4 sets of colored points, calculate and plot least-squares quadratic
lmlist = lmlist2 = NULL
for(i in 1:4) {
  nsites = Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,i,1]
  lmlist = c(lmlist, list(lm(Xarr[,i,3] ~ nsites + I(nsites^2))))
  lmlist2 = c(lmlist2, list(lm(Xarr[,i,3] ~ nsites)))
  coefs = lmlist[[i]]$coef
  ##  coefs2 = lmlist2[[i]]$coef
  newnsit = sort(nsites)
  lines(newnsit, coefs[1]+coefs[2]*newnsit+coefs[3]*newnsit^2,
        lty=2, col=colrs[i], lwd=1.5)
  ##   lines(newnsit, coefs2[1]+coefs2[2]*newnsit, lwd=1.3)
}

legend(26,3, legend=paste0("alpha=",round(log(c(2.0,2.75,3.5,4.25)),2),
                           ", pval =",sapply(lmlist, function(lmobj)
                             round(summary(lmobj)$coef[3,4],4))),  lty=2, col=colrs)


# ============= Fig 2c =============

plot(c(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,,1]),c(Xarr[,,3]), xlab="Number of sites",
     ylab="Affinity", type="n", ylim=c(0.25, 3) )
for(i in 1:4) points(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,i,1],Xarr[,i,3],
                     pch=18, col=colrs[i])

## For each of these 4 sets of colored points, calculate and plot least-squares quadratic
lmlist = lmlist2 = NULL
for(i in 1:4) {
  nsites = Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,i,1]
  lmlist = c(lmlist, list(lm(Xarr[,i,3] ~ nsites + I(nsites^2))))
  lmlist2 = c(lmlist2, list(lm(Xarr[,i,3] ~ nsites)))
  coefs = lmlist[[i]]$coef
  ##  coefs2 = lmlist2[[i]]$coef
  newnsit = sort(nsites)
  lines(newnsit, coefs[1]+coefs[2]*newnsit+coefs[3]*newnsit^2,
        lty=2, col=colrs[i], lwd=1.5)
  ##   lines(newnsit, coefs2[1]+coefs2[2]*newnsit, lwd=1.3)
}

legend(26,3, legend=paste0(
  "pval=",sapply(lmlist, function(lmobj) round(summary(lmobj)$coef[3,4],4)),
  ", R2=", sapply(lmlist, function(lmobj) round(summary(lmobj)$adj.r.squared,4))),  lty=2, col=colrs)




# -------------------------------------------------------------------------------------
#            SEE THE FOLLOWING TEXT FOR DISCUSSION OF PLOTS OF AFFINITY
#            VERSUS NUMBER OF OCCUPIED SITES AS GIVEN IN USSG FIG 1(C)
# -------------------------------------------------------------------------------------

# First of all, note that the Number of Sites on the x-axes of Fig1(a) and of Fig1(c) of USSG commentary mean somewhat different things.
# In Fig.1(a), the number of sites is simply the table-total N in the 2x2 species-pair contingency table of occurrences,
# i.e., the total number of ecologic sites at which occurrence of species A and B are recorded.
# On the other hand, in Fig.1(c), the number of sites is explained in text to be the total number of distinct sites occupied by at least one species.
# That is, in the 2x2 table notation where X is the number of co-occurrences, mA and mB the total numbers of sites
# at which respectively species A and species B are finite, the number of occupied sites is mA + mB - X.
# The number of co-occurrences is subtracted to avoid double-counting of sites occupied by both species among mA + mB.

# Now we know that when mA, mB and N are fixed, Affinity is an increasing function of X.
# Therefore, if a dataset of 2x2 tables is generated by fixing species prevalences mA, mB,
# and total number N of examined ecologic sites and generating co-occurrence counts X independently
# from the Extended Hypergeometric distribution, we should expect that Affinity and Number of Occupied Sites
# will show a definite negative correlation and no particular curvature.
# On the other hand, real ecologic datasets consist of 2x2 tables for all species-pairs among the species
# examined at the ecologic sites in a specific geographic location.
# These species have joint occurrences governed by complicated biological processes
# and so the (Affinity, Number of Occupied Sites) data-pairs generated from all species-pairs should by no means be assumed independent.
# The patterns of correlation or curvature in plots like Fig.1(c) cannot be validly analyzed using least-squares and p-values
# for coefficients computed as though these points were independent.

# The artificial X-data used to create the figures plotting Affinity versus p_gt or versus Xstd
# was also not a set of points generated according to the assumptions underlying linear-regression software.
# The X's in those figures were generated as Extended-Hypergeometric medians with fixed values of mA,mB,N,
# and 4 different alpha's. As we show below, linear-regression to fit
# quadratic curves to Affinity versus Number of Occupied Sites non-significant and/or meaningless (close to zero) curvature.
# This is still an artifact of the (nonstochastic!) way the points were generated.
# Contrast this with the corresponding plots when the X values are pseudo-randomly drawn from
# Extended Hypergeometric (mA,mB,N, exp(alpha)) distributions with the same set of four alphas.
# That figure shows a clear negatively sloped relationshiop between
# Affinity and Number of Occupied sites, but still shows some artifactual curvature
# because of the very restricted range of Affinity values in the examples generated with
# largest and smallest numbers of occupied sites. But note that all of these figures would not be said
# to exhibit curvatures if the points with largest and smallest numbers of occupied sites were excluded from the dataset.



## Generate a new array of meanX, sdX, meanJ, sdJ from these same data

StdzArr  = array(0, c(231,4), dimnames=list(NULL, c("meanX","sdX","meanJ","sdJ")))
for(i in 1:231) {
  mA=Exmp1bc[i,3]; mB=Exmp1bc[i,4]; N=50
  Xlo = max(0, mA+mB-N); Xhi = min(mA,mB)
  meanX = mA*mB/N
  sdX = sqrt(meanX*(N-mA)*(N-mB)/(N*(N-1)))
  meanJ = sum(dhyper(Xlo:Xhi, mA,N-mA,mB)*((Xlo:Xhi)/(mA+mB-(Xlo:Xhi))))
  sdJ = sqrt( sum(
    dhyper(Xlo:Xhi, mA,N-mA,mB)*((Xlo:Xhi)/(mA+mB-(Xlo:Xhi))-meanJ)^2 ) )
  StdzArr[i,] = c(meanX,sdX,meanJ,sdJ) }



# ============= Fig 2b =============

plot( (Exmp1bc[,"X"]-StdzArr[,1])/StdzArr[,2],
      (Exmp1bc[,"X"]/(Exmp1bc[,"mA"]+Exmp1bc[,"mB"]-Exmp1bc[,"X"])-
         StdzArr[,3])/StdzArr[,4],  xlab="Stdzd X", ylab="Stdzd J",
      main="Standardized J vs Standardized X, Compartmented-Matrix Data")



# ============= Fig supplement =============

plot(0,0, xlab="Null-standardized Jaccard", ylab="Affinity",
     type="n", xlim=c(0.5,3.2), ylim=c(0.5,2), main=paste0(
       "Affinity vs Null-Stdizd J in Compartmented-Matrix Data","\n",
       "Generated as Median Extended-Hypergeometric with 4 alpha's"))
for(j in 1:4)
  points((Xarr[,j,1]/(Exmp1bc[,"mA"]+Exmp1bc[,"mB"]-Xarr[,j,1])-StdzArr[,3])/
           StdzArr[,4], Xarr[,j,3], pch=18, col=colrs[j])
legend(0.5,2, legend=paste0("alpha=log(",c("2","2.75","3.5","4.25"),")"),
       pch=rep(18,4), col=colrs, cex=1)
