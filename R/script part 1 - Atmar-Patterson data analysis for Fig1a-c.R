rm(list = ls())

library(CooccurrenceAffinity)
library(ggplot2)

# the following functions comb(), cooc() and veech_p() were copied from https://github.com/giovannistrona/co_occurrence
comb <- function(n, k){
  if (k > n) return(0)
  if (suppressWarnings(is.infinite(factorial(n))) || suppressWarnings(is.infinite(factorial(k)))){
    x <- gmp::factorialZ(n) %/% gmp::factorialZ(k) %/% gmp::factorialZ(n - k)
    x <- as.numeric(x)
    if(is.infinite(x)) stop("n and/or k is too large for the internal comb function; there are likely too many species", "\n",
                            "in the dataset for calculations to be possible")
  } else{
    x <- factorial(n) / factorial (k) / factorial(n - k)
  }
  return(x)
}

cooc <- function(j, N, N1, N2){
  a <- comb(N1, j)
  b <- comb( (N - N1), (N2 - j))
  d <- comb(N, N2)
  return( (a * b) / d)
}

veech_p<-function(N,N1,N2,obs){
  p<-0
  for (i in obs:min(N1,N2)){
    pc<-cooc(i,N,N1,N2)
    p<-p+pc
  }
  return (list('p'=p))
}



# load Atmar_Patterson data copied from https://github.com/giovannistrona/co_occurrence
load("data/Atmar_Patterson_Matrices.RData")
ap_data<-data.Atmar

res_veech<-c()
# the following for loop and some codes elsewhere were copied from https://github.com/giovannistrona/co_occurrence
for (m in 1:length(ap_data)){
  mat<-ap_data[[m]]
  R<-nrow(mat)
  C<-ncol(mat)
  fill<-sum(mat)/(R*C)
  if ((min(R,C)>5)&&((R*C)<500)){
    site_n<-ncol(mat)
    amle_vals<-c()
    prop<-c()
    veech<-c()
    for (i in 1:nrow(mat)){
      for (j in 1:nrow(mat)){
        if (i<j){

          # compute species pair affinity
          A<-mat[i,]
          B<-mat[j,]
          C <- A + B
          N1<-sum(A)
          N2<-sum(B)
          X <- sum(C==2) ##total number of sites for both species
          if (N1*N2>0){ # && max(N1,N2)/site_n<0.8 && min(N1,N2)/site_n>0.2){
            veech<-c(veech,veech_p(site_n,N1,N2,X)$p)
            prop<-rbind(prop,c(N1,N2,X,site_n))
            alpha_mle_list <- ML.Alpha(x=X, c(sum(A),sum(B),site_n),lev=0.9)
            if (alpha_mle_list != "Degenerate co-occurrence distribution!"){
              alpha_mle <- alpha_mle_list$est
              amle_vals<-c(amle_vals,alpha_mle)} else {amle_vals<-c(amle_vals,NaN)
              }
          }
        }
      }
    }
    res_veech<-rbind(res_veech,cbind(prop,amle_vals,veech))
  }
  print (m)
}

head(res_veech)
nrow(res_veech)
colnames(res_veech)<-c('N1','N2','X','site_n','alphaMLE','Veech_p')
res_veech<-res_veech[is.finite(rowSums(res_veech)),]

head(res_veech)
nrow(res_veech)

# the first species pair in the res_veech data (row #1) has the following:
# N1 = 7, N2 = 16, X = 6, sites_n = 18
#      N1 N2 X site_n   alphaMLE    Veech_p
# [1,]  7 16 6     18 -0.4814054 0.86274510

# for this case:
phyper(6, 7, 18-7, 16) # cumulative prob of up to 6 co-occurrences, including 6
1-phyper(6, 7, 18-7, 16) # 1-cumulative prob computed above OR the prob of seeing more than 6 co-occurrences
# to compute the probability of observed co-occurrences (6) and higher co-occurrences, substract 1 from the co-occurrences
1-phyper((6-1), 7, 18-7, 16)
# this is what p_gt is in veech's cooccur function, reported by Ulrich et al as pv or Veech's probability





# Ulrich et al created a subset of res_veech to show some behaviors that they report in their commentary
# ------------------------------------------------------------------------------------------------------

par(mfrow=c(3,4))

table(res_veech[,4])

res_5<-res_veech[res_veech[,5]>5,]
head(res_5)
nrow(res_5)

plot(res_5[,4],res_5[,5],pch=16,cex=1.5,xlab='site number',ylab='Affinity',cex.lab=1.2,cex.axis=1.2,las=1,cex.main=1.5)

diffvec <- rep(NA, nrow(res_5))
r <- 14
for(r in 1:nrow(res_5)) {
  mydiff <- res_5[r,3] - min(res_5[r,1], res_5[r,2])
  diffvec[r] <- mydiff
}
diffvec
summary(diffvec)
# this proves X is min of N1 and N2 for all cases of res_5

summary(res_5[,5] - log(2*res_5[,4]^2))
# this proves that alpha MLE in res_5 was log(2N^2); let's explore this more

res_5 <- cbind(res_5, capped = NA)
for(r in 1:nrow(res_5)) {
  mycap <- log(2*(res_5[r,4]^2))
  if(res_5[r,5] >= mycap) {
    res_5[r,7] <- mycap
  }
}
head(res_5)
plot(res_5[,4], res_5[,5],pch=16,cex=1.5,xlab='site number',ylab='Affinity',cex.lab=1.2,cex.axis=1.2,las=1,main='(a)',cex.main=1.5)
points(res_5[,4], res_5[,7], pch=1,cex=2, col='red')
# this proves affinity = log(2N^2) in the res_5 dataset

plot(res_veech[,6],res_veech[,5],ylab='Affinity',xlab=expression("p"[V]),cex.lab=1.2,cex.axis=1.2,las=1,cex.main=2,pch=16,ylim=c(-8,8))
points(res_5[,6], res_5[,5],pch=16,col='red')
# this proves affinity=5 used by Ulrich et al used an approximate cutoff to separate truncated affinity values
# when infinite positive values were returned by equation


# let's do it in a more proper way
# --------------------------------
head(res_veech)
dim(res_veech)

# create a subset of res_veech named res_veech_posinf below for species-pairs with
# the habitat occupancy of one species being a subset of another
# min(N1, N2) == X; this is when affinity is positively infinite and capped to log(2N^2)
res_veech_posinf <- res_veech[0,]
# and, another subset named res_veech_neginf for species-pairs with
# the habitat occupancy being larger of 0 or N1+N2-site_n
# this is when affinity is negatively infinite and capped to -log(2N^2)
res_veech_neginf <- res_veech[0,]

# save rest of the points (non-extreme cooccurrences) separately
res_veech_noninf <- res_veech[0,]

for(i in 1:nrow(res_veech)) {
  if(min(res_veech[i,1:2]) == res_veech[i,3]){
    res_veech_posinf <- rbind(res_veech_posinf, res_veech[i,])
  } else if(max(0, (res_veech[i,1]+res_veech[i,2]-res_veech[i,4])) == res_veech[i,3]){
    res_veech_neginf <- rbind(res_veech_neginf, res_veech[i,])
  } else {
    res_veech_noninf <- rbind(res_veech_noninf, res_veech[i,])
  }
}

head(res_veech_posinf)
dim(res_veech_posinf)

head(res_veech_neginf)
dim(res_veech_neginf)

head(res_veech_noninf)
dim(res_veech_noninf)

head(res_veech)
dim(res_veech)

nrow(res_veech) == nrow(res_veech_posinf) + nrow(res_veech_neginf) + nrow(res_veech_noninf)


# pie chart of number of species-pairs in each of the three categories
piedata <- data.frame(
  group=c("noninf", "posinf", "neginf"),
  value=c(nrow(res_veech_noninf), nrow(res_veech_posinf), nrow(res_veech_neginf))
)
piedata$label <- paste0(round(piedata$value/sum(piedata$value)*100), "%")

pie <- ggplot(piedata, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(x=1.1, label = label), position = position_stack(vjust = 0.5), size=3, col='white', fontface='bold') +
  scale_fill_manual(values=c("purple", "blue", "red")) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position = "none")
pie

mypiecols <- c("#00BA38", "#619CFF", "#F8766D")
pie <- ggplot(piedata, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(x=1.1, label = label), position = position_stack(vjust = 0.5), size=2.5, col='black') +
  scale_fill_manual(values=mypiecols) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position = "none")
pie



par(mfrow=c(3,4))

# ============= Fig 1a =============
plot(res_veech_noninf[,6],res_veech_noninf[,5],ylab='Affinity',xlab=expression("p"[V]),cex.lab=1.2,cex.axis=1.2,las=1,cex.main=2,pch=16,ylim=c(-8,8), col='blue')
points(res_veech_posinf[,6], res_veech_posinf[,5],pch=16,col='red')
points(res_veech_neginf[,6], res_veech_neginf[,5],pch=16,col='purple')


# ============= Fig 1b =============
# for all species pairs with min(N1, N2) == X, affinity is equal to log(2N^2),
# exhibiting a deterministic log relation with N, number of sites
res_veech_posinf <- res_veech_posinf[order(res_veech_posinf[,4]),]
plot(res_veech_posinf[,4],res_veech_posinf[,5],ylab='Affinity',xlab='Number of sites (N)',cex.lab=1.2,cex.axis=1.2,las=1,cex.main=2,pch=16, col='red')
lines(res_veech_posinf[,4], log(2*res_veech_posinf[,4]^2), col='green3', lwd=2)
points(res_veech_posinf[,4],res_veech_posinf[,5],pch=16, col='red')
text(23.5, 6, expression(log(2*N^2)), col="green3", cex=1.2)
segments(23.5, 6.15, 23.5, 7)




# USSG missed some points in the separate cluster with their cutoff of 5.
# for this entire cluster res_veech_posinf (positive affinity truncated), compute lower Blaker CI end point
out_mat_posinf = t(apply(res_veech_posinf, 1, function(row) {
  tmp = ML.Alpha(row[3], c(row[1],row[2], row[4]))
  c(tmp$est, tmp$CI.Blaker[1])} ))
head(out_mat_posinf)
colnames(out_mat_posinf) <- c("alpha_mle", "Blaker_CI_lowerpt")
nrow(out_mat_posinf)
## A much more sensible summary of values when X = min(mA,mB) is the lower Blaker CI endpt
plot(res_veech_posinf[,6], out_mat_posinf[,2], xlab="Veech_p", ylab="Lower point of Blaker CI of Alpha MLE")



# for the cluster of points with negative affinity truncated, compute upper Blaker CI end point
out_mat_neginf = t(apply(res_veech_neginf, 1, function(row) {
  tmp = ML.Alpha(row[3], c(row[1],row[2], row[4]))
  c(tmp$est, tmp$CI.Blaker[1])} ))
head(out_mat_neginf)
colnames(out_mat_neginf) <- c("alpha_mle", "Blaker_CI_upperpt")
nrow(out_mat_neginf)



# ============= Fig 1c =============

plot(res_veech_posinf[,6], out_mat_posinf[,2], xlab="Veech_p", ylab="End point of Blaker CI of Alpha MLE", col='red',
     ylim = range(c(out_mat_posinf[,2], out_mat_neginf[,2])), xlim=c(0,1))
points(res_veech_neginf[,6], out_mat_neginf[,2], col='purple')

