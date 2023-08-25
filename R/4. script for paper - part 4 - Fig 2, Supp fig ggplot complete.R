ls()
library(CooccurrenceAffinity)


# Points for revision of J Biogeog Commentary
# ===========================================
#
# lines 48-49: The referee would like to see p-value plots replaced by plots in which
#    p-values are transformed by qnorm(). We need a way to highlight the asymmetries
#    that this produces in  some examples.

#---------------------------------------------------------------

StdIndx = function(x, marg, tol=1e-2) {
          mu = marg[1]*marg[2]/marg[3]
          sigsq = mu*(marg[3]-marg[1])*(marg[3]-marg[2])/(marg[3]*(marg[3]-1))
         (x - mu)/sqrt(ifelse(sigsq>0, sigsq, tol)) }


## Idea:
#  first generate an ensemble of (mA, mB, N) triples, or better: fix (N, mA) for
#    #sites and a reference species; then generate a set of mB  values
#    (for species B always wrt species A, NOT to be compared with one another);
#    next generate 110 alpha values at -4, -3.2, -2.4, -1.6, -0.8, 0, 0.8, ..., 4
#    in clumps of 10 and random X with ExtHyp(mA,mB,N, exp(alpha)).
#    Plot standardized values,  qnorm  p-values,  and alph.hat


#--------------------------------------------------------
# Notes to focus revision
#========================================================

# Suggestion for exhibit: find an array of (4 ?) (mA,mB,N)
# configurations exhibiting a variety of mA < mB, mA > mB and mA and mB approximately equal, for each of two alpha values such as 1 and 2.5,
# and then overlay pictures for each of alph.hat, zstand and qnor.pv.
#
# NB caption will have to explain that these interpolated curves are not
# probability densities: it is the probability spikes (the vertical
# ordinates) and not the areas that sum to 1 for each curve.


margarr = rbind(c(60,20,100), c(55,35,120), c(35,55,90), c(20,60,120))
margarr

colr2 = cbind(c("cyan","blue","gray","brown"),
                c("red","orange","green","darkorange"))
colr2



MetricPMF_twoalph  =  function(marg, alpha, plt=T, cap.qpv = 5.5) {
  ## function to calculate prob mass functions for 3 metrics
  # alph.hat, zstand = null-standardized X, qnor.pv = qnorm-transformed p-value
  # Current version 5/23/23 also has a cap on the absolute qnor.pv
  mA=marg[1]; mB=marg[2]; N=marg[3]
  require(BiasedUrn)
  Xrng = max(mA+mB-N,0):min(mA,mB)
  probs = dFNCHypergeo(Xrng, mA, N-mA, mB, exp(alpha))
  inds = which(probs > (0.001)*(0.1^(N>100)))
  Xrng = Xrng[inds]
  probs = probs[inds]
  zstand = alph.hat = qnor.pv = numeric(length(Xrng))
  for(i in 1:length(Xrng)) {
    zstand[i] = StdIndx(Xrng[i], marg)
    qnor.pv[i] = qnorm(phyper(Xrng[i], mA, N-mA, mB))
    alph.hat[i] = ML.Alpha(Xrng[i], marg)$est }
  ## cap abs qnor.pv at cap.qpv
  qnor.pv = pmax(-cap.qpv, pmin(cap.qpv, qnor.pv))
  if(plt) {
    par(mfrow=c(3,1))
    plot(alph.hat,probs, type="b", xlab="alph.hat",
         ylab="prob", main="ExtHyp Prob Masses, versus alph.hat")
    abline(v=alpha, col="blue")
    plot(zstand,probs, type="b", xlab="zstand",
         ylab="prob", main="ExtHyp Prob Masses, versus zstand")
    plot(qnor.pv,probs, type="b", xlab="qnor.pv",
         ylab="prob", main="ExtHyp Prob Masses, versus qnor.pv") }
  par(mfrow=c(1,1))
  list(alph.hat = alph.hat, zstand=zstand, qnor.pv=qnor.pv,
       probs=probs, Xrng=Xrng, inds=inds) }



avalC = c(1.5,3)
margarrC = rbind(c(80,120,300), c(150,90,300), c(90,150,240), c(120,80,240))
calc.listC = NULL
ym.arrC = array(0,c(4,2))
xint.arrC = array(0, c(4,2,2,3))
for (j in 1:2) for(i in 1:4)  {
  aux = MetricPMF_twoalph(margarrC[i,], avalC[j], plt=F, cap.qpv=10)
  aux$affinity <- avalC[j]
  aux$margscenario <- margarrC[i,]
  calc.listC = c(calc.listC, list(aux))
  ym.arrC[i,j] = max(aux$probs)
  for(k in 1:3) xint.arrC[i,j,,k] = range(aux[[k]]) }

ymC = max(c(ym.arrC)) # was ym.arr

labs = c("alph.hat","zstand", "qnor.pv", "qnor.veech") # line copied from above

par(mfrow=c(3,1))
for(k in 1:3) {
  xint = range(xint.arrC[,,,k])
  plot(0,0.01,type="n", xlab=labs[k], ylim=c(0,ymC), xlim=xint,
       ylab="prob", main=paste0("Extended Hypergeom Prob Masses ",
                                "versus ",labs[k],",  alpha's = ",avalC[1],", ",avalC[2]))
  for(j in 1:2) for(i in 1:4)  {
    newlst = calc.listC[[(j-1)*4+i]]
    lines(newlst[[k]], newlst[[4]], type="b",
          lty=i, col=colr2[i,j], cex=1.5)
    if(k==1) abline(v=avalC[j], col="purple")  }
}


# plot alph.hat vs prob
# -----------------------
# merge a specific element of all 8 examples in a dataframe
length(calc.listC)
mergedf <- data.frame()
for (i in 1:length(calc.listC)) {
  calc.listC[[i]]
  rm(temp, myscenario)
  myscenario <- paste0("(", paste(calc.listC[[i]]$margscenario, collapse = ","), ")"); myscenario
  temp <- data.frame(cbind(alph.hat = calc.listC[[i]]$alph.hat, probs = calc.listC[[i]]$probs,
                           affinity = calc.listC[[i]]$affinity, scenario = myscenario))
  mergedf <- rbind(mergedf, temp)
}
head(mergedf)
table(mergedf$affinity)
table(mergedf$scenario)

mergedf$alph.hat <- as.numeric(mergedf$alph.hat)
mergedf$probs <- as.numeric(mergedf$probs)
mergedf$affinity <- factor(mergedf$affinity)
mergedf$scenario <- factor(mergedf$scenario)
sapply(mergedf, class)
fig2.alpha <- ggplot(mergedf, aes(x=alph.hat, y=probs, color=scenario, shape=affinity)) +
  geom_point() +
  geom_line(aes(group = interaction(scenario, affinity))) +
  # geom_smooth(aes(group = interaction(scenario, affinity)), se=FALSE) +
  xlab("Affinity") + ylab("Probability mass") +
  ylim(c(0,0.23)) +
  geom_vline(xintercept = c(1.5, 3), color = c("red", "blue"), size = 0.25) +
  labs(color = "(mA, mB, N)", shape = "alpha") +
  theme(legend.position = c(0.85, 0.55), legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(1, "lines"))

fig2.alpha <- fig2.alpha +
  geom_text(x = 1.5, y = 0.175, label = "alpha\n= 1.5", color="red", angle = 90, size=3) +
  geom_text(x = 3, y = 0.215, label = "alpha\n= 3", color="blue", angle = 90, size=3)
fig2.alpha


# plot zstand vs prob
# -----------------------
# merge a specific element of all 8 examples in a dataframe
length(calc.listC)
mergedf <- data.frame()
for (i in 1:length(calc.listC)) {
  calc.listC[[i]]
  rm(temp, myscenario)
  myscenario <- paste0("(", paste(calc.listC[[i]]$margscenario, collapse = ","), ")"); myscenario
  temp <- data.frame(cbind(zstand = calc.listC[[i]]$zstand, probs = calc.listC[[i]]$probs,
                           affinity = calc.listC[[i]]$affinity, scenario = myscenario))
  mergedf <- rbind(mergedf, temp)
}
head(mergedf)
table(mergedf$affinity)
table(mergedf$scenario)

mergedf$zstand <- as.numeric(mergedf$zstand)
mergedf$probs <- as.numeric(mergedf$probs)
mergedf$affinity <- factor(mergedf$affinity)
mergedf$scenario <- factor(mergedf$scenario)
sapply(mergedf, class)
fig2.zstand <- ggplot(mergedf, aes(x=zstand, y=probs, color=scenario, shape=affinity)) +
  geom_point() +
  geom_line(aes(group = interaction(scenario, affinity))) +
  xlab("Null-standardized co-occurrence count") + ylab("Probability mass") +
  # geom_vline(xintercept = c(1.5, 3), color = "red", size = 0.25) +
  labs(color = "(mA, mB, N)", shape = "alpha") +
  theme(legend.position = c(0.125, 0.575), legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        legend.key.size = unit(1, "lines"))
  # theme(legend.position = "none")
fig2.zstand




# merge a specific element of all 8 examples in a dataframe
length(calc.listC)
mergedf <- data.frame()
for (i in 1:length(calc.listC)) {
  calc.listC[[i]]
  rm(temp, myscenario)
  myscenario <- paste0("(", paste(calc.listC[[i]]$margscenario, collapse = ","), ")"); myscenario
  temp <- data.frame(cbind(qnor.pv = calc.listC[[i]]$qnor.pv, probs = calc.listC[[i]]$probs,
                           affinity = calc.listC[[i]]$affinity, scenario = myscenario))
  mergedf <- rbind(mergedf, temp)
}
head(mergedf)
table(mergedf$affinity)
table(mergedf$scenario)

mergedf$qnor.pv <- as.numeric(mergedf$qnor.pv)
mergedf$probs <- as.numeric(mergedf$probs)
mergedf$affinity <- factor(mergedf$affinity)
mergedf$scenario <- factor(mergedf$scenario)
sapply(mergedf, class)
fig2.qnor.pv <- ggplot(mergedf, aes(x=qnor.pv, y=probs, color=scenario, shape=affinity)) +
  geom_point() +
  geom_line(aes(group = interaction(scenario, affinity))) +
  xlab("qnor.pv") + ylab("Probability") +
  # geom_vline(xintercept = c(1.5, 3), color = "red", size = 0.25) +
  theme(legend.position = "none")
fig2.qnor.pv




## ==================================================================
## Further code for exhibit:

## Try to plot graph of median statistic value as a function of alpha
##  for a range of alphas, for the three statistics alpha.hat, zstand, qnor.pv

MedCurv.Alph = function(alphvec, marg) {
  ## function to calculate prob mass functions for 3 metrics
  ## alph.hat, zstand = null-standardized X, qnor.pv = qnorm-transformed p-value
  mA=marg[1]; mB=marg[2]; N=marg[3]
  require(BiasedUrn)
  nv = length(alphvec)
  meds = numeric(nv)
  for(i in 1:nv)
    meds[i] = qFNCHypergeo(0.5, mA, N-mA, mB, exp(alphvec[i]))
  zstand = alph.hat = qnor.pv = qnor.veech = numeric(nv)
  for(i in 1:nv) {
    zstand[i] = StdIndx(meds[i], marg)
    qnor.pv[i] = qnorm(phyper(meds[i], mA, N-mA, mB))
    alph.hat[i] = ML.Alpha(meds[i], marg)$est
    qnor.veech[i] = qnorm(veech_p(N, mA, mB, meds[i])$p) # Kumar's addition
  }
  list(alph.hat=alph.hat, zstand=zstand, qnor.pv = qnor.pv, qnor.veech = qnor.veech)  }



# now let's generate Fig 1a by replacing pv with qnorm(pv)
# -------------------------------------------------------------------

# ============= Fig S1 =============

head(res_veech_noninf_df); nrow(res_veech_noninf_df)
head(res_veech_posinf_df); nrow(res_veech_posinf_df)
head(res_veech_neginf_df); nrow(res_veech_neginf_df)

for(i in 1:nrow(res_veech_noninf_df)) {
    temp <- res_veech_noninf_df[i,]; temp
    out <- MedCurv.Alph(temp$alphaMLE, c(temp$N1,temp$N2,temp$site_n)); out
    res_veech_noninf_df$zstand[i] <- out$zstand
    res_veech_noninf_df$qnor.pv[i] <- out$qnor.pv
    res_veech_noninf_df$qnor.veech[i] <- out$qnor.veech
}


head(res_veech_noninf_df); tail(res_veech_noninf_df); nrow(res_veech_noninf_df)


head(res_veech_noninf_df)
figS1 <- ggplot(data=res_veech_noninf_df, aes(x=qnor.veech, y=alphaMLE)) +
  geom_point(col="blue", pch=1) +
  xlab(expression("qnorm-transformed p"[V])) + ylab("Affinity") +
  ylim(range(c(res_veech_neginf_df$alphaMLE, res_veech_posinf_df$alphaMLE)))
figS1


head(res_veech_posinf_df)

nrow(res_veech_noninf_df)
res_veech_noninf_df_unq <- unique(res_veech_noninf_df[c("qnor.veech", "alphaMLE")])
nrow(res_veech_noninf_df_unq)

myxlim <- c(res_veech_neginf_df$qnor.veech, res_veech_posinf_df$qnor.veech, res_veech_noninf_df_unq$qnor.veech)
myxlim <- myxlim[is.finite(myxlim)]

figS1.unq <- ggplot(data=res_veech_noninf_df_unq, aes(x=qnor.veech, y=alphaMLE)) +
  geom_point(col="blue", pch=1) +
  xlab(expression("qnorm-transformed p"[V])) + ylab("Affinity") +
  ylim(range(c(res_veech_neginf_df$alphaMLE, res_veech_posinf_df$alphaMLE))) +
  xlim(range(myxlim))
# fig1a.unq <- fig1a.unq +
#   geom_point(data=subset(res_veech_posinf_df_unq, is.finite(qnor.veech)), aes(x=qnor.veech, y=alphaMLE), col="red", pch=1) +
#   geom_point(data=subset(res_veech_neginf_df_unq, is.finite(qnor.veech)), aes(x=qnor.veech, y=alphaMLE), col="purple", pch=1)
figS1.unq



head(fourXdf)
hist(fourXdf$N)
fourXdf2 <- fourXdf
head(fourXdf2)
sapply(fourXdf2, class)

for(i in 1:nrow(fourXdf2)) {
  temp <- fourXdf2[i,]; temp
  out <- MedCurv.Alph(temp$alphMLE, c(temp$mA,temp$mB,temp$N)); out
  fourXdf2$zstand[i] <- out$zstand
  fourXdf2$qnor.pv[i] <- out$qnor.pv
  fourXdf2$qnor.veech[i] <- out$qnor.veech
}
head(fourXdf2)
tail(fourXdf2)

summary(fourXdf2$qnor.veech)
nrow(fourXdf2)
fourXdf2 <- subset(fourXdf2, is.finite(qnor.veech))
nrow(fourXdf2)

ggplot(fourXdf2, aes(x=qnor.veech, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab(expression("qnorm of p"[V])) + ylab("Affinity") +
  theme(legend.position = c(0.825, 0.8), legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "lines")) +
  guides(color = guide_legend(nrow = 2))
fig2a <- ggplot(fourXdf2, aes(x=qnor.veech, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab(expression("qnorm-transformed p"[V])) + ylab("Affinity") +
  # xlim(-2.5,0.75) + ylim(0.35,2.05) +
  theme(legend.position = "none")
fig2a <- fig2a +
  # geom_point(show.legend = FALSE) +
  directlabels::geom_dl(aes(label = cat2), method = "smart.grid", size=5)
fig2a

head(fourXdf2)
data_summary(data = fourXdf2, varnames = "alphMLE", groupnames = "cat")
log(2); log(2.75); log(3.5); log(4.25)


fig2a_alt <- ggplot(fourXdf2, aes(x=qnor.veech, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab(expression("qnorm-transformed p"[V])) + ylab("Affinity") +
  # xlim(-2.5,0.75) + ylim(0.35,2.05) +
  theme(legend.position = "none")
fig2a_alt


rm(patchwork)
patchwork = fig2a + fig2.zstand + fig2b + fig2.alpha + plot_layout(ncol=2, widths=c(1,2)); patchwork

# pdf("/Users/kpmainali/Dropbox/Documents/RESEARCH/Cooccurrence_SA_Paper_Commentary_Response/first revision/plots/fig2.pdf", width = 16/2.54*1.5, height = 5.35/2.54*1.5*2)
#   patchwork
# dev.off()

rm(patchwork)
patchwork = fig2a_alt + fig2.zstand + fig2b + fig2.alpha + plot_layout(ncol=2, widths=c(1,2)); patchwork

# pdf("/Users/kpmainali/Dropbox/Documents/RESEARCH/Cooccurrence_SA_Paper_Commentary_Response/first revision/plots/fig2_alt.pdf", width = 16/2.54*1.5, height = 5.35/2.54*1.5*2)
#   patchwork
# dev.off()



# ============= Fig S6 =============

head(fourXdf2)
figS6 <- ggplot(fourXdf2, aes(x=cooc_stz, y=qnor.veech, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab("Null-standardized co-occurrence") + ylab(expression("qnorm-transformed p"[V])) +
  theme(legend.position = "none")
figS6 <- figS6 +
  # geom_point(show.legend = FALSE) +
  directlabels::geom_dl(aes(label = cat2), method = "smart.grid", size=5)
figS6


rm(patchwork)
patchwork = figS1.unq + figS2 + figS3 + figS4 + figS5 + figS6 + plot_layout(ncol=3); patchwork

# pdf("/Users/kpmainali/Dropbox/Documents/RESEARCH/Cooccurrence_SA_Paper_Commentary_Response/first revision/plots/fig_supp_sixpanel.pdf",
#     width = 16/2.54*1.5, height = 10.5/2.54*1.5)
#   patchwork
# dev.off()
