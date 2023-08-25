# run the following script after running script for paper - part 1 and part 2
# because regenerating ggplot needs the objects created in those previous two files


# =================================================================================================================================
# ================================================== PUBLICATION QUALITY FIGURES ==================================================
# =================================================================================================================================

library(ggplot2)


# ============= Fig 1a =============

res_veech_noninf_df <- data.frame(res_veech_noninf)
res_veech_posinf_df <- data.frame(cbind(res_veech_posinf, out_mat_posinf[,2, drop=F]))
res_veech_neginf_df <- data.frame(cbind(res_veech_neginf, out_mat_neginf[,2, drop=F]))


head(res_veech_noninf_df)
fig1a <- ggplot(data=res_veech_noninf_df, aes(x=Veech_p, y=alphaMLE)) +
  geom_point(col="blue") +
  xlab(expression("Veech probability, p"[V])) + ylab("Affinity") +
  ylim(range(c(res_veech_neginf_df$alphaMLE, res_veech_posinf_df$alphaMLE)))
  # theme(text = element_text(size = 8))
fig1a <- fig1a +
  geom_point(data=res_veech_posinf_df, aes(x=Veech_p, y=alphaMLE), col="red") +
  geom_point(data=res_veech_neginf_df, aes(x=Veech_p, y=alphaMLE), col="purple")
fig1a


# check how many points are unique in the above plot because it makes a difference in a vector graphics to load
head(res_veech_posinf_df)
nrow(res_veech_posinf_df)
temp <- unique(res_veech_posinf_df[c("Veech_p", "alphaMLE")]); nrow(temp)
nrow(res_veech_neginf_df)
temp <- unique(res_veech_neginf_df[c("Veech_p", "alphaMLE")]); nrow(temp)
nrow(res_veech_noninf_df)
temp <- unique(res_veech_noninf_df[c("Veech_p", "alphaMLE")]); nrow(temp)


# the unique points in the two columns above are: 515/21791 in res_veech_posinf_df, 40/17439 in res_veech_neginf_df and
# 1069/8169 in res_veech_noninf_df. make a plot with the unique points for above plot, avoiding a complete overlay of points
fig1a.unq <- ggplot(data=unique(res_veech_noninf_df[c("Veech_p", "alphaMLE")]), aes(x=Veech_p, y=alphaMLE)) +
  geom_point(col="blue", pch=1) +
  xlab(expression("Veech probability, p"[V])) + ylab("Affinity") +
  ylim(range(c(res_veech_neginf_df$alphaMLE, res_veech_posinf_df$alphaMLE)))
fig1a.unq <- fig1a.unq +
  geom_point(data=unique(res_veech_posinf_df[c("Veech_p", "alphaMLE")]), aes(x=Veech_p, y=alphaMLE), col="red", pch=1) +
  geom_point(data=unique(res_veech_neginf_df[c("Veech_p", "alphaMLE")]), aes(x=Veech_p, y=alphaMLE), col="purple", pch=1)
fig1a.unq





# ============= Fig 1b =============
# for all species pairs with min(N1, N2) == X, affinity is equal to log(2N^2),
# exhibiting a deterministic log relation with N, number of sites

res_veech_posinf_df <- res_veech_posinf_df[order(res_veech_posinf_df$site_n),]

fig1b <- ggplot(data=res_veech_posinf_df, aes(x=site_n, y=alphaMLE)) +
  geom_line(aes(x=site_n, y=log(2*site_n^2)), col="green3", lwd=1) +
  geom_point(col='red') +
  xlab("Number of sites (N)") + ylab("Affinity")
fig1b

head(res_veech_posinf_df)
nrow(res_veech_posinf_df)
temp <- unique(res_veech_posinf_df[c("site_n", "alphaMLE")]); nrow(temp)

# 19/21791 unique

fig1b.unq <- ggplot(data=unique(res_veech_posinf_df[c("site_n", "alphaMLE")]), aes(x=site_n, y=alphaMLE)) +
  geom_line(aes(x=site_n, y=log(2*site_n^2)), col="green3", lwd=1) +
  geom_point(col='red') +
  xlab("Number of sites (N)") + ylab("Affinity")
fig1b.unq



# ============= Fig 1c =============

# fig1c <- ggplot(res_veech_noninf_df, aes(x=Veech_p, y=Blaker_CI_lowerpt)) +
#   geom_point(col='blue') +
#   xlab(expression("Veech probability, p"[V])) + ylab("End point of Blaker CI") +
#   ylim(range(c(res_veech_neginf_df$Blaker_CI_lowerpt, res_veech_posinf_df$Blaker_CI_lowerpt)))
#   # theme(text = element_text(size = 8))
# fig1c <- fig1c +
#   geom_point(data=res_veech_posinf_df, aes(x=Veech_p, y=Blaker_CI_lowerpt), col="red") +
#   geom_point(data=res_veech_neginf_df, aes(x=Veech_p, y=Blaker_CI_lowerpt), col="purple")
# fig1c


fig1c <- ggplot(res_veech_posinf_df, aes(x=Veech_p, y=Blaker_CI_lowerpt)) +
  geom_point(col='red') +
  xlab(expression("Veech probability, p"[V])) + ylab("Lower end point of Blaker CI") +
  xlim(0,1.0015)
fig1c

head(res_veech_posinf_df)
nrow(res_veech_posinf_df)
temp <- unique(res_veech_posinf_df[c("site_n", "Blaker_CI_lowerpt")]); nrow(temp)

# 544/21791 unique

fig1c.unq <- ggplot(unique(res_veech_posinf_df[c("Veech_p", "Blaker_CI_lowerpt")]), aes(x=Veech_p, y=Blaker_CI_lowerpt)) +
  geom_point(col='red', pch=1) +
  xlab(expression("Veech probability, p"[V])) + ylab("Lower end point of Blaker CI") +
  xlim(0,1.0015)
fig1c.unq


# create multipanel figure now
# ----------------------------
library(patchwork)
rm(patchwork)
patchwork = fig1a.unq + fig1b.unq + fig1c.unq + plot_layout(ncol=3); patchwork
patchwork + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# pdf("/Users/kpmainali/Dropbox/Documents/RESEARCH/Cooccurrence_SA_Paper_Commentary_Response/first revision/plots/fig1.pdf",
#     width = 16/2.54*1.5, height = 5.35/2.54*1.5*1.05)
#   patchwork
# dev.off()






# compute standardized X, alpha MLE and p_gt for the four schemes of drawing co-occurrences from different affinity levels
# ------------------------------------------------------------------------------------------------------------------------
head(Xarr)
data.frame(cbind(cat=colnames(Xarr)[1], cooc_stz=((Xarr[,1,"Xmed"]-null.mean)/null.sdev), alphMLE=Xarr[,1,"alphMLE"], p_gt=Xarr[,1,"p_gt"]))
fourXdf <- data.frame()
for(i in 1:4) {
  temp <- data.frame(cbind(cat = colnames(Xarr)[i],
                           cooc_stz=((Xarr[,i,"Xmed"]-null.mean)/null.sdev),
                           alphMLE=Xarr[,i,"alphMLE"],
                           p_gt=Xarr[,i,"p_gt"],
                           j_stz = (Xarr[,i,"Xmed"]/(Exmp1bc[,"mA"]+Exmp1bc[,"mB"]-Xarr[,i,1])-StdzArr[,3])/StdzArr[,4]),
                           mA=Exmp2bc$mA, mB=Exmp2bc$mB, N=50)

  # Exmp2bc$Sites.Occ is where at least one of the species was present, and so not total sites (N)
  # In the case of compartmented matrix, the total N is always 50

  fourXdf <- rbind(fourXdf, temp)
}
head(fourXdf)
summary(fourXdf)

for(i in 1:nrow(fourXdf)) {
  if(fourXdf$cat[i] == "exp.alph2") {
    fourXdf$cat[i] <- "log(2)=0.69"
    fourXdf$cat2[i] <- "log(2)"}
  if(fourXdf$cat[i] == "exp.alph2.75") {
    fourXdf$cat[i] <- "log(2.75)=1.01"
    fourXdf$cat2[i] <- "log(2.75)"}
  if(fourXdf$cat[i] == "exp.alph3.5") {
    fourXdf$cat[i] <- "log(3.5)=1.25"
    fourXdf$cat2[i] <- "log(3.5)"}
  if(fourXdf$cat[i] == "exp.alph4.25") {
    fourXdf$cat[i] <- "log(4.25)=1.45"
    fourXdf$cat2[i] <- "log(4.25)"}
}
head(fourXdf)
table(fourXdf$cat)
table(fourXdf$cat2)

fourXdf$cat <- factor(fourXdf$cat, levels = c("log(2)=0.69", "log(2.75)=1.01",  "log(3.5)=1.25", "log(4.25)=1.45"))
fourXdf$cat2 <- factor(fourXdf$cat2, levels = c("log(2)", "log(2.75)",  "log(3.5)", "log(4.25)"))
fourXdf$cooc_stz <- as.numeric(fourXdf$cooc_stz)
fourXdf$alphMLE <- as.numeric(fourXdf$alphMLE)
fourXdf$p_gt <- as.numeric(fourXdf$p_gt)
fourXdf$j_stz <- as.numeric(fourXdf$j_stz)
fourXdf$alphMLE <- as.numeric(fourXdf$alphMLE)


# define colors for ggplot
mycols <- c("#F8766D", "#7CAE00", "#E69F00", "#C77CFF")


# ============= Fig S2 =============

figS2 <- ggplot(fourXdf, aes(x=p_gt, y=alphMLE, col=cat)) +
  geom_point(aes(shape=cat)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab(expression("Veech probability, p"[V])) + ylab("Affinity") +
  # theme(text = element_text(size = 8)) +
  theme(legend.position = c(0.7225, 0.695), legend.title=element_blank())
figS2


# function to get summary of multiple response variables
data_summary <- function(data, varnames, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]]),
      count = length(x[[col]]),
      sd = sd(x[[col]]),
      se = sd(x[[col]])/sqrt(length(x[[col]])))
  }

  data_summall <- data.frame()
  for(i in varnames) {
    sub <- na.omit(subset(data, select = c(i, groupnames)))
    data_summ <- ddply(sub, groupnames, .fun = summary_func,
                       i)
    data_summ <- cbind(responsevar = i, data_summ)
    data_summ[c("responsevar", groupnames)] <- lapply(data_summ[c("responsevar", groupnames)], as.factor)
    data_summall <- rbind(data_summall, data_summ)
  }
  return(data_summall)
}

head(fourXdf)
data_summary(data = fourXdf, varnames = "alphMLE", groupnames = "cat")

log(2); log(2.75); log(3.5); log(4.25)



# ============= Fig S3 =============

## Let's try one more figure for Affinity vs null-standardized cooccurrence
## in the settting of the four-color picture with X's generated as medians
## according to Extended Hypergeometric with four different alpha's


ggplot(fourXdf, aes(x=cooc_stz, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  xlab("Null-standardized co-occurrence") + ylab("Affinity") +
  theme(legend.position = c(0.1, 0.9), legend.title=element_blank())

figS3 <- ggplot(fourXdf, aes(x=cooc_stz, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab("Null-standardized co-occurrence") + ylab("Affinity") +
  # theme(text = element_text(size = 8)) +
  theme(legend.position = "none")
figS3 <- figS3 +
  # geom_point(show.legend = FALSE) +
  directlabels::geom_dl(aes(label = cat2), method = "smart.grid", size=5)
figS3


# ============= Fig S4 =============

figS4 <- ggplot(fourXdf, aes(x=j_stz, y=alphMLE, col=cat2)) +
  geom_point(aes(shape=cat2)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  scale_color_manual(values = mycols) +
  xlab("Null-standardized Jaccard") + ylab("Affinity") +
  theme(legend.position = "none")
figS4 <- figS4 +
  # geom_point(show.legend = FALSE) +
  directlabels::geom_dl(aes(label = cat2), method = "smart.grid", size=5)
figS4



# ============= Fig S5 =============

head(Exmp1bc)
head(StdzArr)

xjplotdf <- data.frame(cbind(X.stdz = (Exmp1bc[,"X"]-StdzArr[,1])/StdzArr[,2],
                              J.stdz = (Exmp1bc[,"X"]/(Exmp1bc[,"mA"]+Exmp1bc[,"mB"]-Exmp1bc[,"X"])-StdzArr[,3])/StdzArr[,4]))
head(xjplotdf)

figS5 <- ggplot(xjplotdf, aes(x=X.stdz, y=J.stdz)) +
  geom_point(col="blue", pch=1) +
  xlab("Null-standardized co-occurrence") + ylab("Null-standardized Jaccard")
figS5



rm(patchwork)
patchwork = figS2 + figS3 + figS4 + figS5 + plot_layout(ncol=2); patchwork

# pdf("/Users/kpmainali/Dropbox/Documents/RESEARCH/Cooccurrence_SA_Paper_Commentary_Response/first revision/plots/fig_supp_fourpanel.pdf", width = 11/2.54*1.5, height = 10.5/2.54*1.5)
#   patchwork
# dev.off()



# ============= Fig 2b =============
colrs=c("black","blue","red","orange")
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


# make a better version of this figure
pvalvec <- sapply(lmlist, function(lmobj) round(summary(lmobj)$coef[3,4],4))
r2vec <- sapply(lmlist, function(lmobj) round(summary(lmobj)$adj.r.squared,4))

head(Exmp1bc)
head(Xarr)
data.frame(cbind(cat=colnames(Xarr)[1], sites_total=(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,1,1]), alphMLE=Xarr[,1,"alphMLE"]))
fourdf <- data.frame()
for(i in 1:4) {
  temp <- data.frame(cbind(cat=colnames(Xarr)[i], sites_total=(Exmp1bc[,3]+Exmp1bc[,4]-Xarr[,i,"Xmed"]), alphMLE=Xarr[,i,"alphMLE"]))
  fourdf <- rbind(fourdf, temp)
}
head(fourdf)
for(i in 1:nrow(fourdf)) {
  if(fourdf$cat[i] == "exp.alph2") fourdf$cat[i] <- "alpha=log(2)"
  if(fourdf$cat[i] == "exp.alph2.75") fourdf$cat[i] <- "alpha=log(2.75)"
  if(fourdf$cat[i] == "exp.alph3.5") fourdf$cat[i] <- "alpha=log(3.5)"
  if(fourdf$cat[i] == "exp.alph4.25") fourdf$cat[i] <- "alpha=log(4.25)"
}

table(fourdf$cat)
fourdf$cat <- factor(fourdf$cat, levels = c("alpha=log(2)", "alpha=log(2.75)",  "alpha=log(3.5)", "alpha=log(4.25)"))
fourdf$sites_total <- as.numeric(fourdf$sites_total)
fourdf$alphMLE <- as.numeric(fourdf$alphMLE)


ggplot(fourdf, aes(x=sites_total, y=alphMLE, col=cat)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=F) +
  xlab("Total number of sites for both species") + ylab("Affinity") +
  scale_color_manual(values = mycols) +
  theme(legend.position = c(0.1, 0.9), legend.title=element_blank())


fig2b <- ggplot(fourdf, aes(x=sites_total, y=alphMLE, col=cat)) +
  geom_point(aes(shape=cat)) +
  scale_shape_manual(values=c(1, 3, 2, 4))+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.65, se=F) +
  xlab("Total number of sites for both species") + ylab("Affinity") +
  scale_color_manual(values = mycols) +
  theme(legend.position = "none")
fig2b

g <- ggplot_build(fig2b)
unique(g$data[[1]]["colour"])
plotcols <- unique(g$data[[1]]["colour"])$colour

for(i in 1:4) {
  lb1 <- paste("italic(P) == ", pvalvec[i])
  lb2 <- paste("R^2 == ", r2vec[i])
  fig2b <- fig2b +
    annotate("text", x=26, y=1.85+i*0.15, label=lb1, parse=TRUE, col=plotcols[i], size=3.375, hjust=0) +
    annotate("text", x=31.5, y=1.85+i*0.15, label=";", col=plotcols[i], size=3.375, hjust=0) +
    annotate("text", x=32.25, y=1.865+i*0.15, label=lb2, parse=TRUE, col=plotcols[i], size=3.375, hjust=0)
}
fig2b
