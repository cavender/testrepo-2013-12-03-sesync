#Routine to calculate PD and FD from sample communities in FAB
#March 26, 2015, written by J. Cavender-Bares

rm(list = ls(all = TRUE))

library(ape)
library(vegan)
library(picante)
library(plotrix)
library(geiger)
require(xlsx)

setwd("~/Dropbox/Cedar Creek FAB/Analysis/DiversityMetrics/")

#Read a phylogeny in here

phy<-read.tree("~/Dropbox/SESYNC.Macroevolution.ES.Trees/Phylos/Tank(Zanne)Tree/Vascular_Plants_rooted.dated.tre")
#phy<-read.tree("Phylogeny/CC_PD_phylo.new")
#phy2<-read.tree("Phylogeny/tank.fia.tre")

#Check the species in the phylogeny
phy$tiplabels

#Read in plots to analyze
dat<-read.csv("~/Dropbox/NSF_NASA-DoB-remote-sensing/Data/Biomass/2014.e120.comm.abund.csv")
dat<-na.omit(dat)
dat<-data.frame(dat)
com<-sample2matrix(dat)

write.csv(com, "~/Dropbox/NSF_NASA-DoB-remote-sensing/Data/Biomass/2014.e120.comm.abund.t.csv")

#dat<-read.csv("Communities/FAB.com.csv")
#dat.abund<-read.csv("Communities/FAB.com.abund.csv")

#Read in traits
#alldat <- read.csv("Traits/Trait_list3.csv", header = TRUE)
alldat <- read.csv("~/Dropbox/NSF_NASA-DoB-remote-sensing/Data/Traits/JCB.traits.Townsend.E120.csv", header = TRUE)
head(alldat)
sla<-alldat[,c(1,5)]
ht<-alldat[,8]
trait.dat<-na.omit(alldat)
#trait.dat<-na.omit(alldat[,-10])

trait.dat<-sla

#Make plot data into a data frame and remove first column
#For presence/absence
cmat<-data.frame(dat)
cmat<-cmat[,2:13]

#For abundance values
cmat.abund<-data.frame(com)
#cmat.abund<-data.frame(dat.abund)
#cmat.abund<-cmat.abund[,2:13]

#Prune phylogeny to match plots
pphy<-prune.sample(cmat,phy)
pphy<-prune.sample(cmat, phy)
pphy2<-prune.sample(cmat, phy2)

#Calculate phylogenetic diversity metrix
##PSV; Matt Helmus phylogenetic species variability. He has another one that uses abundance, PSE
dist<-cophenetic(pphy)
psvobs <- psv(cmat.abund, pphy) 
pseobs<-pse(cmat.abund, pphy)
	
#Faith's pd
faith <- pd(cmat.abund, pphy, include.root=FALSE)
faith.nr2 <- pd(cmat, pphy, include.root=FALSE)

mat<-cbind(faith, psvobs,pseobs)
write.csv(mat, "~/Dropbox/NSF_NASA-DoB-remote-sensing/Data/Biomass/2014.e120.pd.csv")
##MPD Webb's mean phylogenetic distance
nri <- ses.mpd(cmat.abund, dist, null.model = "sample.pool", abundance.weighted = TRUE, runs = 1000) #could also use "sample.pool"; results are probably identical
mpd<-nri$mpd.obs

#MNTD Webb's mean nearest taxon distance
nti<-ses.mntd(cmat.abund, dist, null.model = "phylogeny.pool", abundance.weighted = TRUE, runs = 1000) 
mntd<-nti$mntd.obs

## now with Tank FIA phylogeny 

##PSV and PSE
dist2<-cophenetic(pphy2)
psvobs2 <- psv(cmat, pphy2) 
pseobs2<-pse(cmat.abund, pphy2)
	
#Faith's pd
faith2 <- pd(cmat, pphy2)
faith.nr2 <- pd(cmat, pphy2, include.root = FALSE)	

##MPD
nri2 <- ses.mpd(cmat, dist2, null.model = "sample.pool", abundance.weighted = FALSE, runs = 10)
mpd2<-nri2$mpd.obs

#####Routine for calculating functional diversity from the random communities above
##Normalize traits
##multivariate
##input phylogeny as the dendrogram of traits

#Prune Community data to match traits
keep <- intersect(trait.dat[,1], colnames(cmat.abund))
cmat.traits <- cmat.abund[,keep]
traitmat<-data.frame(cmat.traits)

#standardize trait data to mean = 0 and var = 1
traits<-scale(trait.dat[,2:ncol(trait.dat)]) 

#hierarchical clustering as per Petchey & Gaston
trait.clust<-hclust(dist(traits), method="average") 
#label cluster objects as species names -plot by "plot(trait.clust)"
trait.clust$labels<-keep[trait.clust$order] 
#transform trait.clust into phylo object
trait.phy <- as.phylo(trait.clust) 

#prune phylogeny so only species with traits in the communities are used
prunedtraitphy<-prune.sample(traitmat,trait.phy)

#create cophenetic trait distance matrix
trait.dist = cophenetic.phylo(prunedtraitphy)

#create comdist matrix
#comdist.traits<-comdist(traitmat,trait.dist)

#psv phylogenetic significant variation for traits
psvtrait <- psv(traitmat, prunedtraitphy) ##PSV
psetrait<-pse(traitmat, prunedtraitphy) ##PSE
nti_trait <- ses.mntd(traitmat, trait.dist, null.model = "sample.pool", abundance.weighted=FALSE, runs=10) ## NTI
nri_trait <- ses.mpd(traitmat, trait.dist, null.model = "sample.pool", abundance.weighted=FALSE, runs=10)  ## NRI
TD_trait<-taxondive(traitmat, trait.dist, match.force=FALSE)

mat<-data.frame(cbind(psvtrait, psetrait))
write.csv(mat, "~/Dropbox/NSF_NASA-DoB-remote-sensing/Data/Traits/2014.e120.all.csv")

dd<-data.frame(cbind(psvobs$PSVs,faith$PD,faith.nr$PD,psvobs2$PSVs,faith2$PD,faith.nr2$PD, psvtrait$PSVs, nri_trait$mpd.obs, nti_trait$mntd.obs,TD_trait$Dplus))

colnames(dd)<-c("psv","FaithPD","FaithPD.nr","psv.bl=1","FaithPD.bl=1","FaithPD.nr.bl=1","PSV-trait","mpd-trait","nti-trait","Dplus")

#D <- data.frame(psvtrait$PSVs, nri_trait$mpd.obs, nti_trait$mntd.obs,TD_trait$Dplus)

#colnames(D) <- c("PSV-trait","mpd-trait","nti-trait","Dplus")

#write(D, file = "Cedar Creek/Traits/Output/Trait_Diversity_no_ferns2.2.11.txt",append = FALSE, sep = "\t", ncolumns = 126)

write.csv(dd, "put path here")

