library(INLA)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
devtools::install_github("becarioprecario/INLAMSM")
inla.upgrade()
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLAMSM)
library(maptools)
library(rgdal)
library(sp)
rm(list=ls())
data<-read.csv("C:/Users/fisol/OneDrive/Pictures/RIsk/data1.csv")
data$stunting
#data
#length(data)
#length(data$stunting1)
#length(data$state[2])
library(R2BayesX)
library(ggplot2)
library(tidyverse)
install.packages("tidyverse")
install.packages("R2BayesX")
install.packages("BayesX")
install.packages("spdep")
library(spdep)
map<-read.bnd("C:/Users/fisol/OneDrive/Pictures/RIsk/nigeria.bnd", sorted=FALSE)
# prepare for INLA Graphlibrary(R2
##temp <-poly2nb(bnd)#### not for this model
##nb2INLA("LDN.graph", temp) ### not for this model
library(INLA)
install.packages("INLA", dependencies = TRUE)
library(INLAMSM)
library(ggplot2)
library(rgdal)
library(maptools)
require(sp) 
require(grid)
require(Matrix)
require(foreach)
library("spdep")
require(spData)
require(sf)
library(fortify)
install.packages("fortify")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("mgcv")
getwd()
temp <-poly2nb(map)
H=inla.read.graph(filename="LDN.graph")
image(inla.graph2matrix(H),xlab="",ylab="")
LDN.adj <-paste("getwd()","/LDN.graph",sep="")
W <- inla.graph2matrix(H)
#### create states index (1-37, 38-74, 75-112)
N_idx<-c(data$state,data$state+37, data$state+74)
n=length(data$wasting)
d<-list(
  N_Urban=c(data$Urban, data$Urban, data$Urban),
  N_prim=c(data$prim, data$prim, data$prim),
  N_sec=c(data$sec,data$sec, data$sec),
  N_high=c(data$high, data$high, data$high),
  N_water=c(data$water,data$water, data$water),
  N_toilet=c(data$toilet,data$toilet, data$toilet),
  N_hhmember56=c(data$hhmember56,data$hhmember56, data$hhmember56),
  N_hhmember7=c(data$hhmember7,data$hhmember7, data$hhmember7),
  N_Poorer=c(data$Poorer,data$Poorer, data$Poorer),
  N_Middle=c(data$Middle, data$Middle, data$Middle),
  N_Richer=c(data$Richer, data$Richer, data$Richer),
  N_Richest=c(data$Richest, data$Richest, data$Richest),
  N_breastfeeding=c(data$breastfeeding, data$breastfeeding, data$breastfeeding),
  N_married=c(data$married, data$married, data$married),
  N_bord_23=c(data$bord_23,data$bord_23, data$bord_23),
  N_bord_4=c(data$bord_4, data$bord_4, data$bord_4),
  N_Male=c(data$Male,data$Male, data$Male),
  N_fever=c(data$fever, data$fever, data$fever),
  N_cough=c(data$cough, data$cough, data$cough),
  N_diarrhoea=c(data$diarrhoea, data$diarrhoea, data$diarrhoea),
  N_contraceptive=c(data$contraceptive, data$contraceptive, data$contraceptive),
  #### adding more fixed effect variables####
  N_child_age=c(data$child_age, data$child_age, data$child_age),
  N_mother_age=c(data$mother_age, data$mother_age, data$mother_age),
  ### Response stunting and underweight
  OBS=c(data$stunting,data$underweight,data$wasting)/100,
  idx=N_idx
)
d.add=1E-4
#HC.prior  = "expression:
# sigma = exp(-theta/2);
# gamma = 25;
# log_dens = log(2) - log(pi) - log(gamma);
# log_dens = log_dens - log(1 + (sigma / gamma)^2);
# log_dens = log_dens - log(2) - theta / 2;
# return(log_dens);
#"
#HT.prior = "expression:
# sigma = exp(-theta/2);
#  nu = 3;
#  log_dens = 0 - 0.5 * log(nu * pi) - (-0.1207822);
#  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
#  log_dens = log_dens - log(2) - theta / 2;
# return(log_dens);
#"
UN.prior = "expression:
log_dens = 0 - log(2) - theta / 2;
 return(log_dens);
"
#shyper= list(prec = list(prior = HC.prior))
 shyper=list(prec = list(prior = UN.prior))
#shyper=list(prec = list(prior = "tau.prior", param = c(1, 0.001)))
#shyper=list(prec = list(prior = "pcprior", param = c(1, 0.001)))
### Number of variables(i.e stunting {\xmlamp} underweight)
k<-3;e=rep(0,k)
A<-kronecker(Diagonal(k,1),Matrix(1,ncol= nrow(W),nrow= 1))
install.packages( c("your","list","of","packages"), type="binary")
A<-kronecker(Diagonal(k,1),Matrix(1,ncol=nrow(W),nrow=1))
# IMCAR model (wrapper)
a.min <-  0; a.max <- 1
b.min <- 0; b.max <- 1
alpha.min<-0;alpha.max<-1
model.imcar <- inla.IMCAR.model(k=k,W=W)
formula<-OBS ~ 1 +N_Urban+N_prim+N_sec+N_high+ #f(idx,model=model.imcar,extraconstr=list(A=as.matrix(A),e=e))
  N_water+N_toilet+N_hhmember56+N_hhmember7+
  N_Poorer+N_Middle+N_Richer+N_Richest+N_breastfeeding+N_married+N_bord_23+N_bord_4+
  N_Male+N_fever+N_cough+N_diarrhoea+N_contraceptive+
  f(N_child_age,model="rw2",scale.model=T,diagonal=d.add,hyper=shyper)+
  f(N_mother_age,model="rw2",scale.model=T,diagonal=d.add,hyper=shyper)+
  f(idx,model=model.imcar,extraconstr=list(A=as.matrix(A),e=e))
#Fit model
IMCAR <- inla(formula= formula,
 data = d,
 family= "gaussian", control.predictor =list(compute =TRUE),
 control.compute =list(dic =TRUE, waic =TRUE, cpo=TRUE),control.inla=list(cmin=0),verbose=TRUE)
     library(xtable)
xtable(IMCAR$summary.fixed[,c("mean","0.025quant","0.975quant")],digits=3)
par(mfrow = c(1, 2))
CH = IMCAR$summary.random$N_child_age
plot(CH$ID,CH$mean,type="l",ylim=c(-1,2.0), xlab="Child age in months", ylab="Posterior mean")
lines(CH$ID,CH$`0.025quant`,col="red")
lines(CH$ID,CH$`0.975quant`,col="red")
MTH=IMCAR$summary.random$N_mother_age
plot(MTH$ID, MTH$mean, type="l", ylim=c(-1,1.0), xlab="Mother age in months", ylab="Posterior Mean")
lines(MTH$ID,MTH$`0.025quant`,col="red")
lines(MTH$ID,MTH$`0.975quant`,col="red")

SP = IMCAR$summary.random$i
SP$ID
path = "C:/Users/fisol/OneDrive/Pictures/RIsk/nigeria-lgas/"
new_lga_nigeria_2003.shp


# plot map
require(rgdal)
require(ggplot2)
library(ggpubr)
library(readxl)

#d2=mean(data$state==2)
#d2

shp <- readOGR(dsn = paste0(path,"new_lga_nigeria_2003.shp"), stringsAsFactors = F)

shp.map <- fortify(shp, region = "STATE")

shp.map$estimate=0
#shpProperty = read.csv("C:/Users/fisol/OneDrive/Pictures/RIsk/poverty/poverty/NIGBND.csv")
state2 <- c("Abia.1",        "Abuja.1" , "Adamawa.1",     "Akwa Ibom.1" ,  "Anambra.1"  ,   "Bauchi.1",      "Bayelsa.1",    
            "Benue.1",       "Borno.1" ,      "Cross River.1", "Delta.1",       "Ebonyi.1",      "Edo.1",         "Ekiti.1",      
            "Enugu.1"   ,    "Gombe.1"  ,     "Imo.1",         "Jigawa.1",      "Kaduna.1",    "Kano.1" ,       "Katsina.1",    
            "Kebbi.1",       "Kogi.1"  ,      "Kwara.1"  ,     "Lagos.1" ,    "Nassarawa.1" , 
            "Niger.1",       "Ogun.1" ,       "Ondo.1" ,       "Osun.1" ,       "Oyo.1"    ,     "Plateau.1"   ,  "Rivers.1",     
            "Sokoto.1",      "Taraba.1",      "Yobe.1",        "Zamfara.1")
cloc2 <- c(22,37,8,1,2,3,31,5,6,7,23,32,4,33,24,34,9,25,10,11,12,26,27,13,14,35,15,16,17,
           28,18,19,20,21,29,30,36)
data.frame(state2,cloc2)
est_data <- data.frame(region=state2,
                       est= IMCAR$summary.random$i$mean[1:37][cloc2])



est_data <- data.frame(region=state2,
                       est=data$statemean[1:37][cloc2])

hist(data$state[12])

plot((data$state[3]), data$stunting)


s <- function(region,est){
  shp.map$estimate[which(shp.map$group==as.character(region))] <- est
  return(shp.map)
}

# Mean 

for (i in 1:37) {
  shp.map <-s(est_data[i,1],est_data[i,2]) 
}
#shp.map$estimate[which(shp.map$group=="Lagos.2")] <- est_data[25,2]
#shp.map$estimate[which(shp.map$group=="Lake.1")] <- est_data[6,2]
#shp.map$estimate[which(shp.map$group=="Borno.2")] <-est_data[6,2]
#shp.map$count=ifelse(shp.map$count>= 0.1,0.1,shp.map$count)
ggplot(shp.map, aes(x = long, y = lat, group = group, fill = estimate)) +
  geom_polygon(colour = "black", size = 0.5, aes(group = group)) +
  theme_bw()+scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    #label.position = "left",
    #direction = "vertical",
    barheight = unit(40, units = "mm"),
    barwidth = unit(1, units = "mm"),
    #title.position = "top",
    #title.vjust=0.5,
    #label.vjust = 0.5
  ))+labs(x="",y="",title =" ")

ggsave("stunting1_corrected.pdf", width = 5, height = 4, units = "in")




# Histogram plot

# Labels
Lab = data.frame(state2,cloc2)

Lab$state2 = substr(Lab$state2,1,nchar(Lab$state2)-2)

Lab = Lab[order(Lab$cloc2),]

par(mfrow=c(3,3))

for(s in 1:9){
  hist(filter(data,state==s)$stunting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
  lines(density(filter(data,state==s)$stunting1))
}
s=12
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=29
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data, state==s)$wasting1), col="red", lty=2)

s=11
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=1
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=7
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=20
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=16
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=14
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)

s=18
hist(filter(data,state==s)$wasting1,xlab = Lab$state2[s],main="",probability = T,ylim=c(0,0.4))
lines(density(filter(data,state==s)$wasting1))
abline(v=mean(filter(data,state==s)$wasting1),col="red",lty=2)
