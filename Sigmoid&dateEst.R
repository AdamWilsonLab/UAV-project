#load packages
library(sp)
library(raster)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(minpack.lm)
library(rasterVis)
library(rgdal)
library(gridExtra)

#SOS and EOS for UAV (NDVI)
##spring
spring_grid<-subset(ndvi_grid_d[complete.cases(ndvi_grid_d),],Ndays<200)
plot(spring_grid$Ndays,spring_grid$ndvi)
parms<-par(125,25,0.65,0.3)
spring_models<-spring_grid %>% group_by(grid) %>% do(fit_model=nlsLM(ndvi ~ c/(1+exp(a+b*Ndays))+d,data=.,start = list(a=parms$a,b=parms$b,c=parms$c,d=parms$d),control = nls.lm.control(ftol = 10^-4)))
spring_date_d<-data.frame(grid=character(),sos=integer(),eos=integer(),stringsAsFactors = F)
for (grid in spring_models$grid) {
  model<-spring_models$fit_model[[which(spring_models$grid==grid)]]
  date<- season_date(model,'Ndays',0,200,'spring')
  spring_date_d[nrow(spring_date_d)+1,]<-list(grid,date[1],date[2])
}
write.csv(spring_date_d,file = 'output_data/spring_date_d.csv')
##fall
fall_grid<-subset(ndvi_grid_d[complete.cases(ndvi_grid_d),],Ndays>250)
plot(fall_grid$Ndays,fall_grid$ndvi)
fall_models<-fall_grid%>% group_by(grid) %>% do(fit_model=sigmoid(dataset = ., mid =310 ,dur = -30, c= 0.65, d=0.3,acc = 10^-7))
fall_date_d<- data.frame(grid=character(),sos=integer(),eos=integer(),stringsAsFactors = F)
for (grid in fall_models$grid) { 
  model<- fall_models$fit_model[[which(fall_models$grid==grid)]]
  date<-season_date(model,'Ndays',250,360,'fall')
  fall_date_d[nrow(fall_date_d)+1,]<-list(grid,date[1],date[2])
}
write.csv(fall_date_d,file = 'output_data/fall_date_d.csv')

#SOS and EOS for Sentinel (NDVI)
##spring
spring_grid<-subset(ndvi_grid_s[complete.cases(ndvi_grid_s),],Ndays<200)
plot(spring_grid$Ndays,spring_grid$ndvi)
spring_models<-spring_grid %>% group_by(grid) %>% do(fit_model=sigmoid(dataset = .,mid =125 ,dur = 25, c= 0.65, d=0.3,acc = 10^-7))
spring_date_s<-data.frame(grid=character(),sos=numeric(),eos=numeric(),stringsAsFactors = F)
for (grid in spring_models$grid) {
  model<-spring_models$fit_model[[which(spring_models$grid==grid)]]
  date<- season_date(model,'Ndays',0,200,'spring')
  spring_date_s[nrow(spring_date_s)+1,]<-list(grid,date[1],date[2])
}
write.csv(spring_date_s,file = 'output_data/spring_date_s.csv')
##fall
fall_grid<-subset(ndvi_grid_s[complete.cases(ndvi_grid_s),],Ndays>260)
plot(fall_grid$Ndays,fall_grid$ndvi)
parms<-par(310,-35,0.6,0.25)
fall_models_2<-data.frame(grid=character(),fit_models=character(),stringsAsFactors = F)
for (grid_n in spring_models$grid) {
  sub_data<-subset(fall_grid,grid==grid_n)
  min_val<-min(sub_data$ndvi)
  fit_model<-nlsLM(ndvi ~ (0.85-min_val)/(1+exp(a+b*Ndays))+min_val, data=sub_data ,start = list(a=parms$a,b=parms$b),control = nls.lm.control(ftol  = 0.01))
  fall_models_2[nrow(fall_models_2)+1,]<-list(grid_n,list(fit_model))
}
fall_models<-fall_grid%>% group_by(grid) %>% do(fit_model=nlsLM(ndvi ~ c/(1+exp(a+b*Ndays))+d, data=. ,start = list(a=parms$a,b=parms$b,c=parms$c,d=parms$d),control = nls.lm.control(maxiter = 8)))
fall_date_s<- data.frame(grid=character(),sos=integer(),eos=integer(),stringsAsFactors = F)
for (grid in fall_models$grid) { 
  model<- fall_models$fit_model[[which(fall_models$grid==grid)]]
  date<-season_date(model,'Ndays',260,360,'fall')
  fall_date_s[nrow(fall_date_d)+1,]<-list(grid,date[1],date[2])
}

#NDRE subset for spring and fall
ndre_2018<-subset(ndre_mean_d,date>as.Date('2018-01-01'))
ndre_2018$Ndays<-yday(ndre_2018$date)
plot(ndre_2018$Ndays,ndre_2018$mean)
spring_ndre<-subset(ndre_2018,Ndays<200)
fit_spring_ndre<-sigmoid(spring_ndre,125,35,0.4,0.2,0.0000001)
lines(c(0:200),predict(fit_spring_ndre,data.frame(Ndays=c(0:200))))
fit_spring_ndre_date<-season_date(fit_spring_ndre,'Ndays',0,200,'spring')
abline(v=fit_spring_ndre_date)
abline(v=fit_spring_date,col='red')
fall_ndre<-subset(ndre_2018,Ndays>200)
fit_fall_ndre<-sigmoid(fall_ndre,300,-30,0.4,0.2,0.0000001)
lines(c(200:360),predict(fit_fall_ndre,data.frame(Ndays=c(200:360))))
fit_fall_ndre_date<-season_date(fit_fall_ndre,'Ndays',200,360,'fall')
abline(v=fit_fall_ndre_date)
abline(v=fit_fall_date,col='red')

##SOS and EOS for NDRE
ndre_grid_d$Ndays<-yday(ndre_grid_d$date)
spring_grid<-subset(ndre_grid_d[complete.cases(ndre_grid_d),],Ndays<200)
plot(spring_grid$Ndays,spring_grid$ndre)
spring_models<-spring_grid %>% group_by(grid) %>% do(fit_model=sigmoid(dataset = .,145,35,0.4,0.2,0.0000001))
spring_date_ndre_d<-data.frame(grid=character(),sos=integer(),eos=integer(),stringsAsFactors = F)
for (grid in spring_models$grid) {
  model<-spring_models$fit_model[[which(spring_models$grid==grid)]]
  date<- season_date(model,'Ndays',0,200,'spring')
  spring_date_ndre_d[nrow(spring_date_ndre_d)+1,]<-list(grid,date[1],date[2])
}
write.csv(spring_date_ndre_d,file = 'output_data/spring_date_ndre_d.csv')

##SOF and EOF for NDRE
fall_grid<-subset(ndre_grid_d[complete.cases(ndre_grid_d),],Ndays>250&date>as.Date('2018-01-01'))
plot(fall_grid$Ndays,fall_grid$ndre)
fall_models<-fall_grid%>% group_by(grid) %>% do(fit_model=sigmoid(dataset = .,305,-20,0.3,0.2,0.01))
fall_date_ndre_d<- data.frame(grid=character(),sos=integer(),eos=integer(),stringsAsFactors = F)
for (grid in fall_models$grid) { 
  model<- fall_models$fit_model[[which(fall_models$grid==grid)]]
  date<-season_date(model,'Ndays',250,360,'fall')
  fall_date_ndre_d[nrow(fall_date_ndre_d)+1,]<-list(grid,date[1],date[2])
}
write.csv(fall_date_d,file = 'output_data/fall_date_ndre_d.csv')


#raster
grid_t<-grid1
rownames(spring_date_d)<-substr(spring_date_d$grid,2,4)
grid_data<-SpatialPolygonsDataFrame(grid_t,spring_date_d,match.ID = T)
spring_sos_d<-rasterize(grid_data,rasterSample,'sos')
spring_eos_d<-rasterize(grid_data,rasterSample,'eos')
rownames(fall_date_d)<-substr(fall_date_d$grid,2,4)
grid_data_fall<-SpatialPolygonsDataFrame(grid_t,fall_date_d,match.ID = T)
fall_sos_d<-rasterize(grid_data_fall,rasterSample,'sos')
fall_eos_d<-rasterize(grid_data_fall,rasterSample,'eos')
rownames(spring_date_s)<-substr(spring_date_s$grid,2,4)
grid_data_spring_s<-SpatialPolygonsDataFrame(grid_t,spring_date_s,match.ID = T)
spring_sos_s<-rasterize(grid_data_spring_s,rasterSample,'sos')
spring_eos_s<-rasterize(grid_data_spring_s,rasterSample,'eos')
rownames(spring_date_ndre_d)<-substr(spring_date_ndre_d$grid,2,4)
grid_data<-SpatialPolygonsDataFrame(grid_t,spring_date_ndre_d,match.ID = T)
spring_sos_ndre_d<-rasterize(grid_data,rasterSample,'sos')
spring_eos_ndre_d<-rasterize(grid_data,rasterSample,'eos')
rownames(fall_date_ndre_d)<-substr(fall_date_ndre_d$grid,2,4)
grid_data_fall<-SpatialPolygonsDataFrame(grid_t,fall_date_ndre_d,match.ID = T)
fall_sos_ndre_d<-rasterize(grid_data_fall,rasterSample,'sos')
fall_eos_ndre_d<-rasterize(grid_data_fall,rasterSample,'eos')

spring_date<-merge(spring_date_d,spring_date_s,'grid')
colnames(spring_date)<-c('grid','sos_d','eos_d','sos_s','eos_s')
plot(spring_date$sos_d,spring_date$sos_s)
spring_date_cor<-lm(sos_s~sos_d,data = spring_date)

#raster visualizing of SOS and EOS for NDVI and NDRE
my.col<-colorRampPalette(c('blue', 'green','yellow', 'red'))
r1<-levelplot(spring_sos_d,main='Start of spring season (Number of days)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r2<-levelplot(spring_eos_d,main='End of spring season (Number of days)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r3<-levelplot(fall_sos_d,main='Start of fall season (Number of days)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r4<-levelplot(fall_eos_d,main='End of fall season (Number of days)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
grid.arrange(r1,r2,r3,r4,ncol=2)

r5<-levelplot(spring_sos_ndre_d,main='Start of spring season (NDRE)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r6<-levelplot(spring_eos_ndre_d,main='End of spring season (NDRE)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r7<-levelplot(fall_sos_ndre_d,main='Start of fall season (NDRE)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r8<-levelplot(fall_eos_ndre_d,main='End of fall season (NDRE)',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
grid.arrange(r5,r6,r7,r8,ncol=2)

dif_spring_sos<-overlay(spring_sos_d,spring_sos_ndre_d,fun=function(r1,r2){return(r1-r2)})
cellStats(dif_spring_sos,stat = 'mean')
dif_spring_eos<-overlay(spring_eos_d,spring_eos_ndre_d,fun=function(r1,r2){return(r1-r2)})
cellStats(dif_spring_eos,stat = 'mean')
dif_fall_sos<-overlay(fall_sos_d,fall_sos_ndre_d,fun=function(r1,r2){return(r1-r2)})
cellStats(dif_fall_sos,stat = 'mean')
dif_fall_eos<-overlay(fall_eos_d,fall_eos_ndre_d,fun=function(r1,r2){return(r1-r2)})
cellStats(dif_fall_eos,stat = 'mean')

my.col <- colorRampPalette(c("seagreen", "yellow","firebrick"))
my.colkey<-list(at=seq(-15,15),labels=list())
r9<-levelplot(dif_spring_sos,main='Difference in SOS (Number of days)\nAverage= -6.917',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r10<-levelplot(dif_spring_eos,main='Difference in EOS (Number of days)\nAverage= -15.651',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r11<-levelplot(dif_fall_sos,main='Difference in SOF (Number of days)\nAverage= 12.071',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
r12<-levelplot(dif_fall_eos,main='Difference in EOF (Number of days)\nAverage= 8.261',margin=F,col.regions=my.col,xlab=NULL, ylab=NULL, scales=list(draw=FALSE))
grid.arrange(r9,r10,r11,r12,ncol=2)


#variance
std<-ndvi_grid_d%>%group_by(date)%>%do(std=sd(.$ndvi,na.rm = T))
std$std<-as.numeric(std$std)
std$Ndays<-yday(std$date)
mean(std$std[std$Ndays>250])

#sigmoid function
sigmoid<-function(dataset, mid,dur,c,d,acc){
  parms<-par(mid,dur,c,d)
  fit_model=nlsLM(ndre ~ c/(1+exp(a+b*Ndays))+d,data=dataset,start = list(a=parms$a,b=parms$b,c=parms$c,d=parms$d),control = nls.lm.control(ftol = acc))
  return(fit_model)
}

#parameter for sigmoid
par<-function(mid,dur,c,d){
  a=(10*mid)/dur
  b=-10/dur
  parms<-list(a=a,b=b,c=c,d=d)
  return(parms)
}

#date estimate function
season_date<-function(model,colname,S_day,E_day,season){
  if(season=='spring'){
    newdata<-data.frame(days=c(S_day:E_day))
    colnames(newdata)<-colname
    spline_spring<-smooth.spline(c(S_day:E_day),predict(model,newdata))
    curve_rate<-data.frame(predict(spline_spring,deriv=3))
    min_day<-curve_rate$x[which.min(curve_rate$y)]
    s1<-subset(curve_rate,curve_rate$x<min_day)
    s2<-subset(curve_rate,curve_rate$x>min_day)
    sos<-s1$x[which.max(s1$y)]
    eos<-s2$x[which.max(s2$y)]
    return(c(sos,eos))
  }
  else if(season=='fall'){
    newdata<-data.frame(days=c(S_day:E_day))
    colnames(newdata)<-colname
    spline_fall<-smooth.spline(c(S_day:E_day),predict(model,newdata))
    curve_rate<-data.frame(predict(spline_fall,deriv=3))
    max_day<-curve_rate$x[which.max(curve_rate$y)]
    s1<-subset(curve_rate,curve_rate$x<max_day)
    s2<-subset(curve_rate,curve_rate$x>max_day)
    sos<-s1$x[which.min(s1$y)]
    eos<-s2$x[which.min(s2$y)]
    return(c(sos,eos))
  }
  else {print('Warning: invalid input') }
}

