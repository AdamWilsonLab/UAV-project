#load packages
library(sp)
library(raster)
library(rgdal)
library(gdalUtils)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(minpack.lm)

#set working directory
setwd('E:/Final')

#load dataset
ndvi_grid_d<-read.csv('output_data/ndvi_grid_d.csv',stringsAsFactors = F)
ndvi_grid_d$date<-as.Date(ndvi_grid_d$date)
ndvi_mean_d<-read.csv('output_data/ndvi_mean_d.csv', stringsAsFactors = F)
ndvi_mean_d$date<-as.Date(ndvi_mean_d$date)
ndvi_grid_s<-read.csv('output_data/ndvi_grid_s.csv', stringsAsFactors = F)
ndvi_grid_s$date<-as.Date(ndvi_grid_s$date)
ndvi_mean_s<-read.csv('output_data/ndvi_mean_s.csv', stringsAsFactors = F)
ndvi_mean_s$date<-as.Date(ndvi_mean_s$date)
ndre_mean_d<-read.csv('output_data/ndre_mean_d.csv',stringsAsFactors = F)
ndre_mean_d$date<-as.Date(ndre_mean_d$date)
ndre_mean_d$Ndays<-yday(ndre_mean_d$date)
ndre_grid_d<-read.csv('output_data/ndre_grid_d.csv',stringsAsFactors = F)
ndre_grid_d$date<-as.Date(ndre_grid_d$date)
evi_mean_d<-read.csv('output_data/evi_mean_d.csv',stringsAsFactors = F)
evi_mean_d$date<-as.Date(evi_mean_d$date)
ndvi_grid_150_d<-read.csv('output_data/ndvi_pixel_150_d.csv',stringsAsFactors = F)
ndvi_grid_150_d$date<-as.Date(ndvi_grid_150_d$date)
ndvi_grid_150_d$ID<-as.character(ndvi_grid_150_d$ID)

#plots 
ggplot(data = ndvi_grid_d,aes(x=date,y=ndvi))+
  geom_boxplot(data=ndvi_grid_d,aes(x=date,y=ndvi,group=date),outlier.alpha = 0.05)+
  geom_line(data = ndvi_mean_d,aes(x=date,y=mean),color='red',size=1,alpha=0.7)+
  geom_line(data = ndvi_mean_s,aes(x=date,y=mean),color='blue',size=1,alpha=0.7)+
  scale_x_date(date_labels = '%b',date_breaks = "1 month",limits = as.Date(c('2018-03-15','2018-12-01')))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  labs(x=NULL,y='NDVI')+
  theme_bw()

ggplot()+
  geom_point(data=ndvi_mean_d,aes(x=date,y=mean),color='red',shape=15,size=3)+
  geom_point(data=ndre_mean_d,aes(x=date,y=mean),color='blue',shape=16,size=3)+
  geom_point(data=evi_mean_d,aes(x=date,y=mean),color='dark green',shape=18,size=3)+
  geom_smooth(data=ndvi_mean_d,aes(x=date,y=mean),color='red',size=1,se=F,span=0.15)+
  geom_smooth(data=ndre_mean_d,aes(x=date,y=mean),color='blue',size=1,se=F,span=0.15)+
  geom_smooth(data = evi_mean_d,aes(x=date,y=mean),color='green', size=1,se=F,span=0.12)+
  scale_x_date(date_labels = '%d.%b',date_breaks = "2 weeks",limits = as.Date(c('2018-04-01','2018-12-01')))+
  labs(x=NULL,y='NDVI / NDRE / EVI')+
  theme_bw()

ggplot()+
  geom_boxplot(data=ndvi_grid_d,aes(x=date,y=ndvi,group=date),outlier.alpha = 0.02)+
  geom_line(data=ndvi_mean_d,aes(x=date,y=mean),color='red',alpha=0.5,size=1)+
  scale_x_date(date_labels = '%b',date_breaks = "1 month",limits = as.Date(c('2018-04-01','2018-12-01')))+
  labs(x=NULL,y='NDVI')+
  theme_bw()

ggplot(data = ndvi_grid_150_d,aes(x=date,y=ndvi))+
  geom_boxplot(aes(group=date),outlier.alpha = 0.05)+
  geom_line(data = subset(ndvi_grid_s,grid=='X150'),aes(x=date,y=ndvi),color='red',size=1,alpha=0.65)+
  scale_x_date(date_labels = '%b',date_breaks = "1 month",limits = as.Date(c('2018-03-15','2018-12-01')))+
  labs(x=NULL,y='NDVI')+
  theme_bw()

#sigmoid model for spring summer and fall of NDVI and NDRE
ndvi_grid_d$Ndays<-yday(ndvi_grid_d$date)
ndvi_mean_d$Ndays<-yday(ndvi_mean_d$date)
evi_mean_d$Ndays<-yday(evi_mean_d$date)
spring_ndvi_stats<-subset(ndvi_mean_d,Ndays<150)
fall_ndvi_stats<-subset(ndvi_mean_d,Ndays>250)
summer_ndvi_stats<-subset(ndvi_mean_d,Ndays>150&Ndays<280)
summer_ndre_stats<-subset(ndre_mean_d,Ndays>160&Ndays<250)
summer_evi_stats<-subset(evi_mean_d,Ndays>150&Ndays<250)

fit_spring<-nlsLM(mean ~ Asym/(1+exp((xmid-Ndays)/scal))+0.3,data=spring_ndvi_stats,start = list(Asym=0.60,xmid=110,scal=4))
fit_fall<-nlsLM(mean ~ Asym/(1+exp((xmid-Ndays)/scal))+0.3,data=fall_ndvi_stats,start = list(Asym=0.6,xmid=300,scal=-15))
fit_summer<-lm(mean ~ Ndays,data =summer_ndvi_stats)
fit_summer_ndre<-lm(mean~Ndays,data = summer_ndre_stats)
fit_spring_date<-season_date(fit_spring,'Ndays',50,200,'spring')
fit_spring_ndre_date
fit_fall_date<-season_date(fit_fall,'Ndays',250,360,'fall')
fit_fall_ndre_date

newdata<-data.frame(Ndays=c(250:360))
spline_fall<-smooth.spline(c(250:360),predict(fit_fall,newdata))
curve_rate<-data.frame(predict(spline_fall,deriv=3))
plot(curve_rate$x,curve_rate$y)
curve_rate$x[which.max(curve_rate$y)]
curve_rate$x[which.min(curve_rate$y)]
abline(v=fit_fall_date)

#plots for sigmoid models
ggplot()+
  geom_line(aes(x=c(50:150),y=predict(fit_spring,newdata=data.frame(Ndays=c(50:150)))),size=1,col='dark green')+
  geom_line(aes(x=c(50:160),y=predict(fit_spring_ndre,newdata=data.frame(Ndays=c(50:160)))),size=1,col='blue')+
  geom_line(aes(x=c(250:360),y=predict(fit_fall,newdata=data.frame(Ndays=c(250:360)))),size=1,col='dark green')+
  geom_line(aes(x=c(230:360),y=predict(fit_fall_ndre,newdata=data.frame(Ndays=c(230:360)))),size=1,col='blue')+
  geom_line(aes(x=c(150:250),y=predict(fit_summer,newdata=data.frame(Ndays=c(150:250)))),size=1,col='dark green')+
  geom_line(aes(x=c(160:230),y=predict(fit_summer_ndre,newdata=data.frame(Ndays=c(160:230)))),size=1,col='blue')+
  #geom_point(data = ndvi_mean_d,aes(x=Ndays,y=mean),color='red',size=2,shape=4)+
  geom_vline(xintercept = c(fit_spring_date,fit_fall_date),linetype = "dashed",col='green',size=1)+
  geom_vline(xintercept = c(fit_spring_ndre_date,fit_fall_ndre_date),linetype = "dashed",col='blue',size=1)+
  geom_label(aes(x=c(fit_spring_date,fit_fall_date),y=0.89,label=c('SOS\n123','EOS\n138','SOF\n290','EOF\n331')),size=3, angle=0, vjust=-0.50, hjust=0.5,label.size = 0)+
  geom_label(aes(x=c(fit_spring_ndre_date,fit_fall_ndre_date),y=0.27,label=c('SOS\n131','EOS\n153','SOF\n268','EOF\n313')),size=3, angle=0, vjust=-0.50, hjust=0.5,label.size = 0)+
  labs(x='Number of days',y='NDVI')+
  scale_x_continuous(breaks = seq(50,350,50))+
  scale_y_continuous(breaks = seq(0,1,0.05))+
  theme_bw()


#self-start function
##fit_spring <- nls(mean ~ SSlogis(Ndays, Asym, xmid, scal), data =spring_ndvi_stats)
##fit_fall <- nls(mean ~ SSlogis(Ndays, Asym, xmid, scal), data =fall_ndvi_stats)

