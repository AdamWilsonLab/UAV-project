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
library(RColorBrewer)

#set working directory
setwd('E:/Final')

#load all flight folders
flight<-list.files('data',full=T)
flight

#create raster stack in GeoTiff format for each flight
for (path in flight){
  file_dir<-list.files(path)
  r1<-raster(paste(path,file_dir[str_detect(file_dir,'blue.tif')],sep = '/'))
  r2<-raster(paste(path,file_dir[str_detect(file_dir,'green.tif')],sep = '/'))
  r3<-raster(paste(path,file_dir[str_detect(file_dir,'red.tif')],sep = '/'))
  r4<-raster(paste(path,file_dir[str_detect(file_dir,'red edge.tif')],sep = '/'))
  r5<-raster(paste(path,file_dir[str_detect(file_dir,'nir.tif')],sep = '/'))
  layer_stack<-stack(r1,r2,r3,r4,r5)
  writeRaster(layer_stack,paste('image',paste(substr(path,6,13),'tif',sep = '.'),sep = '/'),options='INTERLEAVE=BAND', overwrite=TRUE,progress='text')
}

#resample 
raster_master<-brick('image/06282018.tif')[[1]]
for (image in list.files('image', full=T)) {
  image1<-brick(image)
  image1<-resample(image1,raster_master,method="bilinear",progress='text')
  writeRaster(image1,sub('image','resample',image),options='INTERLEAVE=BAND', overwrite=TRUE,progress='text')
  print(paste('completed',image,sep=': '))
}

#create a info table for image path
info<-data.frame(path=list.files('resample',full=T))
info$date=as.Date(substr(info$path,10,17),'%m%d%Y')

#load sub-region
region1<-shapefile('regions/region1.shp')

#crop image, calculate NDVI and save result
for (image_path in info$path){
  raster1<-stack(image_path)
  raster2<-mask(crop(raster1,region1),region1)
  ndvi<-overlay(raster2,fun=function(x){return((x[5]-x[3])/(x[5]+x[3]))})
  writeRaster(ndvi, sub('.tif','_ndvi.tif',sub('resample','ndvi',image_path)), format="GTiff", overwrite=TRUE, progress='text')
}
info$ndvi<-sub('.tif','_ndvi.tif',sub('resample','ndvi', info$path))

#EVI
for (image_path in info$path){
  raster1<-stack(image_path)
  raster2<-mask(crop(raster1,region1),region1)
  evi<-overlay(raster2,fun=function(x){return(2.5*((x[5]-x[3])/(x[5]+6*x[3]-7.5*x[1]+1)))})
  writeRaster(evi, sub('.tif','_evi.tif',sub('resample','evi',image_path)), format="GTiff", progress='text',overwrite=T)
}
info$evi<-sub('.tif','_evi.tif',sub('resample','evi', info$path))

#NDRE
for (image_path in info$path){
  raster1<-stack(image_path)
  raster2<-mask(crop(raster1,region1),region1)
  ndvi<-overlay(raster2,fun=function(x){return((x[5]-x[4])/(x[5]+x[4]))})
  writeRaster(ndvi, sub('.tif','_ndre.tif',sub('resample','ndre',image_path)), format="GTiff", overwrite=TRUE, progress='text')
}
info$ndre<-sub('.tif','_ndre.tif',sub('resample','ndre', info$path))

#shadow removal
sample_raster<-brick('resample/10122018.tif')
sample_raster<-mask(crop(sample_raster,region1),region1)
raster1<-sample_raster
raster1[raster1[[1]]<0.01&raster1[[2]]<0.01&raster1[[3]]<0.01]<-NA
plotRGB(sample_raster,r=3,g=2,b=1,stretch='lin',axes=TRUE)
lines(grid1,col='green')
sample2<-raster('ndvi_shadowRM/10302018_ndvi.tif')
my.col<-colorRampPalette(c('white','orange','yellow','dark green'))
plot(levelplot(sample2,col.regions=my.col,margin=F))
plot(sample2,col=rev(terrain.colors(7)))
lines(grid1)

#shadow removal for EVI
for (path in info$path) {
  raster1<-brick(path)
  raster1<-mask(crop(raster1,region1),region1)
  raster1[raster1[[1]]<0.01&raster1[[2]]<0.01&raster1[[3]]<0.01]<-NA
  evi<-raster(info$evi[info$path==path])
  evi<-mask(evi,raster1[[1]])
  writeRaster(evi,sub('evi','evi_shadowRM',info$evi[info$path==path]),format="GTiff",overwrite=T,progress='text')
}

#shadow removal for NDVI
for (path in info$path) {
  raster1<-brick(path)
  raster1<-mask(crop(raster1,region1),region1)
  raster1[raster1[[1]]<0.01&raster1[[2]]<0.01&raster1[[3]]<0.01]<-NA
  ndvi<-raster(info$ndvi[info$path==path])
  ndvi<-mask(ndvi,raster1[[1]])
  writeRaster(ndvi,sub('ndvi','ndvi_shadowRM',info$ndvi[info$path==path]),format="GTiff",overwrite=T,progress='text')
}

#shadow removal for NDRE
for (path in info$path) {
  raster1<-brick(path)
  raster1<-mask(crop(raster1,region1),region1)
  raster1[raster1[[1]]<0.01&raster1[[2]]<0.01&raster1[[3]]<0.01]<-NA
  ndre<-raster(info$ndre[info$path==path])
  ndre<-mask(ndre,raster1[[1]])
  writeRaster(ndre,sub('ndre','ndre_shadowRM',info$ndre[info$path==path]),format="GTiff",overwrite=T,progress='text')
}

#raster statistics
#ndvi
stack1<-stack(list.files('ndvi_shadowRM',full=T))
ndvi_stats<-data.frame(date=names(stack1),mean=cellStats(stack1,mean))
ndvi_stats$date<-substr(ndvi_stats$date,2,9)
ndvi_stats$date<-as.Date(ndvi_stats$date,'%m%d%Y')
write.csv(ndvi_stats,file = 'output_data/ndvi_mean_d.csv')
#ndre
stack1<-stack(list.files('ndre_shadowRM',full=T))
ndre_stats<-data.frame(date=names(stack1),mean=cellStats(stack1,mean))
ndre_stats$date<-substr(ndre_stats$date,2,9)
ndre_stats$date<-as.Date(ndre_stats$date,'%m%d%Y')
write.csv(ndre_stats,file = 'output_data/ndre_mean_d.csv')
#EVI
stack1<-stack(list.files('evi_shadowRM',full=T))
evi_stats<-data.frame(date=names(stack1),mean=cellStats(stack1,mean))
evi_stats$date<-substr(evi_stats$date,2,9)
evi_stats$date<-as.Date(evi_stats$date,'%m%d%Y')
write.csv(evi_stats,file = 'output_data/evi_mean_d.csv')
#Aggregate ndvi from drone data to the grid
ex_ndvi_d<-raster::extract(stack1,grid1,fun=mean,na.rm=T)
ndvi_d<-t(ex_ndvi_d)
colnames(ndvi_d)<-c(1:length(grid1))
rownames(ndvi_d)<-substr(rownames(ndvi_d),2,9)
ndvi_d<-data.frame(ndvi_d,stringsAsFactors = F)
ndvi_d$date<-as.Date(rownames(ndvi_d),'%m%d%Y')
ndvi_d<-gather(ndvi_d,grid,ndvi,X1:X338)
write.csv(ndvi_d,file='output_data/ndvi_grid_d.csv')
#Aggregate UAV data under one sentinel pixel
stack_150<-mask(crop(stack1,grid1[150]),grid1[150])
grid2<-rasterToPolygons(stack_150[[1]])
ex_ndvi_d_2<-raster::extract(stack_150,grid2,weights=F,na.rm=T,df=T)
ndvi_d_2<-gather(ex_ndvi_d_2,date,ndvi,2:30)
ndvi_d_2$date<-as.Date(substr(ndvi_d_2$date,2,9),'%m%d%Y')
ndvi_d_2$ID<-as.character(ndvi_d_2$ID)
write.csv(ndvi_d_2,file = 'output_data/ndvi_pixel_150_d.csv')
#Aggregate ndre from drone data to the grid
stack_ndre<<-stack(list.files('ndre_shadowRM',full = T))
ex_ndre<-raster::extract(stack_ndre,grid1,fun=mean,na.rm=T)
ndre_d<-t(ex_ndre)
colnames(ndre_d)<-c(1:length(grid1))
rownames(ndre_d)<-substr(rownames(ndre_d),2,9)
ndre_d<-data.frame(ndre_d,stringsAsFactors = F)
ndre_d$date<-as.Date(rownames(ndre_d),'%m%d%Y')
ndre_d<-gather(ndre_d,grid,ndre,X1:X338)
write.csv(ndre_d,file='output_data/ndre_grid_d.csv')

#Sentinel data process
#create stack for all ndvi maps and resize to study region
sentinel<-raster('Sentinel_NDVI/S2A_MSIL2A_20180325T155911_N0206_R097_T17TPH_20180325T211852.tif')
sentinel_list<-list.files('Sentinel_NDVI',full=T)
stack_sentinel<-stack(sentinel_list)
stack_sentinel<-mask(crop(stack_sentinel,region1),region1)
#cerate outlines of Sentinel pixels in the region
rasterSample<-raster(sentinel_list[1])
rasterSample<-mask(crop(rasterSample,region1),region1)
plot(rasterSample)
grid1<-polygons(rasterToPolygons(rasterSample))
lines(grid1)
#raster statistics
stats_sentinel<-data.frame(date=names(stack_sentinel),mean=cellStats(stack_sentinel,mean),stringsAsFactors = F)
stats_sentinel$date<-as.Date(substr(stats_sentinel$date,12,19),'%Y%m%d')
write.csv(stats_sentinel,file = 'output_data/ndvi_mean_s.csv')
ggplot(data = stats_sentinel,aes(x=date,y=mean))+geom_line()
#extract ndvi values and data wrangling
ex_S<-raster::extract(stack_sentinel,grid1,fun=mean,na.rm=T)
ndvi_S<-t(ex_S)
colnames(ndvi_S)<-c(1:length(grid1))
rownames(ndvi_S)<-substr(rownames(ndvi_S),12,19)
ndvi_S<-data.frame(ndvi_S,stringsAsFactors = F)
ndvi_S$date<-as.Date(rownames(ndvi_S),'%Y%m%d')
ndvi_S<-gather(ndvi_S,grid,ndvi,X1:X338)
write.csv(ndvi_S,file = 'output_data/ndvi_grid_s.csv')
ggplot(data = ndvi_S,aes(x=date,y=ndvi))+geom_line(aes(fill=grid),alpha=0.05,size=0.1,show.legend = F)+
  geom_line(data = stats_sentinel,aes(x=date,y=mean),color='red',size=1)+
  theme_bw()

#plots
ggplot(data = ndvi_d,aes(x=date,y=ndvi))+
  geom_line(aes(fill=grid),alpha=0.05,size=0.1,show.legend = F)+
  geom_line(data = ndvi_stats,aes(x=date,y=mean),color='red',size=1)+
  geom_line(data = stats_sentinel,aes(x=date,y=mean),color='blue',size=1)+
  scale_x_date(limits = as.Date(c('2018-03-01','2018-12-01')))+
  theme_bw()

ggplot()+
  geom_line(data=ndvi_stats,aes(x=date,y=mean),color='red',size=1)+
  geom_line(data=ndre_stats,aes(x=date,y=mean),color='blue',size=1)+
  geom_line(data = evi_stats,aes(x=date,y=mean),color='green', size=1)+
  scale_x_date(date_labels = '%d-%b',date_breaks = "2 weeks",limits = as.Date(c('2018-04-01','2018-12-01')))+
  theme_bw()

ggplot()+
  geom_boxplot(data=ndvi_d,aes(x=date,y=ndvi,group=date),outlier.alpha = 0.02)+
  geom_line(data=ndvi_stats,aes(x=date,y=mean),color='red',alpha=0.5)+
  scale_x_date(date_labels = '%b',date_breaks = "1 month",limits = as.Date(c('2018-04-01','2018-12-01')))+
  theme_bw()

ggplot(data = ndvi_d_2,aes(x=date,y=ndvi))+
  geom_line(aes(fill=ID),alpha=0.01,size=0.01,show.legend = F)+
  geom_line(data = subset(ndvi_S,grid=='X150'),aes(x=date,y=ndvi),color='red',size=1)+
  scale_x_date(date_labels = '%b',date_breaks = "1 month",limits = as.Date(c('2018-03-01','2018-12-01')))+
  theme_bw()

#sigmoid regression for mean ndvi of spring summer and fall
ndvi_d$Ndays<-yday(ndvi_d$date)
ndvi_stats$Ndays<-yday(ndvi_stats$date)
spring_ndvi_stats<-subset(ndvi_stats,Ndays<150)
fall_ndvi_stats<-subset(ndvi_stats,Ndays>250)
summer_ndvi_stats<-subset(ndvi_stats,Ndays>150&Ndays<280)

ggplot()+
  geom_line(aes(x=c(50:150),y=predict(fit_spring_2,newdata=data.frame(Ndays=c(50:150)))))+
  geom_line(aes(x=c(250:360),y=predict(fit_fall_2,newdata=data.frame(Ndays=c(250:360)))))+
  geom_line(aes(x=c(150:250),y=predict(fit_summer,newdata=data.frame(Ndays=c(150:250)))))+
  geom_point(data = ndvi_stats,aes(x=Ndays,y=mean))+
  labs(x='Number of days',y='NDVI')+
  scale_x_continuous(breaks = seq(50,350,50))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme_bw()

#fit_spring_2<-nlsLM(mean ~ Asym/(1+exp((xmid-Ndays)/scal))+0.3,data=spring_ndvi_stats,start = list(Asym=0.60,xmid=110,scal=4))
#fit_fall_2<-nlsLM(mean ~ Asym/(1+exp((xmid-Ndays)/scal))+0.3,data=fall_ndvi_stats,start = list(Asym=0.6,xmid=300,scal=-15))
fit_summer<-lm(mean ~ Ndays,data =summer_ndvi_stats)



#correlation between UAV and Sentinel data
fall_drone<-subset(ndvi_grid_d, date==as.Date('2018-10-09'))
fall_sentinel<-subset(ndvi_grid_s,date==as.Date('2018-10-09'))
fall<-left_join(fall_drone,fall_sentinel,by=c('grid'='grid'))
cor(fall$ndvi.x,fall$ndvi.y)
fall_cor<-lm(ndvi.x~ndvi.y , data = fall)
fall$diff<-fall$ndvi.x-fall$ndvi.y
mean(abs(fall$diff))
p1<-ggplot(data = fall, aes(x=ndvi.y,y=ndvi.x))+
  geom_point()+
  stat_smooth(method = "lm", col = "red", se=F)+
  geom_abline(intercept = 0,slope = 1)+
  labs(x='NDVI from Sentinel',y='NDVI from UAV',title = 'Fall')

summer_drone<-subset(ndvi_grid_d, date==as.Date('2018-07-09'))
summer_sentinel<-subset(ndvi_grid_s,date==as.Date('2018-07-08'))
summer<-left_join(summer_drone,summer_sentinel,by=c('grid'='grid'))
cor(summer$ndvi.x,summer$ndvi.y)
summer_cor<-lm(ndvi.x~ndvi.y,data = summer)
summer$diff<-summer$ndvi.x-summer$ndvi.y
mean(abs(summer$diff))
p2<-ggplot(data = summer, aes(x=ndvi.y,y=ndvi.x))+
  geom_point()+
  stat_smooth(method = "lm", col = "red", se=F)+
  geom_abline(intercept = 0,slope = 1)+
  labs(x='NDVI from Sentinel',y='NDVI from UAV',title = 'Summer')

spring_drone<-subset(ndvi_grid_d, date==as.Date('2018-05-16'))
spring_sentinel<-subset(ndvi_grid_s,date==as.Date('2018-05-17'))
spring<-left_join(spring_drone,spring_sentinel,by=c('grid'='grid'))
spring<-subset(spring, ndvi.x&ndvi.y)
cor(spring$ndvi.x,spring$ndvi.y)
spring_cor<-lm(ndvi.x~ndvi.y,data = spring)
spring$diff<-spring$ndvi.x-spring$ndvi.y
mean(abs(spring$diff))
p3<-ggplot(data = spring, aes(x=ndvi.y,y=ndvi.x))+
  geom_point()+
  stat_smooth(method = "lm", col = "red", se=F)+
  geom_abline(intercept = 0,slope = 1)+
  labs(x='NDVI from Sentinel',y='NDVI from UAV',title = 'Spring')

grid.arrange(p3,p2,p1,ncol=2)

