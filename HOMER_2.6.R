############################
#######   HOMER 2.6  #######
############################
#### Version 18/03/2013 ####
############################
#### * Estimates urban trends
#### * Correction for rainfall - and cumulative parameters enabled
#### * Graphic modifications in comparison series
####    - all comparisons on a single sheet
####    - standard deviation of noise plotted
#### * Improved visualization of series (coloured polygons option)
#### * Allows moving geographic and correlation neighbourhoods
#### * Neighbourhood options may be changed in the menu
#### 1.06 NEWS
#### * Correction factors put in file in "meta" directory (Brigitte Dubuisson)
#### * Synthetic breaks points figures (in fig, PNG file) (Brigitte Dubuisson)
#### * Fast Quality Control implemented
#### * Automatic removal of validated outliers
#### 1.07 NEWS
#### * Computation of Hannart&Naveau posterior distribution (not used for the moment)
#### * Seasonal pairwise comparisons enabled
#### * New graphics for pairwise (synth_a_...
####
#### 1.08 NEWS
#### * Correction of minor bugs in FQC (neighborhood) and main menu (when quitting 
####   directly
####
#### 1.09 NEWS
#### * Modification of break points menu : option "b" now 
####   allows modify existing xxxxdetected.txt file
####
#### 1.10 NEWS
#### * Already validated breaks in the detected file are 
####   shown as ligth grey lines in the comparison series
####
#### 1.20 NEWS
#### * Minor modification : seasons are DJF MAM JJA SON
#### * Seasonal averages or sums can now be visualized in "v" option
####
#### 1.30 NEWS
#### * Minor modification : change-point dates are automatically sorted
####   prior to correction - since users could make mistakes entering
####   the dates manually
####
#### 1.40 NEWS
#### * Bug correction in gaps reconstitution: cumulative parameter (ratio or log-ratio)
####
#### 1.50 NEWS
#### * Detection on monthly series enabled
#### * Correction of cumulative parameters maybe estimated with varying monthly coeff 
#### * Checks if data are duplicated in different files
####
#### 1.51 NEWS
#### * Correction of a minor bug in typing outliers function
#### * All graphic outputs in ".eps", standardized names
#### 
#### 1.60 NEWS
#### * Integration of multivariate detection procedure
####
#### 1.61 NEWS
#### * Menus modifications
#### * Additional graphic outputs (pdf,png...)
####
#### 1.7 NEWS
#### * Urban trends estimation enabled
####
#### 1.8 NEWS
#### * Multivariate detection enabled
#### * CLIMATOL checks enabled
####
#### 1.9 NEWS
#### * ACMANT detection enabled
####
#### 2.0 NEWS
#### * Neighborhoods enabled for multivariate detection
#### * Monthly coefficients now smoothed
####
#### 2.1 - 2.4 various bug fixes; nicer visualizations ( well, we hope so :-P )
#### 
#### 2.5 NEWS
#### * Station file check (number of columns) : original file save in xxxxxxstations.txt.bak
#### * Checks presence of data files and length of series (removed from the list is <15)
#### * Detects the "missing period between two breaks case" (program no longer stops when correcting)
#### * Detect the "missing data for all series" case (joint detection no longer stops)
#### * ACMANT and Month of change assessment now forced to run on raw data
#### * Successive versions of "detected" file are now stored in "tmp" directory
####
#### 2.6 NEWS
####
#### * Various bugs corrected in "correction" procedure for cumulative parameters and zero series
#### * No longer requires "segclust" - uniseg (cghseg package) used instead of segmeans
#### 
#### 2.6.1
#### * Corrected bug: joint detection could not appear, since range of plots took only into
####   account pairwise detection
####
#### 2.6.2
#### * Corrected small bug: distances computation in CLIMATOL correlogram
#### 
#### 2.6.3
#### * Corrected small (but tedious) bug, when zero series on the whole period
####
#### Known bug : when dealing with homogenized data,cghseg may produce an error and stop the program
####
####
#### Requires : GSL installed on computer (GNU math library)
#### Requires : R library (cghseg); note that cghseg fails if GSL not installed
#### Requires : R library (mapproj,maps)
#### File format : same as benchmark dataset
####
#### Launch : under R session type:  source("home.2.6.R")

###########################################
#### Report bugs to olivier.mestre@meteo.fr
###########################################

#### Copyright (C) 2012 Olivier Mestre, Peter Domonkos, Jose Guijarro, Enric Aguilar
####
#### This program is free software: you can redistribute it and/or modify
#### it under the terms of the GNU General Public License as published by
#### the Free Software Foundation, either version 3 of the License, or
#### (at your option) any later version.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.

#### You should have received a copy of the GNU General Public License
#### along with this program.  If not, see <http://www.gnu.org/licenses/>
####
#########################################################################
#########################################################################

rm(list=ls())
options(show.error.messages=T)


#pepe=require(cghseg)
# if (pepe == F){install.packages("cghseg",repos="http://cran.cict.fr/");library(cghseg)}
# Cannot use the automatic install of cghseg, as the latest version does not work properly. 
#pepe=require(mapproj)
#if (pepe == F){install.packages("mapproj",repos="http://cran.cict.fr/");library(mapproj)}
#pepe=require(maps)
#if (pepe == F){install.packages("maps",repos="http://cran.cict.fr/");library(maps)}
library(cghseg)
require(mapproj)

############################### Dir Create ####################################
if(file.exists('qc')==FALSE){dir.create('qc')}
if(file.exists('fig')==FALSE){dir.create('fig')}
if(file.exists('meta')==FALSE){dir.create('meta')}
if(file.exists('tmp')==FALSE){dir.create('tmp')}
if(file.exists('raw')==FALSE){dir.create('ra')}
if(file.exists('ho')==FALSE){dir.create('ho')}



###################################################################################
read.serie=function(file,month.year.flag,j.flag=F){
  ###################################################################################
  
  ## Reading files
  ## month.year.flag = 1,2,...12 (month), 13 for year,14 for NDJFMA,15 for MJJASO
  ## 16 for DJF, 17 for MAM, 18 for JJA, 19 for SON
  ## Annual AVERAGES are used for additive models
  ## Annual SUMS otherwise
  ######################################################
  
  ## To be improved, loading the full frame once and for all
  ## instead of reading once at a time
  ##########################################################
  
  ## Seasons
  if (month.year.flag ==14){season = c(11,12,1:4)+1}
  if (month.year.flag ==15){season = (5:10)+1}
  
  if (month.year.flag ==16){season = c(12,1:2)+1}
  if (month.year.flag ==17){season = (3:5)+1}
  if (month.year.flag ==18){season = c(6:8)+1}
  if (month.year.flag ==19){season = (9:11)+1}
  if (month.year.flag ==20){season = (c(5,6,7,8,11,12,1,2))+1;
  w.acma =  c(rep(1,3),.5,rep(-1,3),.5)}
  
  data      = read.table(file)
  n         = nrow(data)
  tt        = as.vector(data[,1])
  mean.data = 0
  cpt.data  = 0
  if (month.year.flag == 13){
    y    = rep(miss.flag,n)
    for (i in 1:n) {
      if (sum(data[i,2:13]==miss.flag)==0){
        if (comp.option=="a"){
          y[i]=mean(t(data[i,2:13]))
        } else {
          y[i]      = sum(data[i,2:13])
          mean.data = mean.data+y[i]
          cpt.data  = cpt.data+1}
      }
    }
    if (comp.option!="a" & j.flag==T & cpt.data>0){
      mean.data=mean.data/cpt.data
      if (mean.data!=0) {
        y[y != miss.flag]=y[y != miss.flag]/mean.data
      } else {y[y != miss.flag]=0}
    }
  }
  
  if (month.year.flag > 13 & month.year.flag !=16){
    y    = rep(miss.flag,n)
    for (i in 1:n) {
      if (sum(data[i,season]==miss.flag)==0){
        if (comp.option=="a"){
          y[i]=mean(t(data[i,season]))
        } else {y[i]=sum(data[i,season])}}}
  }
  
  
  if (month.year.flag ==16){
    y    = rep(miss.flag,n)
    for (i in 2:n) {
      if (data[i-1,13]!=miss.flag &
          sum(data[i,2:3]==miss.flag)==0) {
        y[i]=data[i-1,13]+data[i,2]+data[i,3] # beware 1st column == year
        if (comp.option=="a"){y[i]=y[i]/3}}}
  }
  
  if (month.year.flag ==20){
    y    = rep(miss.flag,n)
    for (i in 1:n) {
      if (sum(data[i,]==miss.flag)==0){ ### ensures period compatibility with annual mean
        if (comp.option=="a"){
          y[i]=sum(t(data[i,season])*w.acma)/3.5
        } else {y[i]=NA;options(show.error.messages=F)
        write("ATTEMPT USING ACMANT WITH MULTIPLICATIVE PARAMETER",file="")
        stop()}}}
  }
  
  if (month.year.flag < 13){y=data[,month.year.flag+1]}
  return(cbind(tt,y))}
####################

########################################################
fig.file=function(file,he=6,wi=8,unit="in"){
  ########################################################
  
  ## Opens specified devices for figures
  ######################################
  ######################################
  
  if (dev.str == "ps")  { file=paste(file,".eps",sep="") 
  postscript(file,horizontal=F,paper="special",
             height=he,width=wi)}
  
  if (dev.str == "svg")  { file=paste(file,".svg",sep="") 
  svg(file,paper="special",height=he,width=wi)}
  
  if (dev.str == "pdf") { file=paste(file,".pdf",sep="")
  pdf(file,paper="special",height=he,width=wi)}
  
  if (dev.str == "png") { file=paste(file,".png",sep="")
  png(file,units=unit,height=he,width=wi)}
  return()}

##############################################################################################################
read.zone=function(list.of.files,month.year.flag,list.file){
  ##############################################################################################################
  
  ## Reading files, producing a single frame as output
  ## format following seglm requirements
  ## month.year.flag = 1,2,...12 (month), 13 for year,14 for NDJFMA,15 for MJJASO
  ## 16 for DJF, 17 for MAM, 18 for JJA, 19 for SON
  ## Annual AVERAGES are used for additive models
  ## Annual SUMS otherwise
  ######################################################
  
  ## Seasons
  if (month.year.flag ==14){season = c(11,12,1:4)+1}
  if (month.year.flag ==15){season = (5:10)+1}
  
  if (month.year.flag ==16){season = c(12,1:2)+1}
  if (month.year.flag ==17){season = (3:5)+1}
  if (month.year.flag ==18){season = c(6:8)+1}
  if (month.year.flag ==19){season = (9:11)+1}
  
  an.deb =  1000000
  an.fin = -1000000
  i.ref  = numeric()  ### moving neighbourhood option
  
  for (files in list.of.files){
    i.ref = c(i.ref,match(files,list.file))
    data  = read.table(files)
    n     = nrow(data)
    tt    = as.vector(data[,1])
    test  = diff(tt)
    if (sum(test !=1)>0){options(show.error.messages=F)
      write(paste("Non consecutive years in file",files),file="")
      stop()}
    if (tt[1]<an.deb) {an.deb=tt[1]}
    if (tt[n]>an.fin) {an.fin=tt[n]}
  }
  
  y.frame  = an.deb:an.fin
  n.frame  = length(y.frame)
  p.frame  = length(list.of.files)
  
  logratios = numeric()
  year      = numeric()
  position  = numeric()
  patient   = numeric()
  station   = character()
  index     = character()
  
  group   = rep(1,n.frame*p.frame)
  for (i in 1:p.frame) {
    patient   = c(patient,rep(i,n.frame))
    station   = c(station,rep(station.name[i.ref[i]],n.frame))
    index     = c(index,rep(station.index[i.ref[i]],n.frame))
    position  = c(position,1:n.frame)
    year      = c(year,y.frame)
  }
  
  for (files in list.of.files){
    y    = rep(miss.flag,n.frame)
    data = read.table(files)
    n    = nrow(data)
    tt   = as.vector(data[,1])
    i.sh = tt[1]-an.deb      # shift index for storage in d.frame
    
    if (month.year.flag == 13){
      for (i in 1:n) {
        if (sum(data[i,2:13]==miss.flag)==0){
          if (comp.option=="a"){
            y[i+i.sh]=mean(t(data[i,2:13]))
          } else {y[i+i.sh]=sum(data[i,2:13])}}}}
    
    if (month.year.flag > 13 & month.year.flag !=16){
      for (i in 1:n) {
        if (sum(data[i,season]==miss.flag)==0){
          if (comp.option=="a"){
            y[i+i.sh]=mean(t(data[i,season]))
          } else {y[i+i.sh]=sum(data[i,season])}}}
    }
    
    if (month.year.flag ==16){
      for (i in 2:n) {
        if (data[i-1,13]!=miss.flag &
            sum(data[i,2:3]==miss.flag)==0) {
          y[i+i.sh]=data[i-1,13]+data[i,2]+data[i,3] # beware 1st column == year
          if (comp.option=="a"){y[i+i.sh]=y[i+i.sh]/3}}}
    }
    if (month.year.flag < 13){
      y[(i.sh+1):(i.sh+1+nrows(data))]=data[,month.year.flag+1]}
    logratios=c(logratios,y)
  }
  i.miss            = which(logratios==miss.flag)
  logratios[i.miss] = NA
  d.frame   = data.frame(group,patient,position,logratios,station,index,year)
  return(d.frame)}
################


##############################################################################################################
read.cghseg=function(list.of.files,month.year.flag,list.file){
  ##############################################################################################################
  
  ## Reads data into a data frame meeting CGHSEG requirements
  
  an.deb      =  1000000
  an.fin      = -1000000
  name.cghseg = character()
  
  for (files in list.of.files){
    i.ref       = match(files,list.file)
    name.cghseg = c(name.cghseg,station.index[i.ref])
    data        = read.table(files)
    n           = nrow(data)
    tt          = as.vector(data[,1])
    test        = diff(tt)
    if (sum(test !=1)>0){options(show.error.messages=F)
      write(paste("Non consecutive years in file",files),file="")
      stop()}
    if (tt[1]<an.deb) {an.deb=tt[1]}
    if (tt[n]>an.fin) {an.fin=tt[n]}
  }
  
  n           = an.fin-an.deb+1
  p           = length(list.of.files)
  
  cghseg.frame = matrix(NA,nr=n,nc=p)
  
  k=0
  for (files in list.of.files){
    k                                     = k+1
    data                                  = read.serie(files,month.year.flag,T)
    tt                                    = as.vector(data[,1])
    n.data                                = length(tt)
    i.sh                                  = tt[1]-an.deb+1  # shift index for storage in check.frame
    data                                  = data[,-1]
    data[data==miss.flag]                 = NA
    data                                  = as.vector(t(as.matrix(data)))
    cghseg.frame[i.sh:(i.sh+n.data-1),k]  = data
    cghseg.frame                          = as.data.frame(cghseg.frame)
  }
  ### Checking whether some years are missing for every station in the neighborhood
  ### Since cghseg produces an error in that case
  
  cghseg.frame        = cbind(an.deb:an.fin,cghseg.frame)
  names(cghseg.frame) = c("year",name.cghseg)
  
  i.remove=numeric()
  for (i in 1:nrow(cghseg.frame)){
    i.warning=sum(is.na(cghseg.frame[i,]))
    if (i.warning == k) {i.remove=c(i.remove,i)}
  }
  if (length(i.remove)>0) {
    write(" ! Warning: no data for the following years when computing neigbourhood",file="")
    write(paste(" ",cghseg.frame$year[i.remove],collapse=" "),file="")
    cghseg.frame=cghseg.frame[-i.remove,]
  }
  Y = cghseg.frame
  save(Y,file="tmp/Y.sav")
  
  if (nrow(cghseg.frame)<10){cghseg.frame=F;write("Warning: common period too short (<10 years)",file="")}
  return(cghseg.frame)}
#####################

##############################################################################################################
home.input.data.check=function(list.of.files){
  ##############################################################################################################
  
  ## Reads data into a data frame meeting CLIMATOL
  ## Quality Checks requirements (Jose Guijarro, AEMET)
  #####################################################
  
  an.deb =  1000000
  an.fin = -1000000
  
  for (files in list.of.files){
    data = read.table(files)
    n    = nrow(data)
    tt   = as.vector(data[,1])
    test = diff(tt)
    if (sum(test !=1)>0){options(show.error.messages=F)
      write(paste("Non consecutive years in file",files),file="")
      stop()}
    if (tt[1]<an.deb) {an.deb=tt[1]}
    if (tt[n]>an.fin) {an.fin=tt[n]}
  }
  
  n           = (an.fin-an.deb+1)*12
  p           = length(list.of.files)
  
  check.frame = matrix(miss.flag,nr=n,nc=p)
  
  k=0
  for (files in list.of.files){
    k         = k+1
    data      = read.table(files)
    tt        = as.vector(data[,1])
    n.data    = length(tt)*12
    i.sh      = (tt[1]-an.deb)*12+1         # shift index for storage in check.frame
    data      = data[,-1]
    data      = as.vector(t(as.matrix(data)))
    check.frame[i.sh:(i.sh+n.data-1),k]=data
  }
  return(list(check.frame,an.deb,an.fin))}
########################################

###########################################################
data.check <- function(dat,ne,anyi,anyf,nm=12,nclust=100) {
  ###########################################################
  ### Basic network checks, based on CLIMATOL software
  ### Author: J. Guijarro
  ###########################################################
  
  #total number of data:
  na=anyf-anyi+1
  nd=na*nm
  x=(anyi*nm):(anyf*nm+nm-1)/12
  dat[dat==-999.9]=NA
  numdat  <- apply(!is.na(dat),1,sum)
  
  file.check=paste("meta/",head.str,net.str,"diagnostics",sep="")
  fig.file(file.check,10,12)
  plot(x,numdat,type="l",ylim=c(0,ne),col="blue",xlab="Year",ylab="Nr. of data",main=paste("Nr of",par.str,"data in network",net.str))
  grid()
  #boxplots of data:
  if(nm>1) { #boxplots for every nm (months, seasons, etc)
    dim(dat) <- c(nm,na,ne) #provisional redimensioning
    for(me in 1:nm) { #for every 'month' (or season, or ...)
      z <- data.frame(dat[me,,])
      names(z) <- 1:ne
      #lable of the month (if nm!=12, we will use only the number):
      if(nm==12) labm <- name.month[me] else labm <- me
      labm <- paste(" (",labm,")",sep="")
      if(ne>nclust) hist(as.matrix(z),xlab=par.str,main=paste("Data values of",par.str,labm,sep=" "),col="wheat")
      else {
        boxplot(z,ylab=paste(par.str,"(",unit.str,")"),main=paste("Data values of ",par.str,labm,sep=""),
                names=station.index,col="wheat",border=hsv(.7,1,.9),las=2,cex=.7)
        grid(col=grey(.4))
        abline(h=0)
      }
    }
    dim(dat) <- c(na*nm,ne) #restore the working dimensions
  }
  else { #only one graphic with all station boxplots:
    z <- data.frame(dat)
    names(z) <- 1:ne
    if(ne>nclust) hist(as.matrix(z),xlab=paste(par.str,"(",unit.str,")"),main=paste("Data values of",par.str),col="wheat")
    else {
      boxplot(z,ylab=paste(par.str,"(",unit.str,")"),main=paste("Data values of",par.str),
              col="wheat",border=hsv(.7,1,.9),las=2,cex=0.7)
      grid(col=grey(.4))
      abline(h=0)
    }
  }
  #histogram of all data together (to look at the distribution shape)
  main=paste("Histogram of",par.str,"for network",net.str)
  hist(dat,xlab=paste(par.str,"(",unit.str,")"),main=main,col=hsv(.4,1,.8))
  #correlogram of the first differenced monthly series (r <-> distancia)
  #(if more than nclust stations, only a random sample of nclust elements)
  if(ne>nclust) { splc <- sample(1:ne,nclust); nec <- nclust }
  else { splc <- 1:ne; nec <- ne }
  est.d <<- matrix(NA,nec,nec) #distance matrix
  for(i in 1:(nec-1)) {
    for(j in (i+1):nec) {
      dx <- lon.station[splc[i]]-lon.station[splc[j]]
      dy <- lat.station[splc[i]]-lat.station[splc[j]]
      #convert to km
      dx <- dx*111*cos((lat.station[splc[i]]+lat.station[splc[j]])*pi/360)
      dy <- dy*111
      dz <- (alt.station[splc[i]]-alt.station[splc[j]])*0.001
      d2 <- dx*dx+dy*dy+dz*dz #squared distance
      est.d[i,j] <<- sqrt(d2) #distance
      est.d[j,i] <<- est.d[i,j]  #(the matrix is symmetric)
    }
  }
  cat("\n")
  data <- dat[,splc] #copy of the data
  if(nm>1) { #compute the first differences month by month
    dim(data) <- c(nm,na,nec) #redimensioning
    difd <- array(NA,c(nm,na-1,nec)) #matriz of the first differencies
    for(i in 1:nm) { #for every month...
      difd[i,,] <- diff(data[i,,])
    }
    dim(difd) <- c(nd-nm,nec) #restoring dimensions
  }
  else difd <- diff(data) #global differentiation
  corm <- cor(difd,use="p") #correlation matrix
  #remove |r|==1 (probably due to stations with only 2 data in common):
  corm[corm==1] <- NA; corm[corm==-1] <- NA
  if(ne>nclust) main <- paste("Correlogram of first difference",nclust,"sampled series")
  else main <- "Correlogram of first difference series"
  xd <- as.vector(est.d); y <- as.vector(corm)
  plot(xd,y,xlim=c(0,max(est.d,na.rm=TRUE)),xlab="Distance (km)",
       ylab="Correlation coefficient",main=main,col="blue")
  grid(col=gray(.4))
  
  if(ne>2) {  #dendrogram of the stations
    dism <- dist(corm) #dissimilarity matrix
    #do not attempt clustering with NA's in the dissimilarity matrix:
    if(!sum(is.na(dism))) {
      hc <- hclust(dism,method="ward")
      if(ne>nclust) main <- paste("Dendrogram of",nclust,"sampled stations")
      else main <- "Dendrogram of station clusters"
      plot(hc,xlab="Stations",sub="",ylab="Dissimilarity",main=main)
      #we will cut the tree by the middle of the dissimilarity range
      #(but if the number of groups is greater than 9, we will incremet the
      #cut level by .1):
      if(ne<=nclust) { #(do not make clusters of a sample)
        cutlev <- mean(range(hc$height))
        repeat {
          ct <- cutree(hc,h=cutlev)
          nc <- length(levels(factor(ct)))
          if(nc<10) break
          cutlev <- cutlev + .1
        }
        if(nc>1) abline(h=cutlev,col="red",lty=2)
      }
    } else { nc <- 1; ct <- 1 } #just one cluster if dism has NA's
    if(ne>nclust) { nc <- 1; ct <- 1 } #just one cluster if >nclust stations
    #map of the stations:
    if(nc==1) { col="blue"; main=paste(par.str,"station locations") }
    else {
      col=rainbow(nc,1,.6)[ct]
      main=paste(par.str," station locations (",nc," clusters)",sep="")
    }
    #   if(ne>nclust) { #plot symbols if more than nclust stations
    map.x.lim = range(lon.station)
    map.x.ext = abs(diff(map.x.lim))/8
    map.x.lim[1]=map.x.lim[1]-map.x.ext
    map.x.lim[2]=map.x.lim[2]+map.x.ext
    map.y.lim = range(lat.station)
    map.y.ext = abs(diff(map.y.lim))/8
    map.y.lim[1]=map.y.lim[1]-map.y.ext
    map.y.lim[2]=map.y.lim[2]+map.y.ext
    
    map(xlim=map.x.lim,ylim=map.y.lim,
        xlab="Longitude (degree)",ylab="Latitude (degree)",fill=T,
        col="grey",main=paste("NETWORK",net.str))
    text(lon.station,lat.station,paste("+",station.index),col=col)
    map.axes()
    title(paste("NETWORK",net.str))
    map.grid()
    
    #      plot(lon.station,lat.station,asp=1/(cos(mean(lat.station)*pi/180)),
    #           xlab="Longitude (deg)",
    #           ylab="Latitude (deg)",col=col,pch=ct,main=main)
    #     else plot(est.c[,1:2],asp=1,xlab="X (km)",ylab="Y (km)",col=col,pch=ct,main=main)
    #    }
    #    else { #lable with numbers up to nclust stations
    #      plot(lon.station,lat.station,type="n",asp=1/(cos(mean(lat.station)*pi/180)),xlab="Longitude (deg)",ylab="Latitude (deg)",main=main)
    #     else plot(est.c[,1:2],type="n",asp=1,xlab="X (km)",ylab="Y (km)",main=main)
    #      text(lon.station,lat.station,labels=1:ne,pos=4,col=col)
    #    }
    grid(col=gray(.4))
  }
  rm(data,difd,corm) #remove auxiliary objects
  dev.off()
}

############################################################
qc.serie=function(y1,t1,y.comp=NA,y2=NA,t2=NA,flag.first=T){
  ############################################################
  
  ###############################################
  ## Data matrix computation for Fast QC        #
  ###############################################
  
  if (flag.first == T){
    is.na(t1) = y1==miss.flag
    is.na(y1) = y1==miss.flag
    y.comp    = cbind(t1,y1)} else {
      t1 = y.comp[,1]
      y1 = y.comp[,2]
      is.na(y2) = y2==miss.flag
      is.na(t2) = y2==miss.flag    
      y.comp    = cbind(y.comp,rep(NA,nrow(y.comp)))
      i.col     = ncol(y.comp)
      t.comp    = t1[which(!is.na(t2[match(t1,t2)]))]
      if (!is.na(t.comp[1])) {
        y.comp[match(t.comp,t1),i.col]=y2[match(t.comp,t2)]}
    }
  return(y.comp)}
###############

###################################################
comparison.serie=function(y1,t1,y2,t2,flag.corr=F){
  ###################################################
  
  ###############################################
  ## Computation of relative comparison series  #
  ## Options : difference "a"                   #
  ##           ratio "r"                        #
  ##           logratio "log"                   #
  ## Handles missing values                     #
  ## If flag.corr=TRUE, just computes the       #
  ## correlation of the first order difference  #
  ## series                                     #
  ###############################################
  
  t1.begin  = t1[1]
  t1.end    = t1[length(t1)]
  t2.begin  = t2[1]
  t2.end    = t2[length(t2)]
  
  is.na(t1) = y1==miss.flag
  is.na(t2) = y2==miss.flag
  is.na(y1) = y1==miss.flag
  is.na(y2) = y2==miss.flag
  
  ### Removes null or negative values
  ### in order to avoid +/-Inf in log or ratio comparison
  if (comp.option != "a"){
    is.na(t1) = y1<=0
    is.na(t2) = y2<=0
    is.na(y1) = y1<=0
    is.na(y2) = y2<=0
  }
  
  t.comp    = t1[which(!is.na(t2[match(t1,t2)]))]
  
  if (!is.na(t.comp[1])) {
    if (flag.corr == FALSE) {
      if (comp.option=="a"  ) {y.comp = y1[match(t.comp,t1)]-y2[match(t.comp,t2)]}
      if (comp.option=="r"  ) {y.comp = y1[match(t.comp,t1)]/y2[match(t.comp,t2)]}
      if (comp.option=="log") {y.comp = log(y1[match(t.comp,t1)]/y2[match(t.comp,t2)])}
    } else {y.comp = cor(diff(y1[match(t.comp,t1)]),diff(y2[match(t.comp,t2)]));t.comp = NA}
  } else {t.comp = NA; y.comp=NA}
  return(cbind(t.comp,y.comp))}
#############################


#####################################
distance.serie=function(list.file,j){
  #####################################
  
  ###############################################
  ## Given a set of series and their lat lon    #
  ## computes the list of neighbours whose      #
  ## geographic distance to the candidate is    #
  ## smaller or equal than inter.number (km)    #
  ## At least the n.min closest series are      #
  ## returned, even if d>dmax                   #
  ###############################################
  
  lst           = list.file
  n.serie       = length(lst)
  radius        = 6400
  distance      = sqrt((x.sta-x.sta[j])^2+(y.sta-y.sta[j])^2)
  distance[j]   = 0
  
  i.distance    = sort(distance,decreasing=F,index.return=T)$ix
  distance      = distance[i.distance]
  station.i     = station.index[i.distance]
  station.n     = station.name[i.distance]
  lst           = lst[i.distance]
  n.select      = length(which(distance < inter.number))
  
  if (n.select-1 >= n.min){lst=lst[2:n.select];i.select=2:n.select
  } else {i.select=2:min(n.min+1,n.serie);lst=lst[i.select]}
  
  if(file.exists("tmp/file_neighbours.txt")) { file.remove("tmp/file_neighbours.txt") }
  for (i in i.select){ 
    write(sprintf("%8s %8.2f %s",station.i[i],distance[i],paste("km  ",station.n[i])),file="")
    write(sprintf("%8s %8.2f %s",station.i[i],distance[i],paste("km ",station.n[i])),file="tmp/file_neighbours.txt",append=T) }
  write(" ",file="")
  
  return(lst)}


###################################################
correlogram.serie=function(work.out,month.index){
  ###################################################
  
  ###############################################
  ## Given a set of series and their lat lon    #
  ## computes the list of neighbours whose      #
  ## geographic distance to the candidate is    #
  ## smaller or equal than inter.number (km)    #
  ## At least the n.min closest series are      #
  ## returned, even if d>dmax                   #
  ###############################################
  
  n          = nrow(work.out)
  n.sta      = nlevels(work.out$g)
  n.year.sta = numeric(n.sta)
  g.sta      = levels(work.out$g)
  cov.matrix = matrix(0,nr=n,nc=n)
  geo.dist   = matrix(0,nr=n.sta,nc=n.sta)
  geo.cov    = geo.dist
  v.dist     = numeric()
  v.cov      = numeric()
  x.sta      = unique(work.out$x.lon)
  y.sta      = unique(work.out$y.lat)
  
  
  for (i in 1:n.sta) {
    g.i           = subset(work.out,g==g.sta[i])[,c(1,5)]
    n.year.sta[i] = nrow(g.i)
    t1            = g.i$fmu
    y1            = g.i$e
    for (j in i:n.sta) {
      g.j          = subset(work.out,g==g.sta[j])[,c(1,5)]
      t2           = g.j$fmu
      y2           = g.j$e          
      t.comp       = t1[which(!is.na(t2[match(t1,t2)]))]
      geo.cov[i,j] = cov(y1[match(t.comp,t1)],y2[match(t.comp,t2)])
      geo.cov[j,i] = geo.cov[i,j]
    }}
  
  for (i in 1:(n.sta-1)) {
    for (j in (i+1):n.sta) {
      geo.dist[i,j]=sqrt((x.sta[i]-x.sta[j])^2+(y.sta[i]-y.sta[j])^2)
      geo.dist[j,i]=geo.dist[i,j]
      v.dist = c(v.dist,geo.dist[i,j])
      v.cov  = c(v.cov,geo.cov[i,j])
    }}
  
  plot(v.dist,v.cov,main=paste(month.index))
  points(v.dist,fitted(loess(v.cov~v.dist)),pch="+",col="red")
  
  j = 1
  diagos = numeric()
  for (i.1 in 1:(n.sta-1)){
    k=j-1
    diagos=c(diagos,rep(geo.cov[i.1,i.1],n.year.sta[i.1]))
    for (i.2 in (i.1+1):n.sta){
      k=k+n.year.sta[i.2-1]
      print(paste(i.1,n.year.sta[i.1],i.2,n.year.sta[i.2]))
      print(j:(j+n.year.sta[i.1]-1))
      print((k+1):(k+n.year.sta[i.2]))
      cov.matrix[j:(j+n.year.sta[i.1]-1),(k+1):(k+n.year.sta[i.2])]=geo.cov[i.1,i.2]
      cov.matrix[(k+1):(k+n.year.sta[i.2]),j:(j+n.year.sta[i.1]-1)]=geo.cov[i.1,i.2]
    }
    j=j+n.year.sta[i.1]
  }
  i.1=n.sta
  diagos=c(diagos,rep(geo.cov[i.1,i.1],n.year.sta[i.1]))
  diag(cov.matrix)=diagos
  solve(cov.matrix)
  #print(format(round(cov.matrix,2)))
  cov.matrix = as.data.frame(cov.matrix)
  write.table(format(round(cov.matrix,2)),file="toto.txt",quote=F,row=F,col=F)
  return(cov.matrix)}
###################

###########################################
correlation.distance=function(list.file,j){
  ###########################################
  
  ###############################################
  ## Given a set of series                      #
  ## computes the list of neighbours whose      #
  ## correlation distance to the candidate is   #
  ## greater or equal to cor.min                #
  ## At least the n.min most correlated series  #
  ## are returned, even if r<rmin               #
  ###############################################
  lst        = list.file
  c.min      = inter.number
  n.serie    = length(lst)
  distance   = rep(0,n.serie)
  list.index = 1:n.serie
  data.tmp   = read.serie(lst[j],13)
  tj         = data.tmp[,1]
  yj         = data.tmp[,2]
  for (i in list.index){
    data.tmp    = read.serie(lst[i],13)
    ti          = data.tmp[,1]
    yi          = data.tmp[,2]
    distance[i] = comparison.serie(yj,tj,yi,ti,TRUE)[2]
  }
  
  i.distance    = sort(distance,decreasing=T,index.return=T)$ix
  distance      = distance[i.distance]
  station.index = station.index[i.distance]
  station.name  = station.name[i.distance]
  lst           = lst[i.distance]
  n.select      = length(which(distance > c.min))
  
  if (n.select-1 >= n.min){lst=lst[2:n.select];i.select=2:n.select
  } else {i.select=2:min(n.min+1,n.serie);lst=lst[i.select]}
  
  if(file.exists("tmp/file_neighbours.txt")) { file.remove("tmp/file_neighbours.txt") }
  for (i in i.select){
    write(sprintf("%8s %5.3f %s",station.index[i],distance[i],station.name[i]),file="")
    write(sprintf("%8s %5.3f %s",station.index[i],distance[i],station.name[i]),file="tmp/file_neighbours.txt",append=T)}
  write(" ",file="")
  
  return(lst)}
############

##########################
Cmat=function(Ya,Ys=NULL){
  ##########################
  
  ##################################################
  ### Generates cost matrix used in dynprog function
  ### Allows joint segmentation of Ya(nnual) and Ys 
  ### (if Ys(easonal) is provided)
  ### Uses ACMANT parametrization
  ##################################################
  
  n    = length(Ya)
  matD = matrix(0,n,n)
  
  if (is.null(Ys)){
    for (i in 0:(n-1)){
      for (j in (i+1):n){
        matD[i+1,j]=sum((Ya[(i+1):j]-mean(Ya[(i+1):j]))^2)}}
  } else {
    for (i in 0:(n-1)){
      for (j in (i+1):n){
        matD[i+1,j] = sum((Ya[(i+1):j]-mean(Ya[(i+1):j]))^2)+
          0.5*sum((Ys[(i+1):j]-mean(Ys[(i+1):j]))^2)}}
  }
  return(matD=matD)}
##################

###############################
dynprog=function(matD,Kmax=50){
  ###############################
  
  ################################################################
  ### Dynamic Programming algorithm for segmentation in the mean
  ### of a gaussian process.
  ###
  ### package. Implemented to allow ACMANT joint detection.
  ### Input:    matD cost matrix produced by Cmat function
  ###           Kmax max number of segments
  ###
  ### Translated from original Marc Lavielle MATLAB code 18-08-03
  ################################################################
  N=nrow(matD)
  I   = matrix(Inf,Kmax,N)
  t   = matrix(0,Kmax-1,N)
  I[1,]=matD[1,]
  if (Kmax>2){
    for (k in 2:(Kmax-1)){
      for (L in k:N){
        t[k-1,L] = which.min(I[k-1,1:(L-1)]+t(matD[2:L,L]))
        I[k,L]   = (I[k-1,1:(L-1)]+t(matD[2:L,L]))[t[k-1,L]]
      }
    }
  }
  t[Kmax-1,N] = which.min(I[Kmax-1,1:(N-1)]+t(matD[2:N,N]))
  I[Kmax,N]   = (I[Kmax-1,1:(N-1)]+t(matD[2:N,N]))[t[Kmax-1,N]]
  
  J.est=I[,N]
  
  ### Compute the change-points instants
  t.est = diag(N,Kmax)
  for (K in 2:Kmax){
    for (k in seq(K-1,1,-1)){           
      t.est[K,k] = t[k,t.est[K,k+1]]
    }
  }
  
  ### Removes close change-points
  for (i in 2:Kmax){
    for (j in 1:i){
      for (k in 1:(i-1)){
        if (((abs(t.est[i,j]-t.est[i-1,k])<2) & (t.est[i,j]!=t.est[i-1,k]))
            |  (t.est[i,j]-t.est[i,k]==1)){
          t.est[i,j]=0}}}}
  
  return(list(J.est=J.est,t.est=t.est))}
######################################

###########################
polyg.fill=function(tm,an,i.miss){
  #############################################
  ### Ancillary function for visu             #
  ### Fills with col.str1 and col.str2 colors #
  #############################################
  
  n   = length(tm)
  tmm = mean(tm,na.rm=T)
  
  tm[i.miss]=tmm
  an.cross=an[1]
  tm.cross=tm[1]
  
  for (i in 2:n){
    if ((tm[i]>tmm & tm[i-1]<tmm) | (tm[i]<tmm & tm[i-1]>tmm)){
      pente=tm[i]-tm[i-1]
      orig =tm[i-1]
      an.cross=c(an.cross,an[i-1]+(tmm-orig)/pente,an[i])
      tm.cross=c(tm.cross,tmm,tm[i])
    } else {an.cross=c(an.cross,an[i]);tm.cross=c(tm.cross,tm[i])}
  }
  tm.cross=c(tmm,tm.cross,tmm)
  tm.cross.sup=tm.cross
  tm.cross.inf=tm.cross
  an.cross=c(an[1],an.cross,an[n])
  k=which(tm.cross < tmm)
  tm.cross.sup[k]=tmm
  polygon(an.cross,tm.cross.sup,col=col1.str)
  k=which(tm.cross > tmm)
  tm.cross.inf[k]=tmm
  polygon(an.cross,tm.cross.inf,col=col2.str)
  lines(an,tm,lwd=2)
  
  decal.sup=0
  decal.inf=0
  if (length(i.miss)>0){
    sup.tm=max(tm)
    inf.tm=min(tm)
    for (k in i.miss){
      if (an[k]==an[1]) {decal.inf=-0.15} else {decal.inf=.5}
      if (an[k]==an[n]) {decal.sup=-0.15} else {decal.sup=.5}
      x.poly=c(an[k]-decal.inf,an[k]-decal.inf,an[k]+decal.sup,an[k]+decal.sup,an[k]-decal.inf)
      y.poly=c(inf.tm,sup.tm,sup.tm,inf.tm,inf.tm)
      polygon(x.poly,y.poly,col="white",border="white")
    }
    points(an[i.miss],rep(tmm,length(i.miss)),pch=".",cex=3)
  }}
##

##############################
uni.seg = function(y,kmax=80){
  ########################################################
  #### Change-point detection ############################
  ########################################################
  #### Replaces function "critpen.lik"                   #
  #### Uses Zhang & Siegmund criterion                   #
  #### Uses cghseg package (Franck Picard et al.)        #
  #### and function "uniseg"                             #
  ########################################################
  sink(type="output",file=file.sink)
  CGHd        = new("CGHdata",y)
  CGHo        = new("CGHoptions")
  uni.out     = uniseg(CGHd,CGHo)
  brk         = getbp(uni.out)
  n.brk       = sum(brk$Y[,2])-1
  sink()
  if (n.brk == 0){brk=0} else {brk=which(brk$Y[,2]==1)[1:n.brk]}
  return(brk)}
############

###################################
jump.posterior = function(x,sig_a){
  ###################################
  
  ########################################################
  #### Posterior distribution of a unique change-point   #
  #### in sequence x using bayesian model                #
  #### Hannart & Naveau (WRR,2009)                       #
  #### sig_a  = prior jump amplitude                     #
  ########################################################
  #### Based on original MATLAB code by A. Hannart       #
  ########################################################
  
  nt      = length(x)
  tau1    = 1:(nt-1)
  zz      = seq(nt-1,1,-1)
  lambda1 = tau1*zz/nt^2
  eta1    = sqrt(lambda1)
  s2      = var(x)
  
  ### Difference in partial means
  delta_t = (tau1*sum(x) - nt*cumsum(x[tau1]))/(tau1*zz)
  r_t     = (1-lambda1*(delta_t^2)/s2)
  rn_t    = (r_t^(-0.5*(nt-2)))/eta1
  s_t     = s2*r_t/(nt*lambda1)
  norm_t  = dnorm(0,delta_t,sqrt(s_t+sig_a^2))
  pn_t    = rn_t*norm_t
  
  ### Posterior of change point location and bayes factor
  bf = sum(pn_t)/nt
  p_t = pn_t/(nt*bf) 
  return(p_t)}
############
############################
joint.detection=function(Y){
  ############################
  
  #######################################################
  ####         PROCEDURE FOR JOINT DETECTION            #  
  ####            CALLED FROM j.detect PROCEDURE        #   
  #######################################################
  #### Uses cghseg (Picard, et al.)                     #  
  #######################################################
  
  cgh.year      = Y[,1]
  Y             = Y[,-1];save(Y,file="tmp/Y.sav")
  M             = ncol(Y)
  sink(type="output",file=file.sink)
  CGHd     = new("CGHdata",Y)
  CGHo     = new("CGHoptions")
  wavenorm(CGHo)= "position"
  CGHo           = new("CGHoptions")     
  out      = multiseg(CGHd,CGHo)
  sink()
  return(out)}
############

########################################
##  MAIN PROCEDURES - CALLED BY MENU  ##
########################################

###############################################################
###############################################################
###############################################################
fqc = function(list.file,n.file,head.str){
  ###############################################################
  ###############################################################
  ###############################################################
  
  ########################################################
  ####       MAIN PROCEDURE FOR FAST QUALITY CONTROL     #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  ####         Plots Monthly anomalies of candidate      #
  ####                 vs Neighbours                     #
  ########################################################
  
  if (comp.option == "a"){
    if (unit.str=="c") {y.label=expression(paste(delta,"T (",degree,"C )"))
    } else {y.label=paste(par.str," DIFF. (",unit.str,")",sep="")}
  } else {y.label=""}
  
  ## Loop over files
  ##################
  for (j in 1:n.file){
    
    file = list.file[j]
    write("                      ",file="")
    write(paste(station.index[j],station.name[j]),file="")
    write("=====================================",file="")
    
    list.neighbour = list.file[-j]
    if (inter.option == "g") {list.neighbour=distance.serie(list.file,j)}
    if (inter.option == "c") {list.neighbour=correlation.distance(list.file,j)}
    n.neighbour = length(list.neighbour)
    
    if (comp.option == "a"){
      if (unit.str=="c") {text.label=substitute(paste(delta,"par.str (",degree,"C)",sep=""))
      } else {text.label=substitute(paste(delta,par.str," (",unit.str,")",sep=""))}
    } else {text.label=substitute(paste("./. (",percentage,")"))}
    
    file.control=paste("meta/control_",head.str,station.index[j],sep="")
    fig.file(file.control,he=20,wi=12)
    par(mfrow=c(13,1),mar=c(3,4,2,1))
    
    for (k in 1:13) {
      
      data.tmp     = read.serie(file,k)
      y.candidate  = data.tmp[,2]
      t.candidate  = data.tmp[,1]
      t.lim        = c(t.candidate[1],t.candidate[nrow(data.tmp)])
      y.comp       = qc.serie(y.candidate,t.candidate)
      
      ## Loop over neighboors
      #######################
      
      for (file.ref in list.neighbour){
        brk      = 0
        n.brk    = 0
        i.ref    = match(file.ref,list.file)
        
        data.tmp = read.serie(file.ref,k)
        y.ref    = data.tmp[,2]
        t.ref    = data.tmp[,1]
        y.comp   = qc.serie(y.candidate,t.candidate,y.comp,y.ref,t.ref,flag.first=F)}
      
      #####################################################
      ## Graphic outputs, once detection has been performed
      #####################################################
      
      gt.lim      = c(10*floor(t.lim[1]/10),10*round((t.lim[2]+5)/10))
      v.ticks     = seq(gt.lim[1],gt.lim[2],5)
      month.str   = as.character(k)
      if (k <  10) { month.str = paste("0",month.str,sep="")}
      if (k == 13) { month.str = "ANNUAL"}
      t.comp      = y.comp[,1]
      y.comp      = y.comp[,-1]
      means.comp  = colMeans(y.comp,na.rm=T)
      
      for (i in 1:length(means.comp)){
        if (comp.option == "a"){
          y.comp[,i]=y.comp[,i]-means.comp[i]} else {
            if (means.comp[i]!=0){y.comp[,i]=y.comp[,i]/means.comp[i]}}}
      means.comp  = rowMeans(y.comp,na.rm=T)
      y.comp = y.comp-means.comp
      #improvement: type="o",		   
      plot(t.comp,y.comp[,1],xlab="",ylab=y.label,pch="o",
           xlim=gt.lim,xaxs="i",ylim=range(y.comp,na.rm=T),
           lab=c(10,5,7),col="red", type="o",
           main=paste(name.str[j],month.str))
      for (i in 2:(n.neighbour+1)) {points(t.comp,y.comp[,i],pch="+")}
    }
    dev.off()
  }}
##


###############################################################
rem.out=function(head.str){
  ###############################################################
  ###############################################################
  ###############################################################
  
  ########################################################
  ####       MAIN PROCEDURE FOR REMOVING OUTLIERS IN     #
  ####                   DATA FILES                      #
  ########################################################
  ####         Validated outliers are put "missing"      #
  ####                 in the files                      #
  ########################################################
  
  file.outlier = paste(head.str,net.str,"outliers.txt",sep="")
  if (file.exists(file.outlier)){
    outlier.frame = read.table(file.outlier,
                               colClasses=c("character","numeric","numeric","character"),quote="",sep="\t")
    
    if (nrow(outlier.frame)>0){
      file.list=levels(factor(outlier.frame$V1))
      for (file in file.list){
        data=read.table(file)
        outlier.tmp=subset(outlier.frame,V1==file)
        for (j in 1:nrow(outlier.tmp)){
          data[!is.na(match(data$V1,outlier.tmp$V2[j])),outlier.tmp$V3[j]+1]=miss.flag
        }
        
        write.table(format(round(data,1)),file=file,quote=F,row=F,col=F,sep="\t")
      }}
    write(" ",file="")
    write(paste("Number of outliers removed:",nrow(outlier.frame)),file="")
    write(" ",file="")
  } else {write(paste("Warning:",file.outlier,"does not exists"),file="")}}
### End function rem.out
########################

#############################################################################
#############################################################################
#############################################################################
j.detect = function(list.file,head.str,season,crit="CAU",flag.inter=F){
  #############################################################################
  #############################################################################
  #############################################################################
  
  ########################################################
  ####      MAIN PROCEDURE FOR JOINT DETECTION           #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  #### Uses cghseg (Picard et al.)                       #
  #### cghseg : climate factor treated as fixed effect   #
  ####          LS estimator                             #
  ########################################################
  
  graphics.off()
  flag.prev     = F
  file.detect   = paste(net.str,"detected.txt",sep="")
  
  index.sta     = character()
  type.brk      = character()
  year.sta      = numeric()
  month.sta     = numeric()
  name.sta      = character()
  meta.sta      = character()
  
  modif.year    = numeric()
  modif.index   = character()
  modif.station = character()
  modif.type    = character()
  modif.month   = numeric()
  modif.meta    = character()
  
  if (file.exists(file.detect)){
    prev.frame=read.table(file.detect,head=F,
                          colClasses=c("character","character","numeric","numeric","character","character"),sep="\t")
    flag.prev=T}
  
  if (inter.option == ""){
    Y=read.cghseg(list.file,13,list.file);cgh.year=Y[,1]
    if (!is.logical(Y)){
      out=joint.detection(Y)
    }
  }
  
  status.str = substr(head.str,1,1)
  if (status.str == "h"){status.str="(H)"} else {status.str = ""}
  
  
  #################################
  #################################
  ### COMPARISON PAIRWISE/JOINT ###
  #################################
  #################################
  
  flag.modif=F
  flag.action=T
  for (j in 1:length(list.file)){
    archive.detect=paste("tmp/detect_",head.str,station.index[j],".sav",sep="")
    archive.detect.djf=paste("tmp/detect_",head.str,station.index[j],"_DJF.sav",sep="")
    archive.detect.jja=paste("tmp/detect_",head.str,station.index[j],"_JJA.sav",sep="")
    archive.jdetect   = paste("tmp/jdetect_",head.str,station.index[j],".sav",sep="")
    
    write("                      ",file="")
    write(paste(station.index[j],station.name[j]),file="")
    write("=====================================",file="")
    
    if (inter.option != ""){
      list.neighbour = list.file[-j]
      if (inter.option == "g") {list.neighbour=c(list.file[j],distance.serie(list.file,j))}
      if (inter.option == "c") {list.neighbour=c(list.file[j],correlation.distance(list.file,j))}
      Y=read.cghseg(list.neighbour,13,list.file);cgh.year= Y[,1]
      out=joint.detection(Y)
    }
    
    if (inter.option == ""){patient=j} else {patient=1}
    j.brk    = getbp(out)[[patient]]
    j.nb.brk = sum(j.brk$bp)-1
    if (j.nb.brk >0){
      j.amp.brk = diff(getsegprofiles(out))[,patient]
      j.pos.brk = which(j.amp.brk !=0)
      j.amp.brk = j.amp.brk[j.pos.brk]
      j.pos.brk = cgh.year[j.pos.brk]
      save(j.pos.brk,j.amp.brk,file=archive.jdetect)
    }
    
    if (file.exists(archive.detect)) {
      load(archive.detect)
      if (j.nb.brk>0){range.amp=range(range.amp,j.amp.brk)}
      if (file.exists(archive.detect.jja) & file.exists(archive.detect.djf)) {
        load(archive.detect.djf)
        load(archive.detect.jja)
        gt.lim=range(gt.lim,gt.lim.jja,gt.lim.djf)
        range.amp = range(range.amp,range.amp.jja,range.amp.djf)}
      
      if (flag.inter==F){
        file=paste("fig/detect_",head.str,station.index[j],"_b",s.title,sep="")
        fig.file(file,he=6,wi=8)
      } else  {x11(wi=15,he=6)}     
      plot(gt.lim[1]:gt.lim[2],rep(0,diff(range(gt.lim))+1),
           xlab="",ylab=y.label,pch=".",
           xlim=gt.lim,xaxs="i",ylim=range.amp,lab=c(10,5,7),
           main=paste(name.str[j],status.str,s.title))
      if (flag.prev) {
        prev.brk = subset(prev.frame,V1==station.index[j])$V3
        abline(v=prev.brk,lwd=5,col="lightcyan2")}
      
      if (length(amp.brk)>0){
        sym.triangle = rep(25,length(amp.brk))
        points(pos.brk,amp.brk,pch=25,cex=2)}
      abline(v=v.ticks,lty=3,lwd=.75)
      
      if (file.exists(archive.detect.jja) & file.exists(archive.detect.djf)) {
        if (length(amp.brk.jja)>0){
          points(pos.brk.jja,amp.brk.jja,pch=24,bg="red1")}
        if (length(amp.brk.djf)>0){
          points(pos.brk.djf,amp.brk.djf,pch=24,bg="blue1")}
      }
      
      ### Adding joint detection ("o" symbols)
      if (j.nb.brk >0){
        points(j.pos.brk,j.amp.brk,cex=2,pch="o",col="green")
        points(j.pos.brk,j.amp.brk,cex=2,pch="+",col="green")
      }
      
      if (j.nb.brk >0){
        index.sta = c(index.sta,rep(station.index[j],j.nb.brk))
        name.sta  = c(name.sta,rep(station.name[j],j.nb.brk))
        year.sta  = c(year.sta,j.pos.brk)
        type.brk  = c(type.brk,rep("BREAK",j.nb.brk))
        month.sta = c(month.sta,rep(12,j.nb.brk))
        meta.sta  = c(meta.sta,rep("n",j.nb.brk))
      }
      
      if (flag.inter==T){
        xy=locator(40,type="p",ps=50,pch="+",col="red")
        if (!is.null(xy)){xy=data.frame(xy)
        xy=subset(xy,(xy$x>gt.lim[1]) & (xy$x<gt.lim[2]))}
        ### Removes duplicates
        if (!is.null(xy)>0){
          xy$x=round(xy$x,0)
          dup.flag=!(duplicated(xy$x,fromLast=T)|duplicated(xy$x))
          xy=xy[dup.flag,]}
        if (!is.null(xy)){if (nrow(xy)>0){
          flag.modif    = T   ### At least one modification in the list
          n.modif       = length(xy$x)
          modif.year    = c(modif.year,round(xy$x,0))
          modif.index   = c(modif.index,rep(station.index[j],n.modif))
          modif.station = c(modif.station,rep(station.name[j],n.modif))
          modif.type    = c(modif.type,rep("BREAK",n.modif))
          modif.meta    = c(modif.meta,rep("n",n.modif))
          modif.month   = c(modif.month,rep(12,n.modif))
        }}}
      dev.off()
    } else {write(paste("Warning: run pairwise detection on",station.index[j]),file="");flag.action=F}
  }
  
  ### Automatic constitution of detection file
  detect.frame = cbind(index.sta,type.brk,year.sta,month.sta,meta.sta,name.sta)
  detect.frame = data.frame(detect.frame)
  names(detect.frame)=c("V1","V2","V3","V4","V5","V6")
  
  if (flag.prev == T){detect.frame = rbind(prev.frame,detect.frame)}
  
  ### Takes into account interactively added/removed change-points
  ### Warning : unelegant code
  
  if (flag.inter==T & flag.modif==T){  ### Has to be an interactive session+effective modifications
    i.remove = numeric()
    j.remove = numeric()
    modif.frame = data.frame(V1=modif.index,V2=modif.type,V3=modif.year,V4=modif.month,V5=modif.meta,V6=modif.station)
    for (i in 1:nrow(detect.frame)){
      for (j in 1:nrow(modif.frame)){
        if (as.character(detect.frame[i,1])   == as.character(modif.frame[j,1])
            & as.character(detect.frame[i,2]) == as.character(modif.frame[j,2])
            & as.numeric(as.character(detect.frame[i,3]))==as.numeric(as.character(modif.frame[j,3]))){
          i.remove = c(i.remove,i)
          j.remove = c(j.remove,j)
        }}}
    ### Removes deselected detections
    if (length(i.remove)>0){detect.frame = detect.frame[-i.remove,]
    modif.frame  = modif.frame[-j.remove,]}
    ### Adds pairwise validated change-points
    names(modif.frame)=c("V1","V2","V3","V4","V5","V6")
    detect.frame=rbind(modif.frame,detect.frame)
  }
  
  ### Reorders by station and change-point dates
  detect.frame = detect.frame[order(as.character(detect.frame[,1]),
                                    as.character(detect.frame[,3]),
                                    as.numeric(detect.frame[,4])),]
  i.remove=numeric()
  if (nrow(detect.frame)>1){
    for (i in 2:nrow(detect.frame)){
      if (as.character(detect.frame[i-1,1]) == as.character(detect.frame[i,1]) &
          as.character(detect.frame[i-1,2]) == as.character(detect.frame[i,2]) &
          as.numeric(as.character(detect.frame[i-1,3])) == as.numeric(as.character(detect.frame[i,3]))){
        i.remove=c(i.remove,i)}}
    if (length(i.remove)>0){detect.frame=detect.frame[-i.remove,]}}
  
  if (flag.action){
    write.table(detect.frame,file=file.detect,quote=F,row=F,col=F,sep="\t")}
  return(detect.frame)
}

################################################################
################################################################
################################################################
change.month = function(list.file,list.file.ref,n.file,head.str){
  ################################################################
  ################################################################
  ################################################################
  
  ################################################################
  #### MAIN PROCEDURE FOR ASSESSING PRECISE MONTH OF CHANGES     #
  ####            CALLED FROM MAIN PROGRAM                       #
  ################################################################
  ####   Based on Peter Domonkos ACMANT idea             #
  ####  Requires Prehomogenized ref series               #
  ########################################################
  
  flag.prev   = F
  file.detect = paste(net.str,"detected.txt",sep="")
  if (file.exists(file.detect)){
    prev.frame=read.table(file.detect,head=F,
                          colClasses=c("character","character","numeric","numeric","character","character"),sep="\t")
    flag.prev=T
  }
  
  status.str = substr(head.str,1,1)
  if (status.str == "h"){status.str="(H)"} else {status.str = ""}
  
  
  new.frame   = data.frame()
  flag.modif  = F
  flag.action = T
  
  write(" ",file="")
  for (j in 1:length(list.file)){
    write("                      ",file="")
    write(paste(station.index[j],station.name[j]),file="")
    write("=====================================",file="")
    work.frame = subset(prev.frame,prev.frame$V1==station.index[j])
    
    if (nrow(work.frame) > 0){
      if (file.exists(list.file.ref[j])) {
        
        file        = list.file[j]
        file.ref    = list.file.ref[j]
        data.mo     = read.table(file)         ### Monthly data (candidate)
        n.mo        = nrow(data.mo)
        t.mo        = rep(data.mo[1,1]:data.mo[n.mo,1],each=12)
        month.mo    = rep(1:12,n.mo)
        
        data.mo   = as.matrix(data.mo)
        i.miss    = which(data.mo == miss.flag)
        data.mo[i.miss] = NA
        
        ### remove season cycle - on complete cases only
        ### not optimal, but in presence of missing values, still better than filters or MA
        
        month.means     = colMeans(data.mo[complete.cases(data.mo),2:13],na=T)
        for (i in 1:n.mo){data.mo[i,2:13]= data.mo[i,2:13]-month.means}
        y.mo      = numeric()
        for (i in 1:n.mo) {y.mo=c(y.mo,as.vector(data.mo[i,2:13]))}
        
        data.ref    = read.table(file.ref)     ### Monthly data (reference)
        n.ref       = nrow(data.ref)
        t.ref       = rep(data.ref[1,1]:data.ref[n.ref,1],each=12)
        month.ref   = rep(1:12,n.ref)
        
        data.ref  = as.matrix(data.ref)      
        i.miss    = which(data.ref == miss.flag)
        data.ref[i.miss] = NA
        month.means     = colMeans(data.ref[complete.cases(data.ref),2:13],na=T)
        for (i in 1:n.ref){data.ref[i,2:13]= data.ref[i,2:13]-month.means}
        
        y.ref     = numeric()
        for (i in 1:n.ref) {y.ref=c(y.ref,as.vector(data.ref[i,2:13]))}
        
        i.comp      = which(!is.na(t.ref[match(t.mo,t.ref)]))
        t.comp      = t.mo[i.comp]
        month.comp  = month.mo[i.comp]
        cp.list     = rbind(c(t.comp[1],month.comp[1]),
                            work.frame[,c(3,4)],
                            c(t.comp[length(t.comp)],month.comp[length(t.comp)]))
        meta.cp     = c("n",work.frame[,5],"n")
        cp.list     = cbind(cp.list,meta.cp)
        n.cp        = nrow(cp.list)
        
        if (!is.na(t.comp[1] & comp.option=="a")) {
          
          y.comp = y.mo[(t.mo %in% t.comp) & (month.mo %in% month.comp)]-y.ref[(t.ref %in% t.comp) & (month.ref %in% month.comp)]
          
          for (i.cp in 2:(n.cp-1)){
            if (cp.list[i.cp,3] =="n"){ 
              target.cp   = cp.list[i.cp,c(1,2)]
              i.work      = which(t.comp==cp.list[i.cp-1,1] & month.comp >cp.list[i.cp-1,2])
              i.work      = c(i.work,which(t.comp>cp.list[i.cp-1,1] & t.comp<cp.list[i.cp+1,1]-1))
              y.work      = y.comp[i.work]
              t.work      = t.comp[i.work]
              month.work  = month.comp[i.work]
              
              i.miss.work = is.na(y.work)
              y.work      = y.work[!i.miss.work]      
              t.work      = t.work[!i.miss.work]
              month.work  = month.work[!i.miss.work]
              if (length(y.work)>12){
                matD                 = Cmat(y.work)
                seg.out              = dynprog(matD,2)
                j.brk                = seg.out$t.est[2,1]
                if (j.brk == 0) {
                  write(paste(work.frame[i.cp-1,3],"/",work.frame[i.cp-1,4],"->",work.frame[i.cp-1,3],"/",work.frame[i.cp-1,4],sep=""),file="")
                } else {
                  if (abs(work.frame[i.cp-1,3]-t.work[j.brk])<=2){
                    write(paste(work.frame[i.cp-1,3],"/",work.frame[i.cp-1,4],"->",
                                t.work[j.brk],"/",month.work[j.brk],sep=""),file="")
                    work.frame[i.cp-1,3] = t.work[j.brk]
                    work.frame[i.cp-1,4] = month.work[j.brk]
                    cp.list[i.cp,1]      = t.work[j.brk]
                    cp.list[i.cp,2]      = month.work[j.brk]}
                }
              }
            } else {write(paste(work.frame[i.cp-1,3],"/",
                                work.frame[i.cp-1,4]," validated by metadata",sep=""))}
          }
          new.frame=rbind(new.frame,work.frame)
        } else {write(paste("Warning: run a correction round before",station.index[j]),file="");flag.action=F}
      }
    }}
  
  write.table(new.frame,file=file.detect,quote=F,row=F,col=F,sep="\t")
  return(new.frame)
}
#######################################
#######################################
#######################################
##########################################################################################
##########################################################################################
##########################################################################################
acmant.detect = function(list.file,list.file.ref,n.file,head.str,flag.inter=F){
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  
  ########################################################
  ####    MAIN PROCEDURE FOR ACMANT-STYLE DETECTION      #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  ####   Based on Peter Domonkos ACMANT principle        #
  ####  Prehomogenized ref series + joint detection on   #
  ####  annual+seasonal series                           #
  ########################################################
  
  graphics.off()
  flag.prev   = F
  file.detect = paste(net.str,"detected.txt",sep="")
  if (file.exists(file.detect)){
    prev.frame=read.table(file.detect,head=F,
                          colClasses=c("character","character","numeric","numeric","character","character"),sep="\t")
    flag.prev=T}
  
  status.str = substr(head.str,1,1)
  if (status.str == "h"){status.str="(H)"} else {status.str = ""}
  # if (file.exists(file.detect)){file.remove(file.detect)}
  
  modif.year    = numeric()
  modif.index   = character()
  modif.station = character()
  modif.type    = character()
  modif.meta    = character()
  modif.month   = numeric()
  
  index.sta     = character()
  type.brk      = character()
  year.sta      = numeric()
  month.sta     = numeric()
  name.sta      = character()
  meta.sta      = character()
  
  ########################################
  ########################################
  ### COMPARISON PAIRWISE/JOINT/ACMANT ###
  ########################################
  ########################################
  
  flag.modif  = F
  flag.action = T
  
  write(" ",file="")
  for (j in 1:length(list.file)){
    write("                      ",file="")
    write(paste(station.index[j],station.name[j]),file="")
    write("=====================================",file="")
    archive.detect     = paste("tmp/detect_",head.str,station.index[j],".sav",sep="")
    archive.detect.djf = paste("tmp/detect_",head.str,station.index[j],"_DJF.sav",sep="")
    archive.detect.jja = paste("tmp/detect_",head.str,station.index[j],"_JJA.sav",sep="")
    archive.jdetect    = paste("tmp/jdetect_",head.str,station.index[j],".sav",sep="")
    if (file.exists(archive.detect)) {
      load(archive.detect)
      if (file.exists(archive.detect.jja) & file.exists(archive.detect.djf)) {
        load(archive.detect.djf)
        load(archive.detect.jja)
        gt.lim=range(gt.lim,gt.lim.jja,gt.lim.djf)
        range.amp = range(range.amp,range.amp.jja,range.amp.djf)}
      
      if (flag.inter==F){
        file=paste("fig/detect_",head.str,station.index[j],"_b",sep="")
        fig.file(file,he=6,wi=8)
      } else  {x11(wi=15,he=6)}     
      plot(gt.lim[1]:gt.lim[2],rep(0,diff(range(gt.lim))+1),
           xlab="",ylab=y.label,pch=".",xlim=gt.lim,xaxs="i",
           ylim=range.amp,lab=c(10,5,7),main=paste(name.str[j],status.str))
      abline(v=v.ticks,lty=3,lwd=.75)
      
      if (flag.prev) {
        prev.brk = subset(prev.frame,V1==station.index[j])$V3
        #improvement:
        #             abline(v=prev.brk,lwd=4,col="lightgrey")}
        abline(v=prev.brk,lwd=5,col="lightcyan2")}
      
      if (length(amp.brk)>0){
        sym.triangle = rep(25,length(amp.brk))
        points(pos.brk,amp.brk,pch=25,cex=2)}
      
      if (file.exists(archive.detect.jja) & file.exists(archive.detect.djf)) {
        if (length(amp.brk.jja)>0){
          points(pos.brk.jja,amp.brk.jja,pch=24,bg="red1")}
        if (length(amp.brk.djf)>0){
          points(pos.brk.djf,amp.brk.djf,pch=24,bg="blue1")}
      }
      
      if (file.exists(archive.jdetect)){
        load(archive.jdetect)
        points(j.pos.brk,j.amp.brk,cex=2,pch="o",col="lightgrey")
        points(j.pos.brk,j.amp.brk,cex=2,pch="+",col="lightgrey")
      }
      
      if (file.exists(list.file.ref[j])) {
        
        ### Adding ACMANT detection 
        ###########################
        ###########################
        ###########################
        file     = list.file[j]
        file.ref = list.file.ref[j]
        data.tmp = read.serie(file,13)      ### Annual data
        y.a      = data.tmp[,2]
        t.a      = data.tmp[,1]
        t.lim    = c(t.a[1],t.a[nrow(data.tmp)])
        data.tmp = read.serie(file,20)      ### Seasonal data
        y.s      = data.tmp[,2]
        t.s      = data.tmp[,1]
        
        data.tmp = read.serie(file.ref,13)
        y.a.ref  = data.tmp[,2]
        t.a.ref  = data.tmp[,1]
        data.tmp = read.serie(file.ref,20)
        y.s.ref  = data.tmp[,2]
        t.s.ref  = data.tmp[,1]
        
        comp.a   = comparison.serie(y.a,t.a,y.a.ref,t.a.ref,FALSE)
        comp.s   = comparison.serie(y.s,t.s,y.s.ref,t.s.ref,FALSE)
        
        flag.ok  = F
        if (nrow(comp.a)==nrow(comp.s)){
          if (sum(comp.a[,1]==comp.s[,1])==nrow(comp.a)){
            flag.ok = T
            Ya      = comp.a[,2]
            Ys      = comp.s[,2]
            t.a     = comp.a[,1]
            kmax    = floor(nrow(comp.a)/4)
          }}
        if (flag.ok == T){
          matD       = Cmat(Ya,Ys)
          seg.out    = dynprog(matD,kmax+1)
          pen.lik    = log(seg.out$J.est[2:(kmax+1)]/seg.out$J.est[1])
          
          for (k in 1:kmax) {pen.lik[k]=pen.lik[k]+2*k*log(nrow(comp.a))/(nrow(comp.a)-1)}
          if (min(pen.lik) >= 0) {j.brk=0;kopt=0} else {
            kopt=which.min(pen.lik);j.brk=seg.out$t.est[kopt+1,1:kopt]
            j.brk=j.brk[j.brk !=0]    
            j.nb.brk = length(j.brk)
            abline(v=t.a[j.brk],lwd=2,lty=3)
            index.sta = c(index.sta,rep(station.index[j],j.nb.brk))
            name.sta  = c(name.sta,rep(station.name[j],j.nb.brk))
            year.sta  = c(year.sta,t.a[j.brk])
            meta.sta  = c(meta.sta,rep("n",j.nb.brk))
            type.brk  = c(type.brk,rep("BREAK",j.nb.brk))
            month.sta = c(month.sta,rep(12,j.nb.brk))
          }}
        
        if (flag.inter==T){
          xy=locator(40,type="p",ps=50,pch="+",col="red")
          if (!is.null(xy)){xy=data.frame(xy)
          xy=subset(xy,(xy$x>t.a[1]) & (xy$x<t.a[length(t.a)]))}
          ### Removes duplicates
          if (!is.null(xy)>0){
            xy$x=round(xy$x,0)
            dup.flag=!(duplicated(xy$x,fromLast=T)|duplicated(xy$x))
            xy=xy[dup.flag,]}
          if (!is.null(xy)){if (nrow(xy)>0){
            flag.modif    = T   ### At least one modification in the list
            n.modif       = length(xy$x)
            modif.year    = c(modif.year,xy$x)
            modif.index   = c(modif.index,rep(station.index[j],n.modif))
            modif.station = c(modif.station,rep(station.name[j],n.modif))
            modif.type    = c(modif.type,rep("BREAK",n.modif))
            modif.meta    = c(modif.meta,rep("n",n.modif))
            modif.month   = c(modif.month,rep(12,n.modif))
          }}}} else {write(paste("Warning: prior ACMANT run a correction round",station.index[j]),file="");flag.action=F}
      
      dev.off()} else {write(paste("Warning: prior ACMANT run pairwise detection on",station.index[j]),file="");flag.action=F}
  }
  
  ### Automatic constitution of detection file
  
  detect.frame = cbind(index.sta,type.brk,year.sta,month.sta,meta.sta,name.sta)
  detect.frame = data.frame(detect.frame)
  names(detect.frame)=c("V1","V2","V3","V4","V5","V6")
  
  if (flag.prev == T){detect.frame = rbind(prev.frame,detect.frame)}
  
  ### Takes into account interactively added/removed change-points
  ### Warning : unelegant code
  
  if (flag.inter==T & flag.modif==T){  ### Has to be an interactive session+effective modifications
    i.remove = numeric()
    j.remove = numeric()
    modif.frame = data.frame(V1=modif.index,V2=modif.type,V3=modif.year,V4=modif.month,V5=modif.meta,V6=modif.station)
    for (i in 1:nrow(detect.frame)){
      for (j in 1:nrow(modif.frame)){
        if (as.character(detect.frame[i,1])   == as.character(modif.frame[j,1])
            & as.character(detect.frame[i,2]) == as.character(modif.frame[j,2])
            & as.numeric(as.character(detect.frame[i,3]))==as.numeric(as.character(modif.frame[j,3]))){
          i.remove = c(i.remove,i)
          j.remove = c(j.remove,j)
        }}}
    ### Removes deselected detections
    if (length(i.remove)>0){detect.frame = detect.frame[-i.remove,]
    modif.frame  = modif.frame[-j.remove,]}
    ### Adds pairwise validated change-points
    names(modif.frame)=c("V1","V2","V3","V4","V5","V6")
    detect.frame=rbind(modif.frame,detect.frame)
  }
  
  ### Reorders by station and change-point dates
  detect.frame = detect.frame[order(as.character(detect.frame[,1]),
                                    as.character(detect.frame[,3]),
                                    as.numeric(detect.frame[,4])),]
  i.remove=numeric()
  if (nrow(detect.frame)>1){
    for (i in 2:nrow(detect.frame)){
      if (as.character(detect.frame[i-1,1]) == as.character(detect.frame[i,1]) &
          as.character(detect.frame[i-1,2]) == as.character(detect.frame[i,2]) &
          as.numeric(as.character(detect.frame[i-1,3])) == as.numeric(as.character(detect.frame[i,3]))){
        i.remove=c(i.remove,i)}}
    if (length(i.remove)>0){detect.frame=detect.frame[-i.remove,]}}
  
  if (flag.action){ ### only if ACMANT computation has occured
    write.table(detect.frame,file=file.detect,quote=F,row=F,col=F,sep="\t")}
  return(detect.frame)
}

###############################################################
###############################################################
###############################################################
detect = function(list.file,n.file,head.str,crit="CAU",season){
  ###############################################################
  ###############################################################
  ###############################################################
  
  ########################################################
  ####            MAIN PROCEDURE FOR DETECTION           #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  #### Uses cghseg package (Franck Picard et al.)        #
  #### Uses uni.seg (Olivier Mestre)                     #
  #### + other ancillary functions                       #
  ########################################################
  
  graphics.off()
  
  ## Storage matrices for results
  ## sd : standard dev of noise
  ## n.brk : number of breaks when comparing i to j
  ## brk.arr : storage array for change-points
  #################################################
  
  sd.mat     = matrix(0,n.file,n.file)
  n.brk.mat  = matrix(NA,n.file,n.file)
  brk.arr    = array(NA,c(n.file,n.file,max.brk))
  flag.done  = matrix(FALSE,n.file,n.file)
  if (comp.option == "a"){
    if (unit.str=="c") {y.label=expression(paste(delta,"T (",degree,"C )"))
    } else {y.label=paste(par.str," DIFF. (",unit.str,")",sep="")}
  } else {y.label=""}
  
  status.str = substr(head.str,1,1)
  if (status.str == "h"){status.str="(H)"} else {status.str = ""}
  
  season.str  = ""
  s.title     = ""
  
  
  if (season == 16) {season.str = "_season";s.title="DJF"}
  if (season == 17) {season.str = "_season";s.title="MAM"}
  if (season == 18) {season.str = "_season";s.title="JJA"}
  if (season == 19) {season.str = "_season";s.title="SON"}
  if (season == 20) {season.str = "_season";s.title="ACMANT"}
  
  if (season == 1) {season.str = "_season";s.title="JAN"}
  if (season == 2) {season.str = "_season";s.title="FEB"}
  if (season == 3) {season.str = "_season";s.title="MAR"}
  if (season == 4) {season.str = "_season";s.title="APR"}
  if (season == 5) {season.str = "_season";s.title="MAY"}
  if (season == 6) {season.str = "_season";s.title="JUN"}
  if (season == 7) {season.str = "_season";s.title="JUL"}
  if (season == 8) {season.str = "_season";s.title="AUG"}
  if (season == 9) {season.str = "_season";s.title="SEP"}
  if (season ==10) {season.str = "_season";s.title="OCT"}
  if (season ==11) {season.str = "_season";s.title="NOV"}
  if (season ==12) {season.str = "_season";s.title="DEC"}
  
  
  ## Reads previously detected and validated breaks
  ## in the xxxxxxdetected.txt file
  flag.prev   = F
  file.detect = paste(net.str,"detected.txt",sep="")
  if (file.exists(file.detect)){prev.frame=read.table(file.detect,head=F,
                                                      colClasses=c("character","character","numeric","numeric","character","character"),sep="\t");flag.prev=T}
  
  ## Loop over files
  ##################
  for (j in 1:n.file){
    
    file = list.file[j]
    if (flag.prev) {prev.brk = subset(prev.frame,V1==station.index[j])$V3}
    
    ## Blah blah
    ############
    write("                      ",file="")
    write(paste("Detection: ",station.index[j],station.name[j],s.title),file="")
    write("===========================================================",file="")
    
    ## Reading files
    ## Detection is performed on ANNUAL or SEASONAL or Monthly values
    ## AVERAGES are used for additive models
    ## SUMS otherwise
    ######################################################
    
    data.tmp     = read.serie(file,season)
    y.candidate  = data.tmp[,2]
    t.candidate  = data.tmp[,1]
    t.lim        = c(t.candidate[1],t.candidate[nrow(data.tmp)])
    
    ## Loop over neighboors
    ####################################################
    
    list.neighbour = list.file[-j]
    if (inter.option == "g") {list.neighbour=distance.serie(list.file,j)}
    if (inter.option == "c") {list.neighbour=correlation.distance(list.file,j)}
    
    for (file.ref in list.neighbour){
      brk      = 0
      n.brk    = 0
      i.ref    = match(file.ref,list.file)
      if ((flag.done[i.ref,j]==FALSE) & (flag.done[j,i.ref]==FALSE)){
        write(paste("vs",file.ref),file="")
        
        data.tmp = read.serie(file.ref,season)
        y.ref    = data.tmp[,2]
        t.ref    = data.tmp[,1]
        
        y.comp   = comparison.serie(y.candidate,t.candidate,y.ref,t.ref,FALSE)
        y.model  = numeric(nrow(y.comp))
        
        ## Computation is performed if comparison length >=15
        ## And series are not identical which may occur when
        ## using series form fucking xXXXXxxbiiipp website
        #####################################################
        if (nrow(y.comp) >= 15 ){
          if (sum(abs(y.comp[,2]))>0)
          {#brk = critpen.lik(y.comp[,2],kmax=min(100,nrow(y.comp)/3),crit=crit)
            brk  = uni.seg(y.comp[,2])
            if (brk[1] !=0) {n.brk=length(brk);brk=c(0,brk)}
            brk = c(brk,nrow(y.comp))             
            seg.means  = numeric(n.brk+1)
            for (i in 1:(n.brk+1))
            {seg.means[i]=mean(y.comp[(brk[i]+1):brk[i+1],2])
            y.model[(brk[i]+1):brk[i+1]]=seg.means[i]}
            
            sd.mat[j,i.ref]  = sqrt(sum((y.comp[,2]-y.model)^2)/(nrow(y.comp)-n.brk))
            sd.mat[i.ref,j]  = sd.mat[j,i.ref]
            brk.arr[j,i.ref,1:(n.brk+2)]=brk
            brk.arr[i.ref,j,1:(n.brk+2)]=brk
            n.brk.mat[j,i.ref]=n.brk
            n.brk.mat[i.ref,j]=n.brk
            flag.done[j,i.ref]=TRUE
            flag.done[i.ref,j]=TRUE
          } else {write(paste("Warning, same data in",file,"and",file.ref),file="")}               
        } else {write("Common period too short",file="")}
      }
    }
    
    #####################################################
    ## Graphic outputs, once detection has been performed
    #####################################################
    pos.brk     = numeric()
    amp.brk     = numeric()
    n.neighbour = length(list.neighbour)
    gt.lim      = c(10*floor(t.lim[1]/10),10*round((t.lim[2]+5)/10))
    v.ticks     = seq(gt.lim[1],gt.lim[2],5)
    
    file=paste("fig/detect_",head.str,station.index[j],"_a",s.title,sep="")
    fig.file(file,he=2*n.neighbour,wi=12)
    par(mfrow=c(n.neighbour,1),mar=c(3,4,2,1))
    
    ## Sorting comparisons by increasing standard deviation of residuals
    ## Note, for previous PRODIGE users, now SMALLER SD ARE ON TOP
    ####################################################################
    i.order = sort(sd.mat[j,],index=T)$ix[(sum(is.na(n.brk.mat[j,]))+1):n.file]
    sd.comp  = round(sort(sd.mat[j,],index=T)$x[(sum(is.na(n.brk.mat[j,]))+1):n.file],2)
    
    k       = 0
    for (file.ref in list.file[i.order]){
      k=k+1
      sd.str=sd.comp[k]
      if (comp.option == "a"){
        if (unit.str=="c") {text.label=substitute(paste(sigma,"=",sd.str,degree,"C"))
        } else {text.label=substitute(paste(sigma,"=",sd.str))}
      } else {text.label=substitute(paste(sigma,"=",sd.str))}
      
      data.tmp  = read.serie(file.ref,season)
      y.ref     = data.tmp[,2]
      t.ref     = data.tmp[,1]
      y.comp    = comparison.serie(y.candidate,t.candidate,y.ref,t.ref,FALSE)
      n.brk     = n.brk.mat[j,i.order[k]]
      brk       = brk.arr[j,i.order[k],1:(n.brk+2)]
      seg.means = numeric(n.brk+1)
      y.model   = numeric(nrow(y.comp))
      
      ### Preparing vectors for synthesis a
      for (i in 1:(n.brk+1))
      {seg.means[i]=mean(y.comp[(brk[i]+1):brk[i+1],2])
      y.model[(brk[i]+1):brk[i+1]]=seg.means[i]}
      if (n.brk>0){pos.brk   = c(pos.brk,y.comp[brk[2:(n.brk+1)],1])
      amp.brk   = c(amp.brk,diff(seg.means))}
      
      ### Plotting differences
      plot(y.comp[,1],y.comp[,2],xlab="",ylab=y.label,pch="+",
           xlim=gt.lim,xaxs="i",lab=c(10,5,7),
           main=paste(name.str[j],"vs",station.index[i.order[k]],
                      station.name[i.order[k]],status.str,s.title))
      
      y.pos.text=max(y.comp[,2])-min(y.comp[,2])
      y.pos.text=min(y.comp[,2])+y.pos.text*.95
      #improvement:				
      if (flag.prev) {abline(v=prev.brk,lwd=6,col="lightcyan2")}  #abline(v=prev.brk,lwd=5,col="lightgrey")
      text(gt.lim[1]+.5,y.pos.text,labels=text.label,pos=4,cex=1.5)
      points(y.comp[,1],y.comp[,2],pch="+")
      abline(v=y.comp[brk[1:(length(brk)-1)],1],lty=1,lwd=2)
      abline(v=v.ticks,lty=3,lwd=.75)
      #               lines(y.comp[,1]-1,y.model,type="s",col="red")
      for (i in 1:(n.brk+1)){
        lines(c(y.comp[brk[i]+1,1],y.comp[brk[i+1],1]),rep(seg.means[i],2),
              col="blue",lty=1,lwd=1.5)}
    }
    dev.off()
    
    if (flag.inter == T & season == 13){
      x11(wi=10,he=10,title=name.str[j])
      n.max.plot = min(6,n.neighbour)
      par(mfrow=c(n.max.plot,1),mar=c(3,4,2,1))
      
      ## Sorting comparisons by increasing standard deviation of residuals
      ## Note, for previous PRODIGE users, now SMALLER SD ARE ON TOP
      ####################################################################
      i.order = sort(sd.mat[j,],index=T)$ix[(sum(is.na(n.brk.mat[j,]))+1):n.file]
      sd.comp  = round(sort(sd.mat[j,],index=T)$x[(sum(is.na(n.brk.mat[j,]))+1):n.file],2)
      k       = 0
      for (file.ref in list.file[i.order][1:n.max.plot]){
        k=k+1
        sd.str=sd.comp[k]
        if (comp.option == "a"){
          if (unit.str=="c") {text.label=substitute(paste(sigma,"=",sd.str,degree,"C"))
          } else {text.label=substitute(paste(sigma,"=",sd.str))}
        } else {text.label=substitute(paste(sigma,"=",sd.str))}
        
        data.tmp  = read.serie(file.ref,season)
        y.ref     = data.tmp[,2]
        t.ref     = data.tmp[,1]
        y.comp    = comparison.serie(y.candidate,t.candidate,y.ref,t.ref,FALSE)
        n.brk     = n.brk.mat[j,i.order[k]]
        brk       = brk.arr[j,i.order[k],1:(n.brk+2)]
        seg.means = numeric(n.brk+1)
        y.model   = numeric(nrow(y.comp))
        
        ### Preparing vectors for synthesis a
        for (i in 1:(n.brk+1))
        {seg.means[i]=mean(y.comp[(brk[i]+1):brk[i+1],2])
        y.model[(brk[i]+1):brk[i+1]]=seg.means[i]}
        if (n.brk>0){pos.brk   = c(pos.brk,y.comp[brk[2:(n.brk+1)],1])
        amp.brk   = c(amp.brk,diff(seg.means))}
        
        ### Plotting differences
        plot(y.comp[,1],y.comp[,2],xlab="",ylab=y.label,pch="+",
             xlim=gt.lim,xaxs="i",lab=c(10,5,7),
             main=paste(name.str[j],"vs",station.index[i.order[k]],
                        station.name[i.order[k]],status.str,s.title))
        
        y.pos.text=max(y.comp[,2])-min(y.comp[,2])
        y.pos.text=min(y.comp[,2])+y.pos.text*.95
        #improvement:				
        if (flag.prev) {abline(v=prev.brk,lwd=6,col="lightcyan2")}  #abline(v=prev.brk,lwd=5,col="lightgrey")
        text(gt.lim[1]+.5,y.pos.text,labels=text.label,pos=4,cex=1.5)
        points(y.comp[,1],y.comp[,2],pch="+")
        abline(v=y.comp[brk[1:(length(brk)-1)],1],lty=1,lwd=2)
        abline(v=v.ticks,lty=3,lwd=.75)
        #               lines(y.comp[,1]-1,y.model,type="s",col="red")
        for (i in 1:(n.brk+1)){
          lines(c(y.comp[brk[i]+1,1],y.comp[brk[i+1],1]),rep(seg.means[i],2),
                col="blue",lty=1,lwd=1.5)}
      }
      write("Click on left button to see next detection",file="")
      xy.dummy=locator(1)
      dev.off()
    }
    
    ## Synthesis of detected breaks (author : Olivier, original idea : Alexis HANNART)
    
    file=paste("fig/detect_",head.str,station.index[j],"_b",s.title,sep="")
    fig.file(file,he=6,wi=8)
    if (length(amp.brk)>0){range.amp=range(amp.brk)} else {range.amp=c(-1,1)}
    plot(gt.lim[1]:gt.lim[2],rep(0,diff(range(gt.lim))+1),
         xlab="",ylab=y.label,pch=".",
         xlim=gt.lim,xaxs="i",ylim=range.amp,lab=c(10,5,7),
         main=paste(name.str[j],status.str,s.title))
    #improvement: 
    if (flag.prev) {abline(v=prev.brk,lwd=5,col="lightcyan2")}  #abline(v=prev.brk,lwd=3,col="lightgrey")
    
    if (length(amp.brk)>0){
      col.triangle = rep("white",length(amp.brk))
      if (season == 16) {col.triangle = rep("blue",length(amp.brk))}
      if (season == 17) {col.triangle = rep("green",length(amp.brk))}
      if (season == 18) {col.triangle = rep("red",length(amp.brk))}
      if (season == 19) {col.triangle = rep("orange",length(amp.brk))}
      #improvement: outline colour
      col.outline = rep("black",length(amp.brk))
      if (season == 16) {col.outline = rep("blue",length(amp.brk))}
      if (season == 17) {col.outline = rep("green",length(amp.brk))}
      if (season == 18) {col.outline = rep("red",length(amp.brk))}
      if (season == 19) {col.outline = rep("orange",length(amp.brk))}
      sym.triangle = rep(25,length(amp.brk))
      
      #improvement: 
      #original          points(pos.brk,amp.brk,pch=sym.triangle,bg=col.triangle,cex=2)}
      # sd.comp is sorted sd.mat[j,],  the size according to STD IS NOT WORKING WELL
      points(pos.brk,amp.brk,pch=sym.triangle, col=col.outline,lwd=2,bg=ifelse(amp.brk>=0,"yellow2","skyblue2"), cex=2*rep(sd.comp[1],length(sd.comp))/sd.mat[j,])}
    abline(v=v.ticks,lty=3,lwd=.75)
    dev.off()
    
    j.index=j
    if (season == 13){
      save(j.index,range.amp,gt.lim,name.str,status.str,s.title,pos.brk,amp.brk,v.ticks,y.label,
           file=paste("tmp/detect_",head.str,station.index[j],".sav",sep=""))}
    if (season == 16){
      range.amp.djf = range.amp
      gt.lim.djf    = gt.lim
      pos.brk.djf   = pos.brk
      amp.brk.djf   = amp.brk
      save(range.amp.djf,gt.lim.djf,pos.brk.djf,amp.brk.djf,
           file=paste("tmp/detect_",head.str,station.index[j],"_DJF.sav",sep=""))}
    if (season == 18){
      range.amp.jja = range.amp
      gt.lim.jja    = gt.lim
      pos.brk.jja   = pos.brk
      amp.brk.jja   = amp.brk
      save(range.amp.jja,gt.lim.jja,pos.brk.jja,amp.brk.jja,
           file=paste("tmp/detect_",head.str,station.index[j],"_JJA.sav",sep=""))}
    
    
    
    
    
    
    ## Synthesis of detected breaks (author : Brigitte Dubuisson)   
    file=paste("fig/detect_",head.str,station.index[j],"_c",s.title,sep="")
    fig.file(file,he=6,wi=8)
    k = 0
    for (file.ref in list.file[i.order]){
      k=k+1
      data.tmp = read.serie(file.ref,season)
      y.ref    = data.tmp[,2]
      t.ref    = data.tmp[,1]
      y.comp   = comparison.serie(y.candidate,t.candidate,y.ref,t.ref,FALSE)
      n.brk    = n.brk.mat[j,i.order[k]]
      brk      = brk.arr[j,i.order[k],1:(n.brk+2)]
      
      y.vect=k
      if(n.brk>1) 
      {    { for (ibrk in 2:n.brk)  y.vect=c(y.vect,k) }
      }
      
      if(k==1) 
      { if(n.brk>0)
      { 
        #improvement: adding bars: needs to plot first, then to put bars, and only then the points
        plot(y.comp[brk[2:(n.brk+1)],1],y.vect,pch=".", xlim=gt.lim,ylim=c(0,length(i.order)+1),xlab="",
             ylab="",yaxt="n",xaxs="i",yaxs="i",main=paste(name.str[j],status.str,s.title) )
        #improvement: adding bars: 
        if (flag.prev) {abline(v=prev.brk,lwd=5,col="lightcyan2")}  #abline(v=prev.brk,lwd=5,col="lightgrey")
        
        #improvement:  size of triangles, colour according to season and amplitude - positive or negative
        #PROBLEM:  AMPLITUDE - DOES NOT WORK -   amp.brk is not SORTED
        points(y.comp[brk[2:(n.brk+1)],1],y.vect,pch=25, col=col.outline,lwd=2,bg=ifelse(amp.brk>=0,"yellow2","skyblue2"),  cex=2*sd.comp[1]/sd.comp[k])
      } else { 
        plot(y.comp[brk[2:(n.brk+1)],1],k,pch=".",xlim=gt.lim,ylim=c(0,length(i.order)+1),xlab="",
             ylab="",yaxt="n",xaxs="i",yaxs="i",main=paste(name.str[j],status.str,s.title)) 
      }
        axis(side=2,at=1:length(i.order),labels=station.index[i.order],las=1,hadj=0.8,cex.axis=0.8)
        axis(side=4,at=1:length(i.order),labels=sd.comp,las=1,hadj=0.3,cex.axis=0.8)
        abline(v=v.ticks,lty=3,lwd=.75)
      }
      if(k>1)
        #improvement: 
        #PROBLEM:  AMPLITUDE - DOES NOT WORK -   amp.brk is not SORTED
        #               { if(n.brk>0) {points(y.comp[brk[2:(n.brk+1)],1],y.vect,pch=25,bg="black")}
      { if(n.brk>0) {points(y.comp[brk[2:(n.brk+1)],1],y.vect,pch=25, col=col.outline,lwd=2,bg=ifelse(amp.brk>=0,"yellow2","skyblue2"), cex=2*sd.comp[1]/sd.comp[k])}
        else {points(y.comp[brk[2:(n.brk+1)],1],k,pch=".")}
      }
      
      abline(h=k)  
      
    }
    dev.off()
  }}
##

######################################################
station.factor = function(year,brk.frame,index,month){
  ######################################################
  
  ###############################
  #### Funtion called by correc #
  #### Computes station factor  #
  #### in correction model      #
  ###############################
  
  brk.vector   = c(0,subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "BREAK")$V3-year[1]+1,year[length(year)]-year[1]+1)
  
  ### Taking into account month of change for MONTHLY series 
  if (month < 13) {
    month.vector = c(12,subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "BREAK")$V4,12)
    brk.vector[month>month.vector]=brk.vector[month>month.vector]-1
  }
  
  nu = character(length(year))
  for (i in 1:(length(brk.vector)-1)){
    nu[(brk.vector[i]+1):(brk.vector[i+1])] = paste(index,year[brk.vector[i]+1],"-",year[brk.vector[i+1]])
  }
  #write.table(data.frame(index,month,nu=nu[1:6]),file="diag.txt",append=T,quote=F,row=F,col=F)
  return(nu)}

################################################
station.trend = function(year,brk.frame,index){
  ################################################
  
  ###############################
  #### Funtion called by correc #
  #### Computes station trend   #
  #### in correction model      #
  ###############################
  
  urb.trend   = rep(0,length(year))
  
  if (nrow(subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "BEGTR"))   > 0 &
      nrow(subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "ENDTR"))   > 0){
    
    start.trend = max(subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "BEGTR")$V3-year[1]+1,1)
    end.trend   = min(subset(brk.frame,brk.frame$V1==index & brk.frame$V2 == "ENDTR")$V3-year[1]+1,length(year))
    urb.trend[start.trend:end.trend] = 1:(end.trend-start.trend+1)
  }
  return(urb.trend)}
##################

##############################################################
##############################################################
##############################################################
correc = function(list.file,n.file,head.str,list.file.ref){
  ##############################################################
  ##############################################################
  ##############################################################
  
  ########################################################
  ####            MAIN PROCEDURE FOR CORRECTION          #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  
  flag.tr = 0
  corr.cycle <<- corr.cycle+1
  file.copy(paste(net.str,"detected.txt",sep=""),paste("tmp/",net.str,"detected.txt.",as.integer(corr.cycle),sep=""))
  if (season.corr =="m"){list.month=1:13} else {list.month=c(13,1:12)}
  
  brk.frame = read.table(paste(net.str,"detected.txt",sep=""),
                         colClasses=c("character","character","numeric","numeric","character","character"),quote="",sep="\t")
  
  ### Sorts change-point dates, in case not in chronological order
  brk.frame    = brk.frame[order(brk.frame[,1],
                                 brk.frame[,3],
                                 as.numeric(brk.frame[,4])),]
  summary.def = data.frame()
  name.coeff  = character(0)
  
  ## Loop over files
  ##################
  for (j in 1:n.file) {
    
    coef.corr   = numeric(0)
    coeff_m     = numeric(0)
    file = list.file[j]
    write("                      ",file="")
    write(paste("Correction of: ",station.index[j],station.name[j]),file="")
    write("==============================================",file="")
    str.tmp     = rep(station.index[j],13)
    name.coeff  = c(name.coeff,paste(str.tmp,".",name.month,sep=""))
    file.corr   = paste("ho/ho",substr(head.str,3,4),"m",station.index[j],"d.txt",sep="")
    file.corr.y = paste("ho/ho",substr(head.str,3,4),"y",station.index[j],"d.txt",sep="")
    f.reference = list.file.ref[j,1]
    f.ref.ho    = list.file.ref[j,2]
    file.meta   = paste("meta/meta",substr(head.str,3,4),station.index[j],"d.txt",sep="")
    write(paste(station.index[j],station.name[j]),file=file.meta)
    write("=====================================",file=file.meta,append=TRUE)
    
    ### Neighbour list constitution      
    list.neighbour = list.file[-j]
    if (inter.option == "g") {list.neighbour=distance.serie(list.file,j)}
    if (inter.option == "c") {list.neighbour=correlation.distance(list.file,j)}
    write("CORRECTION NEIGHBORHOOD",file=file.meta,append=TRUE)
    file.append(file.meta,"tmp/file_neighbours.txt")
    write("=====================================",file=file.meta,append=TRUE)
    if(file.exists("tmp/file_neighbours.txt")) { file.remove("tmp/file_neighbours.txt") }
    
    for (month in list.month) {
      y               = numeric()
      mu              = numeric()
      u.trend         = numeric()
      x.lon           = numeric()
      y.lat           = numeric()
      nu              = character()
      group.station   = character()
      
      data.tmp        = read.serie(file,month)
      y.candidate     = data.tmp[,2]
      t.candidate     = data.tmp[,1]
      
      t.lim           = c(t.candidate[1],t.candidate[nrow(data.tmp)])
      year            = t.lim[1]:t.lim[2]
      missing.values  = which(y.candidate == miss.flag)
      no.miss.values  = which(y.candidate != miss.flag)
      x.lon.candidate = rep(x.sta[j],length(no.miss.values))
      y.lat.candidate = rep(y.sta[j],length(no.miss.values))
      group.candidate = rep(station.index[j],length(no.miss.values))
      
      if (comp.option == "a") {
        mult.coeff = 1
      } else {mult.coeff = mean(y.candidate[no.miss.values])
      if (mult.coeff != 0){
        y.candidate[no.miss.values]=y.candidate[no.miss.values]/mult.coeff
      } else {y.candidate[no.miss.values]=0}  ### in case of zero sums on whole series 
      }
      
      if (month == 1)  {data.corr    = t.candidate
      data.ref     = t.candidate}
      if (month == 13) {data.corr.y  = t.candidate}
      
      ## nu : station factor
      ######################
      if (comp.option == "a" | (comp.option != "a" & month == 13)) {
        
        flag.bad.frame = T
        ### Checks inconsistencies in break frame: common error is to have only missing values between two change-points
        while (flag.bad.frame){
          nu.bak              = station.factor(year,brk.frame,station.index[j],month)
          n.period.candidate  = nlevels(factor(nu.bak))
          
          coeff=numeric(n.period.candidate)
          if (nlevels(factor(nu.bak[y.candidate != miss.flag])) != n.period.candidate) {
            write(paste(" Warning! No data for series",station.index[j],"during period",
                        setdiff(levels(factor(nu.bak)),levels(factor(nu.bak[y.candidate != miss.flag]))),
                        "and month",month),file="")
            grouik = readline(" Check inconsistencies in change-point dates (RETURN) or quit HOMER (q): ")
            if (grouik == "") {
              write("",file="")
              break.frame.index = which(brk.frame$V1==station.index[j] & brk.frame$V2 == "BREAK")
              break.frame.tmp   = brk.frame[break.frame.index,]
              break.frame       = brk.frame[-break.frame.index,]
              i.remove          = numeric()
              for (brk.i in 1:nrow(break.frame.tmp)){
                str.tmp = paste(break.frame.tmp[brk.i,c(1,3,4,6)],collapse=" ")
                str.tmp = paste(" keep:",str.tmp,"     y (RETURN) or no (n) :")
                keep.i=readline(str.tmp)
                if (keep.i == "n") {i.remove=c(i.remove,brk.i)}
              }
              if (length(i.remove)>0){break.frame.tmp=break.frame.tmp[-i.remove,]}
              brk.frame    = rbind(break.frame,break.frame.tmp)
              brk.frame    = brk.frame[order(brk.frame[,1],
                                             brk.frame[,3],
                                             as.numeric(brk.frame[,4])),]
              
            } else {options(show.error.messages=F)
              write(" See you soon!",file="")
              stop()}
          } else {flag.bad.frame = F}
        }}
      
      ## u.trend : urban trend
      ########################
      if (comp.option == "a") {u.trend.bak = station.trend(year,brk.frame,station.index[j])}
      
      ## Loop over neighboors
      #######################
      for (file.ref in list.neighbour){
        data.tmp  = read.serie(file.ref,month)
        i.ref     = match(file.ref,list.file)
        y.ref     = data.tmp[,2]
        t.ref     = data.tmp[,1]
        nu.ref    = station.factor(t.ref,brk.frame,station.index[i.ref],month)
        n.control = nlevels(nu.ref)
        
        ### Removing Reference missing values
        i.m.ref   = y.ref  != miss.flag
        t.ref     = t.ref[i.m.ref]
        nu.ref    = nu.ref[i.m.ref]
        y.ref     = y.ref[i.m.ref]
        
        if (comp.option != "a") {
          ref.mult.coeff=mean(y.ref)
          if (ref.mult.coeff != 0){
            y.ref=y.ref/ref.mult.coeff
          } else {y.ref[]=0}  ### in case of zero sums on whole series 
        }
        ### Keeping common values with Candidate
        j.m.ref       = !is.na(match(t.ref,year))
        mu            = c(mu,t.ref[j.m.ref])
        y             = c( y,y.ref[j.m.ref])
        nu.ref        = nu.ref[j.m.ref]
        x.lon         = c(x.lon,rep(x.sta[i.ref],length(nu.ref)))
        y.lat         = c(y.lat,rep(y.sta[i.ref],length(nu.ref)))
        group.station = c(group.station,rep(station.index[i.ref],length(nu.ref)))
        
        if (comp.option == "a") {
          u.trend.ref = station.trend(year,brk.frame,station.index[i.ref])
          u.trend.ref = u.trend.ref[!is.na(match(year,t.ref[j.m.ref]))]
          if (max(u.trend.ref)>0) {
            i.t.max=which.max(u.trend.ref)
            u.trend.ref[1:i.t.max]=u.trend.ref[1:i.t.max]-max(u.trend.ref)
          }
          u.trend     = c(u.trend,u.trend.ref)}
        
        if (nlevels(nu.ref) != n.control){
          print(paste("CAUTION!! No data for series",station.index[i.ref],"during period",setdiff(nu.ref,nu)),quote=F)
          options(show.error.messages=F)
          write("Check inconsistencies in change-point dates",file="")
          stop()
        } else {nu = c(nu,nu.ref)} 
      }
      i.m.candidate = y.candidate != miss.flag
      y             = c(y,y.candidate[i.m.candidate])
      mu            = c(mu,t.candidate[i.m.candidate])
      nu            = c(nu,nu.bak[i.m.candidate])
      mu            = as.character(mu)
      n.coeff.nu    = nlevels(factor(nu))
      n.coeff.mu    = nlevels(factor(mu))
      x.lon         = c(x.lon,x.lon.candidate)
      y.lat         = c(y.lat,y.lat.candidate)
      group.station = c(group.station,group.candidate)
      
      if (comp.option != "a") {
        data.all   = data.frame(dy=y,fmu=factor(mu),fnu=factor(nu),
                                x.lon=x.lon,y.lat=y.lat,g=factor(group.station))
      } else {
        if (max(u.trend.bak)>0) {
          i.t.max=which.max(u.trend.bak)
          u.trend.bak[1:i.t.max]=u.trend.bak[1:i.t.max]-max(u.trend.bak)
        }
        u.trend.candidate = u.trend.bak[i.m.candidate]
        u.trend    = c(u.trend,u.trend.candidate)
        data.all   = data.frame(dy=y,fmu=factor(mu),fnu=factor(nu),
                                u.trend=u.trend,s=factor(substr(nu,1,8)),
                                x.lon=x.lon,y.lat=y.lat,g=factor(group.station))
        write.table(format(data.all),file="tmp/data.txt",row=F)}
      
      if (comp.option !="a") {
        if (mult.coeff !=0) {
          lm.out     = lm(dy~fnu+fmu-1,data=data.all)}
      } else {lm.out = lm(dy~fnu+fmu+u.trend:s-1,data=data.all)
      save(lm.out,file="tmp/lm.sav")
      work.out=cbind(data.all[,c(2,6,7,8)],residuals(lm.out))
      names(work.out)=c("fmu","x.lon","y.lat","g","e")
      #                      covmat=correlogram.serie(work.out,month)
      
      #gls.out = gls(dy~fnu+fmu+u.trend:s-1,data=data.all,correlation=corExp(form =~x.lon+y.lat|g))
      #print(summary(gls.out))
      miss.coeff=which(is.na(coefficients(lm.out)))
      lm.out$coefficients[miss.coeff]=0
      #                      print(paste(station.name[j],month))
      #                      print("==========================")
      #                      print("==========================")
      #                      print(summary(lm.out))
      summary.frame=summary(lm.out)$coefficients
      summary.frame=as.data.frame(summary.frame)
      summary.selec=which(row.names(summary.frame)==paste("u.trend:s",station.index[j],sep=""))
      if (length(summary.selec)>0){
        summary.frame=round(summary.frame[summary.selec,],5)
        summary.frame=cbind(summary.frame,name.month[month],station.index[j])
        summary.def  = rbind(summary.def,summary.frame)}
      }
      lm.coeff.i = paste("fnu",as.character(levels(factor(nu.bak))),sep="")
      lm.coeff   = coefficients(lm.out)[!is.na(match(names(coefficients(lm.out)),lm.coeff.i))]
      y.corr     = y.candidate
      
      mu.ref.i   = paste("fmu",as.character(levels(factor(mu))),sep="")
      mu.ref     = c(0,coefficients(lm.out)[!is.na(match(names(coefficients(lm.out)),mu.ref.i))])
      if (month < 13) {data.ref   = cbind(data.ref,mu.ref)}
      
      if (month == 13) {lm.coeff.def=lm.coeff}
      
      for (k in 1:n.period.candidate){    
        select.period=nu.bak==as.character(levels(factor(nu.bak))[k]) & y.corr != miss.flag
        if (season.corr == "m"){
          correction=lm.coeff[n.period.candidate]-lm.coeff[k]
          coeff[k]=correction
          
          if (comp.option == "a") {
            y.corr[select.period]=(y.corr[select.period]+correction)*mult.coeff
          } else {y.corr[select.period]=(y.corr[select.period])*(1+correction)*mult.coeff}
        } else {
          correction=lm.coeff.def[n.period.candidate]-lm.coeff.def[k]
          y.corr[select.period]=(y.corr[select.period])*(1+correction)*mult.coeff
          if(month == 13) { 
            write(sprintf("\nPeriode %d",k),file=file.meta,append=TRUE) 
            write(levels(factor(nu.bak))[k],file=file.meta,append=TRUE) 
            if(k<n.period.candidate) {   
              write(paste(month,round(1+correction,3),round(mult.coeff,3)),file=file.meta,append=TRUE)
              write(paste("Amplitude : ",round(lm.coeff.def[k+1]-lm.coeff.def[k],2)),file=file.meta,append=TRUE)
            } else {write("UNCHANGED REFERENCE PERIOD",file=file.meta,append=TRUE)}
          }}
      }
      if (month==list.month[1]){coeff_m=coeff} else {coeff_m=rbind(coeff_m,coeff)}
      
      ### Filling missing values
      ##########################
      ##########################
      
      if (month <= 12) {
        if (length(missing.values) >=1){
          missing.values=missing.values[!is.na(match(as.character(year[missing.values]),mu))]
          if (comp.option =="a") {
            data.miss = data.frame(fmu=factor(t.candidate[missing.values]),
                                   fnu=factor(rep(as.character(levels(factor(nu.bak))[n.period.candidate]),
                                                  length(missing.values))),
                                   s=factor(rep(station.index[j],length(missing.values))),
                                   u.trend=u.trend.bak[missing.values])
            y.corr[missing.values]=predict(lm.out,newdata=data.miss)*mult.coeff
          } else {
            data.miss = data.frame(fmu=factor(t.candidate[missing.values]),
                                   fnu=factor(rep(as.character(levels(factor(nu.bak))[n.period.candidate]),
                                                  length(missing.values))))
            if (mult.coeff !=0){
              y.corr[missing.values]=predict(lm.out,newdata=data.miss)*mult.coeff
              y.corr[y.corr<0]=0} else {y.corr[missing.values]=0}}
        }}
      
      if (month <= 12) {
        
        ### Removing trends
        ###################
        ### Trends are removed on the MONTHLY data, not on annual (special hoxxy files)
        ### so that comparisons can be made with/without urban trends
        
        if (comp.option =="a") {
          data.all=subset(data.all,as.character(data.all$s)==station.index[j])
          index.trend = which(diff(data.all$u.trend)>0)+1
          if (length(index.trend)>0){
            flag.tr=1
            if (month == 1) {year.trend  = data.corr} else {year.trend = data.corr[,1]}
            
            start.trend = year.trend[data.all$fmu[index.trend[1]]]-year.trend[1]+1
            end.trend   = year.trend[data.all$fmu[index.trend[length(index.trend)]]]-year.trend[1]+1
            t.trend     = rep(0,length(year.trend))
            n.trend     = end.trend-start.trend+1
            t.trend[1:start.trend]=-n.trend
            t.trend[start.trend:end.trend]=-((n.trend-1):0)
            b.trend=lm.out$coefficients[paste("u.trend:s",station.index[j],sep="")]
            y.corr=-b.trend*t.trend+y.corr
          }} ### End removing trends
        
        if (comp.option != "a") {y.corr[y.corr<0]=0}
        data.corr=cbind(data.corr,y.corr)
      } else {
        if (comp.option == "a") {
          ### Annual correction without removing urban trends
          lm.coeff.yy   = coefficients(lm.out)
          name.coef.str  = names(lm.coeff.yy)
          name.coef.str  = c(name.coef.str,paste("fmu",t.lim[1],sep=""))
          lm.coeff.yy    = c(lm.coeff.yy,0)
          names(lm.coeff.yy)=name.coef.str
          climate.signal=lm.coeff.yy[paste("fmu",data.corr.y,sep="")]
          mean.mu=mean(climate.signal)
          climate.signal=climate.signal-mean.mu
          
          if (length(missing.values) >=1){
            for (k in 1:length(missing.values)){
              y.corr[missing.values[k]] = (lm.coeff.yy[n.period.candidate]+
                                             lm.coeff.yy[paste("fmu",t.candidate[missing.values[k]],sep="")]+
                                             lm.coeff.yy[paste("u.trend:s",station.index[j],sep="")]*
                                             u.trend.bak[missing.values[k]]) * mult.coeff
            }}
          data.corr.y=cbind(data.corr.y,y.corr,climate.signal)
        } else {data.corr.y=cbind(data.corr.y,y.corr)}
      }
    } ### End loop months
    
    ### Smoothing monthly coefficients
    ### Hanning 3 points filter for monthly coefficients ### unelegant code
    if (season.corr == "m"){
      coeff.tmp = coeff_m
      for (k in 1:n.period.candidate){
        for (month in 1:12){
          coeff_m[1,k]=(coeff.tmp[1,k]+coeff.tmp[12,k]+coeff.tmp[2,k])/3
          coeff_m[2,k]=(coeff.tmp[2,k]+coeff.tmp[1,k]+coeff.tmp[3,k])/3
          coeff_m[3,k]=(coeff.tmp[3,k]+coeff.tmp[2,k]+coeff.tmp[4,k])/3
          coeff_m[4,k]=(coeff.tmp[4,k]+coeff.tmp[3,k]+coeff.tmp[5,k])/3
          coeff_m[5,k]=(coeff.tmp[5,k]+coeff.tmp[4,k]+coeff.tmp[6,k])/3
          coeff_m[6,k]=(coeff.tmp[6,k]+coeff.tmp[5,k]+coeff.tmp[7,k])/3
          coeff_m[7,k]=(coeff.tmp[7,k]+coeff.tmp[6,k]+coeff.tmp[8,k])/3
          coeff_m[8,k]=(coeff.tmp[8,k]+coeff.tmp[7,k]+coeff.tmp[9,k])/3
          coeff_m[9,k]=(coeff.tmp[9,k]+coeff.tmp[8,k]+coeff.tmp[10,k])/3
          coeff_m[10,k]=(coeff.tmp[10,k]+coeff.tmp[9,k]+coeff.tmp[11,k])/3
          coeff_m[11,k]=(coeff.tmp[11,k]+coeff.tmp[10,k]+coeff.tmp[12,k])/3
          coeff_m[12,k]=(coeff.tmp[12,k]+coeff.tmp[11,k]+coeff.tmp[1,k])/3
        }}
      
      for (month in 1:12){
        y.corr = as.vector(data.corr[,month+1])
        nu.bak = station.factor(year,brk.frame,station.index[j],month)
        n.period.candidate  = nlevels(factor(nu.bak))
        for (k in 1:n.period.candidate){
          select.period=nu.bak==as.character(levels(factor(nu.bak))[k]) & y.corr != miss.flag
          correction=coeff_m[month,k]-coeff.tmp[month,k]
          if (comp.option == "a") {
            y.corr[select.period]=(y.corr[select.period]+correction)*mult.coeff
          } else {y.corr[select.period]=(y.corr[select.period])*(1+correction)*mult.coeff}
        }
        data.corr[,month+1]=y.corr
      }
    }
    #########################
    #########################   
    
    data.corr=as.data.frame(data.corr)
    save(data.corr.y,file="tmp/out.sav")
    data.corr.y=as.data.frame(data.corr.y,row=F)
    write.table(format(round(data.corr,1)),file=file.corr,quote=F,row=F,col=F,sep="\t")
    write.table(format(round(data.ref,2)),file=f.reference,quote=F,row=F,col=F,sep="\t")
    write.table(format(round(data.ref,2)),file=f.ref.ho,quote=F,row=F,col=F,sep="\t")
    write.table(format(round(data.corr.y,2)),file=file.corr.y,quote=F,row=F,col=F,sep="\t")
    
    if(season.corr == "m"){
      for (k in 1:n.period.candidate){
        write(sprintf("\nPeriode %d",k),file=file.meta,append=TRUE) 
        write(levels(factor(nu.bak))[k],file=file.meta,append=TRUE) 
        for (month in list.month) {
          if(k<n.period.candidate) {
            write(sprintf("%2.2d : %6.2f",month,coeff_m[month,k]),file=file.meta,append=TRUE)} 
          if((k<n.period.candidate)&(month == 13)) {
            write(sprintf("BREAK AMPLITUDE : %6.2f",coeff_m[13,k]-coeff_m[13,k+1]),
                  file=file.meta,append=TRUE)
          }}
        if(k==n.period.candidate) {  write("UNCHANGED REFERENCE PERIOD",file=file.meta,append=TRUE) }
      }}
    
  } ### End loop over files
  if (flag.tr==1){
    names(summary.def)=c("  Estim.","st. dev"," t value","p.value"," mm"," index")
    write.table(format(summary.def),file=paste(net.str,"trends.txt",sep=""),quote=F,col=T,row=F)
  }
  brk.frame    = brk.frame[order(brk.frame[,1],
                                 brk.frame[,3],
                                 as.numeric(brk.frame[,4])),]
  write.table(brk.frame,file=paste(net.str,"detected.txt",sep=""),quote=F,row=F,col=F,sep="\t")
  return(summary.def)}
####################
####################


###############################################################
###############################################################
###############################################################
visu = function(list.file,n.file,head.str){
  ###############################################################
  ###############################################################
  ###############################################################
  
  ########################################################
  ####            MAIN PROCEDURE FOR VISUALISATION       #
  ####            CALLED FROM MAIN PROGRAM               #
  ########################################################
  
  if (season.option == "a") {season = 13 }
  if (season.option ==  "") {season = c(13,16:19)}
  if (season.option == "m") {season = 1:12}
  
  season.str = c(rep("",15),"DJF","MAM","JJA","SON")
  
  status.str = substr(head.str,1,1)
  if (status.str == "h"){status.str="(H)"} else {status.str = ""}
  
  ## Loop over seasons
  ####################
  for (i.season in season) {
    
    ## Loop over files
    ##################
    for (j in 1:n.file){
      
      file = list.file[j]
      
      ## Blah blah
      ############
      write("                      ",file="")
      write(file,file="")
      write("----------------------",file="")
      
      ## Reading files
      ## Annual AVERAGES are used for additive parameters
      ## Annual SUMS otherwise
      ###################################################
      
      data.tmp            = read.serie(file,i.season)
      y.candidate         = data.tmp[,2]
      t.candidate         = data.tmp[,1]
      i.miss              = which(y.candidate == miss.flag)
      y.candidate[i.miss] = NA  
      t.lim               = c(t.candidate[1],t.candidate[nrow(data.tmp)])
      
      ## Graphic output
      #################
      
      file=paste("fig/visu_",head.str,station.index[j],season.str[i.season],sep="")
      fig.file(file,he=6,wi=12)
      
      if (unit.str=="c") {y.label=expression(paste("TEMPERATURE (",degree,"C )"))
      } else {y.label=paste(par.str," (",unit.str,")",sep="")}
      
      
      plot(t.candidate,y.candidate,xlab="",
           ylab=y.label,type="l",xlim=t.lim,xaxs="i",main=paste(name.str[j],season.str[i.season],status.str),lwd=2)
      
      if (col.str=="y") {polyg.fill(y.candidate,t.candidate,i.miss)}
      
      if (trend.str  == "y") {abline(lm(y.candidate~t.candidate),lty=5,lwd=2)}
      write(paste(name.str[j],":",t.lim[1],"-",t.lim[2]),file="")
      trend.stat=summary(lm(y.candidate~t.candidate))
      write(paste("LS trend estimate ",season.str[i.season],":",round(trend.stat$coefficients[2],6)),file="")
      write(paste("Two sided Kendall test p-value : ",
                  round(cor.test(t.candidate,y.candidate+rnorm(length(y.candidate),0,0.00001),method="kendall")$p.value,5)),file="")
      if (smooth.str == "y") {lines(t.candidate,predict(loess(y.candidate~t.candidate),t.candidate),lty=3,lwd=2.5)}
      
      dev.off()
      
      if (flag.inter==T & i.season==13){
        x11(wi=10,he=4,title=paste(name.str[j]))
        plot(t.candidate,y.candidate,xlab="",
             ylab=y.label,type="l",xlim=t.lim,xaxs="i",
             main=paste(name.str[j],season.str[i.season],status.str),lwd=2)
        
        if (col.str=="y") {polyg.fill(y.candidate,t.candidate,i.miss)}
        if (trend.str  == "y") {abline(lm(y.candidate~t.candidate),lty=5,lwd=2)}
        if (smooth.str == "y") {lines(t.candidate,predict(loess(y.candidate~t.candidate),t.candidate),lty=3,lwd=2.5)}
        write("Click on left button to see next plot",file="")
        xy.dummy=locator(1)
        dev.off()
      }
      
    }}} ################### END FUNCTION "visu"


##################################################################
##################################################################
##################################################################
create.brk = function(file,n.file){
  ##################################################################
  ##################################################################
  ##################################################################
  
  ####################################################
  ####         FUNCTION FOR TYPING                   #
  ####          CHANGE-POINT DATES                   #
  ####################################################
  
  write("",file="")
  write("",file="")
  write("",file="")
  write("          BREAK : Abrupt shift",file="")
  write("          BEGTR : Beginning of a progressive shift (linear)",file="")
  write("          ENDTR : End of a progressive shift (linear)",file="")
  write("          ! hint: entering b when date is asked forces BEGTR to the start of the series",file="")
  write("          ! hint: entering e when date is asked forces ENDTR to the end of the series",file="")
  write("          OUTLI : Outlier (Natural : caused by local thunderstorm for example)",file="")
  write("                  Note that erroneous data have to be put as missing in files",file="")
  write("",file="")
  write("          New dates are taken until an empty date is typed (return)",file="")
  write("          Dates have to be typed from most ancient to most recent",file="")
  write("          Different shift types may have same date",file="")
  
  # IMPORTANT : For the moment, only Breaks are taken into account in correction
  
  write("",file="")
  write("",file="")
  write("",file="")
  
  index.sta      = character(0)
  name.sta       = character(0)
  meta.sta       = character(0)
  year.sta       = numeric(0)
  month.sta      = numeric(0)
  type.brk       = character(0)
  
  ## Loop over stations
  #####################
  flag.prev   = F
  file.detect = paste(net.str,"detected.txt",sep="")
  if (file.exists(file.detect)){prev.frame=read.table(file.detect,head=F,
                                                      colClasses=c("character","character","numeric","numeric","character","character"),sep="\t");flag.prev=T}
  
  for (j in 1:n.file){
    
    data.tmp        = read.serie(file[j],1)
    year.b          = data.tmp[1,1]
    year.e          = data.tmp[nrow(data.tmp),1]
    rm(data.tmp)
    
    write("",file="")
    write("",file="")
    write(paste("     ",station.index[j],station.name[j]),file="")
    write("      =========================================",file="")
    write("",file="")
    
    year.str  = "-1"
    year.prev = -1
    while (year.str !="") {
      # While year not empty, continue asking for more shifts
      #######################################################
      flag.trend = 0
      
      write("",file="")
      year.str = readline("          New date (return to quit)            : ")
      if (year.str != "")  {
        if (year.str == "b") {
          flag.trend = 1
          brk.str    = "BEGTR";month.str=12
          index.sta  = c(index.sta,station.index[j])     # station index
          name.sta   = c(name.sta,station.name[j])       # station name
          year.sta   = c(year.sta,year.b)                # year
          month.sta  = c(month.sta,as.numeric(month.str))
          type.brk   = c(type.brk,brk.str)
          meta.sta = c(meta.sta,"n")
        }
        if (year.str == "e") {
          flag.trend = 1
          brk.str    = "ENDTR";month.str=12
          index.sta  = c(index.sta,station.index[j])     # station index
          name.sta   = c(name.sta,station.name[j])       # station name
          year.sta   = c(year.sta,year.e)                # year
          month.sta = c(month.sta,as.numeric(month.str))
          type.brk = c(type.brk,brk.str)
          meta.sta = c(meta.sta,"n")
        }
        
        if (flag.trend==0) {
          year.prev = as.numeric(year.str)
          index.sta = c(index.sta,station.index[j])     # station index
          name.sta  = c(name.sta,station.name[j])       # station name
          year.sta  = c(year.sta,as.numeric(year.str))  # year
          
          month.str = readline("          Month (return for end of year)       : ")
          if (month.str == "") {month.str = 12}
          month.sta = c(month.sta,as.numeric(month.str))
          
          meta.str  = readline("          Metadata? (return for no, y for yes) : ")
          if (meta.str == "") {meta.str="n"} else {meta.str="v"}
          meta.sta = c(meta.sta,meta.str)
          
          write("          BREAK : return",file="")
          write("          BEGTR : b     ",file="")
          write("          ENDTR : e     ",file="")
          write("          OUTLI : o     ",file="")
          brk.str = readline("          Your choice              : ")
          if (brk.str == "")  {brk.str = "BREAK"}
          if (brk.str == "b") {brk.str = "BEGTR"}
          if (brk.str == "e") {brk.str = "ENDTR"}
          if (brk.str == "o") {brk.str = "OUTLI"}
          type.brk = c(type.brk,brk.str)
        }}
      flag.trend = 0
    }
  }
  
  detect.frame = cbind(index.sta,type.brk,year.sta,month.sta,meta.sta,name.sta)
  detect.frame = data.frame(detect.frame)
  names(detect.frame)=c("V1","V2","V3","V4","V5","V6")
  
  if (flag.prev == T){detect.frame = rbind(prev.frame,detect.frame)}
  
  ### Reorders by station and change-point dates
  detect.frame = detect.frame[order(detect.frame[,1],
                                    detect.frame[,3],
                                    as.numeric(detect.frame[,4])),]
  
  write.table(detect.frame,file=file.detect,quote=F,row=F,col=F,sep="\t")
  return(detect.frame)} 

########## END FUNCTION "create.brk"

#######################################
#######################################
#######################################
create.outlier = function(file,n.file){
  #######################################
  #######################################
  #######################################
  
  ####################################################
  ####         FUNCTION FOR TYPING                   #
  ####            OUTLIER DATES                      #
  ####################################################
  
  write("",file="")
  write("",file="")
  write("",file="")
  write("          ! hint: entering 13 when month is asked forces all",file="")
  write("                  months of corresponding year as outliers ",file="")
  
  write("",file="")
  write("",file="")
  write("",file="")
  
  flag.prev      = F
  index.sta      = character(0)
  name.sta       = character(0)
  year.sta       = numeric(0)
  month.sta      = numeric(0)
  type.brk       = character(0)
  
  file.detect = paste(head.str,net.str,"outliers.txt",sep="")
  if (file.exists(file.detect)){prev.frame=read.table(file.detect,head=F,
                                                      colClasses=c("character","numeric","numeric","character"),"\t");flag.prev=T}
  
  
  ## Loop over stations
  #####################
  
  for (j in 1:n.file){
    
    data.tmp        = read.serie(file[j],1)
    year.b          = data.tmp[1,1]
    year.e          = data.tmp[nrow(data.tmp),1]
    rm(data.tmp)
    
    write("",file="")
    write("",file="")
    write(paste("     ",station.index[j],station.name[j]),file="")
    write("      =========================================",file="")
    write("",file="")
    
    year.str  = "-1"
    while (year.str !="") {
      # While year not empty, continue asking for more outliers
      #########################################################
      
      write("",file="")
      year.str = readline("          New date (return to quit)      : ")
      
      if (year.str != "")  {
        if (as.numeric(year.str)>=year.b & as.numeric(year.str)<=year.e) {
          month.str = readline("          Month (13 for whole year) : ")
          if (month.str != "13") {
            index.sta = c(index.sta,file[j])              # station index
            name.sta  = c(name.sta,station.name[j])       # station name
            year.sta  = c(year.sta,as.numeric(year.str))  # year
            month.sta = c(month.sta,as.numeric(month.str))
          } else {
            index.sta = c(index.sta,rep(file[j],12))              # station index
            name.sta  = c(name.sta,rep(station.name[j],12))       # station name
            year.sta  = c(year.sta,rep(as.numeric(year.str),12))  # year
            month.sta = c(month.sta,1:12)
          }}}}
  }
  
  detect.frame = cbind(index.sta,type.brk,year.sta,month.sta,name.sta)
  detect.frame = data.frame(detect.frame)
  names(detect.frame)=c("V1","V2","V3","V4")
  
  if (flag.prev == T){detect.frame = rbind(prev.frame,detect.frame)}
  
  ### Reorders by station and change-point dates
  detect.frame = detect.frame[order(detect.frame[,1],
                                    detect.frame[,2],
                                    as.numeric(detect.frame[,3])),]
  write.table(detect.frame,file=paste(head.str,net.str,"outliers.txt",sep=""),quote=F,row=F,col=F,sep="\t")
} 
##### END FUNCTION "create.outlier"


###################################
###################################
###################################
####     HOME MAIN PROGRAM     ####
####   CAUTION : MASTERPIECE   ####
###################################
###################################
###################################
## Parameters
#############

miss.flag  <<- -999.9
max.brk    <<- 50
flag.corr  = F
corr.cycle <<- 0

net.num    = "000000"
head.str   = "       "
name.month <<- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec","ann")
options(warn=-1)
options(show.error.messages=T)

## Menu
#######
write("          ____________________",file="")
write("          ",file="")
write("              ___  _____    ",file="")
write("            . /,-Y       ~-.  ",file="")
write("            l.Y             ^.   ",file="")        
write("            /\\               _\\_  ",file="")             
write("           i            ___/     \\ ",file="")
write("           |          /     \\   o ! ",file="")  
write("           l         ]     o !__./   ",file="")
write("            \ _  _    \\.___./     ~\\  ",file="")
write("             X \\/ \\            ___./  ",file="")
write("            ( \\ ___.   _..--~~/   ~`-.  DOH!",file="")
write("             ` Z,--   /               \\  ",file="")  
write("               \\__.  (   /       ______) ",file="")
write("                 \\   l  /-----~~  /   ",file="")   
write("                  Y   \\          / ",file="")
write("                  |     x______.^ ",file="")
write("                  |           \\    ",file="")
write("                  j            Y  ",file="")
write("",file="")
write("",file="")
write("               HOMER  2.6",file="")
write("",file="")
write("          ____________________",file="")
write("          ",file="")
write("          ",file="")
write("          Dataset parameters",file="")
write("          ",file="")

dir.flag=T
file.sink<<-"tmp/sink.txt"
if (.Platform$OS.type == "unix") {
  file.sink<<-"/dev/null"
} else {if (.Platform$OS.type == "windows") {
  file.sink<<-"NUL"
}}

### Checking directories - and creates them
### If none of qc or ra exists, HOMER creates missing directories and exits

if ((!file.exists("qc")) & (!file.exists("ra"))) {dir.flag=F}

if (!file.exists("fig")){dir.create("fig");write("fig directory created",file="")}
if (!file.exists("ho")) {dir.create("ho");write("ho directory created",file="")}
if (!file.exists("meta")) {dir.create("meta");write("meta directory created",file="")}
if (!file.exists("qc")) {dir.create("qc");write("qc directory created",file="")}
if (!file.exists("ra")) {dir.create("ra");write("ra directory created",file="")}
if (!file.exists("tmp")) {dir.create("tmp");write("tmp directory created",file="")}

if (dir.flag ==T){
  write("          ",file="")
  net.str     <<- readline("      Network number (ref station file)        : ")
  
  net.str                              = as.character(net.str)
  substr(net.num,6-nchar(net.str)+1,6) = net.str
  net.str                              <<- net.num
  action.str                           = ""
  if (file.exists(paste(net.num,"stations.txt",sep=""))){
    ncol.stations = count.fields(paste(net.num,"stations.txt",sep=""),sep="\t")
    
    if (sum(ncol.stations==9)==length(ncol.stations)) {
      station=read.table(paste(net.num,"stations.txt",sep=""),sep="\t",quote="")
    } else {ncol.stations = max(count.fields(paste(net.num,"stations.txt",sep="")))
    
    ### Checks the most common error in station format
    if (ncol.stations > 9) {write("",file="");
      station.col.names = paste("V",1:ncol.stations,sep="")
      station           = read.table(paste(net.num,"stations.txt",sep=""),col.names=station.col.names,quote="",fill=T)
      write(paste(" More than 9 columns in ",net.num,"stations.txt",sep=""),file="")
      write(" Most probably spaces in station names: HOMER will remove them",file="")
      station$V9=paste(station$V9,station[,10:ncol(station)],sep="")
      station=station[,1:9]
      write(paste(" Check format of ",net.num,"stations.txt",sep=""),file="")
      write(paste(" Original station file saved in ",net.num,"stations.txt.bak",sep=""),file="")
      station=edit(station);write("",file="")
      file.copy(paste(net.num,"stations.txt",sep=""),paste(net.num,"stations.txt.bak",sep=""))
      write.table(station,file=paste(net.num,"stations.txt",sep=""),quote=F,row=F,col=F,sep="\t")
    } else { station =      read.table(paste(net.num,"stations.txt",sep=""),quote="")}
    }
    station.index <<- substr(as.character(station$V1),6,13)
    station.name  <<- as.character(station$V9)
    name.str      <<- paste(station.index,station.name)
    
    raw.str     <<- readline("      Header of input files (ex: ratx, qcrr)   : ")
    check.file    = paste(substr(raw.str,1,2),"/",raw.str,"m",station.index,"d.txt",sep="")
    
    i.remove = numeric()
    k        = 0
    for (i in check.file){
      k = k+1
      if (!file.exists(i)) {
        write("",file="")
        write(paste(i," does not exist: create it"),file="")
        grouik=readline("RETURN when done, or q to exit HOMER: ")
        if (grouik=="q"){
          options(show.error.messages=F)
          write("Inconsistencies between input (station file or header) and datafiles",file="")
          write("",file="")
          stop()}
      } else {check.frame=read.table(i,head=F)
      ### removes from the list short series
      if (nrow(check.frame) < 15 ){i.remove = c(i.remove,k)
      write(paste(" ",i,"removed from the list: less than 15 years of data"))}
      }}
    if (length(i.remove)>0){
      station = station[-i.remove,]
      if (!file.exists(paste(net.num,"stations.txt.bak",sep=""))){
        write(paste(" Original station file saved in ",net.num,"stations.txt.bak",sep=""),file="")
        file.copy(paste(net.num,"stations.txt",sep=""),paste(net.num,"stations.txt.bak",sep=""))}
      write.table(station,file=paste(net.num,"stations.txt",sep=""),quote=F,row=F,col=F,sep="\t")
      station.index <<- substr(as.character(station$V1),6,13)
      station.name  <<- as.character(station$V9)
      name.str      <<- paste(station.index,station.name)
    }    
    
    
    hom.str     <<- raw.str
    substr(hom.str,1,2)="ho"
    par.str     <<- readline("      Parameter name (for graphic outputs)     : ")
    name.str    <<- paste(par.str,name.str)
    
    unit.str    <<- readline("      Unit for graphic outputs (c for celsius) : ")
    if (unit.str == ""){unit.str<<-" "}
    write("",file="")
    write("      Parameter type",file="")
    write("      Physical parameters (Temperature, Pressure, ...)",file="")
    write("      => Additive correction       : additive (return)",file="")
    write("      Cumulative parameters (Rainfall, Sunshine Duration, ...)",file="")
    write("      => Multiplicative correction : log ratio (log) or ratio (r) comparisons",file="")
    comp.option <<- readline("      Type        : ")
    if (comp.option == "") {comp.option<<-"a"}
    if ((comp.option != "log") & (comp.option !="r")) {comp.option<<-"a"}
    
    write("",file="")
    write("      Graphic outputs",file="")
    write("      pdf (return), postscript (ps), svg (svg), png (png)",file="")
    dev.str     <<- readline("      Output option          : ")
    if (dev.str !="ps" & dev.str != "svg" & dev.str != "png") {dev.str<<-"pdf"}
    write("",file="")
    write("      Interactive option",file="")
    write("      Yes (return) or no (n)",file="")
    flag.inter  <<- readline("      Interactive option   : ")
    if (flag.inter=="n"){flag.inter=F} else {flag.inter=T}
    
    
    # Computation of lat/lon in degree and 1/100th of degree
    ########################################################
    ########################################################
    lat.station   <<- station$V2+station$V3/60.+station$V4/3600
    lon.station   <<- station$V5+station$V6/60.+station$V7/3600
    alt.station   <<- station$V8
    xy.station    <-  mapproject(x=lon.station,y=lat.station,proj="lambert",param=range(lat.station))
    x.sta         <<- xy.station$x*6378
    y.sta         <<- xy.station$y*6378
    
    
    write("",file="")
    write("      Intercomparison Neighbourhood",file="")
    write("      All series (return), geographic (g) or correlation (c) distance",file="")
    inter.option <<- readline("      Intercomparison type          : ")
    if (inter.option == "g") {inter.number <<- as.numeric(readline("      Maximum distance (km)         : "))}
    if (inter.option == "c") {inter.number <<- as.numeric(readline("      Minimum correlation r         : "))}
    if (inter.option == "c" | inter.option == "g") {
      write("      !! Warning, next parameter",file="")
      write("      !! superseeds r.min or d.max",file="")
      n.min <<- as.numeric(readline("      Minimum number of neighbours  : "))}
    if (inter.option == "")  {inter.number <<- -1; n.min <<- -1}
    write("",file="")
    write("     Season comparison option for pairwise detection",file="")
    write("     Annual+seasons (return), annual (a) or monthly (m)",file="")
    season.option <<- readline("     Season option        : ")
    
    season.corr <<- "m"
    if (comp.option != "a"){
      write("     Correction option for cumulative parameters",file="")
      write("     Annual coefficient estimation (return, default) or monthly (m) warning: not recommanded",file="")
      season.corr <<- readline("     Correction option        : ")}
    if (season.corr != "m") {season.corr<<-""}
    
    write("",file="")
    write("     Options for series visualization",file="")
    trend.str      <<-  readline("     Linear trend?          yes (return)/n : ")
    smooth.str     <<-  readline("     Smoothing option?      yes (return)/n : ")
    if (trend.str  !="n"){ trend.str  <<- "y" }
    if (smooth.str !="n"){ smooth.str <<- "y" }
    col.str        <<-  readline("     Polygon fill?          yes (return)/n : ")
    if (col.str!="n"){col.fill=readline("     return for red/blue gy=green/yellow.. : ")
    if (col.fill==""){col.fill="rb"}
    col1.str=substr(col.fill,1,1)
    col2.str=substr(col.fill,2,2)
    if (col1.str=="r"){col1.str <<- "red"}
    if (col1.str=="b"){col1.str <<- "blue"}
    if (col1.str=="g"){col1.str <<- "green"}
    if (col1.str=="y"){col1.str <<- "yellow"}
    if (col2.str=="r"){col2.str <<- "red"}
    if (col2.str=="b"){col2.str <<- "blue"}
    if (col2.str=="g"){col2.str <<- "green"}
    if (col2.str=="y"){col2.str <<- "yellow"}
    col.str = "y"
    }
    
    write("",file="")
    write("          ____________________",file="")
    
    if (file.exists(paste("tmp/",raw.str,net.str,"option.sav",sep=""))) {load(paste("tmp/",raw.str,net.str,"option.sav",sep=""))}
    
    ### if station file not found, HOMER exits
  } else {options(show.error.messages=F)
    write("",file="")
    write(paste("!  ",net.num,"stations.txt does not exist: create it",sep=""),file="")
    stop()}
  
  while (action.str != "q") {
    write("",file="")
    write("",file="")
    write("     What do you wish, Master/Mistress?",file="")
    
    # Automatic User gender detection to be implemented soon
    ########################################################
    write("",file="")
    write("      FAST QUALITY CONTROL",file="")
    write("      -> Fast CLiMATOL checks               type i",file="")
    write("      -> Fast QC                            type f",file="")
    write("      -> Outlier file creation?             type o",file="")
    write("      -> Removal of outliers?               type r",file="")
    write("",file="")
    write("      HOMOGENISATION",file="") 
    write("      -> Pairwise detection?                type d",file="")
    write("      -> Joint detection?                   type j",file="")
    if (flag.corr==T & comp.option=="a"){
      write("      -> ACMANT detection?                  type a",file="")}
    if (flag.corr==T){
      write("      -> Assess Month of change             type m",file="")}
    write("      -> Correction?                        type c",file="")
    write("      -> Visualization?                     type v",file="")
    write("      -> New neighbourhood                  type n",file="")
    write("      -> Change hinteraction hoption :-)    type h",file="")
    write("      -> Break file creation/modification?  type b",file="")
    write("      -> Break file edition?                type e",file="")
    write("      -> Quit?                              type q",file="")
    write(" ",file="")
    
    header.actions    = c("i","f","o","r","d","j","v")
    no.header.actions = c("a","c","n","b","e","h","m","q")
    action.str        = readline("     Your choice    : ")
    if (!is.na(match(action.str,header.actions))) {
      head.str     = readline("     raw/qc (return) or corrected (c) files    : ")
      if (head.str != "c") {head.str = raw.str} else {head.str = hom.str}
    } else {head.str = raw.str}
    
    # File list constitution
    ########################
    list.file     = paste(substr(head.str,1,2),"/",head.str,"m",station.index,"d.txt",sep="")
    list.file.ref = cbind(paste("tmp/",head.str,"m",station.index,"d.txt",sep=""),
                          paste("tmp/ho",substr(head.str,3,4),"m",station.index,"d.txt",sep=""))
    n.file        = length(list.file)
    
    # Possible actions
    ##################
    
    if (action.str == "i") {
      cli.dat=home.input.data.check(list.file)
      save(cli.dat,file="tmp/cli.sav")
      data.check(cli.dat[[1]],length(list.file),cli.dat[[2]],cli.dat[[3]])
    }
    
    if (action.str == "f") {fqc(list.file,n.file,head.str)}
    if (action.str == "h") {
      write("",file="")
      write("      Interactive option",file="")
      write("      Yes (return) or no (n)",file="")
      flag.inter  <<- readline("      Interactive option   : ")
      if (flag.inter=="n"){flag.inter=F} else {flag.inter=T}
    }
    
    if (action.str == "r") {rem.out(head.str)}
    
    if (action.str == "d") {
      detect(list.file,n.file,head.str,crit="CAU",13)
      
      if (season.option == "") {
        for (season in 16:19){
          detect(list.file,n.file,head.str,crit="CAU",season)}}
      
      if (season.option == "m") {
        for (season in 1:12){
          detect(list.file,n.file,head.str,crit="CAU",season)}}
    }
    
    if (action.str == "j") {
      toto=j.detect(list.file,head.str,13,crit="CAU",flag.inter)
    }
    
    if (action.str == "c"){ 
      
      correction.output=correc(list.file,n.file,head.str,list.file.ref)
      flag.corr = T
      save(flag.corr,corr.cycle,file=paste("tmp/",head.str,net.str,"option.sav",sep=""))
      head.str  = hom.str 
      write("",file="")
      write("",file="")
      write("",file="")
      
      write("Running Pairwise detection on corrected series",file="")
      write("==============================================",file="")
      list.file = paste(substr(head.str,1,2),"/",head.str,"m",station.index,"d.txt",sep="")
      detect(list.file,n.file,head.str,crit="CAU",13)
      if (season.option == "") {
        for (season in 16:19){
          detect(list.file,n.file,head.str,crit="CAU",season)}}
      if (season.option == "m") {
        for (season in 1:12){
          detect(list.file,n.file,head.str,crit="CAU",season)}}
    }
    
    if (flag.corr==T & comp.option=="a" & action.str=="a"){
      toto=acmant.detect(list.file,list.file.ref,n.file,head.str,flag.inter)
    }
    
    if (flag.corr==T & comp.option=="a" & action.str=="m"){
      toto=change.month(list.file,list.file.ref,n.file,head.str)
    }
    
    if (action.str == "v") {visu(list.file,n.file,head.str)}
    if (action.str == "o") {create.outlier(list.file,n.file)}
    
    if (action.str == "b") {create.brk(list.file,n.file)}
    if (action.str == "e") {
      file.detect = paste(net.str,"detected.txt",sep="")
      if (file.exists(file.detect)) {
        detect.frame = read.table(file.detect,head=F,
                                  colClasses=c("character","character","numeric","numeric","character",
                                               "character"),sep="\t")
        detect.frame = edit(detect.frame)
        write.table(detect.frame,file=file.detect,quote=F,row=F,col=F,sep="\t")
      } else {
        write("",file="")
        write("     Warning: create break file first!",file="")
        write("",file="")
      }
    }
    
    if (action.str == "n") {
      
      write("     Intercomparison option",file="")
      write("     All series (return), geographic (g) or correlation (c) distance",file="")
      write(" ",file="")
      inter.option <<- readline("     Intercomparison type        : ")
      if (inter.option == "g") {inter.number <<- as.numeric(readline("     Maximum distance (km)       : "))}
      if (inter.option == "c") {inter.number <<- as.numeric(readline("     Minimum correlation r       : "))}
      if (inter.option == "c" | inter.option == "g") {
        write("     !! Warning, next parameter",file="")
        write("     !! superseeds r.min or d.max",file="")
        
        n.min <<- as.numeric(readline("     Minimum number of neighbours: "))
      }
      if (inter.option == "")  {inter.number <<- -1; n.min <<- -1}
      write("",file="")}
  }
  
  write(" ",file="")
  write(" ",file="")
  write(" ",file="")
  graphics.off()
  
  write("                      ,;-. ",file="")
  write("                    ,((--\\). ",file="")
  write("                   /        \\ ",file="")
  write("                  |          | ",file="")
  write("                  |          | ",file="")
  write("                 (, `.  ,  `.) ",file="")
  write("                 :     \\/     ;      BYE!",file="")
  write("                 `.o  , `.  o /  ",file="")
  write("                 (|`> `-- `< |)            ,-, ",file="")
  write("     ,.           |/        \\|         ,-./ / ",file="")
  write("   _ | \\,-.       (          )        | `- `--. ",file="")
  write("  ( `  (_/|__      \\   (o   /       ,-      ,-  ",file="")
  write("   ;         )    ,|`.  - , |.       -.   ) \\ ",file="")
  write("   | (    ,-    _/ `-.`   ,-  \\---.   /      ; ",file="")
  write("   |     |   ,-   \\  /\\  / \\  |   |--/       | ",file="")
  write("   |     |_,|    / \\/  \\/   \\/\\   |          | ",file="")
  write("   |     `  \\   |              \\  /        , ",file="")
  write("   |         \\  |              | /      _, ",file="")
  write("   :          \\ ,              `/------' ",file="")
  write("    `-.___,--- )                `. ",file="")
  write("             ,                    \\ ",file="")
  write("            /                      \\ ",file="")
  write("           :                        : ",file="")
  write("           |                      _,| ",file="")
  write("            \\--.___         __,--  ; ",file="")
  write("             `.    `========      , ",file="")
  write("               |                  | ",file="")
  write("               |      .____,      | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |         |        | ",file="")
  write("               |-._____,-|-.____,-| ",file="")
  write("               |_        |_       | ",file="")
  write("             ,   `------ | `----- \\ ",file="")
  write("            /           _|_         \\ ",file="")
  write(",            `--._____,-    `-.___,-  ",file="")
  write("",file="")} else {
    write("",file="")
    write("!  No data directories found",file="")
    write("!  HOMER directories created",file="")
    write("!  Provide data files in ra or qc directories",file="")}


## BYE BYE
##########