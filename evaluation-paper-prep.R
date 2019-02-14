# Script for Zealandia's data preparation:
# - clock adjustment
# - map cropping

options(stringsAsFactors = F)
library(hms)
library(lubridate)
library(ggplot2)
library(dplyr)
library(rjson)
library(rgdal)
library(tidyr)
library(raster)
library(rasterVis)


annotdir = "~/Documents/kiwis/results/manual1007/"
outdir = "~/Documents/kiwis/results/evaluationpaper/"
# Map file downloaded from LINZ data service, data.linz.govt.nz,
# Wellington 0.3m resolution rural aerial photos, BQ31_0506 tile
MAPLOCATION = "~/Documents/kiwis/scripts/zealandia/BQ31_0506.tif"


## 0. FUNCTION INIT

## Converts segments to T/F presence per 1 s blocks
# Needs: input dataframe dfin (in segment format)
#    starttime, a ymd_hms parseable string for initial time at t=0
#    trueTimes, a clock adjustment table with rec & adj columns (or NULL)
convertto01 = function(dfin, starttime, clockadjs){
  dfout = data.frame(secs = seq(hms::as.hms("06:00:00"), hms::as.hms("18:00:00")),
                     ZA = F, ZB = F, ZE = F, ZG = F, ZH = F, ZI = F, ZJ = F)
  for(r in 1:nrow(dfin)){
    # extract call start-end times in seconds, rounded
    callstart = as.numeric(dfin$startAdj[r] - ymd_hms(starttime), units="secs")
    callend = ceiling(callstart + dfin$calllength[r])
    rec = dfin$rec[r]
    
    if(is.data.frame(clockadjs)){
      # actually apply clock adjustments
      adj = clockadjs$adj[clockadjs$rec==rec]
      callstart = callstart+adj
      callend = callend+adj
    }
    
    # mark that as 1 in main df
    dfout[floor(callstart):callend, rec] = TRUE
  }
  return(dfout)
}

# given a nx2 matrix, shifts first column by t \in (-10;10) min
# calculates overlap
binaryPointMatch = function(pair){
  matches = data.frame(ts = -600:600, sharedSeconds=0)
  for(d in seq_along(matches$ts)){
    dd = matches$ts[d]
    if(dd >= 0){
      shifted = c(rep(0, dd), pair[1:(nrow(pair)-dd), 1])  
    } else {
      dd = -dd
      shifted = c(pair[(1+dd):nrow(pair), 1], rep(0, dd))
    }
    
    matches$sharedSeconds[d] = sum(pair[,2] & shifted)
  }
  bestShift = matches[which.max(matches$sharedSeconds),]
  bestShift$names = paste(colnames(pair), collapse="_")
  return(bestShift)
}

plotOverlaps = function(pair, df){
  bestShift = binaryPointMatch(df[,pair])
  print(bestShift)
  
  if(bestShift$ts >= 0){
    df$lagged = lag(df[,pair[1]], n=bestShift$ts, default=0)
  } else {
    df$lagged = lead(df[,pair[1]], n=-bestShift$ts, default=0)
  }
  df$both = as.numeric(df[,pair[1]] & df[,pair[2]])
  df$bothl = as.numeric(df$lagged & df[,pair[2]])
  df = gather(df, key="rec", value="present", ZA:bothl) %>%
    mutate(hour = sprintf("%02d", secs %/% 3600), minsec = secs %% 3600) %>%
    filter(present==1)
  
  return(bestShift)
}

#####

## 1a. MANUALLY ANNOTATED DATA
recorders = c("ZA", "ZB", "ZC", "ZD", "ZE", "ZF", "ZG", "ZH", "ZI", "ZJ", "ZK")

# parse annotations
annot = data.frame()
for(rec in recorders){
  gooddata = list.files(annotdir, pattern=paste0(rec, "_.*wav.data"))
  for(f in gooddata){
    a = fromJSON(file=paste0(annotdir, f))
    if(length(a)>1){
      a = a[-1] # drop metadata
      a = data.frame(t(sapply(a, c))) # to dataframe  
      a$time = parse_date_time(substr(f, 4, 18), "Ymd_HMS")
      a$rec = rec
      annot = rbind(annot, a)
    }
  }
}
# cleanup and change data types
annot$start = annot$time + seconds(annot$X1)
annot$end = annot$time + seconds(annot$X2)
annot$species = unlist(annot$X5)
annot = annot[,6:10]

annot$calllength = annot$end - annot$start
annot$startAdj = annot$start-dhours(12)
annot$dateAdj = date(annot$startAdj)

# sanity checks
table(annot$species)
table(annot$dateAdj)
qplot(annot$calllength)


## 1b. CLOCK ADJUSTMENT

# convert to 0/1 per second blocks
presence6 = convertto01(annot, "2018-10-06 06:00:00", NULL)

bestShifts = combn(colnames(presence6)[-1], 2, function(x) plotOverlaps(x, presence6))
bestShifts = matrix(unlist(bestShifts), ncol=dim(bestShifts)[2], byrow=T)
bestShifts = data.frame(bestShifts)

## OVERLAPS after reviewing on AviaNZ:
### clear matches:
# ZH/ZJ - close recs, +88 on ZH confirmed
# ZA/ZJ - great, +105 on ZA confirmed
# ZA/ZH - +18 on ZA without doubt
# ZB/ZI - best example, -210 on ZB confirmed
# ZG/ZI - -213 on ZG without doubt
# ZB/ZG - +3 on ZB looks legit
# ZE/ZG, ZE/ZJ, ZA/ZB, ZA/ZE, ZA/ZI - no good overlaps found. ZA/ZG - maybe, but mostly not

# infer clock shifts from known
matchedRecs = c("ZH_ZJ", "ZA_ZJ", "ZA_ZH", "ZB_ZI", "ZG_ZI", "ZB_ZG")
bestShifts$confirmed = bestShifts$X3 %in% matchedRecs
bestShifts$rec1 = substr(bestShifts$X3, 1, 2)
bestShifts$rec2 = substr(bestShifts$X3, 4, 5)


### Dijkstra's algorithm to fill time shifts
propagateDistance = function(currRec){
  print(paste("current", currRec))
  visited <<- c(visited, currRec)
  
  distFromCurrent = bestShifts[bestShifts$confirmed & bestShifts$rec1==currRec,]
  print(paste("gonna test 1", distFromCurrent$rec2))
  for(r in distFromCurrent$rec2){
    print(r)
    currPosition = trueTimes$adj[trueTimes$rec == currRec]
    currAdj = -as.numeric(distFromCurrent$X1[distFromCurrent$rec2==r]) + currPosition
    
    if(!r %in% visited){
      if (trueTimes$evidence[trueTimes$rec==r] == "confirmed" &
          trueTimes$adj[trueTimes$rec==r] != currAdj){
        print(paste("conflicting distances ", trueTimes$adj[trueTimes$rec==r], currAdj))
      } else {
        trueTimes$adj[trueTimes$rec==r] <<- currAdj
        trueTimes$evidence[trueTimes$rec==r] <<- "confirmed"
        propagateDistance(r)
      }  
    } else {
      print(paste("skipping", r, "distance would be changed from", trueTimes$adj[trueTimes$rec==r], "to", currAdj))
    }
  }
  distFromCurrent = bestShifts[bestShifts$confirmed & bestShifts$rec2==currRec,]
  print(paste("gonna test 2", distFromCurrent$rec1))
  for(r in distFromCurrent$rec1){
    print(r)
    currPosition = trueTimes$adj[trueTimes$rec == currRec]
    currAdj = as.numeric(distFromCurrent$X1[distFromCurrent$rec1==r]) + currPosition
    
    if(!r %in% visited){
      if (trueTimes$evidence[trueTimes$rec==r] == "confirmed" &
          trueTimes$adj[trueTimes$rec==r] != currAdj){
        print(paste("conflicting distances ", trueTimes$adj[trueTimes$rec==r], currAdj))
      } else {
        trueTimes$adj[trueTimes$rec==r] <<- currAdj
        trueTimes$evidence[trueTimes$rec==r] <<- "confirmed"
        propagateDistance(r)
      }
    } else {
      print(paste("skipping", r, "distance would be changed from", trueTimes$adj[trueTimes$rec==r], "to", currAdj))
    }
  }  
}

# one node is assumed to be true and at 0
trueTimes <<- data.frame(rec=colnames(presence6)[-1], adj=0, evidence="none")
visited <<- "ZA"
trueTimes$evidence[trueTimes$rec == visited] = "confirmed"

# now take all recorders that have confirmed distances to it and propagate
propagateDistance(visited)
trueTimes

# repeat with another starting node because graph is separated
visited <<- "ZG"
trueTimes$evidence[trueTimes$rec == visited] = "confirmed"
propagateDistance(visited)
trueTimes

write.table(trueTimes, paste0(outdir, "clockadjs.tsv"), sep="\t", quote=F, row.names=F, col.names=T)


## 2. GPS PREP
# create recorder grid, using fixed distances for now
gpspos = read.table("ZealandiaPointsOld.txt", sep=" ", nrows=10, h=F)
gpspos$V1 = paste0("Z", substr(gpspos$V1, 2, 2))
ggplot(gpspos, aes(y=V2, x=V3)) + geom_point() + geom_text(aes(label=V1), nudge_y=0.0005) +
  theme(aspect.ratio = 1)
coordinates(gpspos) = c("V3", "V2")
proj4string(gpspos) = CRS("+proj=longlat +datum=WGS84")

# crop main huge map tiff:
gpsmap = readGDAL(MAPLOCATION)
trapbounds = bbox(spTransform(gpspos, CRS(proj4string(gpsmap))))
trapbounds[,2] - trapbounds[,1] # 800 x 1500 m section
# so add 200 m buffer left and down:
trapbounds[,1] = trapbounds[,1] - 200
# and small buffer up and left:
trapbounds[,2] = trapbounds[,2] + 20
# get offset in m from start
trapbounds = trapbounds - gridparameters(gpsmap)[,1]
# get offset in pixels from start
trapbounds = round(trapbounds/0.3)
# extract and save subset
writeGDAL(gpsmap[trapbounds[2,1]:trapbounds[2,2], trapbounds[1,1]:trapbounds[1,2]],
          paste0(outdir, "BQ31_subset.tif"), drivername="GTiff", type="Byte")


## 3. PLOT THE RASTER + RECORDER MAP
# read in cropped map
gpsmap = readGDAL(paste0(outdir, "BQ31_subset.tif"))
bbox(gpsmap)
bbox(spTransform(gpspos, CRS(proj4string(gpsmap))))

# make sure projections match
gpspos = spTransform(gpspos, CRS(proj4string(gpsmap)))
gpsposM = data.frame(gpspos)
colnames(gpsposM) = c("rec", "east", "north", "optional")
gpsposM = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK"))

# for experimenting:
#submap = gpsmap[1:500, 1:500]

# convert RGB to a colortable and assign it
gpsmapcols = SGDF2PCT(gpsmap[,,1:3])
gpsmap = raster(gpsmap)
gpsmap = setValues(gpsmap, gpsmapcols$idx)
colortable(gpsmap) = gpsmapcols$ct

# save it - note the maxpixels arg!!
p = gplot(gpsmap, maxpixels=8e6) +
  geom_raster(aes(fill=value)) +
  geom_point(aes(x=east, y=north), color="white", pch=4, data=gpsposM) +
  geom_text(aes(x=east, y=north, label=rec), color="white", data=gpsposM, vjust=0, nudge_y=20) +
  coord_cartesian(ylim=c(5425000, 5426300), xlim=c(1745150, 1746000)) +
  xlab("Easting, m") + ylab("Northing, m") + 
  scale_fill_gradientn(colours=colortable(gpsmap), guide="none") +
  theme_bw() + theme(aspect.ratio = 1)

ggsave(plot=p, filename=paste0(outdir, "recorder_map.tiff"), device = "tiff", dpi=300)

