## Script for Zealandia's data analysis

options(stringsAsFactors = F)
library(lubridate)
library(ggplot2)
library(dplyr)
library(rjson)
library(tidyr)
library(ascr)
library(rgdal)
library(cowplot)

outdir = "~/Documents/kiwis/results/evaluationpaper/"

## Read in results from prep script
trueTimes = read.table(paste0(outdir, "clockadjs.tsv"), h=T)


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

# Converts to make.capt style capture history.
# define separate calls as having >5 s breaks,
# extract only the start row for each call,
# and gather to get separate row for each recorder.
converttoCapt = function(dfin){
  # everything done in three copies
  # to keep detections from non-overlapping recorder groups separate
  o1 = filter(dfin, occ=="o1") %>%
    mutate(pause=secs-lag(secs, default=0), callID=cumsum(pause>5)) %>%
    group_by(callID) %>%
    summarize(start = paste0(min(secs), unique(occ)), ZA=any(ZA), ZH=any(ZH), ZJ=any(ZJ)) %>%
    gather(key="rec", value="pres", ZA:ZJ) %>%
    filter(pres) %>%
    select(-pres, -callID)
  
  o2 = filter(dfin, occ=="o2") %>%
    mutate(pause=secs-lag(secs, default=0), callID=cumsum(pause>5)) %>%
    group_by(callID) %>%
    summarize(start = paste0(min(secs), unique(occ)), ZE=any(ZE)) %>%
    gather(key="rec", value="pres", ZE:ZE) %>%
    filter(pres) %>%
    select(-pres, -callID)
  
  o3 = filter(dfin, occ=="o3") %>%
    mutate(pause=secs-lag(secs, default=0), callID=cumsum(pause>5)) %>%
    group_by(callID) %>%
    summarize(start = paste0(min(secs), unique(occ)), ZB=any(ZB), ZG=any(ZG), ZI=any(ZI)) %>%
    gather(key="rec", value="pres", ZB:ZI) %>%
    filter(pres) %>%
    select(-pres, -callID)
  
  # merge, add fake variables, return
  bind_rows(o1, o2, o3) %>%
    mutate(session=1, occ=1) %>%
    return()
}

# Reads all AviaNZ style annotations from directory dir,
# over recorders recs (string vector).
# Does some basic date conversions and returns a df.
readAnnots = function(dir, recs){
  annot = data.frame()
  for(rec in recs){
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
    if(length(gooddata)==0) print(paste("Warning: no files found for recorder", rec))
  }
  if(nrow(annot)==0) print("Warning: no files read!")
  
  # actual start and end of call
  annot$start = annot$time + seconds(annot$X1)
  annot$end = annot$time + seconds(annot$X2)
  
  # JSON annotations are read in as list columns, so convert/drop those
  annot$species = unlist(annot$X5)
  annot = annot[,6:10]
  
  annot$calllength = annot$end - annot$start
  # 12 hour shift makes working with nights easier:
  # all recordings from one night start on same day and run 6 to 18,
  # instead of 18 to 00 and then 00 to 06 on another day
  annot$startAdj = annot$start-dhours(12)
  annot$dateAdj = date(annot$startAdj)
  
  return(annot)
}


## 1. MANUALLY ANNOTATED DATA
recorders = c("ZA", "ZB", "ZC", "ZD", "ZE", "ZF", "ZG", "ZH", "ZI", "ZJ", "ZK")
annotdir = "~/Documents/kiwis/results/manual1007/"

# parse annotations
annot = readAnnots(annotdir, recorders)

# sanity checks
table(annot$species)
table(annot$dateAdj)
qplot(annot$calllength)

# apply clock adjustments and re-convert
# (clock adjustments calculated in PREP script)
presence6 = convertto01(annot, "2018-10-06 06:00:00", trueTimes)


## 2. AUTO-DETECTED DATA
annotdir = "~/Documents/kiwis/results/detect1007/"
annot2 = readAnnots(annotdir, recorders)

# sanity checks
table(annot2$species)
table(annot2$dateAdj)
qplot(annot2$calllength)

# convert to 0/1 per second blocks
presence6A = convertto01(annot2, "2018-10-06 06:00:00", trueTimes)


## SAME FOR THE OTHER DAY
annotdir = "~/Documents/kiwis/results/detect1008F/"
annot3 = readAnnots(annotdir, recorders)

# sanity checks
table(annot3$species)
table(annot3$dateAdj)
qplot(annot3$calllength)

# convert to T/F per second blocks
presence7A = convertto01(annot3, "2018-10-07 06:00:00", trueTimes)


## LOAD RAW AUTO ANNOTATIONS
annotdir = "~/Documents/kiwis/results/detect1007/raw/"
annot4 = readAnnots(annotdir, recorders)

# sanity checks
table(annot4$species)
table(annot4$dateAdj)
qplot(annot4$calllength)

annotdir = "~/Documents/kiwis/results/detect1008F/raw/"
annot5 = readAnnots(annotdir, recorders)

# sanity checks
table(annot5$species)
table(annot5$dateAdj)
qplot(annot5$calllength)


## 3. GET TP/FP MEASURES
## (AFTER HUMAN REVIEW)

# plot of all annotations for Oct 6
bind_rows(man=annot, auto_rev=annot2, auto_raw=annot4, .id="detection") %>%
  ggplot(aes(x=start)) +
  #geom_point(aes(y=detection), col="red", pch="|", size=4) +
  geom_segment(aes(y=detection, yend=detection, xend=end), col="red", lwd=6) +
  facet_wrap(~rec, nrow=4) +
  scale_x_datetime(date_minor_breaks = "1 hour", date_breaks = "3 hours", date_labels = "%H:%M") + 
  theme_bw() + theme(axis.title = element_blank(), legend.title=element_blank(), legend.position = "bottom")

# marked fractions
bind_rows(man=annot, auto_rev=annot2, auto_rev=annot3, auto_raw=annot4, auto_raw=annot5, .id="detection") %>%
  group_by(detection) %>%
  summarize(sum(calllength) / 60)

# marked fractions by filter
bind_rows(man=annot, auto_rev=annot2, auto_rev=annot3, auto_raw=annot4, auto_raw=annot5, .id="detection") %>%
  group_by(detection, species) %>%
  summarize(sum(calllength) / 60)

# total num of calls
nrow(annot)
table(annot$calllength>5)

# if we allow overlaps, how many calls were detected by both methods?
annot$tpsOverlapping = FALSE
for(a in 1:nrow(annot)){
  callstart = floor(as.numeric(ymd_hms(annot$startAdj[a])-ymd_hms("2018-10-06 06:00:00"), units="secs"))
  callend = ceiling(callstart + annot$calllength[a])
  rec = annot$rec[a]
  callstart = callstart + trueTimes$adj[trueTimes$rec==rec]
  callend = callend + trueTimes$adj[trueTimes$rec==rec]
  annot$tpsOverlapping[a] = any(presence6A[callstart:callend, rec])
}
table(annot$tpsOverlapping[annot$calllength>10])

# long df of results for each second w/ each method
concordanceF = full_join(gather(presence6, key="rec", value="man", ZA:ZJ),
          gather(presence6A, key="rec", value="auto", ZA:ZJ),
          by=c("rec", "secs"))
concordanceF$annot = "TN"
concordanceF$annot[concordanceF$man & concordanceF$auto] = "TP"
concordanceF$annot[concordanceF$man & !concordanceF$auto] = "FN"
concordanceF$annot[!concordanceF$man & concordanceF$auto] = "FP"
ggplot(concordanceF[concordanceF$annot!="TN",], aes(x=secs, y=rec)) +
  geom_segment(aes(xend=max(secs), x=min(secs), yend=rec), col="grey") +
  geom_point(aes(color=annot), pch="|", size=4) +
  scale_color_manual(values=c("TP"="chartreuse3", "TN"="grey", "FP"="darkorange", "FN"="red")) +
  theme_bw()


# OVERALL CONCORDANCE MEASURES, 1 second resolution
concordance = group_by(concordanceF, man, auto) %>%
  summarize(duration=n(), rate=duration /7 /12/3600) # total time = 7 recs * 12 hr
concordance
summarize(concordanceF, sens = mean(man&auto)/mean(man),
          fpr = mean(auto& !man)/mean(!man),
          prec = mean(man&auto)/mean(auto),
          spec = mean(!man & !auto)/mean(!man),
          acc = mean(man & auto) + mean(!man & !auto), ncall=mean(man))

## side things: export concordance for 6 hr of recorders
filter(concordanceF, secs/3600 >= 10, secs/3600 < 16) %>%
  group_by(man, auto) %>%
  summarize(duration=n(), rate=duration /7 /6/3600) # total time = 7 recs * 12 hr
filter(concordanceF, secs/3600 >= 10, secs/3600 < 16) %>%
  summarize(sens = mean(man&auto)/mean(man),
          fpr = mean(auto& !man)/mean(!man),
          prec = mean(man&auto)/mean(auto),
          spec = mean(!man & !auto)/mean(!man),
          acc = mean(man & auto) + mean(!man & !auto), ncall=mean(man))
filter(annot4, startAdj >= ymd_hms("2018-10-06 10:00:00"),
        startAdj < ymd_hms("2018-10-06 16:00:00")) %>% 
  group_by(species) %>%
  summarize(sum(calllength)/60)
filter(annot, startAdj >= ymd_hms("2018-10-06 10:00:00"),
       startAdj < ymd_hms("2018-10-06 16:00:00")) %>% 
  summarize(n(), sum(calllength)/60)
filter(annot, startAdj >= ymd_hms("2018-10-06 10:00:00"),
       startAdj < ymd_hms("2018-10-06 16:00:00")) %>% 
  filter(calllength >= 5) %>%
  summarize(n(), sum(calllength)/60)
filter(annot2, startAdj >= ymd_hms("2018-10-06 10:00:00"),
       startAdj < ymd_hms("2018-10-06 16:00:00")) %>% 
  group_by(species) %>%
  summarize(sum(calllength)/60)


# MEASURES BY RECORDER
group_by(concordanceF, rec, man, auto) %>%
  summarize(duration=n(), rate=duration /12/3==600*100) %>%
  print.data.frame

specbyrec = group_by(concordanceF, rec) %>%
  summarize(sens = mean(man&auto)/mean(man),
            fpr = mean(auto& !man)/mean(!man),
            prec = mean(man&auto)/mean(auto),
            spec = mean(!man & !auto)/mean(!man),
            acc = mean(man & auto) + mean(!man & !auto), ncall=mean(man))
specbyrec
write.table(format(data.frame(specbyrec), digits=4, scientific=F),
            paste0(outdir, "concordance_by_recs.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

group_by(concordanceF, rec) %>%
  summarize(sens = mean(man&auto)/mean(man)-0.689,
          fpr = mean(auto& !man)/mean(!man),
          prec = mean(man&auto)/mean(auto)-0.8348,
          spec = mean(!man & !auto)/mean(!man)-0.99264,
          acc = mean(man & auto) + mean(!man & !auto)-0.9771, ncall=mean(man)-0.0512)

# DROP COUPLE RECORDERS AND SEE HOW IT CHANGES
recstodrop = list(c("none"), c("ZA", "ZI"), c("ZA", "ZB"), c("ZB", "ZI"), c("ZB", "ZH"))
concordance2 = data.frame()
for(d in recstodrop){
  res = filter(concordanceF, ! rec %in% d) %>%
    summarize(sens = mean(man&auto)/mean(man),
              fpr = mean(auto& !man)/mean(!man),
              prec = mean(man&auto)/mean(auto),
              spec = mean(!man & !auto)/mean(!man),
              acc = mean(man & auto) + mean(!man & !auto), ncall=mean(man)) %>%
    mutate(droppedRecs = paste(d, collapse=", "))
  
  concordance2 = bind_rows(concordance2, res)
  
  if(d!="none"){
    # subtract full data values
    res[1, 1:6]  = res[1, 1:6] - concordance2[1, 1:6]
    concordance2 = bind_rows(concordance2, res)  
  }
}
concordance2
# ZA, ZI - improves everything a bit, mostly acc
# ZB, ZH - improves sens more, keeps acc and ncall unchanged
write.table(format(data.frame(concordance2), digits=4, scientific=F),
            paste0(outdir, "concordance_dropped.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

## OVERALL CONCORDANCE, 15 second resolution
concordanceF15 = group_by(concordanceF, t=secs %/% 15, rec) %>%
  summarize(anyTman = any(man), anyTauto = any(auto)) %>%
  ungroup()
concordance15 = group_by(concordanceF15, anyTman, anyTauto) %>%
  summarize(duration=n(), rate=duration /7 /12/60/4) # total time = 7 recs * 12 hr * 60 min * 4 windows/min

concordance15
summarize(concordanceF15, sens = mean(anyTman&anyTauto)/mean(anyTman),
          fpr = mean(anyTauto& !anyTman)/mean(!anyTman),
          prec = mean(anyTman&anyTauto)/mean(anyTauto),
          spec = mean(!anyTman & !anyTauto)/mean(!anyTman),
          acc = mean(anyTman & anyTauto) + mean(!anyTman & !anyTauto), ncall=mean(anyTman))

# DROP COUPLE RECORDERS, 15 second resolution
recstodrop = list(c("none"), c("ZA", "ZI"), c("ZB", "ZH"))
concordance2 = data.frame()
for(d in recstodrop){
  res = filter(concordanceF15, ! rec %in% d) %>%
    summarize(sens = mean(anyTman&anyTauto)/mean(anyTman),
              fpr = mean(anyTauto& !anyTman)/mean(!anyTman),
              prec = mean(anyTman&anyTauto)/mean(anyTauto),
              spec = mean(!anyTman & !anyTauto)/mean(!anyTman),
              acc = mean(anyTman & anyTauto) + mean(!anyTman & !anyTauto), ncall=mean(anyTman)) %>%
    mutate(droppedRecs = paste(d, collapse=", "))
  
  concordance2 = bind_rows(concordance2, res)
  
  if(d!="none"){
    # subtract full data values
    res[1, 1:6]  = res[1, 1:6] - concordance2[1, 1:6]
    concordance2 = bind_rows(concordance2, res)  
  }
}
concordance2
write.table(format(data.frame(concordance2), digits=4, scientific=F),
            paste0(outdir, "concordance_dropped_15.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

## 4. PREP DATA FOR ASCR

## separate out ZE, ZA-ZH-ZJ, ZB-ZG-ZI into occasions
presence6 = bind_rows(o1=filter(presence6, ZA | ZH | ZJ),
          o2=filter(presence6, ZE),
          o3=filter(presence6, ZB | ZG | ZI), .id="occ")
presence6A = bind_rows(o1=filter(presence6A, ZA | ZH | ZJ),
          o2=filter(presence6A, ZE),
          o3=filter(presence6A, ZB | ZG | ZI), .id="occ")
presence7A = bind_rows(o1=filter(presence7A, ZA | ZH | ZJ),
          o2=filter(presence7A, ZE),
          o3=filter(presence7A, ZB | ZG | ZI), .id="occ")

# GPS DATA:
# create recorder grid, using fixed distances for now
gpspos = read.table("ZealandiaPointsOld.txt", sep=" ", nrows=10, h=F)
gpspos$V1 = paste0("Z", substr(gpspos$V1, 2, 2))
ggplot(gpspos, aes(y=V2, x=V3)) + geom_point() + geom_text(aes(label=V1), nudge_y=0.0005) +
  theme(aspect.ratio = 1)

# project lat / long to easting / northing in meters 
coordinates(gpspos) = c("V3", "V2")
proj4string(gpspos) = CRS("+proj=longlat +datum=WGS84")
gpsposM = spTransform(gpspos, CRS("+proj=utm +zone=60G"))
gpsposM = data.frame(gpsposM)
colnames(gpsposM) = c("rec", "east", "north", "optional")
gpsposM = mutate(gpsposM, east=east-mean(east), north=north-mean(north), optional=NULL)
ggplot(gpsposM, aes(y=north, x=east)) +
  geom_point() + geom_text(aes(label=rec), nudge_y=50) +
  theme_minimal() + theme(aspect.ratio = 1)

# 7 recs only
filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK")) %>%
  ggplot(aes(y=north, x=east)) +
  geom_point() + geom_text(aes(label=rec), nudge_y=50) +
  theme_minimal() + theme(aspect.ratio = 1)

trapgrid = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK")) %>%
  mutate(recid=1:n())
traps = as.matrix(trapgrid[,c("east", "north")])

capt6 = converttoCapt(presence6)
capt6A = converttoCapt(presence6A)
capt7A = converttoCapt(presence7A)

# check detections over recorder combos
bind_rows(.id="detection", man6=capt6, auto6=capt6A, auto7=capt7A) %>%
  group_by(detection, start) %>%
  summarize(recs=paste(unique(rec), collapse=" ")) %>%
  group_by(recs, detection) %>%
  summarize(n=n()) %>%
  print.data.frame

# convert recNames to recIDs & make sure columns are in right order
capt6 = left_join(capt6, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
capt6A = left_join(capt6A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
capt7A = left_join(capt7A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]

calls6 = create.capt(as.data.frame(capt6), n.traps = 7, n.sessions = NULL)
calls6A = create.capt(as.data.frame(capt6A), n.traps = 7, n.sessions = NULL)
calls7A = create.capt(as.data.frame(capt7A), n.traps = 7, n.sessions = NULL)

calls67A = bind_rows(mutate(capt6A, start=paste0("6_", start)),
             mutate(capt7A, start=paste0("7_", start))) %>%
  as.data.frame() %>%
  create.capt(., n.traps = 7, n.sessions = NULL)

# create mask for all 7 traps
mask = create.mask(traps, buffer=700)
plot(mask)


## 5. FIT ASCR

# set this to TRUE for the first run (will take some hours)
bootSEs = FALSE

## full 2 nights, automatic annotations
capt6A = converttoCapt(presence6A)
capt7A = converttoCapt(presence7A)

# bind, convert to capt:
capt67A = bind_rows(mutate(capt6A, start=paste0("6_", start)),
                    mutate(capt7A, start=paste0("7_", start)))
capt67A = inner_join(capt67A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
calls67A = create.capt(as.data.frame(capt67A), n.traps=7, n.sessions = NULL)

cr = fit.ascr(calls67A, traps, mask)
summary(cr)

show.detfn(cr, xlim=c(0, 700))
show.detsurf(cr)

# example locations
par(mfrow=c(2,3), mar=c(0,0,1,1))
for(i in c(1, 50, 100, 150, 200, 250)){
  locations(cr, i, xlim=c(-500,500), ylim=c(-600, 400), levels=c(0.1, 0.3, 0.5, 0.7, 0.9), plot.estlocs = TRUE)  
}
par(mfrow=c(1,1), mar=c(5,5,4,4))

# bootstrap SEs or read saved ones?
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 4)
  summary(bb)
  save(bb, file=paste0(outdir, "/data/boot-67a.RData"))
} else {
  load(paste0(outdir, "/data/boot-67a.RData"), verbose=T)
  bb67A = bb
}

# we can confirm that the number of boot iterations is sufficient by
# get.mce(bb, "se")

speedplot = data.frame()
speedplot = rbind(speedplot, c(2, coef.ascr(cr), stdEr(cr), stdEr(bb)))
names(speedplot) = c("duration", "D", "g0", "sigma", "SE_D", "SE_g0", "SE_sigma", "SE_D_b", "SE_g0_b", "SE_sigma_b")


## 6. SUBSAMPLE DATA (speed and power)

for(fr in c(0.25, 0.5, 1, 1.5)){
  # take first X min of each hour:
  capt6A = filter(presence6A, secs %% 3600 <= fr/2*60*60) %>%
    converttoCapt()
  capt7A = filter(presence7A, secs %% 3600 <= fr/2*60*60) %>%
    converttoCapt()
  
  # bind, convert to capt:
  capt67A = bind_rows(mutate(capt6A, start=paste0("6_", start)),
                      mutate(capt7A, start=paste0("7_", start)))
  capt67A = inner_join(capt67A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
  calls67A = create.capt(as.data.frame(capt67A), n.traps=7, n.sessions = NULL)
  
  # fit ascr.
  # Note: we're not using the survey.length parameter,
  # because bootstrap generates surveys of length 1
  # which gives bad SEs (i.e. corresponding to surveys of length 1, not 0.25 or whatever)
  cr = fit.ascr(calls67A, traps, mask)
  summary(cr)
  if(bootSEs){
    bb = boot.ascr(cr, 200, n.cores = 4)
    summary(bb)
    save(bb, file=paste0(outdir, "/data/boot-67a-", fr, ".RData"))
  } else {
    load(paste0(outdir, "/data/boot-67a-", fr, ".RData"), verbose=T)
  }
  
  # save result:
  speedplot = rbind(speedplot, c(fr, coef.ascr(cr), stdEr(cr), stdEr(bb)))
}
speedplot


# Figure 3
# single frame version
# (note the division by survey length to account for the issue discussed above)
gather(speedplot, key="method", value="SE", c(SE_D, SE_D_b)) %>%
  ggplot(aes(x=duration, y=SE/duration, col=method, group=method)) + geom_point() +
  scale_color_discrete(labels=c("asymptotic", "bootstrapped"), name=NULL) +
  geom_smooth(method="nls", formula = y ~ a/sqrt(x), method.args=list(start=c(a=1)), se=F) +
  theme_bw() + xlab("duration, nights") + ylab(expression(SE(hat(D)))) +
  theme(legend.position = c(0.8,0.8))


## 7. COMPARE WITH MANUAL ANNOTATIONS

# manual:
cr = fit.ascr(calls6, traps, mask)
summary(cr)
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6.RData"))
} else {
  load(paste0(outdir, "/data/boot-6.RData"), verbose=T)
}

# inspect detection function g0*(exp(-1/2*d^2/s^2))
show.detfn(cr, xlim=c(0, 700))
1*exp(-1/2*200^2/200.46^2)
abline(h=0.608, col="red")
abline(v=200, col="red")

# automatic, one night:
cr = fit.ascr(calls6A, traps, mask)
summary(cr)
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6a.RData"))
} else {
  load(paste0(outdir, "/data/boot-6a.RData"), verbose=T)
}

# can check the other automatic night as well
# cr = fit.ascr(calls7A, traps, mask)
# summary(cr)

# inspect detection function
show.detfn(cr, xlim=c(0, 700))
abline(h=0.543, col="red")
abline(v=200, col="red")

# EDR
detectionfn = function(x){
  s = 180.934
  1 * exp(-1/2 * x^2 / s^2)
}
# equiv area of certain detection
eda = integrate(function(r) 2*pi*r*detectionfn(r), 0, Inf)$value
edr = sqrt(eda/pi)
detectionfn(edr)
# example: missed counts within EDR
pi*edr^2 - integrate(function(r) 2*pi*r*detectionfn(r), 0, edr)$value
# detected outside EDR
integrate(function(r) 2*pi*r*detectionfn(r), edr, Inf)

## 8. DROP SOME RECORDERS

# dropping recorders ZA, ZI
trapgrid = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK", "ZA", "ZI")) %>%
  mutate(recid=1:n())
traps = as.matrix(trapgrid[,c("east", "north")])

mask = create.mask(traps, buffer=700)

capt6 = converttoCapt(presence6)
capt6 = inner_join(capt6, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
calls6 = create.capt(as.data.frame(capt6), n.traps=5, n.sessions = NULL)

capt6A = converttoCapt(presence6A)
capt6A = inner_join(capt6A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
calls6A = create.capt(as.data.frame(capt6A), n.traps=5, n.sessions = NULL)

# fit without ZA, ZI
cr = fit.ascr(calls6, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6-zazi.RData"))
} else {
  load(paste0(outdir, "/data/boot-6-zazi.RData"), verbose=T)
}

cr = fit.ascr(calls6A, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6a-zazi.RData"))
} else {
  load(paste0(outdir, "/data/boot-6a-zazi.RData"), verbose=T)
}

# dropping recorders ZB, ZH
trapgrid = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK", "ZH", "ZB")) %>%
  mutate(recid=1:n())
traps = as.matrix(trapgrid[,c("east", "north")])
mask = create.mask(traps, buffer=700)

capt6 = converttoCapt(presence6)
capt6 = inner_join(capt6, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
calls6 = create.capt(as.data.frame(capt6), n.traps=5, n.sessions = NULL)

capt6A = converttoCapt(presence6A)
capt6A = inner_join(capt6A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
calls6A = create.capt(as.data.frame(capt6A), n.traps=5, n.sessions = NULL)

# fit without ZB, ZH
cr = fit.ascr(calls6, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6-zbzh.RData"))
} else {
  load(paste0(outdir, "/data/boot-6-zbzh.RData"), verbose=T)
}

cr = fit.ascr(calls6A, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))
if(bootSEs){
  bb = boot.ascr(cr, 200, n.cores = 6)
  print(summary(bb))
  save(bb, file=paste0(outdir, "/data/boot-6a-zbzh.RData"))
} else {
  load(paste0(outdir, "/data/boot-6a-zbzh.RData"), verbose=T)
}


