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
bind_rows(man=annot, auto_rev=annot2, auto_raw=annot4, .id="detection") %>%
  group_by(detection) %>%
  summarize(sum(calllength) / 60)

# marked fractions by filter
bind_rows(man=annot, auto_rev=annot2, auto_raw=annot4, .id="detection") %>%
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


# OVERALL CONCORDANCE MEASURES
concordance = group_by(concordanceF, man, auto) %>%
  summarize(duration=n(), rate=duration /7 /12/3600) # total time = 7 recs * 12 hr
concordance
summarize(concordanceF, sens = mean(man&auto)/mean(man),
          fpr = mean(auto& !man)/mean(!man),
          prec = mean(man&auto)/mean(auto),
          spec = mean(!man & !auto)/mean(!man),
          acc = mean(man & auto) + mean(!man & !auto), ncall=mean(man))

# MEASURES BY RECORDER
group_by(concordanceF, rec, man, auto) %>%
  summarize(duration=n(), rate=duration /12/3600*100) %>%
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
cr = fit.ascr(calls6, traps, mask)
summary(cr)

# BOOTSTRAP FOR COMPARISON
#bb = boot.ascr(cr, 400, n.cores = 7)
#summary(bb)
# SE D 0.2073 (was 0.2153)
# SE g0 0.0412 (was 1e-06)
# SE sigma 7.7121 (was 7.1693)
#get.mce(bb, "se")
# g0 0.0019, sigma 0.281, D 0.00255

show.detfn(cr, xlim=c(0, 700))
show.detsurf(cr)

# example locations
par(mfrow=c(2,3), mar=c(0,0,1,1))
for(i in c(1, 50, 100, 150, 200, 250)){
  locations(cr, i, xlim=c(-500,500), ylim=c(-600, 400), levels=c(0.1, 0.3, 0.5, 0.7, 0.9), plot.estlocs = TRUE)  
}
par(mfrow=c(1,1), mar=c(5,5,4,4))


## 6. DROP SOME RECORDERS

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

cr = fit.ascr(calls6A, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))


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

cr = fit.ascr(calls6A, traps, mask)
summary(cr)
show.detfn(cr, xlim=c(0, 700))


## 7. SPEED AND POWER
# idea: more data should have smaller SEs, by 1/sqrt(n).
# test: take various amounts of auto-processed data, fit ASCR, measure SE(density).

trapgrid = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK")) %>%
  mutate(recid=1:n())
traps = as.matrix(trapgrid[,c("east", "north")])

mask = create.mask(traps, buffer=700)

speedplot = data.frame()

## FULL 2 NIGHTS
cr = fit.ascr(calls67A, traps, mask)
summary(cr)
speedplot = rbind(speedplot, c(2, coef.ascr(cr), stdEr.ascr(cr)))
names(speedplot) = c("duration", "D", "g0", "sigma", "SE_D", "SE_g0", "SE_sigma")

## REALISTIC SUBSETS (night 1, night 2)
cr = fit.ascr(calls6A, traps, mask)
summary(cr)
speedplot = rbind(speedplot, c(1, coef.ascr(cr), stdEr.ascr(cr)))

# store these as reference
bestD = coef.ascr(cr)[["D"]]
bestSE = stdEr.ascr(cr)[["D"]]

cr = fit.ascr(calls7A, traps, mask)
summary(cr)
speedplot = rbind(speedplot, c(1, coef.ascr(cr), stdEr.ascr(cr)))

# store these as reference
bestD = (coef.ascr(cr)[["D"]] + bestD)/2
bestSE = (stdEr.ascr(cr)[["D"]] + bestSE)/2

## SUBSAMPLING
for(fr in c(0.25, 0.5, 1, 1.5)){
  # take first 30 min of each hour:
  capt6A = filter(presence6A, secs %% 3600 >= fr/2*60*60) %>%
    converttoCapt()
  capt7A = filter(presence7A, secs %% 3600 >= fr/2*60*60) %>%
    converttoCapt()
  
  # bind, convert to capt:
  capt67A = bind_rows(mutate(capt6A, start=paste0("6_", start)),
                      mutate(capt7A, start=paste0("7_", start)))
  capt67A = inner_join(capt67A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
  calls67A = create.capt(as.data.frame(capt67A), n.traps=7, n.sessions = NULL)
  
  # fit ascr:
  cr = fit.ascr(calls67A, traps, mask)
  summary(cr)
  # save result:
  speedplot = rbind(speedplot, c(fr, coef.ascr(cr), stdEr.ascr(cr)))
}
speedplot

# Figure 3
p1 = mutate(speedplot, D=D/duration, SE_D=SE_D/duration) %>%
  mutate(theor=bestSE / sqrt(duration)) %>% 
  ggplot() + geom_point(aes(x=duration, y=SE_D)) +
  geom_line(aes(x, y=bestSE/sqrt(x)), data=data.frame(x=seq(0.2, 2, 0.1)), col="grey50") +
  theme_bw() + xlab("duration, nights")

mutate(speedplot, D=D/duration) %>%
  ggplot(aes(x=duration)) + geom_point(aes(y=D)) +
  geom_ribbon(aes(ymin=D-1.96*SE_D/duration, ymax=D+1.96*SE_D/duration), alpha=0.2, fill="grey") +
  theme_bw()


# given effect size, what is the power as f(number of nights)?
muDiff = 2.89*0.2
mutate(speedplot, power = 1-pnorm(zCr, mean=muDiff/sqrt(2*(SE_D/duration)^2), sd=1))

p2 = data.frame(x=seq(0.1, 5, 0.05)) %>%
  mutate(power=1-pnorm(zCr, mean=muDiff/sqrt(2*bestSE^2/x))) %>%
  ggplot() + geom_line(aes(x=x, y=power)) + 
  ylim(c(0,1)) +
  theme_bw() + xlab("duration, nights")

plot_grid(p1, p2, nrow=2, labels="AUTO")


## 8. CASE STUDY: TIME OF SAMPLING
# as for speed, but instead of taking various size at fixed time (start of hour),
# take various times of fixed time (3 hr?)

trapgrid = filter(gpsposM, !rec %in% c("ZC", "ZF", "ZK")) %>%
  mutate(recid=1:n())
traps = as.matrix(trapgrid[,c("east", "north")])
mask = create.mask(traps, buffer=700)

timeplot = data.frame(time=c("full night", "18-21", "21-00", "00-03", "03-06",
                             "18-20", "20-22", "22-00", "00-02", "02-04", "04-06"),
                      timemin=c(6, 6, 9, 12, 15, 6, 8, 10, 12, 14, 16),
                      timemax=c(18, 9, 12, 15, 18, 8, 10, 12, 14, 16, 18),
                      D=0, g0=0, sigma=0, SE_D=0, SE_g0=0, SE_sigma=0)

for(t in 1:nrow(timeplot)){
  # take first 30 min of each hour:
  capt6A = filter(presence6A, secs / 3600 >= timeplot$timemin[t],
                  secs / 3600 < timeplot$timemax[t]) %>%
    converttoCapt()
  capt7A = filter(presence7A, secs / 3600 >= timeplot$timemin[t],
                  secs / 3600 < timeplot$timemax[t]) %>%
    converttoCapt()
  
  # bind, convert to capt:
  capt67A = bind_rows(mutate(capt6A, start=paste0("6_", start)),
                      mutate(capt7A, start=paste0("7_", start)))
  capt67A = inner_join(capt67A, trapgrid, by="rec")[,c("session", "start", "occ", "recid")]
  calls67A = create.capt(as.data.frame(capt67A), n.traps=7, n.sessions = NULL)
  
  # fit ascr:
  cr = fit.ascr(calls67A, traps, mask)
  summary(cr)
  # save result:
  timeplot[t,4:6] = coef.ascr(cr)
  timeplot[t,7:9] = stdEr.ascr(cr)
}
timeplot

p1 = mutate(timeplot, D=D/(timemax-timemin)*12, SE_D=SE_D/(timemax-timemin)*12) %>%
  ggplot(aes(x=reorder(time, 1:nrow(timeplot)))) + geom_point(aes(y=D)) +
  geom_errorbar(aes(ymin=D-1.96*SE_D, ymax=D+1.96*SE_D)) +
  ylim(c(0, 12)) + theme_bw() + xlab("recording time")
p2 = ggplot(timeplot, aes(x=reorder(time, 1:nrow(timeplot)))) +
  geom_point(aes(y=sigma)) +
  geom_errorbar(aes(ymin=sigma-1.96*SE_sigma, ymax=sigma+1.96*SE_sigma)) +
  ylim(c(120, 250)) + theme_bw() + xlab("recording time")

plot_grid(p1, p2, nrow=2, labels="AUTO")
