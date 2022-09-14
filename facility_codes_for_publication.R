####################
### Codes for the project that applied generalized synthetic control on facility closure
### in California 
### This document included main analysis and Figure 3/4/S1
### Prepared by Chen Chen (chc048@ucsd.edu) on 9/8/2022
####################
## set up the stage
library(data.table)

indir1 <- "data"
indir2 <- "data/hospitalization_exposure_data"
outdir2 <- "results"
outdir3 <- "figures"

dataset <- "hyads_report_40km_by_site" ##main analysis

## define exposed and unexposed zipcodes using hyad data (median in the month of 
## closure during the year before within 40km as cutpoint) and reported retirement date used
## created one analysis for each facility closure that contains all case zctas and potential control zctas 
library(raster) # for spatial analysis
library(fst)
####################
ca_impacts_all <- read.fst(file.path(indir1, "zips_exposures_byfacility_CA_coal_oil.fst"), as.data.table = T) #Hyads data
zip <-fread(file.path(indir1, "zctas_combo_2010_and_2000_CZ.csv")) #California zip data
hyads_ca <- ca_impacts_all[ca_impacts_all$ZIP %in% zip$ZCTA5, ] #subset hyads data to CA only

##US census bureau datum
crs.us <- CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80")
## create a spatial variable for zcta population center
pop.cen <- SpatialPointsDataFrame(zip[, .(FinalLon, FinalLat)], 
                                  zip[, "ZCTA5"], 
                                  proj4string = crs.us)

## Date of closure for oil and coal facilities in California
closure <- fread(file.path(indir1, "coal_oil_ca.csv"), colClasses = c("EIA_ID"="character"))[, c(1:2, 5:12)]
cl <- SpatialPointsDataFrame(closure[, .(LONG, LAT)], 
                             closure, 
                             proj4string = crs.us)

## Identify case and control zctas
buf <- 40000
## create case buffers for facility locations
fac.buf <- buffer(cl, width=buf, dissolve=FALSE)
## identify zctas within case buffers
buff.zp <- over(fac.buf, pop.cen, returnList = TRUE)
names(buff.zp) <- cl$EIA_ID
closure$withinbuf <- sapply(buff.zp, nrow)

## identify the reported closure dates
closure$pr_month <- month(closure$Retire)
closure$pr_yr <- year(closure$Retire)- 1
closure$hyads_av <- "no available"

## list of variables of interest
vb <- c("date", "zcta5", "GEO_ID", "NAME", "pop", #basic info
        "mean_pm25", "tmpt", "dptp", #exposures
        "resp_disease", "cardio_disease" #health outcomes
)
record.case <- data.frame(zcta5=zip$ZCTA5,
                          setNames(replicate(length(closure$EIA_ID), NA, simplify = FALSE), 
                                   paste0("included_", closure$EIA_ID)))
closure$ncase <- closure$ncontrol <- closure$ncase.av <- closure$ncontrol.av <- 0
closure$ncase.inside <- closure$ncase.outside <- 0
closure$m.outside <- closure$threshold <- NA


for (j in 1:length(closure)) {
  dt <- data.table()
  temp <- hyads_ca[uID==closure$EIA_ID[j], ]
  if (nrow(temp[temp$year==closure$pr_yr[j] & temp$month==closure$pr_month[j]]) > 0) {
    temp.w <- reshape(temp[, 1:4], v.names = "hyads", timevar = "yearmonth", idvar = "ZIP", direction = "wide")
    temp.w[is.na(temp.w)] <- 0 ##replace missing value with 0--some zips don't have hyads value
    temp <- temp.w[, c("ZIP", paste0("hyads.", closure$pr_yr[j], closure$pr_month[j])), with=FALSE]
    names(temp)[2] <- "hyads"
    
    closure$hyads_av[j] <- "year prior"
    closure$threshold[j] <- median(temp$hyads[temp$ZIP %in% unlist(buff.zp[[closure$EIA_ID[j]]])], na.rm = TRUE)
    closure$m.outside[j] <- median(temp$hyads[!(temp$ZIP %in% unlist(buff.zp[[closure$EIA_ID[j]]]))], na.rm = TRUE)
    cases <- temp$ZIP[temp$hyads > closure$threshold[j] & temp$ZIP %in% unlist(buff.zp[[closure$EIA_ID[j]]])]
    controls <- temp$ZIP[temp$hyads <= closure$threshold[j] & !(temp$ZIP %in% unlist(buff.zp[[closure$EIA_ID[j]]]))]
    closure$ncase[j] <- closure$ncase.inside[j] <- length(cases)
    closure$ncase.outside[j] <- length(temp$ZIP[temp$hyads > closure$threshold[j] & !(temp$ZIP %in% unlist(buff.zp[[closure$EIA_ID[j]]]))])
    closure$ncontrol[j]  <- length(controls)
    day <- closure$Retire[j]
    range <- c(day-28, day+28, day-182, day+182) #time range around retirement date
    
    for (i in cases){
      file.name <- file.path(indir2, paste0(i, ".rds"))
      if (file.exists(file.name)) {
        bar <- readRDS(file.name)
        bar <- bar[zcta5==i, vb, with=FALSE] ## only keep related variables and corresponding zcta
        bar <- bar[bar$date>=range[3] & bar$date<=range[4], ] ## only keep 26 weeks before and after retirement
        ## only keep zctas with at least 3 measurements in 4 weeks before and 4 weeks after retirement
        if (nrow(bar[!is.na(bar$mean_pm25) & bar$date>=range[1] & bar$date<day]) > 3 &
            nrow(bar[!is.na(bar$mean_pm25) & bar$date<=range[2] & bar$date>day]) > 3) {
          bar$case <- 1
          dt <- rbind(dt, bar)
          record.case[record.case$zcta5==i, paste0("included_", closure$EIA_ID[j])] <- 1
          closure$ncase.av[j] <- closure$ncase.av[j] + 1
        }
        # check for exposure missingness
        bar <- NULL
      } # existence of data file
    } # loop for potential case zctas for each facility
    
    for (i in controls){
      file.name <- file.path(indir2, paste0(info$FIPS[info$ZCTA5==i], ".rds"))
      if (file.exists(file.name)) {
        bar <- readRDS(file.name)
        bar <- bar[zcta5==i, vb, with=FALSE] ## only keep related variables and corresponding zcta
        bar <- bar[bar$date>=range[3] & bar$date<=range[4], ] ## only keep 26 weeks before and after retirement
        ## only keep zctas with at least 3 measurements in 4 weeks before and 4 weeks after retirement
        if (nrow(bar[!is.na(bar$mean_pm25) & bar$date>=range[1] & bar$date<day]) > 3 &
            nrow(bar[!is.na(bar$mean_pm25) & bar$date<=range[2] & bar$date>day]) > 3) {
          bar$case <- 0
          dt <- rbind(dt, bar)
          record.case[record.case$zcta5==i, paste0("included_", closure$EIA_ID[j])] <- 0
          closure$ncontrol.av[j] <- closure$ncontrol.av[j] + 1
        }
        # check for exposure missingness
        bar <- NULL
      } # existence of data file
    } # loop for potential control zctas for each facility
    
    dt$EIA_ID <- closure$EIA_ID[j]
    dt$closure.date <- day
    dt$after <- ifelse(dt$date<day, 0, 1)
    dt$after[dt$date==day] <- 2 #set the day of retirement to 2
    
    ## each file contains all case zctas and potential control zctas for a facility
    ## after controling on data missingness
    saveRDS(dt, file.path(indir1, dataset, paste0(closure$EIA_ID[j], ".rds")))
  }
} # loop for facility
write.csv(record.case, file.path(indir1, dataset, "summary_zctas.csv"), row.names = FALSE)
write.csv(closure, file.path(indir1, dataset, "summary_closures.csv"), row.names = FALSE)
####################

## Create weekly sum datast for all closed facility, apply synthetic control analysis
## and make gsc plots (Figure 3 and supplementary Figure 1)
# remotes::install_github("xuyiqing/gsynth")
library(gsynth)
####################
record.case <- fread(file.path(indir1, dataset, "summary_closures.csv"))
sites <- unique(record.case[hyads_av != "no available" & ncase.av > 0, EIA_ID])

## create dataset
bar.dt <- numeric()
for (site in sites) {
  dt <- readRDS(file.path(indir1, dataset, paste0(site, ".rds")))
  dt$temp <- as.numeric(difftime(dt$date, dt$closure.date, units = "week"))
  dt$difweek <- floor(dt$temp) # treat the day of closure as exposed
  dt$after <- ifelse(dt$after==2, 1, dt$after)
  
  foo1 <- aggregate(cbind(mean_pm25, tmpt, dptp) ~ zcta5 + EIA_ID + GEO_ID + after + case + difweek,
                    FUN= mean, data = dt)
  foo2 <- aggregate(cbind(resp_disease, cardio_disease) ~ zcta5 + closure.date + pop + difweek,
                    FUN= sum, data = dt)
  baz <- merge(foo1, foo2, by = c("zcta5", "difweek"))
  
  baz$treatedpost <- baz$case * baz$after
  baz$cvdresp_rate <- (baz$resp_disease + baz$cardio_disease)/baz$pop * 10000
  
  bar.dt <- rbind(bar.dt, baz)
}
saveRDS(bar.dt, file.path(indir1, dataset, "all_sites_by_week.rds"))

## run gsc
nboots <- 1000
vbs <- c("mean_pm25",
         "cvdresp_rate")
tls <- c(expression(Weekly~PM[2.5]~concentration~(mu*g/m^3)),
         expression("Weekly rate of cardiorespiratory hosp. ("*10^3~"PT)"))

bar.dt <- readRDS(file.path(indir1, dataset, "all_sites_by_week.rds"))

bar <- numeric()
for (site in sites) {
  dt <- bar.dt[bar.dt$EIA_ID == site, ]
  setDT(dt)
  dt <- dt[difweek <= 3, ] # included all 26 weeks before closure and 4 weeks after closure
  
  baz <- numeric()
  for (i in 1:length(vbs)){
    model <- reformulate(c("treatedpost", "tmpt", "dptp"), response = vbs[i])
    
    f <- tryCatch({
      gsynth(model, data = dt, EM = F, index = c("zcta5", "difweek"), na.rm=TRUE,
             inference = "parametric", se = TRUE, nboots = nboots, seed = 1234, r = c( 0, 5),
             CV = TRUE, force = "two-way", parallel = TRUE, cores = 4)
    }, condition = function(cond) {
      cat("\t", site, vbs[i], as.character(cond))
      cond$call <- NULL
      cond
    })
    
    if (class(f)=="gsynth") {
      outfile <- file.path(outdir3, dataset, paste0("Weekly_facility_", site,"_", vbs[i], "_GSC.png"))
      png(outfile,
          width=6, height=4, units="in", res = 600, bg="white", family="sans")
      print(plot(f, type = "counterfactual", raw = "band", main = "", legendOff = TRUE,
                 xlab = "Time to closure (week)", ylab = " ", theme.bw = TRUE, shade.post = TRUE))
      grid.text(tls[i], 0.02, .55, rot=90)
      
      dev.off()
    }
    
    baz <- c(baz, list(f))
  }
  names(baz) <- vbs
  bar <- c(bar, list(baz))
}
names(bar) <- sites
saveRDS(bar, file.path(outdir2, dataset, "Weekly_GSC_all_facility.rds"))
####################

## pooled results across retirements and make forest plot (Figure 4)
library(meta)
library(grid)
####################
vbs <- c("mean_pm25",
         "cvdresp_rate")
tls <- c(expression(Weekly~PM[2.5]~concentration~(mu*g/m^3)),
         expression("Weekly rate of cardiorespiratory hosp. ("*10^3~"PT)"))

bar <- readRDS(file.path(outdir2, dataset, "Weekly_GSC_all_facility.rds"))
for (i in 1:length(vbs)) {
  vb <- vbs[i]
  foo <- lapply(bar, "[[", vb)
  fal <- foo[sapply(foo, inherits, what = "condition")]
  cat(vb, "reasons for failing:", unique(unlist(fal)), "\n")
  cat(vb, "number of failed facilities:", length(fal), "\n")
  foo <- foo[!sapply(foo, inherits, what = "condition")]
  
  dif <- as.data.frame(do.call("cbind", lapply(foo, function(a){
    data.frame(est=paste0(round(a$est.att[,"ATT"], 1), " (", round(a$est.att[,"S.E."], 2), ")"))
  })))
  names(dif) <- names(foo)
  dif$difweek <- -26:3
  print(dif)
  
  d <- as.data.frame(do.call("rbind", lapply(foo, function(a){
    a$est.avg
  })))
  d$facility <- names(foo)
  out <- metagen(TE = d$Estimate, seTE = d$S.E., studlab = d$facility)
  print(summary(out))
  
  png(file.path(outdir3, dataset, paste0("forestplot of", vb, ".png")),
      width=8, height=4.5, units="in", res = 600, bg="white", family="sans")
  forest.meta(out, layout = "JAMA", xlab=" ",
              leftcols = c("studlab", "effect.ci"), leftlabs = c("Facility", "Estimate [95% CI]"))
  if (vb=="mean_pm25") {
    grid.text(tls[i], .7, .1, gp=gpar(cex=1)) ## for PM
  } else {
    grid.text(tls[i], .7, .05, gp=gpar(cex=1))  ## for ha
  }
  dev.off()
}
####################
