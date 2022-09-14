library(tidyverse)
library(ggplot2)
library(sf)
library(ggthemes)
library(gridExtra)
library(mapproj)
library(lubridate)
library(ggrepel)

#load data
ca <- subset(map_data("state"), region == "california")
county <- subset(map_data("county"), region == "california")
fac <-read.csv("fac_id.csv")
em<-read.csv("fac_emissions_so2_for_map2.csv")
ca<-st_read("exp_maps.shp")

#prep data
#create new facility numbers
fac2<-filter(fac, exclude == 0) %>% 
      arrange(fac_id) %>%
      select(-exclude) %>%
      mutate(fac_id2= 1:11)

dat<-merge(em, fac2, by.x="EIA_ID", by.y="fac_id")

dat2<-mutate(dat, so2_cat=cut(so2_year, breaks=c(0, 10, 14.9, 100, 150), labels=c("<10", "10-14", "15-100", ">100")))

#create map
setwd("figures")
map<-ggplot(ca)+
  geom_polygon(data=county, aes(x=long, y=lat, group=group), fill="white")+
  geom_polygon(data=county, aes(x=long, y=lat, group=group), color="#a8a7a7", fill="NA", alpha=0.8, size=0.5)+
  
  geom_text_repel(data=dat, aes(x=long, y=lat, label=fac_id2), segment.size=0.25, size=2.5, min.segment.length = 0,
                  max.overlaps = Inf, box.padding=unit(.45, "lines"))+

  geom_point(data=dat2, aes(x=long, y=lat, fill=so2_cat), shape=21, stroke=.05, size=2)+
  scale_fill_manual(values=c("#ffffd9", "#c7e9b4", "#1d91c0", "#081d58"))+
 
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  theme_map(base_size=9, base_family="")

ggsave(filename="fac_map.png", 
       plot=map, 
       dpi=600, 
       width=7.7,
         height=5.2,
         units="in")
  





