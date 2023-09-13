# Alberio Lab Image analysis script - generalized

library(tidyverse)
library(ggplot2)
library(BiocManager)
library(readr)
library(dplyr)
library(data.table)
library(plotly)
library(factoextra)
library(viridis)
library(reshape2)
library(gplots)
library(ggpubr)
library(remotes)
library(gatepoints)
library(scSHC)

channelid<-"x=dapi y=T z=sox2 t=foxa2"
idX<-"DAPI"
idY<-"TBXT"
idZ<-"TFAP2C"
idT<-"SOX17"

# Set the path to the directory where your csv files are
setwd("D:/One drive/OneDrive - The University of Nottingham/General/AlberioLab/Andrew/Zeiss and LSM900 images/ssegmentation test/csv")
getwd()->wd
# Creates an empty list to populate your csvs into - call it whatever you want - Gastruloids/Embryo/Unicorns
Gastruloids <- list()
# Creates the list of all the csv files in the directory
listcsv <- dir(pattern = "*") 

# Populate you list with all the csvs from the directory
for (k in 1:length(listcsv)){
  
  Gastruloids[[k]] <- read.csv(listcsv[k])
  
}

# Add file name as a column - th path will be the same as your working directory
all_paths <-
  list.files(path = wd,
             pattern = "*.csv",
             full.names = TRUE)

all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")

all_filenames <- all_paths %>%
  basename() %>%
  as.list()

all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

all_result <- rbindlist(all_lists, fill = T)

# Change column name
names(all_result)[3] <- "File.Path"

Gastruloid_data <- rbindlist(Gastruloids, fill = TRUE)

# Create a data frame with the file names (all_results$V1) appended onto the csv data as the first column
Gastruloid_data <- cbind(all_result$V1, Gastruloid_data)

# Check the column names you have specified and whether the transform command has worked
colnames(Gastruloid_data)
head(Gastruloid_data)


# Split by for 2D and 3D analysis the Fiji macro will output 5 csvs per image (based on a 4 channel image) with each nuclei/cell represented by a row in each of the 5 csvs. We now need to combine those rows but first we need to split the data into each of the channels (and the morphology measurements) for each object.
Gastruloid_morphology <- subset(Gastruloid_data, V1 %like% "Morphological")
Gastruloid_channel_1 <- subset(Gastruloid_data, V1 %like% "C1")
Gastruloid_channel_2 <- subset(Gastruloid_data, V1 %like% "C2")
Gastruloid_channel_3 <- subset(Gastruloid_data, V1 %like% "C3")
Gastruloid_channel_4 <- subset(Gastruloid_data, V1 %like% "C4")


# Combine the channels across rows so that all the information for each object is in one row
Gastruloid_data_all_channels <- cbind(Gastruloid_morphology, x = Gastruloid_channel_1, y = Gastruloid_channel_2, z = Gastruloid_channel_3, t = Gastruloid_channel_4)

# There will now be many duplicate and empty columns - don't worry these can be ignored or you can remove them if you wish
colnames(Gastruloid_data_all_channels)

# Make a new column called Gastruloid_number_and_stage or any combine the two columns with these individual pieces of data. Also useful for combining stage and conditions e.g. Gastruloid_data_combined$condition <-  paste(Gastruloid_data_combined$Ratio_of_PT_to_SC,Gastruloid_data_combined$PT_cell_line)
Gastruloid_data_all_channels$Gastruloid_number_and_stage <- paste(Gastruloid_data_all_channels$Gastruloid_number, Gastruloid_data_all_channels$Gastruloid_stage, sep = "_")

# Check that this has worked 
head(Gastruloid_data_all_channels$Gastruloid_number_and_stage)

# Make a three new columns with the centre position of the xyz coordinates for each object - useful for spatially plotting data later (for 2D data you don't need to have the z coordinate)
Gastruloid_data_all_channels$X_centre <- (Gastruloid_data_all_channels$Xmin..pix. + Gastruloid_data_all_channels$Xmax..pix.)/2
Gastruloid_data_all_channels$Y_centre <- (Gastruloid_data_all_channels$Ymin..pix. + Gastruloid_data_all_channels$Ymax..pix.)/2
Gastruloid_data_all_channels$Z_centre <- (Gastruloid_data_all_channels$Zmin..pix. + Gastruloid_data_all_channels$Zmax..pix.)/2


# Filter for objects which have a high enough DAPI signal - make a histogram of the mean DAPI expression in each object. This could be x/y/z/t.Mean depending on which channel of your image DAPI was.
ggplot(data = Gastruloid_data_all_channels, aes(x=x.Mean)) + 
  geom_histogram(binwidth=1) 

# Set limits on the x-axis to make it easier to interpret the histogram
#ggplot(data = Gastruloid_data_all_channels, aes(x=x.Mean)) + 
#  geom_histogram(binwidth=20) +
#  coord_cartesian(xlim = c(0, 19000))+
#  scale_x_continuous(breaks = seq(0, 12000, 500))





# Threshold from histogram suggests any object with a mean DAPI signal below X is not a real nuclei.
# Set threshold on the DAPI mean with an upper and lower limit which excludes doublets and erroneous objects. use this line to test
# Filter for nuclei which are too small to be real
ggplot(data = Gastruloid_data_all_channels, aes(x=Vol..unit.)) + 
  geom_histogram(binwidth=20) 



ggplot(data = Gastruloid_data_all_channels, aes(x=x.Mean)) + 
  geom_histogram(binwidth=200) 
 # coord_cartesian(xlim = c(0, 12000))+
#  scale_x_continuous(breaks = seq(0, 12000, 500))





Gastruloid_data_all_channelsfilt <- subset(Gastruloid_data_all_channels, x.Mean>400)
Gastruloid_data_all_channelsfilt <- subset(Gastruloid_data_all_channelsfilt, x.Mean<17000)

ggplot(data = Gastruloid_data_all_channelsfilt, aes(x=x.Mean)) + 
  geom_histogram(binwidth=20) +
  coord_cartesian(xlim = c(0, 12000))+
  scale_x_continuous(breaks = seq(0, 12000, 500))

#use this line to do the actual filtering
Gastruloid_data_all_channels <- subset(Gastruloid_data_all_channels, x.Mean>400)
Gastruloid_data_all_channels <- subset(Gastruloid_data_all_channels, x.Mean<7000)


# Filter for nuclei which are too small to be real
ggplot(data = Gastruloid_data_all_channelsfilt, aes(x=Vol..unit.)) + 
  geom_histogram(binwidth=20) 




#flip z axis data
Gastruloid_data_all_channels$Z_centre<-max(Gastruloid_data_all_channels$Z_centre)+10-Gastruloid_data_all_channels$Z_centre



#Scatter plot

# Here you can make a scatter plot of the integrated density of channel 2 against channel 3 to look for co-expressing cells for example
ggplot(data = Gastruloid_data_all_channels, mapping = aes(x = (y.IntDen), y = (x.IntDen))) +
  geom_point(size = 0.1, col="DarkOrchid")+
  scale_colour_viridis_d()+
  scale_fill_viridis_d()+
  labs(title = "Co-expression scatter plot of bracy vs foxa2", y = "IntDen of FOXA2 protein level", x = "IntDen of T protein level")  +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=0.5), 
        axis.title.y = element_text(size=12, vjust = 0.5), 
        strip.text=element_text(size=10, angle=45),
        plot.title = element_text(hjust = 0.5, size = 15)) +
  theme(legend.position = "right")






# Spatial plots - here you need to subset your dataframe so you are only working with a single colony/embryo/gastruloid e.g. Gastruloid_data_all_channels_Gastruloid_2 <- subset(Gastruloid_data_all_channels, Gastruloid_number_and_stage %like% "Gastruloid 2")


# Look into plotting the cells back into the embryo - This will make a 2D plot of all the objects in their respective position colour coded by the integrated density of channel 3
ggplot(Gastruloid_data_all_channels, aes(x = -X_centre, y = -Y_centre, colour = z.IntDen )) +
  geom_point(alpha = 10, size = 2)   +
  facet_wrap(~Gastruloid_number_and_stage, nrow = 10, ncol =6) +
  scale_colour_viridis(option = 'F') +
  theme_minimal() +
  coord_fixed()



#####################here for 2d plotly














# Plotly for 3D graphs - for 3D graphs - make sure you are working with the subsetted dataframe with only one object

plot_ly(data = Gastruloid_data_all_channels, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~y.Mean, marker = list(size = 5),
        opacity = 0.3) %>%
  layout(title = 'Bracy expression', plot_bgcolor = "#e5ecf6")


plot_ly(data = Gastruloid_data_all_channels, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~t.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'FOXA2 expression', plot_bgcolor = "#e5ecf6")


plot_ly(data = Gastruloid_data_all_channels, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~z.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'SOX2 expression', plot_bgcolor = "#e5ecf6")

plot_ly(data = Gastruloid_data_all_channels, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~x.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'DAPI expression', plot_bgcolor = "#e5ecf6")

#subsetting only streak and eipbast (by sox2)
subset(Gastruloid_data_all_channels, z.Mean >3000  |  y.Mean>5000 )->epi

plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~y.Mean, marker = list(size = 5), opacity = 0.9)
plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~t.Mean, marker = list(size = 5), opacity = 0.9)

plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~y.Mean, marker = list(size = 5),
        opacity = 0.9) %>%
  layout(title = 'Bracy expression', plot_bgcolor = "#e5ecf6")


plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~t.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'FOXA2 expression', plot_bgcolor = "#e5ecf6")


plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~z.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'SOX2 expression', plot_bgcolor = "#e5ecf6")

plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", mode="markers", color= ~x.Mean, marker = list(size = 5), 
        opacity = 0.9)%>%
  layout(title = 'DAPI expression', plot_bgcolor = "#e5ecf6")


hist(Gastruloid_data_all_channels$y.Mean, breaks=200)
axis(side=1, at=seq(0,30000,500))

hist(Gastruloid_data_all_channels$t.Mean, breaks=200)
hist(Gastruloid_data_all_channels$z.Mean, breaks=200)
axis(side=1, at=seq(0,30000,100))




# You could also give your objects a +/- rating for gene expression based on the range of values of expression found by looking at a histogram e.g. see example below from an embryo script

brc=950
fox=3250
sox=1075



Gastruloid_data_all_channels <- transform(Gastruloid_data_all_channels, 
                                          Brachy = ifelse(y.Mean>5000,
                                                          (ifelse(y.Mean>((max(Gastruloid_data_all_channels$y.Mean)+5000)*0.5), "Brachy_high", "Brachy_Low")),
                                                         "Brachy_no"))



Gastruloid_data_all_channels <- transform(Gastruloid_data_all_channels, 
                                          FOXA2 = ifelse(t.Mean>2000,
                                                         (ifelse(t.Mean>((max(Gastruloid_data_all_channels$t.Mean)+4500)*0.5), "FOXA2_high", "FOXA2_Low")),
                                                           "FOXA2_no"))


Gastruloid_data_all_channels <- transform(Gastruloid_data_all_channels, 
                                          simBrachy = ifelse(y.Mean>brc,"Brachy+","Brachy-"))



Gastruloid_data_all_channels <- transform(Gastruloid_data_all_channels, 
                                          simFOXA2 = ifelse(t.Mean>fox, "FOXA2+","FOXA2-"))











Gastruloid_data_all_channelssox<-Gastruloid_data_all_channels
Gastruloid_data_all_channelssim<-Gastruloid_data_all_channels

Gastruloid_data_all_channelssox <- transform(Gastruloid_data_all_channels, 
                                          SOX2 = ifelse(z.Mean>3000,
                                                         (ifelse(t.Mean>(max(Gastruloid_data_all_channels$t.Mean)*.5), "sox_high", "sox_Low")),
                                                         "sox_no"))


#visualise data based on each channel's value as a histogram


hist(Gastruloid_data_all_channels$y.Mean, breaks=200)
hist(Gastruloid_data_all_channels$t.Mean, breaks=200)
hist(Gastruloid_data_all_channels$z.Mean, breaks=200)
quantile(Gastruloid_data_all_channels$y.Mean, prob=c(.25,.5,.75))
quantile(Gastruloid_data_all_channels$t.Mean, prob=c(.25,.5,.75))
quantile(Gastruloid_data_all_channels$z.Mean, prob=c(.25,.5,.75,.95))


#label each cell by its +/- state of the protein in channels Y and T

Gastruloid_data_all_channels$Gene_expression <- paste(Gastruloid_data_all_channels$Brachy, Gastruloid_data_all_channels$FOXA2,sep = "_")
Gastruloid_data_all_channelssim$Gene_expression <- paste(Gastruloid_data_all_channels$simBrachy, Gastruloid_data_all_channels$simFOXA2,sep = "_")
subset(Gastruloid_data_all_channelssim, z.Mean >sox  |  y.Mean>brc   )->Gastruloid_data_all_channelssim

Gastruloid_data_all_channelssox$Gene_expression <- paste(Gastruloid_data_all_channelssox$SOX2,sep = "_")
subset(Gastruloid_data_all_channels, z.Mean >sox  |  y.Mean>brc   )->epi


# in silico recapitulation of object with poitn colouation based on  Channels Y (brachy) and T (FOXA2) +/- state
plot_ly(data = epi, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", 
        mode="markers", color= ~Gene_expression, size = ~Vol..unit., marker=list(sizemode = 'diameter'), 
        opacity = 0.5, sizes = c(7,15)) %>% 
  layout(xaxis = list(autorange = "reversed"))%>% 
  layout(title = '
         +/-ve state of cells within object',
         scene = list(xaxis = list(title = 'X-axis)',
                                   range = c(-20, (max(epi$X_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwidth = 2),
                      yaxis = list(title = 'Y-axis',
                                   range = c(0, (max(epi$y_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwith = 2),
                      zaxis = list(title = 'Z-axis',  range = c(-20, (max(epi$z_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   
                                   
                                   gridwith = 2)))


# Cell counting object based on Channels Y (brachy) and T (FOXA2)
subset(epi, ( Brachy=="Brachy_no")  &   FOXA2=="FOXA2_Low") ->lowfox  
subset(epi, ( Brachy=="Brachy_no")  &   FOXA2=="FOXA2_high")->highfox

subset(epi,( ( Brachy=="Brachy_high") ))->HighT
subset(epi,( ( Brachy=="Brachy_Low") ))->LowT

subset(epi, ( Brachy=="Brachy_Low")  &   FOXA2=="FOXA2_Low") ->lowfoxlowt
subset(epi, ( Brachy=="Brachy_Low")  &   FOXA2=="FOXA2_high") ->highfoxlowt
subset(epi, ( Brachy=="Brachy_high")  &   FOXA2=="FOXA2_Low") ->lowfoxhight
subset(epi, ( Brachy=="Brachy_high")  &   FOXA2=="FOXA2_high") ->highfoxhight

hist(epi$z.Mean, breaks=200)
subset(Gastruloid_data_all_channels, z.Mean >quantile(Gastruloid_data_all_channels$z.Mean, prob=c(.75)) )->sox2


as.data.frame(c(nrow(HighT), nrow(LowT), nrow(highfox), nrow(lowfox), nrow(lowfoxlowt),nrow(highfoxlowt),nrow(lowfoxhight),nrow(highfoxhight), nrow(sox2)))->cellcounts
rownames(cellcounts)<- c("High T","low T","High fox no T", "low fox no T", "low fox low t", "high fox low t", "low fox high t", "high fox high t", "SOX2")
colnames(cellcounts)<-"no. of cells"
cellcounts







#calculation of number of singel/ double positive cells from channels T and Y


subset(Gastruloid_data_all_channelssim,( ( simBrachy=="Brachy+") ))->Bracypos
subset(Gastruloid_data_all_channelssim,( ( simFOXA2=="FOXA2+") ))->Foxpos
subset(Gastruloid_data_all_channelssim, ( simBrachy=="Brachy+")  &   simFOXA2=="FOXA2+") ->dubpos

hist(epi$z.Mean, breaks=200)
subset(Gastruloid_data_all_channels, z.Mean >1400 )->sox2


as.data.frame(c(nrow(Bracypos), nrow(Foxpos), nrow(dubpos)))->simcellcounts
rownames(simcellcounts)<- c("TBXT","FOXA2","Co-expressing")
colnames(simcellcounts)<-"no. of cells"
simcellcounts





# 3d plot with colouration data based on channel Z thresholding to check epiblast
plot_ly(data = Gastruloid_data_all_channelssox, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", 
        mode="markers", color= ~Gene_expression, size = ~Vol..unit., marker=list(sizemode = 'diameter'), 
        opacity = 1, sizes = c(7,10)) %>% 
  layout(xaxis = list(autorange = "reversed"))%>% 
  layout(title = 'Whole mount embryo',
         scene = list(xaxis = list(title = 'X-axis)',
                                   range = c(-20, (max(Gastruloid_data_all_channels$X_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwidth = 2),
                      yaxis = list(title = 'Y-axis',
                                   range = c(0, (max(Gastruloid_data_all_channels$y_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwith = 2),
                      zaxis = list(title = 'Z-axis',  range = c(-20, (max(Gastruloid_data_all_channels$z_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   
                                   
                                   gridwith = 2)))
Gastruloid_data_all_channelssim$X_centre<-( max(Gastruloid_data_all_channelssim$X_centre)+10-Gastruloid_data_all_channelssim$X_centre)
# Plotly plot simplified fox and t
plot_ly(data = Gastruloid_data_all_channelssim, x= ~X_centre, y= ~Y_centre, z= ~Z_centre, type="scatter3d", 
        mode="markers", color= ~Gene_expression, size = ~Vol..unit., marker=list(sizemode = 'diameter'), 
        colors = c(NULL, "forestgreen", "firebrick3", "Blue"),
        opacity = 1, sizes = c(7,10)) %>% 
  layout(xaxis = list(autorange = "reversed"))%>% 
  layout(title = 'Whole mount embryo',
         scene = list(xaxis = list(title = 'X-axis',
                                   range = c(-20, (max(Gastruloid_data_all_channelssim$X_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwidth = 2),
                      yaxis = list(title = 'Y-axis',
                                   range = c(0, (max(Gastruloid_data_all_channelssim$y_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   gridwith = 2),
                      zaxis = list(title = 'Z-axis',  range = c(-20, (max(Gastruloid_data_all_channelssim$z_centre)+100)),
                                   zerolinewidth = 1,
                                   ticklen = 5,
                                   
                                   
                                   gridwith = 2)))












# 2d-plot (sans z axsis) using 2 variables

plot_ly(data = Gastruloid_data_all_channelssim, x = ~Z_centre, y = ~X_centre, 
        type = "scatter", mode = "markers", color = ~Gene_expression, 
        size = ~Vol..unit., marker = list(size = 7, sizemin = 0, sizemax = 0), 
        colors = c("Grey", "#00FF00", "firebrick3", "Blue"), opacity = 1, 
        sizes = c(7,10)) %>% 
  
  layout(showlegend = FALSE)%>% 
  layout(title = '',
         xaxis = list(title = 'X-axis (A.U) ', range = c(0, 1560),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2),
         yaxis = list(title = 'Y-axis (A.U)', range = c(0, 1270),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2))






#2dplot of your object from above: points colored by channel Y intensity


plot_ly(data = Gastruloid_data_all_channelssim, x = ~Y_centre, y = ~X_centre, 
        type = "scatter", mode = "markers", color = ~y.Mean, colors = "Blues",
        marker = list(size = 15, sizemin = 0, sizemax = 0, showscale = FALSE, 
                      hovertemplate = paste("X: %{x}<br>",
                                            "Y: %{y}<br>",
                                            "z.Mean: %{marker.color}<br>")), 
        opacity = 1) %>% 
  layout(xaxis = list(autorange = "reversed"))%>% 
  layout(showlegend = FALSE)%>% 
  layout(title = '',
         xaxis = list(title = 'X-axis (A.U) ', range = c(0, 2000),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2),
         yaxis = list(title = 'Y-axis (A.U)', range = c(0, 1270),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2))

#subset your data based on Z.mean value, and replot the object with colorization based on Y.mean value.Hover over values displayed are
#X and Y co-ordinates
Gastruloid_data_all_channels_5kzmean <- subset(Gastruloid_data_all_channelssim, z.Mean > 5500)
plot_ly(data = Gastruloid_data_all_channels_5kzmean, x = ~Y_centre, y = ~X_centre, 
        type = "scatter", mode = "markers", color = ~y.Mean,colors = "Greens", 
        marker = list(size = 15, sizemin = 0, sizemax = 0, showscale = TRUE, 
                     
                      hovertemplate = paste("X: %{x}<br>",
                                            "Y: %{y}<br>",
                                            "z.Mean: %{marker.color}<br>")), 
        opacity = 1) %>% 

  layout(showlegend = F)%>%
  
  layout(title = paste(idY,"expression in",idZ, "positve cell", sep=" "),
         xaxis = list(title = 'X-axis (A.U) ', range = c(1400, 0),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2),
         yaxis = list(title = 'Y-axis (A.U)', range = c(0, 1400),
                      zeroline = FALSE, ticklen = 5, gridwidth = 2))


