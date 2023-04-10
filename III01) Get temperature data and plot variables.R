##### I) Download and extract sea surface temperature (SST) data (adapted from code by Pieter Arnold) #####

library(curl)
library(ncdf4)
library(RCurl)
library(R.utils)
library(httr)
library(purrr)
library(dplyr)
library(raster)
library(weathermetrics)
library(tidyr)
library(ggplot2)

setwd("D:/SST data/2018/2018-07")                               # INPUT: change working directory as desired
# Define year and month of dataset downloading and extracting SSTs
year   <- 2018        
month  <- "08"                                                                       ## INPUT

setdays <- function(){
  if (month == "01" || month == "03" || month == "05" || month == "07"
      || month == "08" || month == "10" || month == "12") {
    day <- seq(1, 31, 1)   # JANUARY, MARCH, MAY, JULY, AUGUST, OCTOBER, DECEMBER
  } else if (month == "04" || month == "06" || month == "09" || month == "11") {
    day <- seq(1, 30, 1)   # APRIL, JUNE, SEPTEMBER, NOVEMBER
  } else if (month == "02" && year %% 4 == 0)  {
    day <- seq(1, 29, 1) 
  } else if (month == "02") {
    day <- seq(1, 28, 1)
  } 
  day <- formatC(day, width = 2, format = "d", flag = "0")
  print(day)
  
  # Workaround for leap years
  if (year %% 4 == 0) {
    url_date <<- seq(1, 366, 1)
  } else url_date <<- seq(1, 365, 1)
  print(max(url_date))
  
  days <- day[1:31]
  days <<- days[!is.na(days)]
  days
}

setdays()

# Set the name of the final URL segment and the day with appropriate zeros in front
# This also sets which files and the total number of files to download
# Specify days of the month
# Define number of days in the month based on the month defined above (if less than the whole month)
# NB: sometimes the server denies access to certain files for seemingly no reason
#     if this happens, move on and download the next batch without the problem file
#     then go back to the problem file individually later


# Check ftp server file directory and set this value for where 
# the first of the files for the month to download starts
# For example: seq(1, (1+length(days)-1), 1) where days <- day[1:31] will run through all of January
# and seq(61, (61+length(days)-1), 1) where days <- day[1:30] will run through all of March in a leap year
# i.e. the folders in the ftp server use Julian date rather than calendar date, but we write the files as 
# the calendar date but have to set the download by Julian date.
# ftp://ftp.nodc.noaa.gov/pub/data.nodc/ghrsst/L4/GLOB/JPL_OUROCEAN/G1SST/

## This is to know the day number of any chosen date (INPUT: all dates below)
chosen_date <- as.Date("2018-08-01") ##INPUT
dates <- seq(as.Date("2018-01-01"), as.Date("2018-12-31"), by="days") #change year only
day_numbers <- as.data.frame(url_date)
day_numbers$dates <- dates
day_numbers[day_numbers$dates == chosen_date, 1]
inday <- day_numbers[day_numbers$dates == chosen_date, 1]
##

url_name <- seq(inday, (inday+length(days)-1), 1)                                                      ## INPUT: change initial day
url_name <- formatC(url_name, width = 3, format = "d", flag = "0")                              
url_name

day_name <- formatC(url_name, width = 2, format = "d", flag = "0")
day_name

# do the number of days and number of files to download match? If not, check url_name and day_name above
length(days) : length(url_name) 

Sys.time()
# Set URLs, download files from set URLs and write to working directory
for(i in 1:length(url_name)){
  url <- paste0("ftp://ftp.nodc.noaa.gov/pub/data.nodc/ghrsst/L4/GLOB/JPL_OUROCEAN/G1SST/", 
                year, "/", url_name[i], "/")
  curl = getCurlHandle()
  curlSetOpt(.opts = list(forbid.reuse = 1), curl = curl) 
  filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE, crlf = TRUE, ssl.verifypeer = FALSE, 
                     curl = curl)
  filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
  filenamesNC = filenames[1] # subset only the netcdf files
  download.file(url = filenamesNC, mode = "wb", method = "wininet",
                destfile = paste0(year, month, days[i], 
                                  "-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2"))
  rm(curl)
  gc()	
  Sys.sleep(3)
}

# Unzip the .bz files
bzfolder <- "D:/SST data/2018/2018-08"                        # INPUT: change destination of bzfolder as desired
bzfiles  <- list.files(bzfolder, pattern = "*.bz2", full.names = T)
n_zip    <- length(bzfiles)
for (i in 1:n_zip) {
  bunzip2(bzfiles[i], destname = gsub("[.]bz2$", "", bzfiles[i]), overwrite=FALSE, remove=FALSE)
}


###### II) Create variables ######

# Trim nc files to south east Australia
library(raster)
library(weathermetrics)

for(i in 2011:2018){
  year <- i
  for(j in 1:12){
    month <- j
    month <- formatC(month, width = 2, format = "d", flag = "0")
    
    setwd(paste0("~/Small ncFiles/", year, "/", year, "-", month))
    
    ncfolder <- paste0("D:/SST data/", year, "/", year, "-", month)     
    files <- list.files(ncfolder, pattern = "*.nc$", full.names = T)
    e <- extent(134.076, 156, -44.25, -35.273)
    
    for(i in 1:length(files)){
      aux1 <-  raster(files[i])
      aux1 <- crop(aux1, e)
      values(aux1) <- kelvin.to.celsius(values(aux1))
      start <- sum(regexpr(paste0(year, "-", month, "/"), files[i])[1] + 8)                             
      stop <- start + 7
      date <- substr(files[i], start, stop)
      writeRaster(aux1, paste0("SSTAus_", date, ".nc"), format = "CDF", overwrite = TRUE)
    }
  }
}

# Stack all rasters
library(raster)
library(tidyr)
library(RColorBrewer)
library(grDevices)
library(mapplots)

stack <- stack()
for (year in 2010:2019){
  if(year == 2010){
    ini <- 7
    end <- 12
  } else if (year == 2019) {
    ini <- 1
    end <- 8
  } else {
    ini <- 1
    end <- 12
  }
  for (j in ini:end){
    month <- j
    month <- formatC(month, width = 2, format = "d", flag = "0")
    days <- formatC(1:31, width = 2, format = "d", flag = "0")
    ncfolder <- paste("~/Small ncFiles/", year, "/", year, "-", month, sep = "")
    files <- list.files(ncfolder, pattern = "*.nc$", full.names = T)
    map <- stack(files[1:length(files)])
    #val <- values(map) #This lines check if there are any values less than 0 (an error in the data that sometimes occur)
    #val[val < 0] <- NA #and replaces them with NAs, before adding the layers to the main stack
    #values(map) <- val
    #rm(val)
    stack <- addLayer(stack, map)
    rm(map)
    print(paste("year", year, "month", month, "done"))
  }
}

#writeRaster(stack, "C~/SSTAU_stack.nc", format = "CDF")

## This datemap let me sample the stack by the desire years, months, or days by then creating an index
dates <- seq(as.Date("2010-07-01"), as.Date("2019-08-27"), by="days")
dayn <- 1:length(dates)
datemap <- data.frame(dates, dayn)
datemap <- separate(datemap, col = dates, into = c("year", "month", "day"), sep = "-", remove = FALSE)
#names(stack) <- as.character(datemap$dates)
#

## I) The grand mean (2011-2018) (GM)
gm <- mean(stack[[185:3106]], na.rm = TRUE)
plot(gm)
writeRaster(gm, paste(src, "/Data/gm.nc", sep = ""), format = "CDF", overwrite = TRUE)
min(values(gm), na.rm = TRUE)
max(values(gm), na.rm = TRUE)

## II) FLUCTUATIONS
### BIO7 = Temperature Annual Range (BIO5-BIO6) -> mean annual range (mar)
ar <- stack()
for(i in 2011:2018){
  year <- i
  index <- datemap$dayn[datemap$year == year]
  aux <- max(stack[[index]]) - min(stack[[index]])
  ar <- addLayer(ar, aux)
}
mar <- mean(ar, na.rm = TRUE)
plot(mar)
writeRaster(mar, paste(src, "/Data/", "mar.nc", sep = ""), format = "CDF", overwrite = TRUE)

### III) NOISE
## Noise colour (B)
#In the cluster
detach("package:tidyr", unload = TRUE)
n <- 0
ext.betas <- function(timeseries){
  n <<- n+1
  timeseries <- timeseries[!is.na(timeseries)]
  if(length(timeseries) == 0){
    print(c(n, NA))
    return(NA)
  } else {
    spec <- spectrum(timeseries, plot = FALSE)
    x <- spec$freq
    y <- spec$spec
    b <- NA
    tryCatch({ 
      bfit <- nls(y ~ 1/x^b, start = list(b = 1), control = list(maxiter = 800))
      b <- summary(bfit)$coefficients[1]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if(is.numeric(b)){ 
      print(c(n, b))
      return(b)
    } else {
      print(c(n, NA))
      return(NA)
    }
  }                   
}
index <- datemap$dayn[datemap$year != 2010 & datemap$year != 2019]
noise <- stack[[index]]
noisemap <- calc(noise, fun = ext.betas)
noisemap
writeRaster(noisemap, "/home/cristobg/tn20/noisemap.nc", format = "CDF", overwrite = TRUE)
#

#In local pc
file <- "~/Data/noisemap.nc"
noisemap <- stack(file)
noiseval <- values(noisemap) 
noiseval[noiseval > 1] <- NA  #a few values go up to 2 (probably errors, mess up plot)      
values(noisemap) <- noiseval
plot(noisemap)
writeRaster(noisemap, "~/Data/noisemap.nc", format = "CDF", overwrite = TRUE)
