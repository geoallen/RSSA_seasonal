#!/usr/bin/env Rscript
################################################################################
# ncdf_to_shapfile_attribute_transfer.R
################################################################################
# George Allen, May 2019

################################################################################
# Description:

# Maps discharge from Lin et al. to flowlines from MERIT Hydro:

# NOTE: to make nice figures, sort MERIT flowlines by stream order in ArcMap 
# or using R code

################################################################################
library(ncdf4)
library(foreign)

################################################################################


# Define paths: 
wd = "E:/research/2020_01_01_RSSA_seasonal/git/RSSA_seasonal"
inDir = paste0(wd, "/in")
outDir = paste0(wd, "/out")

ncDirIn = paste0(inDir, "/GRADES_seasonal")
mhDirIn = paste0(inDir, "/MERIT_Hydro/riv")

mhDirOut = paste0(outDir, "/MERIT_RAPID")

# list files:
ncPaths = list.files(ncDirIn, pattern="nc", full.names=T)
mHpaths = list.files(mhDirIn, pattern="dbf", full.names=T)
mHpathsOut = sub(mhDirIn, mhDirOut, mHpaths)

# for each netCDF file: 
for (i in 1:length(ncPaths)){
  
  print(i)
  
  # read in netCDF file:
  ncIn = nc_open(ncPaths[i])
  # print(ncIn)
  # sum(ncIn$var$Qmean$varsize[1]*12*13)
  # ncvar_get(ncIn,"Qmean",verbose=T)
  
  # read in river shapefile dbf file:
  dbf = foreign::read.dbf(mHpaths[i])
  # dim(dbf)
  # head(dbf)
  # tail(dbf)
  
  # match netCDF and DBF:
  # o = order(dbf$COMID)
  # j = match(c(1:nrow(dbf)), o)
  # 
  # concatenate columns to shapefile attribute table:
  dbfOut = data.frame(dbf, 
                      Qmn = rowMeans(ncvar_get(ncIn,"Qmean")[,c(1:12)]),
                      Qstd = rowMeans(ncvar_get(ncIn,"Qstd")[,c(1:12)]),
                      # Qjan = ncvar_get(ncIn,"Q50")[,1],
                      # Qfeb = ncvar_get(ncIn,"Q50")[,2],
                      # Qmar = ncvar_get(ncIn,"Q50")[,3],
                      # Qapr = ncvar_get(ncIn,"Q50")[,4],
                      # Qmay = ncvar_get(ncIn,"Q50")[,5],
                      # Qjun = ncvar_get(ncIn,"Q50")[,6],
                      # Qjul = ncvar_get(ncIn,"Q50")[,7],
                      # Qaug = ncvar_get(ncIn,"Q50")[,8],
                      # Qsep = ncvar_get(ncIn,"Q50")[,9],
                      # Qoct = ncvar_get(ncIn,"Q50")[,10],
                      # Qnov = ncvar_get(ncIn,"Q50")[,11],
                      # Qdec = ncvar_get(ncIn,"Q50")[,12],
                      Q000 = rowMeans(ncvar_get(ncIn,"Q0")[,c(1:12)]),
                      Q010 = rowMeans(ncvar_get(ncIn,"Q10")[,c(1:12)]),
                      Q020 = rowMeans(ncvar_get(ncIn,"Q20")[,c(1:12)]),
                      Q030 = rowMeans(ncvar_get(ncIn,"Q30")[,c(1:12)]),
                      Q040 = rowMeans(ncvar_get(ncIn,"Q40")[,c(1:12)]),
                      Q050 = rowMeans(ncvar_get(ncIn,"Q50")[,c(1:12)]),
                      Q060 = rowMeans(ncvar_get(ncIn,"Q60")[,c(1:12)]),
                      Q070 = rowMeans(ncvar_get(ncIn,"Q70")[,c(1:12)]),
                      Q080 = rowMeans(ncvar_get(ncIn,"Q80")[,c(1:12)]),
                      Q090 = rowMeans(ncvar_get(ncIn,"Q90")[,c(1:12)]),
                      Q100 = rowMeans(ncvar_get(ncIn,"Q100")[,c(1:12)]))
  
  
  # write out shapefile attribute table:
  foreign::write.dbf(dbfOut, mHpathsOut[i])

}

system("rundll32 user32.dll,MessageBeep -1")



# NOTE: to make nice figures, sort MERIT flowlines by stream order in ArcMap 

