################################################################################
# DHG_analysis.R
################################################################################
# George H. Allen, Dec. 2019

# Description:
# Reads in joined flow occurrence and GRADES Q data and runs regressions for 
# global hydraulic geometry relationships.

# 

################################################################################
# load libraies:
################################################################################

packages = c("foreign", "zyp", "ggplot2")

ipak = function(pkg){
  new.pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){ install.packages(new.pkg, dependencies = T) }
  sapply(pkg, require, character.only=T)
}

ipak(packages)


################################################################################
# Define paths:
################################################################################

wDir = "E:/research/2020_01_01_RSSA_seasonal/git/RSSA_seasonal"

# input directory paths:
inGRADESdir = paste0(wDir, "/out/MERIT_RAPID/")
centroidDir = paste0(wDir, "/in/hydroBASINS/MERIT_riv_centroid_hyBas_spJoin/")
inXsecDir = paste0(wDir, "/in/GRWL/GRWL_Xsections_MERIT_spJoin/")
inHyBasDir = paste0(wDir, "/in/hydroBASINS/hydroBASINS/")



# ouptut directory paths:
outGRADESdir = paste0(wDir, "/out/MERIT_RAPID_HYBAS/")
outXsecDir = paste0(wDir, "/out/GRWL/GRWL_Xsections_DHG/")
outHyBasDir = paste0(wDir, "/out/hydroBASINS/")
outGRADES_DHGdir = paste0(wDir, "/out/MERIT_RAPID_HYBAS_DHG/")
outGRADESmodWdir = paste0(wDir, "/out/MERIT_RAPID_HYBAS_modW/")


# specify hardcoded parameters:
minPfaf = 1 # set max basin size to run regression on (min pfaf basin)
maxPfaf = 10 # set min basin size to run regression on (max pfaf basin)
minLength = 100 # N width measurements
minUniqA = 5 # N unique drainage area measurements
minUniqQ = 100 # N unique discharge measurements 

qIntervals = c("Q000","Q010","Q020","Q030","Q040","Q050","Q060","Q070","Q080","Q090","Q100")
wIntervals = c("w100","w090","w080","w070","w060","w050","w040","w030","w020","w010","w000")



################################################################################
# Functions:
################################################################################

# extract the p-value from linear regression:
lmp = function(modelobject){
  if (class(modelobject) != "lm"){ stop("Not an object of class 'lm' ") }
  f = summary(modelobject)$fstatistic
  p = pf(f[1], f[2], f[3], lower.tail=F)
  attributes(p) = NULL
  return(p)
}


################################################################################
# set up paths:
################################################################################

# copy any shapfiles:
inPaths = list.files(inGRADESdir, full.names=T)
if (length(grep("pfaf_09", inPaths))>0){ inPaths = inPaths[-grep("pfaf_09", inPaths)] }
file.copy(inPaths, sub(inGRADESdir, outGRADESdir, inPaths), overwrite=F)

inPaths = list.files(inXsecDir, full.names=T)
if (length(grep("pfaf_09", inPaths))>0){ inPaths = inPaths[-grep("pfaf_09", inPaths)] }
file.copy(inPaths, sub(inXsecDir, outXsecDir, inPaths), overwrite=F)

# get in paths:
inGRADESpaths_raw = list.files(inGRADESdir, ".dbf", full.names=T)
inCentroidPaths_raw = list.files(centroidDir, ".dbf", full.names=T)
inXsecPaths_raw = list.files(inXsecDir, ".dbf", full.names=T)
inHyBasPath = list.files(inHyBasDir , ".dbf", full.names=T)

# rm Greenland bc innaccurate GRADES Q data:
inGRADESpaths = inGRADESpaths_raw[-grep("pfaf_09", inGRADESpaths_raw)] 
inCentroidPaths = inCentroidPaths_raw[-grep("pfaf_09", inCentroidPaths_raw)] 
inXsecPaths = inXsecPaths_raw[-grep("pfaf_09", inXsecPaths_raw)] 

# output file paths:
gradesPfafOutPaths = sub(inGRADESdir, outGRADESdir, inGRADESpaths)
xsecPfafOutPaths = sub(inXsecDir, outXsecDir, inXsecPaths)
hybasOutPath = sub(inHyBasDir, outHyBasDir, inHyBasPath)
hybasOutPaths = paste0(gsub(".dbf", "", hybasOutPath), "_", qIntervals, ".dbf")
outGRADES_DHGpaths = sub(inGRADESdir, outGRADES_DHGdir, inGRADESpaths)
outGRADESmodWpaths = sub(inGRADESdir, outGRADESmodWdir, inGRADESpaths)




################################################################################
# join pfafstetter codes to GRADES riv and xSection shapefiles:
################################################################################


# # for each GRADES shapefile: 
# for (h in 1:length(inGRADESpaths)){
#   print(h)
#   print("reading in dbfs...")
# 
#   grades = foreign::read.dbf(inGRADESpaths[h])
#   centroid = foreign::read.dbf(inCentroidPaths[h])
#   xsec = foreign::read.dbf(inXsecPaths[h])
#   
#   # hydroBASIN pfafstetter code columns:
#   pfafColInd = grep("PFAF_", names(centroid))
#   
#   # join pfafstetter codes to GRADES riv:
#   mInd = match(grades$COMID, centroid$COMID)
#   hybasPfafTab = centroid[mInd, pfafColInd]
#   gradesPfaf = data.frame(grades, hybasPfafTab)
#   
#   # join pfafstetter codes to cross sections 
#   mInd = match(xsec$COMID, centroid$COMID)
#   hybasPfafTab = centroid[mInd, pfafColInd]
#   xsecPfaf = data.frame(xsec, hybasPfafTab)
#   
#   print("writing out...")
#   foreign::write.dbf(gradesPfaf, gradesPfafOutPaths[h])
#   foreign::write.dbf(xsecPfaf, xsecPfafOutPaths[h])
#   
# }




################################################################################
# Pfafstetter optimization process:
################################################################################


for (h in 1:length(qIntervals)){
  print(paste("h =", h))

  oTabNames = c("nSeg", "nUniqQ", "a_coef", "b_exp", "R2", "pVal")
  oTab = as.data.frame(array(NA, c(nrow(hybas), length(oTabNames))))
  names(oTab) = oTabNames
  
  
  # loop through each river shapefile: 
  for (i in 1:length(xsecPfafOutPaths)){
    print(paste("i =", i))
    
    # read in river file:
    xSec = read.dbf(xsecPfafOutPaths[i])
    
    # get continental integer vector:
    # contInt = xSec$PFAF_1[1]
    # uniqContInt = unique(contInt)
    
    # identify reaches to remove from analysis, and set data to NA: 
    qCol = grep(qIntervals[h], names(xSec))
    wCols = grep(wIntervals[h], names(xSec))
    wCol = wCols[grep("flag", names(xSec)[wCols], invert=T)]
    flagCol = wCols[grep("flag", names(xSec)[wCols])]
    
    
    w = xSec[,wCol] * xSec$width_m
    Q = xSec[,qCol]
    flag = xSec[,flagCol]
    
    
    exclude = 
      w < 30 |
      Q <= 0 |
      is.na(w) | 
      is.na(Q) | 
      is.infinite(w) |
      is.infinite(Q) |
      flag == 1 
    
    w[exclude] = NA
    Q[exclude] = NA
    
    # take log of width and discharge:
    logW = log(w)
    logQ = log(Q)
    
    # isolate pfafstetter codes: 
    pfafCols = grep("PFAF", names(xSec))
    pfafTab = xSec[, pfafCols]
    
    # set up pfafstetter data tables: 
    pfafMatchTab = 
      nTab =  
      aTab = 
      bTab = 
      R2Tab = 
      pTab = array(NA, dim(pfafTab))
    
    # for each pfaf level, find unique codes:
    uniqPfafCodesWithZeros = apply(pfafTab, 2, unique)
    uniqPfafCodes = lapply(uniqPfafCodesWithZeros, function(x) {x[x!=0]})
    
    print("Calculating statistics for each Pfafstetter level...")
    for (j in minPfaf:maxPfaf){
      print(paste("pfaf level:", j))
      pfafMatchTab[,j] = match(pfafTab[,j], uniqPfafCodes[[j]])
      
      # for each unique pfaf code: 
      for (k in 1:length(uniqPfafCodes[[j]])){
        # identify basins that match the unique pfaf code: 
        keep = which(pfafMatchTab[,j] == k & !exclude)
        
        # if the basin has less than a minimum length of river or 
        # if the basin has less than a minimum unique discharges, skip:
        if (length(keep) < minLength | 
            length(unique(logQ[keep])) < minUniqQ){ next() }
        
        # run a linear regression in log space: 
        zlogW = logW[keep]
        zlogQ = logQ[keep]
        
        # thiel-sen:
        # must sample vector to save memory:
        if (length(zlogQ) > 1e4){
          int = round(seq(1, length(zlogQ), length.out=1e4))
          zlogW = zlogW[int]
          zlogQ = zlogQ[int]
        }
        sen = zyp.sen(zlogW ~ zlogQ)
        
        # least-squares:
        lreg = lm(zlogW ~ zlogQ)
        
              # occationally when the logQs are very close, b_exp is NA. Skip these cases: 
        if (is.na(exp(lreg$coefficients[[1]])) | is.na(lreg$coefficients[[2]])){ next() }
        
        # calculate model parameters:
        nSeg = length(keep)
        nUniqQ = length(unique(zlogQ))
        a_coef = exp(sen$coefficients[[1]])
        b_exp = sen$coefficients[[2]]
        R2 = cor(zlogW, zlogQ)^2
        pVal = lmp(lreg)
        
        # for least squares linear regression:
        # nSeg = length(keep)
        # nUniqQ = length(unique(logQ[keep]))
        # a_coef = exp(lreg$coefficients[[1]])
        # b_exp = lreg$coefficients[[2]]
        # R2 = summary(lreg)$r.squared
        # pVal = lmp(lreg)
        
        # match given pfaf basin to hybas shapefile attribute tab:
        hybasInd = which(hybasPfafTab[,j] == uniqPfafCodes[[j]][k])
        
        # define merge operation. If there is already an estimate 
        # for the given hydrobasin from a larger basin, select a method
        # e.g.: replace with the higher R2, take mean, weighted by  R2,
        # use smaller basin as long as it has an R2 > 0.5: 
        
        # R2 must be equal or greater than predefined larger basin R2:
        predefined = which(!is.na(oTab$R2[hybasInd]))
        if (length(predefined) > 0){
          gtR2 = oTab$R2[hybasInd][predefined] <= R2
          oTab[hybasInd[predefined[gtR2]], ] = mget(names(oTab))
        }else{
          oTab[hybasInd, ] = mget(names(oTab))
        }
        
        
        # plot:
        # x = Q[keep]
        # y = w[keep]
        # plot(x, y,
        #      log = 'xy',
        #      xlab = "Drainage Area (cms)",
        #      ylab = "GRWL Width (m)",
        #      main = "Basin-wide Width-Area Scaling")
        # xseq = seq(min(x), max(x), length.out=100)
        # yseq = exp(lreg$coefficients[[1]])*xseq^lreg$coefficients[[2]]
        # lines(xseq, yseq, col=4)
        # mtext(paste("j=", j, "   k=", k, "   a=", round(lreg$coefficients[[1]]),
        #             "   b=", round(lreg$coefficients[[2]],3),
        #             "   R2=", round(summary(lreg)$r.squared, 2)))
        
      }
    }
    paste0(print(length(which(!is.na(oTab$b_exp)))/nrow(oTab)), "% done")
  }
  
  
  
  # add new columns to hydrobasins attribute table:
  hybasOut = cbind(hybas, oTab)
  hybasOut[is.na(hybasOut)] = -9999
  #hybasOut[hybasOut==-9999] = NA
  
  write.dbf(dbf_mean, hybasOutPaths[h])

}
system("rundll32 user32.dll,MessageBeep -1")




hybasOutPaths = list.files(outHyBasDir, ".dbf", full.names=T)
dbf_list = as.list(NA)


# create a large list of each recurrence interval: 
for (i in 1:length(hybasOutPaths)){
  dbf = foreign::read.dbf(hybasOutPaths[i])
  dbf[dbf==-9999] = NA
  dbf_list[[i]] = dbf
}

# remove the 0% and 100% recurrence interval due to issues discussed in 
# Allen et al., 2020 and take the mean of each interval:
dbf_list = dbf_list[-c(1,length(dbf_list))]
dbf_mean = Reduce(`+`, dbf_list)/length(dbf_list)

summary(dbf_mean)


# join mean DHG parameters to GRADES riv data: 
for (i in 1:length(gradesPfafOutPaths)){
  print(i)
  grades = foreign::read.dbf(gradesPfafOutPaths[i])
  
  j = match(grades$PFAF_12, dbf_mean$PFAF_12)
  transferCols = dbf_mean[j, grep("nSeg|a_coef|b_exp|R2|pVal", names(dbf_mean))] 
  grades_out = data.frame(grades, transferCols)
  
  foreign::write.dbf(grades_out, outGRADES_DHGpaths[i])
  
}


# apply DHG equations to estimate widths at different Q recurrence intervals:
for (i in 1:length(outGRADES_DHGpaths)){
  print(i)
  grades = foreign::read.dbf(outGRADES_DHGpaths[i])
  
  modWtab = as.data.frame(array(NA, c(nrow(grades), length(qIntervals))))
  names(modWtab) = paste0("modW_", qIntervals)

  for (j in 1:ncol(modWtab)){
    modWtab[,j] = grades$a_coef * grades[, grep(qIntervals[j], names(grades))] ^ grades$b_exp
  }

  grades_out = data.frame(grades, modWtab)
  
  if (i == 1){ 
    grades_bigTab = grades_out
  }else{
    colMatchInd = match(names(grades_bigTab), names(grades_out))
    grades_bigTab = rbind(grades_bigTab, grades_out[,colMatchInd])
  }
  
  foreign::write.dbf(grades_out, outGRADESmodWpaths[i])
  
}

# i sorted these shapefiles by stream order with Arc sort tool and display by width

# column means:
colMeans(grades_bigTab, na.rm=T)

# sum up river area: 
sum(grades_bigTab$modW_Q010 * grades_bigTab$Length, na.rm=T)/1e6
# q010 river area = 4,742,570 km2
# q050 river area = 6,025,868 km2
# q090 river area = 7,978,474 km2












#### RECYCLING BIN:

# ##############
# 
# 
# # read in hybas attribute table and find pfafstetter codes:
# hybas = foreign::read.dbf(hybasInPath)
# pfafCols = grep("PFAF", names(hybas))
# pfafTab = hybas[, pfafCols]
# 
# # setup new output table:
# oTabNames = c("nSeg", "nUniqQ", "a_coef", "b_exp", "R2", "pVal")
# oTab = as.data.frame(array(NA, c(nrow(hybas), length(oTabNames))))
# names(oTab) = oTabNames
# 
# 
# 
# # loop through each continent:
# for (h in 4:length(xsecPfafOutPaths)){
# 
#   #####
#   print(paste("h =", h))
# 
#   # read in river file:
#   xSec = read.dbf(xsecPfafOutPaths[h])
#   xsec = xsec[, order(names(xsec))]
# 
#   QrecCols = grep("Q[[:digit:]]", names(xSec))
#   wCols_rev = grep("w[[:digit:]]", names(xSec))
#   wCols = rev(wCols_rev)
#   wOccCols = wCols[grep("flag", names(xSec)[wCols], invert=T)]
#   wFlagCols = wCols[grep("flag", names(xSec)[wCols])]
#   pfafCols = grep("PFAF", names(xSec))
# 
#   # multiply width occurence levels by cross section length:
#   xSec[, wOccCols] = xSec[, wOccCols] *  xSec$width_m
# 
#   # filter:
#   f = xSec$width_m > 90
#   #f = 1:nrow(xSec)
# 
#   # create smaller, tractable tables:
#   qTab = xSec[f, QrecCols]
#   wTab = xSec[f, wOccCols]
#   fTab = xSec[f, wFlagCols]
#   pTab = xSec[f, pfafCols]
# 
#   # set to NA measurements with values of zero or less and measuremets taken
#   # where river water is located at the end of their cross section segments:
#   rmBoo = wTab<=0 | qTab<=0 | fTab==1
#   qTab[rmBoo] = NA
#   wTab[rmBoo] = NA
# 
#   print("Calculating statistics for each recurrence invterval...")
#   # for each recurrence interval:
#   for (i in 1:ncol(qTab)){
# 
#     print(paste("i =", i))
# 
#     # this is funky and a results of merging old code into new framework:
#     # generate copies of oTab for each recurrence interval:
#     if (h != 1){
#       # get seperate oTab for each recurrence interval:
#       oTab = get(paste0("oTab_", names(qTab)[i]))
#     }
# 
#     # set up vars:
#     Q = qTab[,i]
#     w = wTab[,i]
# 
#     # if there are no positive values, skip to next recurrence interval:
#     boo = Q>0 & w>0
#     if (!(T %in% boo)){next}
# 
#     # take log of Q and w data:
#     logQ_raw = log(Q)
#     logW_raw = log(w)
# 
#     logQ = logQ_raw[!is.na(logQ_raw) | !is.na(logW_raw)]
#     logW = logW_raw[!is.na(logQ_raw) | !is.na(logW_raw)]
# 
#     # set up pfafstetter data tables:
#     pfafMatchTab = array(NA, dim(pTab))
# 
#     # for each pfaf level, find unique codes:
#     uniqPfafCodesWithZeros = apply(pTab, 2, unique)
#     uniqPfafCodes = lapply(uniqPfafCodesWithZeros, function(x) {x[x!=0]})
# 
#     print("Calculating statistics for each Pfafstetter level...")
#     for (j in 1:ncol(pTab)){
#       print(paste("pfaf level:", j))
#       pfafMatchTab[,j] = match(pTab[,j], uniqPfafCodes[[j]])
# 
#       # for each unique pfaf code:
#       for (k in 1:length(uniqPfafCodes[[j]])){
#         # identify basins that match the unique pfaf code:
#         keep = which(pfafMatchTab[,j] == k)
# 
#         # if the basin has less than a minimum length of river or
#         # if the basin has less than a minimum unique discharges, skip:
#         if (length(keep) < minLength |
#             length(unique(logQ[keep])) < minUniqQ){ next() }
# 
#         # run a linear regression in log space:
#         lreg = lm(logW[keep] ~ logQ[keep])
# 
#         # occationally when the logQs are very close, b_exp is NA. Skip these cases:
#         if (is.na(exp(lreg$coefficients[[1]])) | is.na(lreg$coefficients[[2]])){ next() }
# 
#         # calculate and store model parameters:
#         nSeg = length(keep)
#         nUniqQ = length(unique(logQ[keep]))
#         a_coef = exp(lreg$coefficients[[1]])
#         b_exp = lreg$coefficients[[2]]
#         R2 = summary(lreg)$r.squared
#         pVal = lmp(lreg)
# 
#         # match given pfaf basin to hybas shapefile attribute tab:
#         hybasInd = which(pTab[,j] == uniqPfafCodes[[j]][k])
# 
#         # define merge operation. If there is already an estimate
#         # for the given hydrobasin from a larger basin, select a method
#         # e.g.: replace with the higher R2, take mean, weighted by  R2,
#         # use smaller basin as long as it has an R2 > 0.5:
# 
#         # R2 must be equal or greater than predefined larger basin R2:
#         predefined = which(!is.na(oTab$R2[hybasInd]))
#         if (length(predefined) > 0){
#           gtR2 = oTab$R2[hybasInd][predefined] <= R2
#           oTab[hybasInd[predefined[gtR2]], ] = mget(names(oTab))
#         }else{
#           oTab[hybasInd, ] = mget(names(oTab))
#         }
# 
#         # plot:
#         # x = Q[keep]
#         # y = w[keep]
#         # plot(x, y,
#         #      log = 'xy',
#         #      xlab = "Q (cms)",
#         #      ylab = "Width (m)",
#         #      main = "")
#         # xseq = seq(min(x, na.rm=T), max(x, na.rm=T), length.out=100)
#         # yseq = exp(lreg$coefficients[[1]])*xseq^lreg$coefficients[[2]]
#         # lines(xseq, yseq, col=2, lwd=2)
#         # legend("topleft", paste0("h=", h, "   i=", i, "   j=", j, "   k=", k,
#         #                          "   a=", round(lreg$coefficients[[1]]),
#         #                          "   b=", round(lreg$coefficients[[2]],3),
#         #                          "   R2=", round(summary(lreg)$r.squared, 2)), text.col = 2)
# 
#       } # end k
#     } # end j
# 
#     # this is funky and a results of merging old code into new framework:
#     # generate copies of oTab for each recurrence interval:
#     assign(paste0("oTab_", names(qTab)[i]), oTab)
# 
#   } # end i
# 
#   paste0(print(length(which(!is.na(oTab$b_exp)))/nrow(oTab)), "% done of this recurrence interval")
# 
# } # end h
# 
# 
# 
# # write out hydroBASINS with pfafstetter optimized DHG values:
# for (i in 1:ncol(qTab)){
# 
# 
# 
#   oTab = get(paste0("oTab_", names(qTab)[i]))
# 
#   hybasOut = cbind(hybas, oTab)
#   foreign::write.dbf(hybasOut, hyBasOutPath)
# 
# 
# 
# }
# 
# 
# 
# system("rundll32 user32.dll,MessageBeep -1")
# 
# 
# 
# 
# 
# 
# 
# ##############################################  ##############################################  ##############################################
# 
# 
# 
# ################################################################################
# # DHG:
# ################################################################################
# # sort out headers:
# tab = gTab
# tab = tab[, order(names(tab))]
# QrecCols = grep("Q[[:digit:]]", names(tab))
# wCols_rev = grep("w[[:digit:]]", names(tab))
# wCols = rev(wCols_rev)
# wOccCols = wCols[grep("flag", names(tab)[wCols], invert=T)]
# wFlagCols = wCols[grep("flag", names(tab)[wCols])]
# 
# # multiply width occurence levels by cross section length: 
# tab[, wOccCols] = tab[, wOccCols] *  tab$width_m 
# 
# # filter:
# f = tab$strmOrder > 3 & tab$w100 > 0 & tab$width_m > 90
# #f = 1:nrow(tab)
# 
# qTab = tab[f, QrecCols]
# wTab = tab[f, wOccCols]
# fTab = tab[f, wFlagCols]
# 
# # set to NA measurements with values of zero or less and measuremets taken 
# # where river water is located at the end of their cross section segments:
# rmBoo = wTab<=0 | qTab<=0 | fTab==1 
# qTab[rmBoo] = NA
# wTab[rmBoo] = NA
# 
# 
# 
# 
# 
# # plot DHG colored by recurrence interval:
# qCols = rainbow(length(QrecCols), start=0, end=0.8)
# qCols_trans = rainbow(length(QrecCols), alpha=0.1)
# 
# # set up empty tables and vectors to store regression data: 
# aVec = bVec = r2Vec = rep(NA, ncol(qTab))
# xseqTab = as.data.frame(array(NA, c(20, ncol(qTab))))
# yseqTab = as.data.frame(array(NA, c(20, ncol(qTab))))
# 
# # set up empty tables and vectors to store Theil-Sen data: 
# aVec_sen = bVec_sen = rmse = rep(NA, ncol(qTab))
# xseqTab_sen = as.data.frame(array(NA, c(20, ncol(qTab))))
# yseqTab_sen = as.data.frame(array(NA, c(20, ncol(qTab))))
# 
# for (i in 1:ncol(qTab)){
#   
#   print(i)
#   
#   # set up vars:
#   Q = qTab[,i]
#   w = wTab[,i]
#   
#   # if there are no positive values, skip to next recurrence interval:
#   boo = Q>0 & w>0
#   if (!(T %in% boo)){next}
#   
#   # take log of Q and w data: 
#   logQ_raw = log(Q)
#   logW_raw = log(w)
#   
#   logQ = logQ_raw[!is.na(logQ_raw) | !is.na(logW_raw)]
#   logW = logW_raw[!is.na(logQ_raw) | !is.na(logW_raw)]
#   
#   
#   ####
#   # take least squares linear regression:
#   reg = lm(logW ~ logQ)
#   
#   # regression parameters:
#   aVec[i] = reg$coefficients[[1]]
#   bVec[i] = reg$coefficients[[2]]
#   r2Vec[i] = summary(reg)$r.squared
#   
#   # regression line vertices:
#   xseqTab[,i] = seq(min(Q, na.rm=T), max(Q, na.rm=T), length.out = nrow(xseqTab))
#   yseqTab[,i] = exp(aVec[i])*xseqTab[,i]^bVec[i]
#   
#   
#   
#   ####
#   # take Theil-Sen median estimator (this takes a long time to run):
#   ts = zyp.sen(logW ~ logQ)
#   
#   # regression parameters:
#   aVec_sen[i] = ts$coefficients[[1]]
#   bVec_sen[i] = ts$coefficients[[2]]
#   
#   xseqTab_sen[,i] = seq(min(Q, na.rm=T), max(Q, na.rm=T), length.out = nrow(xseqTab))
#   yseqTab_sen[,i] = exp(aVec_sen[i])*xseqTab[,i]^bVec_sen[i]
#   
# }
# 
# 
# 
# # set up plot:
# plot(range(xseqTab, na.rm=T), range(yseqTab, na.rm=T),
#      type="n", log="xy", las=1,
#      main="Downstream Hydraulic Geometry",
#      xlab = "Discharge (cms)",
#      ylab = "Width (m)",  cex.axis=0.8)
# 
# # points(Q, w, pch=16, col=qCols_trans[i], cex=0.8)
# matlines(xseqTab, yseqTab, lwd=2, lty=1, col=qCols)
# 
# legend("bottomright", title = "Flow Recurrence", 
#        paste0(names(tab[QrecCols]),
#               "  a=", round(aVec, 2),
#               "  b=", round(bVec, 2),
#               "  r2=", round(r2Vec, 2)), 
#        lty=1, lwd=2, col=qCols, box.col=NA, bg=NA)
# 
# 
# 
# 
# #pdfOut = "E:/research/2020_01_01_RSSA_seasonal/figs/drafts/DHG_fits.pdf"
# #pdf(pdfOut)
# 
# 
# # for each pfaf shapefile: 
# for (h in 1:length(inXsecPaths)){
#   
#   ##############################################
#   # Read in and process each continental shapefile:
#   ##############################################
#   
#   print(h)
#   print("reading in dbfs...")
#   
#   tab = foreign::read.dbf(inXsecPaths[h])
#   centroid = foreign::read.dbf(inCentroidPaths[h])
#   
#   # hydroBASIN pfafstetter code columns:
#   pfafColInd = grep("PFAF_", names(centroid))
#   
#   # join pfafstetter codes to cross section shapefile and 
#   mInd = match(tab$COMID, centroid$COMID)
#   
#   
#   
#   length(unique(tab$COMID))
#   
#   # sort out headers:
#   tab = tab[, order(names(tab))]
#   QrecCols = grep("Q[[:digit:]]", names(tab))
#   wCols_rev = grep("w[[:digit:]]", names(tab))
#   wCols = rev(wCols_rev)
#   wOccCols = wCols[grep("flag", names(tab)[wCols], invert=T)]
#   wFlagCols = wCols[grep("flag", names(tab)[wCols])]
#   
#   # multiply width occurence levels by cross section length: 
#   tab[, wOccCols] = tab[, wOccCols] *  tab$width_m 
#   
#   # filter:
#   f = tab$width_m > 90
#   #f = 1:nrow(tab)
#   
#   qTab = tab[f, QrecCols]
#   wTab = tab[f, wOccCols]
#   fTab = tab[f, wFlagCols]
#   
#   # set to NA measurements with values of zero or less and measuremets taken 
#   # where river water is located at the end of their cross section segments:
#   rmBoo = wTab<=0 | qTab<=0 | fTab==1 
#   qTab[rmBoo] = NA
#   wTab[rmBoo] = NA
#   
# }
# 
# 
# 
# 
# ################################################################################
# # Summary Stats:
# ################################################################################
# # take a quick look at the regression tables: 
# # lsTab = gTab[,grep("ls_", names(gTab))]
# # summary(lsTab)
# # hist(lsTab$AHG_ls_b[lsTab$AHG_ls_b<5], breaks=100)
# # abline(v=0.26, col=2)
# 
# tsTab = gTab[,grep("ts_", names(gTab))]
# summary(tsTab)
# hist(tsTab$AHG_ts_b[tsTab$AHG_ts_b<5], breaks=100)
# abline(v=0.26, col=2)
# 
# length(which(!is.na(tsTab$AHG_ts_b)))
# 
# 
# 
# 
# 
# R2 = rep(NA, nrow(tsTab))
# R2[tsTab$AHG_ts_r2>=0 & tsTab$AHG_ts_r2<0.5] = "0.00 < R2 < 0.5"
# R2[tsTab$AHG_ts_r2>=0.5 & tsTab$AHG_ts_r2<0.8] = "0.50 < R2 < 0.8"
# R2[tsTab$AHG_ts_r2>=0.8 & tsTab$AHG_ts_r2<0.9] = "0.8 < R2 < 0.90"
# R2[tsTab$AHG_ts_r2>=0.9 & tsTab$AHG_ts_r2<1] = "0.90 < R2 < 1.00"
# 
# keep = tsTab$AHG_ts_b < 5
# b = tsTab$AHG_ts_b[keep]
# r2 = R2[keep]
# bHistTab = data.frame(b, r2)
# 
# ggplot(bHistTab, aes(x=b, fill=r2)) +
#   geom_histogram()
# 
# 
# 
# 
# i = rowInd[1]
# ls = lm(as.numeric(wTabLog[i,]) ~ as.numeric(qTabLog[i,]))
# plot(as.numeric(qTab[i,keepCols]), wTab[i,keepCols])
# xSeq = seq(min(qTab[i,keepCols], na.rm=T), max(qTab[i,keepCols], na.rm=T), length.out=100)
# ySeq = exp(ls$coefficients[[1]])*xSeq^ls$coefficients[[2]]
# lines(xSeq, ySeq)
# 
# plot(as.numeric(wTabLog[i,]) ~ as.numeric(qTabLog[i,]))
# abline(ls)

# plot(c(0,1), c(0,1), type="n")
# for (i in (rowInd[seq(10000)])){
#   plotQ = as.numeric(qTab[i,keepCols])[!is.na(as.numeric(qTab[i,keepCols]))]
#   plotW = as.numeric(wTab[i,keepCols])[!is.na(as.numeric(wTab[i,keepCols]))]
#   normQ = (plotQ - min(plotQ)) / (max(plotQ) - min(plotQ))
#   normW = (plotW - min(plotW)) / (max(plotW) - min(plotW))
#   lines(normQ, normW, col=rgb(0,0,0,0.01))
# }
# 
# 
# 
# 
# 
# # plot AHG exponent and coefficient colored by R2:
# pal = colorRampPalette(c("red", "blue"))
# rCol = pal(10)[as.numeric(cut(lsTabLog$r2, breaks=10))]
# plot(lsTab$b, lsTab$a, log="y", col=rCol)
# legend('topright', levels(cut(lsTabLog$r2, breaks = 10)), col= pal(10), pch=16)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # set up plot:
# plot(range(xseqTab_sen, na.rm=T), range(yseqTab_sen, na.rm=T),
#      type="n", log="xy", las=1,
#      main="Downstream Hydraulic Geometry",
#      xlab = "Discharge (cms)",
#      ylab = "Width (m)",  cex.axis=0.8)
# 
# # points(Q, w, pch=16, col=qCols_trans[i], cex=0.8)
# matlines(xseqTab_sen, yseqTab_sen, lwd=2, lty=1, col=qCols)
# 
# legend("bottomright", title = "Flow Recurrence", 
#        paste0(names(tab[QrecCols]),
#               "  a=", round(aVec_sen, 2),
#               "  b=", round(bVec_sen, 2)),
#               #"  r2=", round(r2Vec, 2)), 
#        lty=1, lwd=2, col=qCols, box.col=NA, bg=NA)
# 
# 
# 
# 
# 
# 
# 
# 
#
# Yellow River
# plot(tab$Q50[f], tab$w050[f], log="xy")
# yel = tab$Q50[f] >5e2 & tab$Q50[f] < 5e3
# points(tab$Q50[f][yel], tab$w050[f][yel], col=2)
# fSeq = which(f)[yel]
# 
# 
# 
# 
# ################################################################################
# # plot Q-w relations over a sample of cross sections:
# pdfOut = pdfOutPaths[h]
# pdf(pdfOut)
# 
# # skip cross sections without at least 3 valid measurements:
# validRows = which(apply(qTab, 1, function(x) length(which(!is.na(x))))>8)
# fSeq = sample(validRows, 100) 
# 
# for (i in fSeq){
#   notNA = !is.na(qTab[i,])
#   Q = as.numeric(qTab[i,notNA])
#   w = as.numeric(wTab[i,notNA])
#   
#   plot(Q, w,
#        main=i, 
#        log="xy",
#        xlab="Q (cms)",
#        ylab="W (m)",
#        pch=16,
#        col=qCols[notNA])
#   
#   # plot power law regression:
#   ls = lm(log(w) ~ log(Q))
#   xSeq = seq(min(Q), max(Q), length.out=100)
#   ySeq = exp(ls$coefficients[[1]])*xSeq^ls$coefficients[[2]]
#   lines(xSeq, ySeq, col=rgb(0,0,0,0.5))
#   text(min(Q), max(w),
#        paste0(
#          "w = ", round(exp(ls$coefficients[[1]]),1), " * Q ^ ", round(ls$coefficients[[2]], 2), "\n",
#          "R2 = ", round(summary(ls)$r.squared, 2)),
#        adj = c(0, 1))
#   legend("bottomright", names(tab[QrecCols]), pch=16, col=qCols)
#   
#   # plot expoential regression:
#   ls = lm(log(w) ~ log(Q))
#   
# }
# 
# 
# dev.off()
# cmd = paste('open', pdfOut)
# system(cmd)
# 
# 
# 
# 
# 