################################################################################
# HG_analysis.R
################################################################################
# George H. Allen, Nov. 2019

# Description:
# Reads in joined flow occurrence and GRADES Q data and runs regressions for 
# global hydraulic geometry relationships.


################################################################################
# load libraies:
################################################################################

if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)

if (!"zyp" %in% rownames(installed.packages())){
  install.packages("zyp")}; require(zyp)


################################################################################
# Define paths:
################################################################################

wDir = "E:/research/2020_01_01_RSSA_seasonal/git/RSSA_seasonal"

inDir = paste0(wDir, "/in/GRWL/GRWL_Xsections_MERIT_spJoin/")
outDir = paste0(wDir, "/out/GRWL/GRWL_Xsections_HG/")


################################################################################
# functions:
################################################################################
DHG <- function( tab, QrecCols, qTab, wTab){
  print("hello world")
}


################################################################################
# set up files:
################################################################################

# copy any shapfiles:
inPaths = list.files(inDir, full.names=T)
file.copy(inPaths, sub(inDir, outDir, inPaths), overwrite=F)

# get in paths:
inDBFpaths_raw = list.files(inDir, ".dbf", full.names=T)

# rm Greenland bc innaccurate GRADES Q data:
inDBFpaths = inDBFpaths_raw[-grep("pfaf_09", inDBFpaths_raw)] 

outDBFpaths = sub(inDir, outDir, inDBFpaths)
pdfOutPaths = sub(".dbf", ".pdf", outDBFpaths)


################################################################################
# Main code
################################################################################

pdfOut = "E:/research/2020_01_01_RSSA_seasonal/figs/drafts/AHG_fits.pdf"
pdf(pdfOut)


# for each pfaf shapefile: 
for (h in 1:length(inDBFpaths)){
  
  inDBFpath = inDBFpaths[h] # "E:/research/2020_01_01_RSSA_seasonal/GIS/out/GRWL/GRWL_Xsections_MERIT_spJoin_HG_merge/GRWL_Xsections_MERIT_spJoin_HG_merge.dbf"# 
  
  print(h)
  print("reading in...")
  
  tab = foreign::read.dbf(inDBFpath)
  
  # sort out headers:
  tab = tab[, order(names(tab))]
  QrecCols = grep("Q[[:digit:]]", names(tab))
  wCols_rev = grep("w[[:digit:]]", names(tab))
  wCols = rev(wCols_rev)
  wOccCols = wCols[grep("flag", names(tab)[wCols], invert=T)]
  wFlagCols = wCols[grep("flag", names(tab)[wCols])]
  
  # multiply width occurence levels by cross section length: 
  tab[, wOccCols] = tab[, wOccCols] *  tab$width_m 
  
  # filter:
  #f = tab$strmOrder > 3 #& tab$w100 > 0 #& tab$width_m > 100
   f = 1:nrow(tab)
  
  qTab = tab[f, QrecCols]
  wTab = tab[f, wOccCols]
  fTab = tab[f, wFlagCols]
  
  # set to NA measurements with values of zero or less and measuremets taken 
  # where river water is located at the end of their cross section segments:
  rmBoo = wTab<=0 | qTab<=0 | fTab==1 
  qTab[rmBoo] = NA
  wTab[rmBoo] = NA
  
  
  
  ################################################################################
  #### AHG:
  ################################################################################
  
  
  # filter out low-quality data:
  
  # first, remove 0th and 100th percentile w-Qs due to low reliability (Allen et al., 2019). 
  # Now max number w-Q pairs is equal 9. 
  keepCols = grep("Q000|Q100", names(qTab), invert=T)
  qTabLog = log(qTab[,keepCols])
  wTabLog = log(wTab[,keepCols])
  
  # second, exlude percentiles without at least 5 unique width measurements:
  nUniqWidthObs = as.numeric(apply(wTabLog, 1, function(x){length(which(!is.na(unique(as.numeric(x)))))}))
  rowInd = which(nUniqWidthObs >= 5)
  
  # generate least squares exponent and coefficient table:
  lsTabLog = as.data.frame(array(NA, c(nrow(qTabLog), 3)))
  names(lsTabLog) = c("a", "b", "r2")
  
  # generate theil-sen exponent and coefficient table:
  tsTabLog = as.data.frame(array(NA, c(nrow(qTabLog), 3)))
  names(tsTabLog) = c("a", "b", "r2")
  
  print(paste(length(rowInd), "cross sections to analyze."))
  print("running regressions...")
  
  # for each cross section, calculate several statistics related to AHG: 
  for (i in rowInd){
    
    w = as.numeric(wTabLog[i,])[!is.na(wTabLog[i,])]
    Q = as.numeric(qTabLog[i,])[!is.na(qTabLog[i,])]
    
    ###
    # least squares linear regression:
    ls = lm(w ~ Q)
    
    # add statistics to table: 
    lsTabLog[i,] = c(ls$coefficients, suppressWarnings(summary(ls)$r.squared), RMSE, MAE, MB)
    
    
    ###
    # Thiel-Sen estimate:
    ts = zyp.sen(w ~ Q)
    
    # R2 using Theil-Sen fit:
    R2 = cor(w, ts$coefficients[[2]]*Q+ts$coefficients[[1]]) ^ 2
    
    # add statistics to table: 
    tsTabLog[i,] = c(ts$coefficients, R2)
    
    
    # plot: 
    plot(w ~ Q, asp=1, main=paste0("h=", h, "  i=", i))
    abline(ls, col=2)
    abline(ts$coefficients, col=4)
    
    legend("bottomright", 
           legend = c("Least Squares", 
                      paste(names(lsTabLog), round(lsTabLog[i,], 2), sep="="),
                      "",
                      "Theil-Sen", 
                      paste(names(tsTabLog), round(tsTabLog[i,], 2), sep="=")),
           text.col=c(rep(2,5), rep(4,5)))
    
    
    
    # # RMSE:
    # RMSE = mean((w - (ts$coefficients[[2]]*Q+ts$coefficients[[1]]))^2)^0.5
    # 
    # # mean absolute error:
    # MAE = mean(abs((ts$coefficients[[2]]*Q+ts$coefficients[[1]]) - w))
    # 
    # # mean bias:
    # MB = mean((ts$coefficients[[2]]*Q+ts$coefficients[[1]]) - w)
    
  }
  
  # process regression tables: 
  lsTab = lsTabLog
  lsTab$a = exp(lsTabLog$a)
  tsTab = tsTabLog
  tsTab$a = exp(tsTabLog$a)
  
  # take a quick look at the regression tables: 
  summary(lsTab)
  summary(tsTab)
  
  summary(lsTab[lsTab$r2 > 0.5, ])
  summary(tsTab[tsTab$r2 > 0.5, ])

  
  # add these best-fit regressions to shapefile dbf:
  tab$AHG_ls_a = tab$AHG_ls_b = tab$AHG_ls_r2 = 
    tab$AHG_ts_a = tab$AHG_ts_b = tab$AHG_ts_r2 = NA
  
  tab$AHG_ls_a[f] = lsTab$a
  tab$AHG_ls_b[f] = lsTab$b
  tab$AHG_ls_r2[f] = lsTab$r2
  tab$AHG_ts_a[f] = tsTab$a
  tab$AHG_ts_b[f] = tsTab$b
  tab$AHG_ts_r2[f] = tsTab$r2
  
  print("writing out...")
  foreign::write.dbf(tab, outDBFpaths[h])
  
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)








