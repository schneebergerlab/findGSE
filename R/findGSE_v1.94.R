# findGSE: estimating genome size by fitting k-mer frequencies in short reads
#          with a skew normal distribution model.
#
# Copyright (C) 2016-2017 Hequan Sun (MPIPZ)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# For any questions on usage or potential bugs, please contact via schneeberger@mpipz.mpg.de or sunhequan@gmail.com .
#
##  Assumption and process
###  assume a skew normal distribution of kmer-freq in
###    the unique regions and subsequent repetitive regions
###  copy the observed values on the right side of the peak to form
###    a sequence (to initialize mean and variance of a skew-normal-distribution)
###  optimize the skew-normal dsitribution curve to fit the observed
###    values from proper 'xleft' to 'xright'
## Result: the estimate is corrected to be larger or smaller (due to
##   the left/right-shift of the fitted peak; estimated counts from 1, 2,...
##   can result in the increased estimates).
## with first fitting if 'heterozygous' known! -- to use fitting for heterozygous,
##   expected Kmer coverage for homozygous region must be provided.
## het_fitting_delta.
##
## output:
##  curve-fit .pdf: *curvefitted.pdf, and a corresponding
##            .txt: listing gs estimates (naive, fitted and final-corrected)
##                  and some other info like repetitive ration (and hetero-rate).

## function for minimizing difference from original observed kmer-freq dr: mean, sd, xi, and scaling factor ##
error_minimize<-function(tooptimize, x, end, xfit, xfit_left, xfit_right, d, min_valid_pos, itr)
{
  # x            : histogram replicated from d, which is formed from dr
  # xfit         : x-values to get the skew normal distribution
  # xfit-left    : a left-side  position to calculate initial mean and sd
  # xfit-right   : a right-side position to calculate initial mean and sd
  # d            : the observed kmer-freq that will be fitted
  # min_valid_pos: a left-side  position, from which the observed kmer-freq will be fitted
  # end          : a rigth-side position, till which the observed kmer-freq will be fitted
  sd_scale_factor    <- tooptimize[1]
  yfit_scale_factor  <- tooptimize[2]*100000/itr
  mean_scale_factor  <- tooptimize[3]
  xi_scale_factor    <- tooptimize[4]
  xfit               <- seq(min(x),max(x),length=end)
  meanfit            <- mean(x[x>=xfit_left & x<=xfit_right])*mean_scale_factor
  sdfit              <- sd(x[  x>=xfit_left & x<=xfit_right])*sd_scale_factor
  xifit              <- 1*xi_scale_factor
  yfit               <- dsnorm(xfit,mean=meanfit,sd=sdfit, xi=xifit)*yfit_scale_factor*length(x)
  diff0              <- sqrt(sum((d[min_valid_pos:end, 2] - yfit[min_valid_pos:end])^2))
  #
  return(diff0)
}
# function for tunning final fitting if heterozgyous genomes
error_minimize2<-function(tooptimize, h_het, h_hom, h_target)
{
  # h_het   : raw fitting for the heterozygous region
  # h_hom   : raw fitting for the homozygous   region
  # h_target: before valley: fitted cnt, after valley: raw cnt
  het_delta <- 1-tooptimize[1]
  hom_delta <- 1 # reserved parameter: might be applied in future
  h_merge   <- het_delta*h_het + hom_delta*h_hom
  diff2     <- sum(abs(h_merge - h_target))
  #
  return(diff2)
}
## function for minimizing difference from remainder kmer-freq dr: mean, sd, and scaling factor ##
error_minimize3<-function(tooptimize, x, end, xfit_left, xfit_right, d, min_valid_pos, itr)
{
  # x            : histogram
  # xfit         : x-values to get the normal distribution
  # xfit-left    : a left-side  position to calculate initial mean and sd
  # xfit-right   : a right-side position to calculate initial mean and sd
  # d            : the observed kmer-freq that will be fitted
  # min_valid_pos: a left-side  position, from which the observed kmer-freq will be fitted
  # end          : a rigth-side position, till which the observed kmer-freq will be fitted
  sd_scale_factor    <- tooptimize[1]
  yfit_scale_factor  <- tooptimize[2]*100000/itr
  mean_scale_factor  <- tooptimize[3]
  xfit               <- seq(min(x),max(x),length=end)
  meanfit            <- mean(x[x>=xfit_left & x<=xfit_right])*mean_scale_factor
  sdfit              <- sd(x[  x>=xfit_left & x<=xfit_right])*sd_scale_factor
  yfit               <- dnorm(xfit,mean=meanfit,sd=sdfit)*yfit_scale_factor*length(x)
  diff0              <- sqrt(sum((d[min_valid_pos:end, 2] - yfit[min_valid_pos:end])^2))
  #
  return(diff0)
}
# recover count 0: in initial count, making kmer freq consecutive.
initial_count_recover <- function(d0)
{
  # dr: the initial kmer count from softwares like jellyfish
  dr <- cbind(1:d0[length(d0[,1]), 1], rep(0, d0[length(d0[,1]), 1]))
  i  <- 0
  while( i < length(d0[,1]))
  {
    i               <- i + 1
    dr[d0[i, 1], 2] <- d0[i, 2]
  }
  return (dr)
}
# modify kmer freq before fitting
kmer_count_modify <- function(start, end, left_right, histx)
{
  # start or end does not include the peak position
  # modify rigth part according to left
  if(left_right == 0)
  {
    for (i in start:end)
    {
      histx[2*end+2-i, 2]   <- histx[i, 2]
    }
  }else
  {
    for (i in start:end)
    {
      histx[2*start-2-i, 2] <- histx[i, 2]
    }
  }
  return (histx)
}
################################################## main #############################################################
#' @title Estimating genome size by fitting k-mer frequencies in short reads
#' with a skew normal distribution model.
#'
#' @description findGSE is a function for (heterozygous diploid or homozygous)
#' genome size estimation by fitting k-mer frequencies iteratively
#' with a skew normal distribution model. (version still under testing)
#'
#' @description To use findGSE, one needs to prepare a histo file,
#' which contains two tab-separated columns.
#' The first column gives frequencies at which k-mers occur in reads,
#' while the second column gives counts of such distinct k-mers.
#' Parameters k and related histo file are required for any estimation.
#'
#' @description For heterozygous genomes, another parameter about
#' the average k-mer coverage for the homozygous regions must be provided.
#'
#' @param histo is the histo file (mandatory).
#' @param sizek is the size of k used to generate the histo file (mandatory).
#' K is involved in calculating heterzygosity if the genome is heterozygous.
#' @param outdir is the path to write output files (optional).
#' If not provided, by default results will be written in the folder
#' where the histo file is.
#' @param exp_hom a rough average k-mer coverage for finding the homozygous regions.
#' In general, one can get peaks in the k-mer frequencies file, but has to
#' determine which one is for the homozygous regions, and which one is for the
#' heterozygous regions. It is optional, however, it must be provided
#' if one wants to estimate size for a heterozygous genome.
#' VALUE for exp_hom must satisfy fp < VALUE < 2*fp, where fp is the freq for homozygous peak.
#' If not provided, 0 by default assumes the genome is homozygous.
#' @param species an optional parameter only applied in calculating heterozygosity for human.
#' This is used to indicate that (lx-ly)*hom_c/2 k-mers should be removed from het-kmers,
#' where lx is length of chromosome X, ly is length of chromosome Y, and hom_c is the
#' average k-mer coverage for the homozygous k-mers.
#' Two estimates will be provied as ORIGINAL_EST CORRECTED_EST:
#' for males,   select the second (CORRECTED_EST);
#' for females, select the first  (ORIGINAL_EST)..
#' @export
findGSE <- function(histo="", sizek=0, outdir="", exp_hom=0, species="")
{
  # initial values
  if(missing(histo))
  {
    stop(paste("Cannot find histo file: ", histo, ". Program exited.\n",sep=""))
  }
  if(missing(sizek))
  {
    stop(paste("Cannot find sizek info: ", sizek, ". Program exited.\n",sep=""))
  }
  # defaults
  if(missing(exp_hom))  exp_hom <- 0
  if(missing(outdir))   outdir  <- getwd()
  ######################################## libararies required ###############################################
  # find all peaks in a time series
  suppressWarnings(suppressMessages(library("pracma")))
  # provide skew normal distrbution fitting (dsnorm(...))
  suppressWarnings(suppressMessages(library("fGarch")))
  ## version id
  vers <- 'v1.94.'
  ##
  cat("\nfindGSE initialized...\n")
  ## running time
  start_t <- Sys.time()
  ## check args
  # expert-tunning het fitting (later it would be automatically adjusted)
  het_fitting_delta <- 0.11
  histof  <- F
  kpro    <- F
  if(histo != "") histof <- T
  if(sizek >  0)  kpro   <- T
  if(histof == F || kpro == F)
  {
    stop("Option -histo or -sizek was provided incorrectly. Program exited.\n")
  }
  ## default args
  countingfile    <- histo
  if(file.exists(countingfile)==F)
  {
    stop(paste("Cannot find histo file: ", countingfile, ". Program exited.\n",sep=""))
  }
  cat(paste("   Info: histo file provided as ", countingfile, "\n", sep=""))
  #
  targetsizek     <- sizek
  if(targetsizek  < 15)
  {
    stop(paste("Unexpected size of k: ", targetsizek, ". Size k>=15 required. Program exited.\n",sep=""))
  }
  cat(paste("   Info: size k set as ", targetsizek, "\n", sep=""))
  #
  path_set        <- F
  path            <- outdir
  if(path==".")
  {
    path <- getwd()
  }else
  if(path == "")
  {
     path         <- dirname(normalizePath(countingfile))
  }
  path_set       <- T
  cat(paste("   Info: output folder set as ", path, "\n", sep=""))
  #
  het_observed    <- F
  # expected maximum value for the hom k-mer peak! set as fpeak~2*fpeak
  exp_hom_cov_max <- exp_hom
  cat(paste("   Info: expected coverage of homozygous k-mers set as ", exp_hom_cov_max, "\n", sep=""))
  # derived  maximum value for the het k-mer peak!
  der_het_cov_max <- 0
  ## if output folder does not exist, create it; if existing, give a warning.
  dir.create(file.path(path), showWarnings = T)
  if(exp_hom_cov_max > 0)
  {
    het_observed      <- T
    # derived  maximum value for the het k-mer peak.
    der_het_cov_max   <- round(0.5*exp_hom_cov_max);
    cat(paste("   Info: het observed set as true ==> heterozygous fitting asked. \n", sep=""))
    if(het_fitting_delta != 0.11)
    {
      cat(paste("   Info: value for tuning het fitting set as ", het_fitting_delta, "\n", sep=""))
    }
  }else
  {
    het_observed      <- F
    only_one_het_peak <- F
    only_one_hom_peak <- T
    cat(paste("   Info: het observed set as false ==> heterozygous fitting not asked. \n", sep=""))
  }
  ##
  prefix <- ''
  samples          <- basename(countingfile)
  for (sample in samples)
  {
    cat(paste("\nIterative fitting process for sample ", sample, " started... \n", sep=""))
    # all accessions -- where i can change and modify
    # initialize variables where outputs are recorded
    size_all       <-  0   # Gsize with               error-k-mer-freq
    size_exl       <-  0   # Gsize without            error-k-mer-freq
    size_fit       <-  0   # Gsize with              fitted-k-mer-freq
    size_cat       <-  0   # Gsize with fitted-cat-observed-k-mer-freq
    size_cor       <-  0   # Gsize corrected with size_cat and skewness; size_cor=size_cat/skewness (no)
    size_cor2      <-  0   # Gsize corrected with fitted-cat-observed-freq and upperbound freq 'e' (yes)
    i              <-  0   # the ith sample/accession
    kmer_peak_cnt  <-  0   # k-mer count at the peak position of k-mer counting
    kmer_peak_pos  <-  0   # genome-wide k-mer coverage
    kmer_peak_pos2 <-  0   # genome-wide k-mer coverage according to fitted curve
    labels         <- NA   # labels according input below
    namebank       <- NA
    border_pos     <-  1
    min_valid_pos  <-  1
    first_peak_pos <-  1
    first_mean     <-  0
    first_sd       <-  1
    first_mean_raw <-  0
    # initialize parameters for several iterations of fitting
    mymeanfit      <-  0
    mysdfit        <-  1
    myxifit        <-  1
    myscalefit     <-  1
    # special fit at itr 1
    norm_fit_itr1  <-  F
    # open pdf to write visualization of fitting
    pdf(paste(path, '/',prefix, vers, 'est.', sample, '.sizek', targetsizek, '.curvefitted.pdf',
              sep=""))
    first_fit  <- TRUE
    for (sizek in targetsizek[1]) {
      normal   <- 0
      if(file.exists(countingfile)==TRUE)
      {
        dhet   <- read.table(countingfile)
        ## recover positions in dr whose initial counts are 0
        dr     <- initial_count_recover(dhet[1:length(dhet[,1]),])
        rm("dhet")
        dhet   <- dr
        # if heterozygosity oberved in k-mer freq,
        #    fit het and remove fitted values from subsequent hom-fitting
        selected    <- min(2000, length(dhet[,1]))
        if(selected < exp_hom_cov_max)
        {
          selected  <- min(round(1.5*exp_hom_cov_max), length(dhet[,1]))
        }
        main_peak_is_hom <- T
        if(het_observed)
        {
          # initialize count as error for fitting
          error    <- dr
          cat(paste("    Size ", sizek , " fitting for het k-mers \n", sep=""));
          # find hom and het peaks (if any).
          peaks    <- findpeaks(error[2:exp_hom_cov_max, 2])
          if(length(peaks) == 0)
          {
            plot(dhet[, 1], dhet[, 2],
                 col="gray",
                 xlab="K-mer freq", ylab="Number of k-mers",
                 xlim=c(1, 4*(which.max(dhet[5:selected, 2])+4)), ylim=c(1, 1.5*max(dhet[5:selected, 2])))
            text(1.5*(which.max(dhet[5:selected, 2])+4), 1.1*max(dhet[5:selected, 2]),
                 paste("Main peak (after freq=5) at k-mer freq: ",
                       which.max(dhet[5:selected, 2])+4,
                       sep=""))
            dev.off();
            stop("No k-mer freq peak was found with the cutoff given by -exp_hom.
                  You may need to increase the cutoff!")
            quit("no")
          }
          hom_peak_pos <- 0
          het_peak_pos <- 0
          hom_peak_pos <- peaks[which.max(peaks[,1]), 2]+1   # hom peak
          peaks[which.max(peaks[,1]), 1] <- 0
          het_peak_pos <- peaks[which.max(peaks[,1]), 2]+1   # het peak
          peaks[which.max(peaks[,1]), 1] <- 0
          # if two peaks found, the larger one is hom; the other is het
          if(het_peak_pos > hom_peak_pos)
          {
            tmp               <- het_peak_pos
            het_peak_pos      <- hom_peak_pos
            hom_peak_pos      <- tmp
            main_peak_is_hom  <- F
            cat("    Warning: two peaks found in k-mer freq with given -exp_hom,
                  peak with lower height as hom-peak!\n")
          }
          if(het_peak_pos>0 & hom_peak_pos>0)
          {
            # in case two peask are near each other (e.g. it may happen with k-mer cov >100x)
            ratio        <- hom_peak_pos/het_peak_pos;
            if (round(ratio)<1.5) # nearly equal
            {
              hom_peak_pos <- max(hom_peak_pos, het_peak_pos)
              het_peak_pos <- max(hom_peak_pos, het_peak_pos)
            }
          }
          # if one peak  found, need to check if it is hom or het accordingly
          only_one_hom_peak <- F
          only_one_het_peak <- F
          if(het_peak_pos == hom_peak_pos)
          {
            unkownpeak    <- het_peak_pos
            if(unkownpeak >  der_het_cov_max)
            {
              het_peak_pos      <- round(0.5*hom_peak_pos)
              cat("    Warning: only one peak in k-mer freq with given -exp_hom,
                  determined as hom-peak!\n")
              cat("            -exp_hom needs to be increased
                  if you think it is a het peak.\n")
              only_one_hom_peak <- T
            }else
              if(unkownpeak  < der_het_cov_max)
              {
                hom_peak_pos      <- 2*het_peak_pos
                cat("    Warning: only one peak in k-mer freq with given -exp_hom,
                    determined as het-peak!\n")
                cat("          -exp_hom needs to be decreased
                    if you think it is a hom peak.\n")
                main_peak_is_hom  <- F
                only_one_het_peak <- T
              }else
              {
                plot(dhet[, 1], dhet[, 2],
                     col="gray",
                     xlab="K-mer freq", ylab="Number of k-mers",
                     xlim=c(1, selected), ylim=c(1, 1.5*max(dhet[3:selected, 2])))
                text(100, max(dhet[3:selected, 2]),
                     paste("Main peak at k-mer freq: ",
                           which.max(dhet[3:selected, 2])+2,
                           sep=""))
                dev.off()
                stop("Unexpected peaks found -
                     please check the raw k-mer curve provided in pdf,
                     and reset the option -exp_hom with a proper value x,
                     which must satisfy 2*hom_peak>x>hom_peak !\n\n");
              }
          }
          cat('    Info: het_peak_pos for het fitting: ', het_peak_pos, '\n')
          cat('    Info: hom_peak_pos for hom fitting: ', hom_peak_pos, '\n')
          het_min_valid_pos <- which.min(error[1:het_peak_pos, 2])
          # fit het distribution
          het_xfit_left     <- max(min(het_min_valid_pos, het_peak_pos-10), 3)
          het_xfit_right   <- het_peak_pos + round(0.5*(hom_peak_pos-het_peak_pos))
          cat('    Info: het_xfit_left  for het fitting: ', het_xfit_left, '\n')
          cat('    Info: het_xfit_right for het fitting: ', het_xfit_right, '\n')
          if(het_xfit_right > het_peak_pos-1 - het_xfit_left + het_peak_pos-1)
          {
            het_xfit_right  <- het_peak_pos-1 - het_xfit_left + het_peak_pos-1
          }
          histx2    <- kmer_count_modify(het_xfit_left, het_peak_pos-1, 0, error)
          x         <- rep(1:het_xfit_right,ceiling(abs(histx2[1:het_xfit_right,2])/(100000/1)))
          xfit      <- seq(min(x),max(x),length=het_xfit_right)
          # optimizing parameters
          v<-optim(par=c(1, 1, 1, 0.8),
                   fn=error_minimize, x=x, end=het_xfit_right, xfit=xfit,
                   xfit_left=het_xfit_left, xfit_right=het_xfit_right,
                   d=histx2, min_valid_pos=het_min_valid_pos, itr=1)
          #
          # scale and correct yfit with optimized values in v$par
          xfit2     <- seq(0.01,selected,0.01)
          meanfit   <- mean(x[x>=het_xfit_left & x<=het_xfit_right])*v$par[3]
          if(meanfit < 0)
          {
            cat(paste("    Warning: data does not follow
                      assumed distribution anymore; fitting stopped.\n", sep=""))
          }
          sdfit     <- sd(x[  x>=het_xfit_left & x<=het_xfit_right])*v$par[1] # raw sd
          #
          xifit     <- 1*v$par[4];
          yfit0     <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/1)*v$par[2]*length(x)
          #
          ## scale fitting according to observation
          hetfit    <- yfit0[seq(100, length(yfit0), 100)]
          hetfit    <- hetfit*dhet[het_peak_pos,2]/max(hetfit[1:selected])
          ## end scaling ##
          #
          # set d0 for further hom-fitting
          d0        <- dr
          if(only_one_hom_peak)
          {
            d0[1:selected , 2] <- d0[1:selected , 2] - (1-het_fitting_delta)*hetfit[1:selected]
          }else
          {
            d0[1:selected , 2] <- d0[1:selected , 2] - hetfit[1:selected]
          }
        }else
        {
          d0                 <- dr
          hetfit             <- 0
          hetfit[1:selected] <- 0
          het_min_valid_pos  <- 100000
          only_one_hom_peak  <- T
        }
        i    <- i + 1
        if(i > length(targetsizek)) { break }
        rm("dr")
        # dr can be different from raw k-mer counting due to het fit
        dr   <- d0
        #
        # initialize count as error for fitting: initial error vector
        error        <- dr
        # set the right-bound for fitting
        maxright4fit <- dr[length(dr[,1])-1, 1]
        maxright4fit <- min(10000,length(dr[,1]))
        # set the number of total iterations corresponding
        #     to raw count given by softwares such as jellyfish
        totalitr     <- 10
        myitr        <-  0
        myyfit       <-  0
        homfit       <-  0
        #
        first_fit    <- TRUE
        itr = 1
        while(itr <= totalitr)
        {
          if(first_fit)
          {
            cat(paste("    Size ", sizek , " at itr ", itr, "\n", sep=""))
            if(itr>=2 && 2*mymeanfit[itr-1] > maxright4fit) { myitr <- itr - 1; break;}
            if(itr == 1) # original data
            {
              # 1st-col    : height,
              # 2nd-col    : index of maximum,
              # 3rd/4th-col: indices that peak begins/ends
              if(het_observed)
              {
                # for  border_pos
                peaks           <- findpeaks(error[3:selected, 2])
                peaks[,2:4]     <- peaks[,2:4]+2
                border_pos      <- peaks[1, 3]          # caution: peaks[peakindex, 3]
                # for others: caution!
                if(main_peak_is_hom==F)
                {
                  peaks           <- findpeaks(error[round(0.5*(het_peak_pos+hom_peak_pos)):selected, 2])
                  peaks[,1:4]     <- peaks[,1:4]+round(0.5*(het_peak_pos+hom_peak_pos))-1
                }else
                {
                  peaks           <- findpeaks(error[het_peak_pos:selected, 2])
                  peaks[,1:4]     <- peaks[,1:4]+het_peak_pos-1
                }
              }
              else
              {
                peaks           <- findpeaks(error[3:selected, 2])
                peaks[,2:4]     <- peaks[,2:4]+2
                border_pos      <- peaks[1, 3]          # caution: peaks[peakindex, 3]
              }
              peakindex         <- which.max(peaks[,1]) # caution: if double peaks
              min_valid_pos     <- peaks[1, 3]          # caution: peaks[peakindex, 3]
              original_peak_pos <- peaks[peakindex, 2]
              original_peal_val <- peaks[peakindex, 1]
              if(het_observed && only_one_het_peak)
              {
                original_peak_pos <- 2*het_peak_pos
                original_peal_val <- peaks[original_peak_pos, 1]
                min_valid_pos     <- het_peak_pos
              }
              datacheck         <- min(300, maxright4fit)
              dci               <- 0
              cat('    Info: min_valid_pos: ', min_valid_pos,'\n')
              while(min_valid_pos == original_peak_pos)
              {
                min_valid_pos     <- min_valid_pos + 1
                original_peak_pos <- which.max(error[min_valid_pos:maxright4fit,2])+(min_valid_pos-1)
                original_peal_val <- max(error[min_valid_pos:maxright4fit,2])
                dci               <- dci + 1
                if(dci > datacheck) { normal=1;break; }
              }
              if(dci   > datacheck)
              {
                cat(paste("Error: weired initial count-data of sample ",
                          sample, ". Please check!!!\n",
                          sep=""))
                normal <- 1
                break
              }
              cat('    Info: signal error border: ', border_pos,'\n')
              first_peak_pos  <- original_peak_pos
              # copy
              d               <- error
              # complement the left side of the peak position
              #    according to observations on the right side of the peak
              for (li in 1:(original_peak_pos-1))
              {
                d[li, 2] <- d[min(2*original_peak_pos-li, length(dr[,1])),2]
              }
              end <- min(original_peak_pos+round((original_peak_pos-min_valid_pos)*5/12),length(dr[,1]))
              if(end - original_peak_pos <= 3 && original_peak_pos>=10)
              {
                end      <- end + 5
              }
            }
            else # error reminder
            {
              # caution: valid kmer count is at least 3
              min_valid_pos     <- round(mymeanfit[itr-1])+3
              original_peak_pos <- which.max(error[min_valid_pos:maxright4fit, 2])+(min_valid_pos-1)
              original_peal_val <- max(error[min_valid_pos:maxright4fit, 2])
              datacheck         <- min(300, maxright4fit)
              dci               <- 0
              while(min_valid_pos == original_peak_pos)
              {
                min_valid_pos     <- min_valid_pos + 1
                original_peak_pos <- which.max(error[min_valid_pos:maxright4fit,2])+(min_valid_pos-1)
                original_peal_val <- max(error[min_valid_pos:maxright4fit,2])
                dci               <- dci + 1
                if(dci > datacheck) { break; }
              }
              if(dci   > datacheck)
              {
                cat(paste("    Warning: weired initial count-data of sample ",
                          sample, ". Please check!!!\n", sep=""))
                normal   <- 0
                xlimmax  <- 100
                break
              }
              # copy
              d          <- abs(error)
              # complement the left side of the peak position according to observations on the right side of the peak
              thisnormal <- 0
              for (li in (original_peak_pos-round(mymeanfit[1])):(original_peak_pos-1)) {
                if(2*original_peak_pos-li > maxright4fit){ thisnormal<-1; break}
                d[li, 2] <- d[2*original_peak_pos-li,2]
              }
              if(thisnormal==1){normal=0; break;}
              end        <- min(original_peak_pos+round(mymeanfit[itr-1]), maxright4fit)
              if(end - original_peak_pos <= 3 && original_peak_pos>=10)
              {
                end      <- end + 5
              }
            }
          }
          else
          {
            # due to the large diff in the first fitting and original curve, increase end for fitting.
            end           <- min(2*original_peak_pos-1, maxright4fit)
            min_valid_pos <- min_valid_pos + 1 ##### 2016-09-16 new
            d             <- error
            first_fit     <- TRUE
          }
          # prepare the histogram data with replication
          forrep     <- ceiling(abs(d[1:end,2])/(100000/itr))
          anyNA      <- which(is.na(forrep))
          if (length(anyNA) >= 1)
          {
            d[anyNA, 2] <- dr[anyNA, 2]; # caution: missing values could happen here!!!
            cat(paste("    Warning: ", length(anyNA),
                      " data-points of d has been replaced as 1 to properly prepare hisgram data x
                      for function optim at itr ", itr, "\n", sep=""))
          }
          x          <- rep(1:end,ceiling(abs(d[1:end,2])/(100000/itr)))
          xfit       <- seq(min(x),max(x),length=end)
          # normal fit
          matchset   <- '0'
          if(!is.na(match(sizek,matchset)))
          {
            xfit_left  <- min_valid_pos
            xfit_right <- 2*original_peak_pos-xfit_left
            if(first_peak_pos>=200)
            {
              xfit_left  <- original_peak_pos - 50
              xfit_right <- original_peak_pos + 50
            }
            v<-optim(par=c(1, 1, 1, 1),
                     fn=error_minimize, x=x, end=end, xfit=xfit,
                     xfit_left=xfit_left, xfit_right=xfit_right,
                     d=error,
                     min_valid_pos=min_valid_pos, itr=itr)
          }
          else
          {
            xfit_left  <- min_valid_pos
            xfit_right <- 2*original_peak_pos-1
            if(first_peak_pos>=100 && itr<2)
            {
              xfit_left  <- original_peak_pos - 30
              xfit_right <- original_peak_pos + 30
            }else
            if(first_peak_pos>=30  && itr<2)
            {
              xfit_left  <- original_peak_pos - 10
              xfit_right <- original_peak_pos + 10
            }
            v<-optim(par=c(1, 1, 1, 0.8),
                     fn=error_minimize, x=x, end=end, xfit=xfit,
                     xfit_left=xfit_left, xfit_right=xfit_right,
                     d=error,
                     min_valid_pos=min_valid_pos, itr=itr)
          }
          cat('    Info: hom_xfit_left  for hom fitting at itr ', itr,  ': ', xfit_left, '\n')
          cat('    Info: hom_xfit_right for hom fitting at itr ', itr,  ': ', xfit_right, '\n')
          #
          # scale and correct yfit with optimized values in v$par
          xfit2     <- seq(0.01,maxright4fit,0.01)
          meanfit   <- mean(x[x>=xfit_left & x<=xfit_right])*v$par[3]
          if(meanfit < 0)
          {
            cat(paste("    Warning: data does not follow assumed distribution anymore at itr ",
                      itr, ", fitting stopped.\n", sep=""))
            normal  <- 0
            if(myyfit==0 && itr==1)
            {
              v<-optim(par=c(1, 1, 1),
                       fn=error_minimize3, x=x, end=end,
                       xfit_left=xfit_left, xfit_right=xfit_right,
                       d=error,
                       min_valid_pos=min_valid_pos, itr=itr)
              xfit2    <- seq(0.01,maxright4fit,0.01)
              yfit0    <- dnorm(xfit2,
                                mean=meanfit,
                                sd=sdfit)*(100000/itr)*v$par[2]*length(x)
              myyfit   <- myyfit + yfit0
              norm_fit_itr1 <- T
              if(het_observed && itr == 1)
              {
                homfit  <- yfit0
              }
            }
            break
          } # stop if the data does not follow normal distribution anymore
          sdfit     <- sd(x[  x>=xfit_left & x<=xfit_right])*v$par[1]
          xifit     <- 1*v$par[4]
          yfit0     <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/itr)*v$par[2]*length(x)
          yfit      <- yfit0[seq(100, length(yfit0), 100)]
          ## repeat first fitting if fitting is not satisfactory
          if(itr==1 && end<min(2*original_peak_pos-1, maxright4fit))
          {
            fitting_check <- error[1:maxright4fit, 2] - yfit
            n_less_than_0 <- 0
            for(fci in c((original_peak_pos+3):min(round(2*original_peak_pos), length(dr[,1]))))
            {
              if(fitting_check[fci] < 0)
              {
                n_less_than_0 <- n_less_than_0 + 1
              }
            }
            if(n_less_than_0 >= 3)
            {
              first_fit <- FALSE
              cat(paste("       Fitting has to be repeated: sizek ", sizek,
                        " at itr ", itr, "\n",
                        sep=""))
              i
              next
            }
          }
          # plot
          if(i == 1)
          {
            # plot(error[1:maxright4fit, 1], abs(error[1:maxright4fit, 2]),
            #      xlim=c(0, min(maxright4fit, 4*meanfit)), ylim=c(0, 1.5*max(error[xfit_left:maxright4fit, 2])),
            #      type='b', col="gray",
            #      xlab="K-mer Coverage", ylab="Genomewide k-mer Frequency")
            # title(main=paste('Example: fitting at itr ', itr, 'of sample ', sample, 'k=', sizek, ')'),
            #      cex.main=0.8)
            # lines(xfit2[100:length(xfit2)], yfit0[100:length(xfit2)], col="blue", lwd=2)
          }
          #### begin
          if(max(yfit) > 1.2*max(error[min_valid_pos:maxright4fit, 2]))
          {
            cat(paste("    Note on hom fitting: fitting stopped at iter ",
                      itr, ", expected: ", totalitr, "\n",sep=""))
            if(myyfit==0 && itr==1)
            {
              v<-optim(par=c(1, 1, 1),
                       fn=error_minimize3, x=x, end=end,
                       xfit_left=xfit_left, xfit_right=xfit_right,
                       d=error,
                       min_valid_pos=min_valid_pos, itr=itr)
              xfit2    <- seq(0.01,maxright4fit,0.01)
              yfit0    <- dnorm(xfit2,
                                mean=meanfit,
                                sd=sdfit)*(100000/itr)*v$par[2]*length(x)
              if(sum(yfit0) > 0)
              {
                myyfit   <- myyfit + yfit0
                norm_fit_itr1 <- T
                if(het_observed && itr == 1)
                {
                  homfit  <- yfit0
                }
              }else # abnormal case: caution!: no fitting performed on hom fitting
              {
                yfit0    <- rep(0, length(xfit2))
                myyfit   <- myyfit + yfit0
                norm_fit_itr1 <- T
                if(het_observed && itr == 1)
                {
                  homfit  <- yfit0
                }
              }
            }
            break;
          }
          #### end
          # estimate error
          error[1:maxright4fit, 2] <- abs(error[1:maxright4fit, 2] - yfit)
          # get more accurate fitted values
          xfit2    <- seq(0.01,maxright4fit,0.01);
          yfit0    <- dsnorm(xfit2,mean=meanfit,sd=sdfit, xi=xifit)*(100000/itr)*v$par[2]*length(x)
          myyfit   <- myyfit + yfit0
          # record fitted parameters
          mymeanfit[itr]  <- meanfit
          mysdfit[itr]    <- sdfit
          myxifit[itr]    <- xifit
          myscalefit[itr] <- (100000/itr)*v$par[2]*length(x)
          if(itr == 1)
          {
            first_mean     <- meanfit
            first_sd       <- sdfit
          }
          #
          ## het est: new feature
          if(het_observed && itr == 1)
          {
            homfit  <- yfit0
          }
          ## end of het rate est: new feature
          #
          itr <- itr + 1
        }
        cat("Iterative fitting done.\n\n")
        if(normal == 0)
        {
          # original count
          xlimmax <- min(3*first_peak_pos, selected)
          plot(dhet[, 1], abs(dhet[, 2]),
               xlim=c(0,xlimmax), ylim=c(0, 1.5*max(dhet[border_pos:xlimmax, 2])),
               type='b', col="gray",
               xlab="K-mer freq", ylab="Number of k-mers")
          title(main=paste('Sample ', sample, ' k=', sizek), cex.main=0.8) #
          # fitted count
          if(sum(hetfit[1:selected] > 0)) # plus het
          {
            # fitted peak
            fitted_peak_value <- which.max(myyfit)/100
            if(norm_fit_itr1 == T)
            {
              fitted_peak_value <- hom_peak_pos
            }
            fittedyvalues <- myyfit[seq(100, length(myyfit), 100)]

            if(only_one_hom_peak)
            {
              # begin find optinum het_fitting_delta 2016-12-20
              h_target <- c(hetfit[1:min(border_pos, het_min_valid_pos)],
                            as.vector(dhet[min(border_pos, het_min_valid_pos)+1:(length(dr[,2])-min(border_pos, het_min_valid_pos)),2]))
              v2 <- optim(par=c(0.11, 1),
                          fn=error_minimize2,
                          h_het=hetfit[1:first_peak_pos],
                          h_hom=fittedyvalues[1:first_peak_pos],
                          h_target=h_target[1:first_peak_pos])
              ## update delta
              het_fitting_delta <- v2$par[1]
              cat(paste("    Info: het_fitting_delta optimized as ", het_fitting_delta, "\n\n", sep=""))
              #cat(paste("    Info: hom_fitting_delta optimized as ", v2$par[2],         "\n\n", sep=""))
              # end   find optinum het_fitting_delta 2016-12-20
            }

            if(only_one_hom_peak)
            {
              fittedyvalues[1:selected] <- fittedyvalues[1:selected] + (1-het_fitting_delta)*hetfit[1:selected]
            }else
            {
              fittedyvalues[1:selected] <- fittedyvalues[1:selected] + hetfit[1:selected]
            }
            lines(1:length(fittedyvalues), fittedyvalues, col="blue", lwd=5)
            # het
            if(only_one_hom_peak)
            {
              lines(1:selected, (1-het_fitting_delta)*hetfit[1:selected], col="cyan", lwd=2, lty="dashed")
            }
            else
            {
              lines(1:selected,  hetfit[1:selected], col="cyan", lwd=2, lty="dashed")
            }
          }
          else
          {
            lines(xfit2[100:length(xfit2)], myyfit[100:length(xfit2)], col="blue", lwd=5)
            # fitted peak
            fitted_peak_value <- which.max(myyfit)/100
            fittedyvalues     <- myyfit[seq(100, length(myyfit), 100)]
          }
          # error rate of missing fitting: from raw k-mer counting dhet which includes het k-mers
          sum(fittedyvalues[border_pos:maxright4fit]-dhet[border_pos:maxright4fit,2])/sum(dhet[border_pos:maxright4fit,2])
          # fitted catting original with border_pos <- min_valid_pos; from raw k-mer counting dhet which includes het k-mers
          yfit2 <- c(fittedyvalues[1:min(border_pos, het_min_valid_pos)],
                     as.vector(dhet[min(border_pos, het_min_valid_pos)+1:(length(dr[,2])-min(border_pos, het_min_valid_pos)),2]))
          #
          lines(dhet[, 1], yfit2, col="red", lwd=1.5)
          #title(main=paste('Sample ', sample, ' k=', sizek, ')'), cex.main=0.8)
          #
          # calculate genome size (Mb) according to original data
          # with all kmers
          genome_size                 <- round(sum(dr[,1]*dhet[,2]/first_peak_pos))
          # with valid kmers: caution: 1.peak-value 2. valid kmer with cov >= min_valid_pos (or border_pos)
          genome_size_filtering_error <- round(sum(dr[border_pos:length(dhet[,1]),1]*dhet[border_pos:length(dr[,1]),2]/first_peak_pos))
          ##
          if(het_observed)
          {
            genome_size                 <- round(sum(dr[,1]*dhet[,2]/hom_peak_pos))
            genome_size_filtering_error <- round(sum(dr[border_pos:length(dhet[,1]),1]*dhet[border_pos:length(dr[,1]),2]/hom_peak_pos))
          }
          # genomic kmers
          # cat(paste("        err-exl genomic ", sizek, "mers: ", sum(dhet[border_pos:length(dhet[,2]),2]), " Mio\n", sep=""));
          # fitted data with fitted peak
          genome_size_fitted <- round(sum(c(1:length(fittedyvalues))*fittedyvalues[1:length(fittedyvalues)]/fitted_peak_value))
          # genomic kmers
          # cat(paste("        fit     genomic ", sizek, "mers: ", round(sum(fittedyvalues)), " Mio\n", sep=""));
          # catted data with fitted peak
          genome_size_catted <- round(sum(c(1:length(yfit2))*yfit2[1:length(yfit2)]/fitted_peak_value))
          # genomic kmers
          # cat(paste("        catted  genomic ", sizek, "mers: ", round(sum(yfit2)), " Mio\n", sep=""));
          # corrected genome size (not meaningful)
          genome_size_corrected  <- round(genome_size_catted/myxifit[1])
          # genome size: catted data + raw mean
          if(sum(hetfit[1:selected] > 0))
          {
            end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*2/4) # >7% error; k should be small
          }else
          if((sizek=="15" | sizek=="17") && myxifit[1]>=0.93 && myxifit[1]<=1.07) # normal case: reduce effect of copy-2 k-mers
          {
            end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*1/4)
          }else
            if(myxifit[1]  > 1.07)
            {
              end_for_mean <- round(fitted_peak_value*2 + fitted_peak_value*2/4)
            }else
              if(myxifit[1]  < 0.93)
              {
                end_for_mean <- round(fitted_peak_value*2 - fitted_peak_value*4/8)
              }else
              {
                end_for_mean <- round(fitted_peak_value*2)
              }
          #### caution: end_for_mean should be not larger than vector size
          end_for_mean <- min(end_for_mean, length(dr[,1]))       ####
          ####
          if(het_observed && only_one_het_peak) ## new ##
          {
            end_for_mean     <- round(het_peak_pos*4 - het_peak_pos*4/4)
            end_for_mean     <- min(end_for_mean, length(dr[,1])) ####
            #cat("end_for_mean (only one het peak): ", end_for_mean, '\n')
            xtmp             <- rep(1:end_for_mean,ceiling(abs(round(hetfit[1:end_for_mean]))/1000))
            first_mean_raw   <- 2*mean(xtmp[xtmp>=1 && xtmp<=end_for_mean])
          }else
            if(het_observed & main_peak_is_hom==F)
            {
              het_end_for_mean <- first_peak_pos
              het_end_for_mean <- min(het_end_for_mean, length(dr[,1])) ####
              xtmp             <- rep(1:het_end_for_mean,ceiling(abs(round(hetfit[1:het_end_for_mean]))/1000))
              first_mean_raw   <- 2*mean(xtmp[xtmp>=1 && xtmp<=het_end_for_mean])
            }else
              if(het_observed && only_one_hom_peak) # caution: expert suppl.
              {
                dtmp           <- ceiling(abs(round(yfit2[1:end_for_mean]-(1-het_fitting_delta)*hetfit[1:end_for_mean]))/1000) # 2016-08-25
                ## start of specific in v1.94: from [0.5*het_peak, 1.5*het_peak]
                offset_count   <- round((fittedyvalues[1:end_for_mean]-yfit2[1:end_for_mean])/1000)
                offl <- round(0.5*which.max(hetfit))
                offr <- round(1.5*which.max(hetfit))
                for (off_i in c(offl:offr))
                {
                  if(offset_count[off_i] > 0)
                  {
                    dtmp[off_i] <- dtmp[off_i] + offset_count[off_i]
                  }
                  else
                  {
                    dtmp[off_i] <- dtmp[off_i] - offset_count[off_i]
                  }
                }
                ## end   of specific in v1.94
                xtmp           <- rep(1:end_for_mean,dtmp)
                first_mean_raw <- mean(xtmp[xtmp>=1 && xtmp<=end_for_mean])
              }else
              {
                xtmp           <- rep(1:end_for_mean,ceiling(abs(round(yfit2[1:end_for_mean]-hetfit[1:end_for_mean]))/1000))
                first_mean_raw <- mean(xtmp[xtmp>=1 && xtmp<=end_for_mean]);
              }
          #
          genome_size_corrected2 <- round(sum(c(1:length(yfit2))*yfit2[1:length(yfit2)]/first_mean_raw))
          #
          #
          ## heterozygosity est: new feature - 2016-11-18
          if(het_observed)
          {
            # alpha = 1 / [ (O/E)*K + K ]
            # het
            het_end <- round(first_mean_raw)
            hom_end <- min(round(2*first_mean_raw), selected)
            if(only_one_hom_peak)
            {
              # only for humans: do the correction for those k-mers from X-Y
              het_correction <- 0
              if(species=="human")
              {
                het_correction <- (155270560 - 59373566)*first_mean_raw/2;
              }
              #sum(c(1:het_end)*(1-het_fitting_delta)*hetfit[1:het_end])
              alpha     <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/  sum(seq(1:het_end)*hetfit[1:het_end])/(1-het_fitting_delta)*as.numeric(sizek) + as.numeric(sizek) )
              alpha_cor <- 0
              alpha_cor <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/ (sum(seq(1:het_end)*hetfit[1:het_end])-het_correction)/(1-het_fitting_delta)*as.numeric(sizek) + as.numeric(sizek) )
            }
            else
            {
              # only for humans: do the correction for those k-mers from X-Y
              het_correction <- 0
              if(species=="human")
              {
                het_correction <- (155270560 - 59373566)*first_mean_raw/2;
              }
              #sum(c(1:het_end)*hetfit[1:het_end])
              alpha     <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/  sum(seq(1:het_end)*hetfit[1:het_end])*as.numeric(sizek) + as.numeric(sizek) )
              alpha_cor <- 0
              alpha_cor <- 1 / ( sum(seq(1:hom_end)*homfit[seq(100, length(homfit), 100)][1:hom_end])/ (sum(seq(1:het_end)*hetfit[1:het_end])-het_correction )*as.numeric(sizek) + as.numeric(sizek) )
            }
          }
          ## end of heterozygosity rate est: new feature - 2016-11-18
          #
          #
          # single copy region of genome
          singlestart = 1
          singleend   = min(round(2*first_mean_raw),length(dr[,1])-1) #### caution: out of bounds
          singlesize  = sum(dhet[singlestart:singleend, 1]*yfit2[singlestart:singleend]/first_mean_raw)
          # begin: new on 2016-08-08
          repsize  = round(sum(dhet[(singleend+1):length(yfit2), 1]*yfit2[(singleend+1):length(yfit2)]/first_mean_raw))
          #legend(45, 1.36*max(dr[border_pos:60, 2]),
          if(het_observed)
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'cyan', 'white', 'blue',
                              'white', 'white', 'white', 'red', 'white', 'white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)', sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', het_min_valid_pos,                 sep=""),
                              paste('Heterozygous fitting (intermediate info)'),
                              paste('* het rate: ',
                                    round(alpha, digits=8),  " (ori) : (cor) ",
                                    round(alpha_cor, digits = 8),                                sep=""),
                              paste('Fitted count with fitted k-mer_cov: ',
                                    round(genome_size_fitted/10^6, digits = 3), ' Mb',           sep=""),
                              paste('* k-mer_cov fit.: ',
                                    round(fitted_peak_value, digits=2),                          sep=""),
                              paste('* 1st-sd:',
                                    round(first_sd, digits=2),                                   sep=""),
                              paste('* 1st-skewness: ',
                                    round(myxifit[1], digits=2),                                 sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',       sep=""),
                              paste('* k-mer_cov cor.: ',
                                    round(first_mean_raw, digits=2),                             sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                      sep="")),
                   border="NA",
                   cex=0.7)
          }
          else
          {
            legend("topright",
                   fill   = c('gray','white', 'white', 'blue', 'white',
                              'white', 'white', 'red', 'white', 'white'),
                   legend = c(paste('Observed (obs.): ',
                                    round(genome_size/10^6, digits = 3), ' Mb (error-excluded: ',
                                    round(genome_size_filtering_error/10^6, digits = 3), ' Mb)',               sep=""),
                              paste('* k-mer_cov obs.: ',     first_peak_pos, sep=""),
                              paste('* signal_error_border ', border_pos,                     sep=""),
                              paste('Fitted count with fitted k-mer_cov: ',
                                    round(genome_size_fitted/10^6, digits = 3), ' Mb',                         sep=""),
                              paste('* k-mer_cov fit.: ',     round(fitted_peak_value, digits=2),              sep=""),
                              paste('* 1st-sd:',              round(first_sd, digits=2),                       sep=""),
                              paste('* 1st-skewness: ',       round(myxifit[1], digits=2),                     sep=""),
                              paste('Fitted+obs. with corrected k-mer_cov: ',
                                    round(genome_size_corrected2/10^6, digits = 3), ' Mb',                     sep=""),
                              paste('* k-mer_cov cor.: ',     round(first_mean_raw, digits=2),                 sep=""),
                              paste('* repetitive_ratio ',
                                    round(repsize/genome_size_corrected2, digits=2), ': ',
                                    round(repsize/10^6, digits = 3), ' Mb',                                    sep="")),
                   border="NA",
                   cex=0.7)
          }
          # end: new on 2016-08-08
          # genome sizes
          size_all[i]  <- genome_size
          size_exl[i]  <- genome_size_filtering_error
          size_fit[i]  <- genome_size_fitted
          size_cat[i]  <- genome_size_catted
          size_cor[i]  <- genome_size_corrected
          size_cor2[i] <- genome_size_corrected2
          # rbind(size_all, size_exl, size_fit, size_cat)
          # genome-wide average kmer coverage related info
          kmer_peak_cnt[i]  <- dr[first_peak_pos, 2]
          kmer_peak_pos[i]  <- first_peak_pos
          kmer_peak_pos2[i] <- fitted_peak_value
        }
        # labels for later plotting purpose
        labels[3*i-1]  <- paste("", sizek, sep="")
        namebank[i]    <- sizek
      }
      else
      {
        dev.off()
        stop(paste("file ", countingfile, " not found!",sep=""))
        quit("no")
      }
    }
    # summarize all kmer counting plots
    ylimmax <- 0
    for (sizek in targetsizek[1]) {
      histx<-read.table(countingfile)
      peaks2             <- findpeaks(histx$V2[3:selected])
      peaks2[,2:4]       <- peaks2[,2:4]+2
      if(1.5*max(peaks2[,1]) > ylimmax)
      {
        ylimmax <- 1.5*max(peaks2[,1])
      }
    }
    # colors
    colors <- rainbow(length(size_all))
    #
    ## plot comparision among gsize-estimations using different info
    genome_size_summary  <-rbind(size_all, size_exl, size_cat, size_fit, size_cor2)
    genome_size_summary2 <-rbind(size_exl, size_fit, size_cor2)
    cat(paste("Genome size estimate for ", sample, ": ", size_cor2, " bp.\n\n", sep=""))
    # caution
    if(length(which(is.na(genome_size_summary))) > 0)
    {
      cat("\n caution: some samples have NA predicted for observed genome size!\n")
      genome_size_summary[is.na(genome_size_summary)]   <- 0
    }
    if(length(which(is.na(genome_size_summary2))) > 0)
    {
      cat("\n caution: some samples have NA predicted for fitted   genome size!\n")
      genome_size_summary2[is.na(genome_size_summary2)] <- 0
    }
    if(genome_size_summary!=0 && genome_size_summary2!=0)
    {
      # record in file
      write.table(genome_size_summary,
                  paste(path, '/',prefix, vers, 'est.', sample, ".genome.size.estimated.k",
                        min(targetsizek), 'to', max(targetsizek), ".fitted.txt", sep=""),
                  quote = FALSE, row.names = TRUE, col.names=FALSE)
      # bar plot
      mingsize <- min(genome_size_summary2[!is.na(genome_size_summary2)])
      maxgsize <- max(genome_size_summary2[!is.na(genome_size_summary2)])
      colors <- grey((2:0)/3)
      mp <- barplot(genome_size_summary2, beside=TRUE, col=colors,
                    xlab=paste("Size k", sep=""), ylab="Genome size (bp)",
                    axes = FALSE, axisnames  = FALSE, border=NA,
                    ylim=c(0.5*max(genome_size_summary2),1.2*max(genome_size_summary2)),
                    xpd=FALSE)
      legend("topright",
             legend = c(paste('EXL: ', round(mean(genome_size_summary2[1,]/10^6), digits=3), 'Mb'),
                        paste('FIT: ', round(mean(genome_size_summary2[2,]/10^6), digits=3), 'Mb'),
                        paste('COR: ', round(mean(genome_size_summary2[3,]/10^6), digits=3), 'Mb')),
             fill = colors, cex=.8, border=NA)
      text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.8)
      axis(2, cex.axis=.8)
      title(main=paste("Genome size estimation by error-excluding, fitting, and correcting\n", sep=""))
      dev.off()
    }
  } # end of sample
  # output est. on het if asked.
  if(het_observed)
  {
    write(paste("Het_rate ", round(alpha, digits = 8), " ", round(alpha_cor, digits = 8), sep=""),
          file = paste(path, '/',prefix, vers, 'est.',
                       sample, ".genome.size.estimated.k",
                       min(targetsizek), 'to', max(targetsizek),
                       ".fitted.txt", sep=""),
          append = T,
          sep = " ")
  }
  # output est. on ratio of repeats
  write(paste("Est. ratio of repeats ", round(repsize/genome_size_corrected2, digits=8), sep=""),
        file = paste(path, '/',prefix, vers, 'est.',
                     sample, ".genome.size.estimated.k",
                     min(targetsizek), 'to', max(targetsizek),
                     ".fitted.txt", sep=""),
        append = T,
        sep = " ")
  # output est. on avg k-mer cov
  write(paste("Final k-mer cov ", round(first_mean_raw, digits=8), sep=""),
        file = paste(path, '/',prefix, vers, 'est.',
                     sample, ".genome.size.estimated.k",
                     min(targetsizek), 'to', max(targetsizek),
                     ".fitted.txt", sep=""),
        append = T,
        sep = " ")
  if(het_observed && only_one_hom_peak)
  {
    write(paste("het_fitting_delta ", round(het_fitting_delta, digits=8), sep=""),
          file = paste(path, '/',prefix, vers, 'est.',
                       sample, ".genome.size.estimated.k",
                       min(targetsizek), 'to', max(targetsizek),
                       ".fitted.txt", sep=""),
          append = T,
          sep = " ")
  }
  #
  end_t <- Sys.time()
  cat('Time consumed: ', format(end_t-start_t), '\n')
}
