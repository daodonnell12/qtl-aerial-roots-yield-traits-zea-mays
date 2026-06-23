#load stats
  setwd("C:/Users/natas/Documents/cl_hub/QTL_overlap/test")
  stats <- read.table("overlap_stats.txt", fill = TRUE)
#load 
  list=c("AR_simulated_overlap.txt","GDM_simulated_overlap.txt","PDM_simulated_overlap.txt","SC_simulated_overlap.txt",
         "DPS_simulated_overlap.txt","GTN_simulated_overlap.txt","PTNP_simulated_overlap.txt")
  #head(list)
  
  for (file in list) {
    QTL <- strsplit(file, "[_.]")[[1]][1]
    assign(QTL, scan(file, what = numeric()))
  }

#AR
  m <- stats$V3[1]
  sstdev <- stats$V4[1]
  observed_AR <- stats$V1[1]

png("QTL_AR_simulated_histogram.png", width = 1200, height = 1000, res = 150)
  {
    hist(AR,
         xlim = c(0, 0.7),
         ylim = c(0, 600),
         xlab = "Simulated overlap",
         ylab = "frequency",
         main = "QTL_AR: simulated values with Mean and Stdev",
         col = "grey")
    
    abline(v = m, col = "red", lwd = 2, lty = 1)
    abline(v = m + sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = observed, col = "blue", lwd = 2, lty = 1)
    
    legend("topright",
           legend = c(paste("mean_simulated:", round(m, 6)),
                      paste("sstdev_simulated:", round(sstdev, 6)),
                      paste("observed_mean:", observed)),
           col = c("red", "red", "blue"),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 1))
  }
  dev.off()
 
 #DPS
  m <- stats$V3[2]
  sstdev <- stats$V4[2]
  observed_DPS <- stats$V1[2]
  
  png("QTL_DPS_simulated_histogram.png", width = 1200, height = 1000, res = 150)
  {
    hist(DPS,
         xlim = c(0, 0.3),
         ylim = c(0, 1000),
         xlab = "Simulated overlap",
         ylab = "frequency",
         main = "QTL_DPS: simulated values with Mean and Stdev",
         col = "grey")
    
    abline(v = m, col = "red", lwd = 2, lty = 1)
    abline(v = m + sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = observed, col = "blue", lwd = 2, lty = 1)
    
    legend("topright",
           legend = c(paste("mean_simulated:", round(m, 6)),
                      paste("sstdev_simulated:", round(sstdev, 6)),
                      paste("observed_mean:", observed)),
           col = c("red", "red", "blue"),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 1))
  }
  dev.off()
  
  sum(DPS > observed)/1000
  
  #PDM
  m <- stats$V3[5]
  sstdev <- stats$V4[5]
  observed_PDM <- stats$V1[5]
  
  png("QTL_PDM_simulated_histogram.png", width = 1200, height = 1000, res = 150)
  {
    hist(PDM,
         xlim = c(0, 0.6),
         ylim = c(0, 300),
         xlab = "Simulated overlap",
         ylab = "frequency",
         main = "QTL_PDM: simulated values with Mean and Stdev",
         col = "grey")
    
    abline(v = m, col = "red", lwd = 2, lty = 1)
    abline(v = m + sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = m - sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = observed, col = "blue", lwd = 2, lty = 1)
    
    legend("topright",
           legend = c(paste("mean_simulated:", round(m, 6)),
                      paste("sstdev_simulated:", round(sstdev, 6)),
                      paste("observed_mean:", observed)),
           col = c("red", "red", "blue"),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 1))
  }
  dev.off()
  
  sum(PDM > observed)/1000
  
  #SC
  m <- stats$V3[7]
  sstdev <- stats$V4[7]
  observed_SC <- stats$V1[7]
  
  png("QTL_SC_simulated_histogram.png", width = 1200, height = 1000, res = 150)
  {
    hist(SC,
         xlim = c(0, 0.8),
         ylim = c(0, 400),
         xlab = "Simulated overlap",
         ylab = "frequency",
         main = "QTL_SC: simulated values with Mean and Stdev",
         col = "grey")
    
    abline(v = m, col = "red", lwd = 2, lty = 1)
    abline(v = m + sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = m - sstdev, col = "red", lwd = 2, lty = 2)
    abline(v = observed, col = "blue", lwd = 2, lty = 1)
    
    legend("topright",
           legend = c(paste("mean_simulated:", round(m, 6)),
                      paste("sstdev_simulated:", round(sstdev, 6)),
                      paste("observed_mean:", observed)),
           col = c("red", "red", "blue"),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 1))
  }
  dev.off()
  
  sum(SC > observed)/1000
  
 
 ###########QTL joint boxplot
  
  QTL <- cbind(AR, DPS, PDM, SC)
  
  group.colors <- c(AR = "green", DPS = "orange", PDM = "darkblue", SC = "purple")
  sim_means <- c(mean(AR), mean(DPS), mean(PDM), mean(SC))
  obs_vals <- c(observed_AR, observed_DPS, observed_PDM, observed_SC)
  groups <- c("AR", "DPS", "PDM", "SC")
  
  png("boxplot.png", width = 1200, height = 1100, res = 150)
  {
  boxplot(QTL, 
          #main = "QTL", 
          ylab = "overlap",
          ylim = c(0,.8),
          col = alpha(group.colors, 0.5))
  points(1:4,obs_vals,
         pch = 18, col = "red", cex = 2)
  
  legend("topright",
         legend = paste(
           c("AR", "DPS", "PDM", "SC"),
           ": ",
           "sim_m=", round(sim_means, 3),",",
           "obs_m=", round(obs_vals, 3)
         ),
         fill = alpha(group.colors, 0.5),
         border = "black",
         bty = "n")
  }
  dev.off()                 
  