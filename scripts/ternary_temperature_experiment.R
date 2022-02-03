library(ggtern)

## Run this script only after run the DE script
library(ggtern)

pm_logcpm5 <- cpm(y)[y$AveLogCPM >5,]


ggtern(data = data.frame(st20 = rowMeans(pm_logcpm5[,4:7]),
                                 st26 = rowMeans(pm_logcpm5[,10:12]),
                                 st30 = rowMeans(pm_logcpm5[,16:18]),
                                 logCPM = y$AveLogCPM[y$AveLogCPM >5]),
               aes(st20, st26, st30, size=logCPM)) + 
  geom_point(alpha=0.4) +
  theme_rgbw(base_size = 15) 
  labs(title = "Expression distribution among 20, 26 and 30 at stationary phase")    +
  guides(color = "none", fill = "none", alpha = "none")
