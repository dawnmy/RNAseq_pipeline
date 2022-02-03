library(ggtern)

## Run this script only after run the DE script
pm_group_logcpm5 <- as.data.frame(cpmByGroup(y)[y$AveLogCPM >5,])

ggtern(data = data.frame(st20 = pm_group_logcpm5$`20_st`,
                         st26 = pm_group_logcpm5$`26_st`,
                         st30 = pm_group_logcpm5$`30_st`,
                         logCPM = y$AveLogCPM[y$AveLogCPM >5]),
               aes(st20, st26, st30, size=logCPM)) + 
  geom_point(alpha=0.4) +
  theme_rgbw(base_size = 15) 
  labs(title = "Expression distribution among 20, 26 and 30 at stationary phase")    +
  guides(color = "none", fill = "none", alpha = "none")
