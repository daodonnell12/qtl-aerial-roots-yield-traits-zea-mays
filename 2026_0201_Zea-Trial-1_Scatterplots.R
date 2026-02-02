#Zea-trial-1_Scatterplots
##31st January 2026, David O'Donnell

setwd("/Users/David_ODonnell/Documents/Dissertation/Publication/Publication_Excel-files/")
data <- read.csv("Zea-Trial-1_Metadata.csv", header=T, na.strings=c("NA"))

install.packages("ggplot2")
library("ggplot2")
install.packages("ggpubr")
library("ggpubr")
install.packages("ggpmisc")
library("ggpmisc")


data$NDFA_T1_B73_ref <- as.numeric(data$NDFA_T1_B73_ref)
data$NODES_w_AR <- as.numeric(data$NODES_w_AR)
data$HEIGHT_cm <- as.numeric(data$HEIGHT_cm)
data$DAYS_TO_POLLEN <- as.numeric(data$DAYS_TO_POLLEN)

#scatterplot AR-nodes vs. NDFA-T1

fit <- lm(data$NODES_w_AR ~ data$NDFA_T1_B73_ref) 
summary(fit)
summary(fit)$r.squared
r_squared_val <-summary(fit)$r.squared

ggplot(data = data, aes(x = NODES_w_AR, y = NDFA_T1_B73_ref, color = Group)) +
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
    axis.ticks = element_line(colour = "black", size = 0.5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(geom = "text", x = 7.5, y = 63,
           label = paste0("R-squared = ", round(r_squared_val, 3))) +
  labs(title = "Aerial Root Nodes vs. NDFA_T1",
       x = "Aerial Root Nodes", 
       y = "NDFA Time-point 1") +
  ggtitle("Aerial Root Nodes vs. NDFA Time-point 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 18)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 14), # Adjust legend item text size
        legend.title = element_text(size = 16)) +  # Adjust legend title size
  theme(legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
#  theme_minimal()


#scatterplot Height vs. AR-nodes

fit <- lm(data$HEIGHT_cm ~ data$NODES_w_AR) 
summary(fit)
summary(fit)$r.squared
r_squared_val <-summary(fit)$r.squared

ggplot(data = data, aes(x = HEIGHT_cm, y = NODES_w_AR, color = Group)) +
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(geom = "text", x = 400, y = 4.4,
           label = paste0("R-squared = ", round(r_squared_val, 3))) +
  labs(title = "Height vs. Aerial Root Nodes",
       x = "Height (cm)", 
       y = "Aerial Root Nodes") +
  ggtitle("Height vs. Aerial Root Nodes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 18)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 14), # Adjust legend item text size
        legend.title = element_text(size = 16)) +  # Adjust legend title size
  theme(legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
#  theme_minimal()


#scatterplot DPS vs. AR-nodes

fit <- lm(data$DAYS_TO_POLLEN ~ data$NODES_w_AR) 
summary(fit)
summary(fit)$r.squared
r_squared_val <-summary(fit)$r.squared

ggplot(data = data, aes(x = DAYS_TO_POLLEN, y = NODES_w_AR, color = Group)) +
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(geom = "text", x = 120, y = 3.5,
           label = paste0("R-squared = ", round(r_squared_val, 3))) +
  labs(title = "Days to Pollen Shed vs. Aerial Root Nodes",
       x = "Days to Pollen Shed", 
       y = "Aerial Root Nodes") +
  ggtitle("Days to Pollen Shed vs. Aerial Root Nodes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 18)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 14), # Adjust legend item text size
        legend.title = element_text(size = 16)) +  # Adjust legend title size
  theme(legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
#  theme_minimal()

#scatterplot DPS vs. Height

fit <- lm(data$DAYS_TO_POLLEN ~ data$HEIGHT_cm) 
summary(fit)
summary(fit)$r.squared
r_squared_val <-summary(fit)$r.squared

ggplot(data = data, aes(x = DAYS_TO_POLLEN, y = HEIGHT_cm, color = Group)) +
  theme(axis.text.x = element_text(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(geom = "text", x = 97, y = 210,
           label = paste0("R-squared = ", round(r_squared_val, 3))) +
  labs(title = "Days to Pollen Shed vs. Height",
       x = "Days to Pollen Shed", 
       y = "Height (cm)") +
  ggtitle("Days to Pollen Shed vs. Height") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size = 18)) +
  theme(plot.title = element_text(face = "bold")) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(face = "bold")) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.y = element_text(face = "bold")) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 14), # Adjust legend item text size
        legend.title = element_text(size = 16)) +  # Adjust legend title size
  theme(legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
  #theme_minimal()
