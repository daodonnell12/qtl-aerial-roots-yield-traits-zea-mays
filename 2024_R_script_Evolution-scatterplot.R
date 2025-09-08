setwd("~/Documents/Dissertation/Publication/Publication_Excel-files/")
Evolution <- read.csv("2023_Evolution_Metadata_Averages_Final_for-R-figure.csv")
Evolution
print(Evolution)

plot(x = Evolution$NAME,NDFA_T1_B73.ref_Mean = Evolution$NDFA_T1_B73.ref_Mean, 
     xlab = "Entry", 
     ylab = "NDFA_T1_B73.ref_Mean", 
     main = "Plot"
)

#Syntax:
  
  #plot(x,y,main,xlab,ylab,sub,asp)

#Parameters:
  
  #x:-the x coordinates of points in the plot
#y:-the y coordinates of points in the plot
#main:-an overall title for the plot
#sub:-a subtitle for the plot
#xlab:-a title for the x-axis
#ylab:-a title for the y-axis
#asp:-the y/x aspect ratio
#Return:
  
  #Scatter plot of the given x and y values.
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("gridExtra")
install.packages("grid")
install.packages("cowplot")
install.packages("gtable")
install.packages("gridtext")
library("tidyverse")
library("ggplot2")
library("gridExtra")
library("cowplot")
library("gtable")
library("gridtext")
library("grid")

Zea_Trial_1 <- read.csv("2023_Evolution_Metadata_Averages_Final_for-R-figure.csv")

#plot_NDFA_T1_B73 <- ggplot(Zea_Trial_1, aes(Zea_Trial_1$Group, Zea_Trial_1$NDFA_T1_B73.ref_Mean), scale="globalminmax") +
  #geom_vline(xintercept = 0, linetype = 2) +
  #geom_hline(yintercept = 0, linetype = 2) +
  #geom_point() +
  #theme_minimal()

#plot_NDFA_T1_B73

element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

#ggplot(Zea_Trial_1,aes(x=Group,y=NDFA_T1_B73.ref_Mean,col=Group))+geom_point()

plot_NDFA_T1_B73_1 <- ggplot(Zea_Trial_1,aes(x=Name,y=NDFA_T1_B73.ref_Mean,col=Group))+geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin=NDFA_T1_B73.ref_Mean-NDFA_T1_B73.ref_SE, ymax=NDFA_T1_B73.ref_Mean+NDFA_T1_B73.ref_SE), width=.1, linetype = 2,
                position=position_dodge(.9)) +
  geom_point(shape = 1) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(y= "Percent NDFA", x = "Group")

plot_NDFA_T1_B73_1
legend = cowplot::get_plot_component(plot_NDFA_T1_B73_1, 'guide-box-right', return_all = TRUE)
legend <- cowplot::ggdraw(legend)
legend

plot_Trait <- ggplot(Zea_Trial_1,aes(x=Name,y=Trait,shape=Trait))+geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  labs(y= "Trait", x = "Group")

plot_Trait
legend2 = cowplot::get_plot_component(plot_Trait, 'guide-box-right', return_all = TRUE)
legend2 <- cowplot::ggdraw(legend2)
legend2

plot_NDFA_T1_B73 <- ggplot(Zea_Trial_1,aes(x=Name,y=NDFA_T1_B73.ref_Mean,col=Group))+geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin=NDFA_T1_B73.ref_Mean-NDFA_T1_B73.ref_SE, ymax=NDFA_T1_B73.ref_Mean+NDFA_T1_B73.ref_SE), width=.1, linetype = 2,
                position=position_dodge(.9)) +
  #coord_fixed(ratio = .1) +
  geom_point(shape = 1) +
  ylim(-36,120) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  #labs(y = NULL, x = "Group") +
  labs(y = NULL, x = NULL) +
  panel_border(color = "black") +
  theme(panel.grid.major.x = element_blank()) +
  theme(strip.background = element_rect(fill = "gray80")) +
  #ggtitle("Percent NDFA") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "% NDFA") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +
  theme(legend.position = "none")

plot_NDFA_T1_B73

#plot_AR_Nodes_Mean <- ggplot(Zea_Trial_1, aes(Group, AR_Nodes_Mean), scale="globalminmax") +
  #geom_vline(xintercept = 0, linetype = 2) +
  #geom_hline(yintercept = 0, linetype = 2) +
  #geom_point() +
  #theme_minimal()

#plot_AR_Nodes_Mean

#ggplot(Zea_Trial_1,aes(x=Group,y=AR_Nodes_Mean,col=Group, jitter(x, factor = 1, amount = NULL)))+geom_point()

plot_AR_Nodes_Mean <- ggplot(Zea_Trial_1,aes(x=Name,y=AR_Nodes_Mean,col=Group))+geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin=AR_Nodes_Mean-AR_Nodes_StErr, ymax=AR_Nodes_Mean+AR_Nodes_StErr), width=.1, linetype = 2,
                position=position_dodge(.9)) +
  #coord_fixed(ratio = 1.95) +
  geom_point(shape = 2) +
  ylim(-1.8,6) +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  #labs(y = NULL, x = "Group") +
  labs(y = NULL, x = NULL) +
  panel_border(color = "black") +
  theme(panel.grid.major.x = element_blank()) +
  theme(strip.background = element_rect(fill = "gray80")) +
  #ggtitle("AR Nodes") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "AR Nodes") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  ) +
  theme(legend.position = "none")

plot_AR_Nodes_Mean

#Group <- textGrob("Group", gp = gpar(fontsize = 14), hjust = 1.5, vjust = -10)
Group <- textGrob("Group", gp = gpar(fontsize = 14), hjust = 2.7, vjust = 0.2)
Mean <- textGrob("Mean", gp = gpar(fontsize = 14), rot=90, vjust = 1)

#plot_grid(plot_NDFA_T1_B73, plot_AR_Nodes_Mean, labels=c("A", "B"), ncol = 2, nrow = 1)
grid.arrange(plot_NDFA_T1_B73, plot_AR_Nodes_Mean, legend, legend2, ncol = 4, widths=c(2.3, 2.3, 0.8, 0.8), nrow = 1, bottom = Group, left = Mean)
#grid.arrange(plot_NDFA_T1_B73, legend, plot_AR_Nodes_Mean, ncol = 2, widths=c(2.3, 0.8), nrow = 2, heights=c(10, 10))
#plots <- grid.arrange(plot_NDFA_T1_B73, plot_AR_Nodes_Mean, legend, ncol = 3, widths=c(2.3, 2.3, 0.8), nrow = 1)
#plots

