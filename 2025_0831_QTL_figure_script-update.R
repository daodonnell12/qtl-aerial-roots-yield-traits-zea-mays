### David O'Donnell
### Figure for QTL and Elevated Mexicana Introgression in Totontepec Maize

dat <- structure(list(chromosome = 1:10, size = c(300927100L, 236965704L, 230537786L, 240523061L, 216206258L, 168818752L, 176024588L, 175698433L, 155721069L, 149097751L)), .Names = c("chromosome","size"), class = "data.frame", row.names = c(NA, -10L))
marks <- structure(list(Chromosome = c(1L, 3L, 9L, 3L, 8L, 10L, 1L, 9L, 9L, 1L, 2L, 4L, 10L, 1L), 
                        Position = c(262937140L, 16566244L, 108296473L, 211685092L, 123844193L, 
                                     108117991L, 294813085L, 102497421L, 139410267L, 61021047L, 
                                     5106835L, 237313637L, 147956241L, 4996244L), 
                        Type = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L), .Label = c("A", "B"), class = "factor")), 
                   .Names = c("Chromosome", "Position", "Type"), class = "data.frame", row.names = c(NA,-10L))

#New Code from https://stackoverflow.com/questions/33727432/how-to-plot-positions-along-a-chromosome-graphic
#Load Data

library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot, you can replace with `geom_text` if you want
library("scales") # for axis labels notation
library("ggpattern") #for introgression
library("png")

# install.packages("remotes")
#remotes::install_github("coolbutuseless/ggpattern")

# insert your steps to load data from tabular files or other sources here; 
# dummy datasets taken directly from files shown in this example

# data with the QTL for the sample
#sample_cns <- structure(list(gene = c("SC_1_LH82", "SC_3_B73", "SC_9_MO17", "DPS_3_MO17", 
#                                      "DPS_8_MO17", "DPS_10_B73", "AR_1_LH82", "AR_9_LH82", 
#                                      "AR_9_MO17", "PTNP_1_MO17", "PTN_2_B73", "PDM_4_MO17", "PDM_10_LH82", 
#                                      "GTN_GDM_1_MO17"), 
#                                      chromosome = c("chr1", "chr3", "chr9", "chr3", "chr8", 
#                                      "chr10", "chr1", "chr9", "chr9", "chr1", "chr2", "chr4", "chr10", "chr1"), 
#                                      start = c(258353240L, 10683505, 22782253L, 209024681L, 
#                                                         118971709L, 100252518L, 293329598L, 33670979L,
#                                                         131569767L, 43756650L, 4990077L, 67143858L, 147811484L, 
#                                                         3408430L), end = c(268813141L, 162587460L, 119839322L, 
#                                                          219411033L, 143492696L, 128967696L, 295309095L, 
#                                                          111431608L, 143804017L, 71127215L, 5565686L, 237313637L, 148581023L, 
#                                                          5813154L), cn = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
#                                                                                      1L, 1L, 1L, 1L, 1L, 1L, 1L), 
#                                                          Trait = c("SC", "SC", "SC", "DPS", "DPS", "DPS", "AR", 
#                                                          "AR", "AR", "PTNP", "PTN", "PDM", "PDM", "GTN_GDM")), 
#                                                          .Names = c("gene", "chromosome", "start", "end", "cn", "Trait"), 
#                                                          row.names = c(NA, 14L), class = "data.frame")

sample_cns <- structure(list(gene = c("5.04, 17.9%", "3.88, 17.1%", "5.65, 15.0%", "4.35, 18.1%", 
                                      "6.23, 44.1%", "5.50, 28.0%", "6.25, 18.4%", "4.27, 18.9%", 
                                      "4.20, 21.6%", "4.44, 29.7%", "6.30, 32.3%", "3.91, 19.0%", "4.07, 14.8%", 
                                      "4.32, 32.80%"), 
                                      chromosome = c("chr1", "chr3", "chr9", "chr3", "chr8", 
                                      "chr10", "chr1", "chr9", "chr9", "chr1", "chr2", "chr4", "chr10", "chr1"), 
                                      start = c(258353240L, 10683505, 22782253L, 209024681L, 
                                       118971709L, 100252518L, 293329598L, 33670979L,
                                       131569767L, 43756650L, 4990077L, 67143858L, 147811484L, 
                                       3408430L), end = c(268813141L, 162587460L, 119839322L, 
                                                          219411033L, 143492696L, 128967696L, 295309095L, 
                                                          111431608L, 143804017L, 71127215L, 5565686L, 237313637L, 148581023L, 
                                                          5813154L), cn = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
                                                                            1L, 1L, 1L, 1L, 1L, 1L, 1L), 
                                      Key = c("QTL_SC", "QTL_SC", "QTL_SC", "QTL_DPS", "QTL_DPS", "QTL_DPS", "QTL_AR", 
                                       "QTL_AR", "QTL_AR", "QTL_PTNP", "QTL_PTN", "QTL_PDM", "QTL_PDM", "QTL_GTN_GDM")), 
                                      .Names = c("gene", "chromosome", "start", "end", "cn", "Key"), 
                                       row.names = c(NA, 14L), class = "data.frame")

#     QTL_Name	chromosome	    start	      end	cn	Key   LOD   %Variance
# 1	  SC_1_LH82	      chr1	258353240	268813141	1	  SC      5.04  17.86
# 2	  SC_3_B73	      chr3	10683505	162587460	1	  SC      3.88  17.06
# 3	  SC_9_MO17	      chr9	22782253	119839322	1	  SC      5.65  15.03
# 4	  DPS_1_MO17	    chr3	209024681	219411033	1	  DPS     4.35  18.06
# 5	  DPS_8_MO17	    chr8	118971709	143492696	1	  DPS     6.23  44.13
# 6	  DPS_10_B73	    chr10	100252518	128967696	1	  DPS     5.50  28.00
# 7	  AR_1_LH82	      chr1	293329598	295309095	1	  AR      6.25  18.44
# 8   AR_9_LH82	      chr9	33670979	111431608	1	  AR      4.27  18.89
# 9	  AR_9_MO17	      chr9	131569767	143804017	1	  AR      4.20  21.61
# 10	PTNP_1_MO17	    chr1	43756650	71127215	1	  PTNP    4.44  29.72
# 11	PTN_2_B73	      chr2	4990077	  5565686	  1	  PTN     6.30  32.32
# 12	PDM_4_MO17	    chr4	67143858	237313637	1	  PDM     3.91  19.00
# 13	PDM_10_LH82	    chr10	147811484	148581023	1	  PDM     4.07  14.80
# 14	GTN_GDM_1_MO17	chr1	3408430	  5813154	  1	  GTN_GDM 4.32  32.80

# hg19 chromosome sizes
chrom_sizes <- structure(list(chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                                             "chr8", "chr9", "chr10"), size = c(300927100L, 236965704L, 230537786L, 
                                              241322382L, 216206258L, 168818752L, 176024588L, 175698433L, 155721069L, 
                                              149097751L)), .Names = c("chromosome", "size"), class = "data.frame", 
                                              row.names = c(NA, -10L))

# > head(chrom_sizes)
#   chromosome      size
# 1       chr1  300927100
# 2       chr2  236965704
# 3       chr3  230537786
# 4       chr4  241322382
# 5       chr5  216206258
# 6       chr6  168818752
# 7       chr7  176024588
# 8       chr8  175698433
# 9       chr9  155721069
# 10      chr10 149097751

# hg19 centromere locations
centromeres <- structure(list(chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                                             "chr9", "chr10"), start = c(133300000L, 89300000L, 94600000L, 
                                              104200000L, 101600000L, 49800000L, 55300000L, 45900000L, 68600000L, 
                                              59300000L), end = c(133900000L, 91100000L, 95400000L, 105000000L, 
                                              108600000L, 50400000L, 55700000L, 48000000L, 69200000L, 60700000L)), 
                                              .Names = c("chromosome", "start", "end"), class = "data.frame", 
                                              row.names = c(NA, -10L))

# > head(centromeres)
#   chromosome     start       end
# 1       chr1 133300000 133900000
# 2       chr2  89300000  91100000
# 3       chr3  94600000  95400000
# 4       chr4 104200000 105000000
# 5       chr5 101600000 108600000
# 6       chr6  49800000  50400000
# 7       chr7  55300000  55700000
# 8       chr8  45900000  48000000
# 9       chr9  68600000  69200000
# 10      chr10 59300000  60700000

#introgression <- structure(list(chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
#                                             "chr9", "chr10"), start = c(150000000L, 80000000L, 50000000L, 
#                                                                         50000000L, 50000000L, 50000000L, 50000000L, 50000000L, 50000000L, 
#                                                                         50000000L), end = c(170000000L, 110000000L, 70000000L, 70000000L, 
#                                                                                             70000000L, 70000000L, 70000000L, 70000000L, 70000000L, 70000000L)), 
#                         .Names = c("chromosome", "start", "end"), class = "data.frame", 
#                         row.names = c(NA, -10L))

#overlapping introgression regions
introgression <- structure(list(chromosome = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3", 
                                               "chr3", "chr3", "chr4", "chr4", "chr4", "chr4", "chr4", "chr4", "chr4", "chr4", "chr4", "chr5", "chr5", "chr5", "chr5", "chr5", "chr6", "chr6", "chr6", "chr6", "chr6", "chr6", 
                                               "chr6", "chr6", "chr7", "chr7", "chr7", "chr7", "chr7", "chr7", "chr7", "chr8", "chr8", "chr8", "chr8", "chr8", "chr9", "chr9", "chr9", "chr9", "chr10", "chr10", "chr10", "chr10"), 
                                                start = c(82702375L,	108621763L,	119003899L, 219231030L, 10232L, 46524503L, 69180911L, 83540212L,	109878087L,	118981561L,	127959743L,	132951761L, 144194202L, 158408759L,	176451542L,	
                                                          35750177L, 38602012L, 44095158L, 49562562L,	63108209L,	112098658L,	129291770L,	136908028L,	162297619L,	41851034L,	49413636L,	62218127L,	67020744L,	124954610L,	131700301L,	
                                                          168739151L,	168753601L, 240524176L,	32919099L,	37077825L, 89317871L, 127604425L,	139234251L,	9490668L,	22664410L,	36090447L,	40875092L,	75398920L,	124531652L,	125187777L,	
                                                          133816320L,	2596362L,	30156469L, 45117728L,	49020436L,	79101130L,	92725739L,	146548837L,	55711353L,	75988259L,	107089606L,	128973724L,	136415654L,	32946981L,	63841178L,	
                                                          81425959L,	153736443L,	27541115L, 33605592L, 108468473L,	121412033L), 
                                                end = c(84633600L,	112639745L,	148628650L,	220785740L,	1268910L,	46524503L,	78500359L,	105782820L,	117919800L,	126763017L,	130956103L,	140661606L,	148193360L,	163727058L,	177874858L,	
                                                        35750177L,	41403869L,	47574787L,	53675536L,	66312613L,	114790935L,	132483085L,	139736412L,	165022441L,	44504375L,	56312802L,	66732444L,	71704923L,	130571862L,	142582826L,	
                                                        168739151L,	168976868L,	241322382L,	34801447L,	39454408L,	126448883L,	136200439L,	142711487L,	18991117L,	25581860L,	38746722L,	74599737L,	78099789L,	124709894L,	128308620L,	
                                                        138187940L,	2596362L,	36558700L,	45117728L,	75749010L,	84124234L,	99955636L,	148488081L,	62665056L,	82560840L,	109391995L,	131113149L,	138213242L,	42289811L,	69060908L,	
                                                        87873575L,	153736443L,	31492794L,	56757623L,	112586972L,	123376961L), 
                           cn = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                  1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
                           Key = c("Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", 
                                   "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana", "Introgression_mexicana")), 
#.Names = c("gene", "chromosome", "start", "end", "cn", "Key"), 
#row.names = c(NA, 14L), class = "data.frame")
                           .Names = c("chromosome", "start", "end", "cn", "Key"), 
                           row.names = c(NA, -66L), class = "data.frame")

#Adjust Data
# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr10", "chr9", "chr8", "chr7", "chr6", "chr5", "chr4", 
                 "chr3", "chr2", "chr1")
chrom_key <- setNames(object = as.character(c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)), 
                      nm = chrom_order)

#if want to reverse chromosome order for figure, use below command
#chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)
sample_cns[["chromosome"]] <- factor(x = sample_cns[["chromosome"]], 
                                     levels = chrom_order)
centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)
introgression[["chromosome"]] <- factor(x = introgression[["chromosome"]], 
                                      levels = chrom_order)
# create a color key for the plot
group.colors <- c(QTL_SC = "purple", QTL_DPS = "orange", QTL_AR = "green", QTL_PTNP = "blue", QTL_PTN = "sky blue", QTL_PDM = "dark blue", QTL_GTN_GDM = "brown", Introgression_mexicana = "red")

#Make Plot
ggplot(data = chrom_sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                xmax = as.numeric(chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
  
  # add bands for centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.35, 
                                    xmax = as.numeric(chromosome) + 0.35, 
                                    ymax = end, ymin = start)) +
  #add bands for introgression regions
  geom_rect(data=introgression,aes(xmin = as.numeric(chromosome) - 0.08, 
                                   xmax = as.numeric(chromosome) + 0.08, 
#                                   ymax = end, ymin = start), fill='red',color='red') +
                                   ymax = end, ymin = start, fill = Key), alpha=0.5) +
  
  # add bands for Trait value
  geom_rect(data = sample_cns, aes(xmin = as.numeric(chromosome) - 0.2, 
                                   xmax = as.numeric(chromosome) + 0.2, 
                                   ymax = end, ymin = start, fill = Key), alpha=0.5) + 
  scale_fill_manual(values = group.colors) +
  
  # add 'SC' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "SC"), 
                  aes(x = chromosome, y = start, label = gene), nudge_x = c(-0.35, -0.35, 0.35),
                  color = "black", show.legend = FALSE) +
  # add 'DPS' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "DPS"), 
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'AR' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "AR"), 
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'PTNP' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "PTNP"), 
#                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(0.35),
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'PTN' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "PTN"), 
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'PDM' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "PDM"), 
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'GTN_GDM' gene markers
  geom_text(data = subset(sample_cns, sample_cns$Key == "GTN_GDM"), 
#                  aes(x = chromosome, y = end, label = gene ), nudge_x = c(-0.35),
                  aes(x = chromosome, y = start, label = gene ), nudge_x = c(-0.35),
                  color = "black", show.legend = FALSE) +
  # add 'Mexicana_introgression' gene markers

#  ggtitle("Quantitative Trait Loci & Introgression Regions") + theme(plot.title = element_text(hjust = 0.5)) + 
#        theme(plot.title = element_text(size=16)) +
  
  # supress scientific notation on the y-axis
  scale_y_continuous(labels = comma) +
  ylab("region (bp)") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12), # Adjust legend item text size
        legend.title = element_text(size = 14)) # Adjust legend title size

#total genes in maize genome
totalgenesmaize <- 1:109871
#total genes in QTL regions
totalgenesqtl <- 1:26474
#sample total genes in Mex to Tpec introgression regions (2455 genes) from maize genome 1000 times
randomsample <- sample(totalgenesmaize, size = 2455)
#determine number of times (in 1000) that overlap between introgression genes and QTL genes is greater than actual overlap (741 genes)
sum(randomsample < 26475)

simfunc <- function(){
  randomsample <- sample(totalgenesmaize, size = 2455)
  total <- sum(randomsample < 26475)
  print(total)
}
vector <- replicate(1000, simfunc())
hist(vector)
