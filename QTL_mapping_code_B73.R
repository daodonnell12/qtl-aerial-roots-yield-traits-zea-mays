#Input library
library(qtl)

#Input file
setwd("~/Documents/Dissertation/Dissertation_Figures/Dissertation_Figures_Ch4/Genotypes_Phenotypes/2019NewFiles")
mapthis <- read.cross("csv", "./", "genotypes_phenotypes_B73_052320.csv", estimate.map=FALSE,genotypes=c("A","H","B"))

#summarize cross data
summary(mapthis)

#plot total missing data
plotMissing(mapthis)

#plot missing by marker, by individual
par(mfrow=c(1,1), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

#drop marker that missing data
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 170])
mapthis <- drop.markers(mapthis, todrop)

#look for duplicate individuals
par(mfrow=c(1,1), las=1)
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

#Remove highly duplicate individuals
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh
g <- pull.geno(mapthis)
for(i in 1:nrow(wh)) {
  tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
  mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA}

#Identify duplicate markers
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))

#Example to delete duplicate markers
mapthis <- drop.markers(mapthis, c("S1_258521836","S1_258353240","S1_260771404","S1_258875033","S1_256514018","S1_202158540","S1_182197270","S1_113972428","S1_121493056","S1_107793132","S1_164951739","S1_105951834","S1_93908120","S1_83523459","S1_83010833","S1_62676068","S1_52161960","S1_23557922","S1_8541287"))
mapthis <- drop.markers(mapthis, c("S5_8814769","S5_17977491","S5_20987360","S5_24609559","S5_25490067","S5_29809787","S5_33980476","S5_51355494","S5_51667851","S5_86824716","S5_132657078","S5_135734277","S5_144176697","S5_156969781","S5_158546122","S5_168261823","S5_169579338","S5_178476407","S5_179270099","S5_186678634","S5_193652783","S5_193734961"))
mapthis <- drop.markers(mapthis, c("S4_40348572","S4_43589697","S4_74696529","S4_158862960","S4_169150576","S4_174011015","S4_199730034","S4_199798785","S4_214783304","S4_225505497","S4_226006473","S4_233442688","S4_235268683"))
mapthis <- drop.markers(mapthis, c("S6_38271909","S6_44824453","S6_56115016","S6_74033234","S6_85588171","S6_86346265","S6_124661162","S6_134469173","S6_138562829","S6_145646257","S6_147897417","S6_148828165","S6_148960226","S6_164366833"))
mapthis <- drop.markers(mapthis, c("S3_40065288","S3_48102165","S3_56900402","S3_123868401","S3_86874194","S3_114432547","S3_127101202","S3_137240207","S3_166690697","S3_197383588","S3_198214924","S3_219373746"))
mapthis <- drop.markers(mapthis, c("S2_12352551","S2_163825099","S2_200456487","S2_207041195","S2_209506613","S2_209536855","S2_217650053"))
mapthis <- drop.markers(mapthis, c("S7_45385946","S7_46479272","S7_85263428","S7_94239335","S7_95431664","S7_100265669","S7_123505354","S7_126448104","S7_135833150","S7_139724303","S7_141010338","S7_141142142","S7_152181399","S7_155745191","S7_165246730"))
mapthis <- drop.markers(mapthis, c("S9_19286687","S9_33441018","S9_43791221","S9_47044003","S9_47214192","S9_47452489","S9_55347766","S9_56143183","S9_83965115","S9_84224710","S9_57833394","S9_60025656","S9_60038259","S9_60085674","S9_61814159","S9_73473889","S9_77309892","S9_82461294","S9_86864186","S9_87812663","S9_88651683","S9_95741461","S9_102497421","S9_135392726","S9_140353844","S9_140773445"))
mapthis <- drop.markers(mapthis, c("S8_33137479","S8_55020828","S8_55934247","S8_83046893","S8_120149919","S8_121895690","S8_132879659"))
mapthis <- drop.markers(mapthis, c("S10_23805682","S10_54170353","S10_30141313","S10_46801060","S10_52223901","S10_60402962","S10_70921408","S10_90243326","S10_100683137","S10_129586074","S10_136961591","S10_137286265","S10_137830927","S10_138027028","S10_138215309"))
  
#check genotype frequencies by individual
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))

#check for segregation distortion
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]

#drop worst segregating markers
todrop <- rownames(gt[gt$P.value < 1e-4,])
mapthis <- drop.markers(mapthis, todrop)
#example mapthis <- drop.markers(mapthis, c("BGSNP-220"))

#estimate recombination fraction
mapthis <- est.rf(mapthis)

#Plot LOD score vs. recombination fraction
par(mfrow=c(1,1), las=1)
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#cleanGeno(cross, maxdist=2.5, maxmark=2, verbose=TRUE)

#Form linkage groups (checking only)
lg <- formLinkageGroups(mapthis, max.rf=0.3, min.lod=5.6)
table(lg[,2])

#Form linkage groups in actual data
mapthis <- formLinkageGroups(mapthis, max.rf=0.3, min.lod=5.6, reorgMarkers=TRUE)

#check linkage map
pull.map(mapthis,chr=1)
pull.map(mapthis,chr=2)
pull.map(mapthis,chr=3)
pull.map(mapthis,chr=4)
pull.map(mapthis,chr=5)
pull.map(mapthis,chr=6)
pull.map(mapthis,chr=7)
pull.map(mapthis,chr=8)
pull.map(mapthis,chr=9)
pull.map(mapthis,chr=10)
pull.map(mapthis,chr=11)
pull.map(mapthis,chr=12)
pull.map(mapthis,chr=13)
pull.map(mapthis,chr=14)
pull.map(mapthis,chr=15)
pull.map(mapthis,chr=16)


mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 9) 

mn <- markernames(mapthis, chr=c(12)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 8)

mn <- markernames(mapthis, chr=c(13,14,15,16)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)

#mapthis <- drop.markers(mapthis, c("S8_834154","S8_2288877"))
#mapthis <- drop.markers(mapthis, c("S10_148513275","S10_148533441","S10_146165305"))

#move markers
mn <- markernames(mapthis, chr=c(2)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(5)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 2)
mn <- markernames(mapthis, chr=c(3)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 12)
mn <- markernames(mapthis, chr=c(4)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 3)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 4)
mn <- markernames(mapthis, chr=c(12)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 5)
mn <- markernames(mapthis, chr=c(6)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 6)
mn <- markernames(mapthis, chr=c(7)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 7)
mn <- markernames(mapthis, chr=c(8)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 8)
mn <- markernames(mapthis, chr=c(9)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 9)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)

#order markers in linkage groups
mapthis <- orderMarkers(mapthis, chr=10)
mapthis <- orderMarkers(mapthis, chr=9)
mapthis <- orderMarkers(mapthis, chr=8)
mapthis <- orderMarkers(mapthis, chr=7)
mapthis <- orderMarkers(mapthis, chr=6)
mapthis <- orderMarkers(mapthis, chr=5)
mapthis <- orderMarkers(mapthis, chr=4)
mapthis <- orderMarkers(mapthis, chr=3)
mapthis <- orderMarkers(mapthis, chr=2)
mapthis <- orderMarkers(mapthis, chr=1)

#check linkage map
pull.map(mapthis,chr=1)
pull.map(mapthis,chr=2)
pull.map(mapthis,chr=3)
pull.map(mapthis,chr=4)
pull.map(mapthis,chr=5)
pull.map(mapthis,chr=6)
pull.map(mapthis,chr=7)
pull.map(mapthis,chr=8)
pull.map(mapthis,chr=9)
pull.map(mapthis,chr=10)

#Switch order of markers 
#mapthis <- switch.order(mapthis, chr=1, c(102,100,99,101,98:93,91,92,90:55,53:51,54,50:10,8,9,7:5,1,3,4,2, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=1, c(1,	4,	3,	5,	2,	102,	101,	6,	7,	100,	8:11,	99,	12:14,	98,	15,	97,	16,	17,	96,	95,	93,	94,	18:20,	92:90,	21,	89,	88,	22,	87,	23,	24,	86,	25:29,	85,	30,	32,	33,	84,	31,	83,	34:36,	82,	37,	39,	38,	81:79,	40,	78:75,	41,	42,	74,	43,	73,	72,	44,	45,	71,	46,	70,	69,	68,	47:50,	67,	51,	66,	52,	65,	53,	54,	64:60,	55,	59,	57,	56,	58, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=1, c(102, 	100, 	99, 	101, 	98, 	97, 	96, 	95, 	94, 	93, 	91, 	92, 	90, 	89, 	88, 	87, 	86, 	85, 	84, 	83, 	82, 	81, 	80, 	79, 	78, 	77, 	76, 	75, 	74, 	73, 	72, 	71, 	70, 	69, 	68, 	67, 	66, 	65, 	64, 	63, 	62, 	61, 	60, 	59, 	58, 	57, 	56, 	55, 	53, 	52, 	51, 	54, 	50, 	49, 	48, 	47, 	46, 	45, 	44, 	43, 	42, 	41, 	40, 	39, 	38, 	37, 	36, 	35, 	34, 	33, 	32, 	31, 	30, 	29, 	28, 	27, 	26, 	25, 	24, 	23, 	22, 	21, 	20, 	19, 	18, 	17, 	16, 	15, 	14, 	13, 	12, 	11, 	10, 	8, 	9, 	7, 	6, 	5, 	1, 	3, 	4, 	2,  error.prob=0.005))
mapthis <- switch.order(mapthis, chr=2, c(2, 	1, 	3, 	5, 	4, 	6, 	61, 	7, 	60, 	8, 	9, 	10, 	59, 	11, 	12, 	58, 	57, 	56, 	55, 	54, 	53, 	52, 	51, 	50, 	13, 	14, 	49, 	15, 	16, 	48, 	17, 	18, 	19, 	20, 	47, 	21, 	46, 	45, 	22, 	44, 	23, 	43, 	42, 	41, 	40, 	39, 	38, 	37, 	36, 	24, 	25, 	26, 	27, 	28, 	29, 	30, 	31, 	33, 	32, 	34, 	35, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=3, c(61,62,1:4,60,59,57,58,55,56,54,5,6,53,50,52,51,49,48,46,47,7:9,45:43,10,42,41,11,12,40,39,13,38,36,37,35,14,15,34,33,16,32,17,18,31:29,19,28,20,27,26,21,22,24,25,23, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=3, c(1:12,	27,	26,	14,	13,	25,	15:22,	24,	23,	28:57,	60,	61,	59,	58,	62, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=3, c(61, 	62, 	1, 	2, 	3, 	4, 	5, 	6, 	7, 	8, 	9, 	10, 	60, 	59, 	12, 	11, 	58, 	13, 	14, 	57, 	56, 	55, 	54, 	53, 	52, 	15, 	51, 	16, 	17, 	18, 	19, 	50, 	49, 	20, 	21, 	22, 	48, 	47, 	45, 	46, 	44, 	43, 	23, 	42, 	41, 	40, 	39, 	38, 	37, 	36, 	35, 	34, 	33, 	32, 	31, 	30, 	29, 	28, 	27, 	25, 	24, 	26, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=4, c(67:64,1,63:57,2,56,3:21,55:52,23,22,51:49,47,48,25,24,46,26:33,36:34,38,37,39:45, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=4, c(1,	67,	2:4,	66:60,	5,	59,	6:8,	58:56,	9,	55,	10,	54:52,	11:17,	51:49,	18,	47,	48,	19,	46,	45,	43,	44,	41,	42,	20,	40:35,	21,	22,	34,	33,	23,	31,	32,	24:26,	30,	29,	27,	28, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=4, c(67, 	66, 	65, 	64, 	1, 	63, 	62, 	61, 	60, 	59, 	58, 	57, 	56, 	2, 	3, 	4, 	5, 	6, 	7, 	55, 	8, 	9, 	10, 	11, 	12, 	13, 	14, 	15, 	16, 	17, 	54, 	53, 	18, 	19, 	20, 	21, 	22, 	52, 	23, 	51, 	24, 	25, 	27, 	26, 	29, 	28, 	50, 	49, 	48, 	30, 	31, 	32, 	47, 	46, 	45, 	33, 	34, 	44, 	36, 	35, 	43, 	42, 	41, 	37, 	38, 	40, 	39, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=5, c(2,3,1,5,4,6:12,67,66,13:18,65,64,19,63,20,62,21,61,60,22:31,59,32,58,33,57,34,35,56:53,36:41,52,42:45,51:49,46,48,47, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=5, c(2,	3,	1,	4:12,	65,	64,	66,	67,	13:15,	61:63,	16,	60,	17,	59,	18,	58,	19,	57,	20,	56,	21,	55,	22,	54,	23,	24,	53,	25,	52:47,	26:29,	46,	30,	45,	44,	31,	32,	43:41,	33,	34,	40:38,	35,	37,	36, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=5, c(2, 	3, 	1, 	5, 	4, 	6, 	7, 	8, 	9, 	10, 	11, 	12, 	15, 	16, 	14, 	13, 	20, 	21, 	22, 	19, 	17, 	18, 	23, 	24, 	25, 	26, 	27, 	28, 	29, 	30, 	31, 	32, 	33, 	34, 	35, 	36, 	37, 	38, 	39, 	40, 	41, 	42, 	43, 	44, 	45, 	46, 	47, 	48, 	49, 	52, 	50, 	51, 	53, 	54, 	55, 	56, 	57, 	58, 	59, 	60, 	61, 	62, 	63, 	64, 	67, 	65, 	66, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=6, c(1,61,2,60,3,59:51,4,50,5,6,49,48,7:10,47:36,11:14,35,15,34:32,16,31:23,17,22,21,18,20,19, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=6, c(2,	1,	3:5,	61:57,	6,	56:52,	7,	51,	8,	9,	50,	10,	49:43,	11,	42,	12,	41,	13,	40,	14,	39:34,	15,	33,	32,	16,	31:29,	17,	28:25,	18,	24,	19,	20,	23,	21,	22, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=6, c(2, 	1, 	3, 	4, 	5, 	61, 	60, 	59, 	58, 	57, 	56, 	55, 	54, 	53, 	6, 	7, 	52, 	51, 	8, 	50, 	49, 	9, 	10, 	11, 	12, 	13, 	14, 	15, 	16, 	48, 	17, 	47, 	46, 	45, 	18, 	44, 	43, 	42, 	41, 	40, 	39, 	19, 	38, 	20, 	21, 	37, 	36, 	35, 	22, 	34, 	33, 	32, 	23, 	31, 	30, 	29, 	24, 	25, 	28, 	26, 	27, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=7, c(1,52,53,3,2,4:7,51,49,50,48,47,8,46:43,9,42,10,41,40,11:15,39:36,16,35:33,17,32,18,19,31:28,20,27,26,21:23,25,24, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=7, c(1,	52,	53,	3,	2,	4:7,	51,	50,	8,	49,	48,	9,	47,	10,	11,	46:35,	12:15,	34,	16,	33,	32,	17:21,	31,	22,	30:25,	23,	24, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=7, c(53, 	52, 	1, 	50, 	51, 	49, 	48, 	47, 	46, 	2, 	45, 	44, 	3, 	4, 	43, 	5, 	42, 	6, 	7, 	41, 	40, 	8, 	9, 	10, 	39, 	38, 	37, 	36, 	35, 	11, 	12, 	34, 	13, 	33, 	14, 	32, 	15, 	16, 	17, 	18, 	31, 	19, 	30, 	20, 	21, 	22, 	23, 	24, 	25, 	26, 	27, 	29, 	28, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=8, c(1:50, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=9, c(4,3,1,2,5:11,15:19,12:14,20:39, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=9, c(36,	37,	39,	38,	35,	33,	34,	32:29,	27:24,	28,	23:1, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=9, c(4, 	3, 	1, 	2, 	5, 	6, 	7, 	8, 	9, 	10, 	11, 	12, 	34, 	35, 	36, 	13, 	33, 	32, 	31, 	30, 	14, 	29, 	28, 	15, 	27, 	16, 	26, 	17, 	18, 	19, 	20, 	21, 	22, 	25, 	23, 	24, 	37, 	38, 	39, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=10, c(1,3,2,4:6,8,9,7,10,35:27,11:15,26,16,25,17,24,18,23:19,36,38,37, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=10, c(4,	6,	5,	7,	38,	37,	35,	36,	8:13,	34,	33,	14,	15,	32,	31,	16:19,	30,	29,	20,	28,	21,	27,	22:26,	3,	1,	2, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=10, c(1, 	3, 	2, 	4, 	5, 	6, 	7, 	9, 	8, 	10, 	35, 	34, 	33, 	32, 	11, 	12, 	13, 	31, 	14, 	30, 	29, 	28, 	27, 	26, 	15, 	24, 	25, 	16, 	17, 	23, 	18, 	19, 	20, 	21, 	22, 	36, 	38, 	37, error.prob=0.005))

#plot linkage map
plotMap(mapthis)
plotMap(mapthis, show.marker.names=F)

#plot heat map
plotRF(mapthis,alternate.chrid=T)
plotRF(mapthis, alternate.chrid=T, chr=1)

#lod threshold via permutation
#'arg' should be one of “em”, “imp”, “hk”, “ehk”, “mr”, “mr-imp”, “mr-argmax”
operm.hk <-scanone(mapthis,method="hk",n.perm=1000)
summary(operm.hk,alpha=0.05)
operm.em <-scanone(mapthis,method="em",n.perm=1000)
summary(operm.em,alpha=0.05)
operm.mr <-scanone(mapthis,method="mr",n.perm=1000)
summary(operm.mr,alpha=0.05)

#single analysis
out.mr <- scanone(mapthis, method="mr",pheno.col=40)
summary(out.mr, threshold=3)
max(out.mr)
out.mr[ out.mr$chr == 2, ]
out.mr[ out.mr$chr == 10, ]

#plot single analysis
plot(out.mr, show.marker.names=F)
plot(out.mr, show.marker.names=F, chr=c("2"))
plot(out.mr, show.marker.names=F, chr=c("10"))

#SIM
out.em <- scanone (mapthis, pheno.col=c(42),method="em")
summary(out.em, threshold=3)
max(out.em)
out.em[ out.em$chr == 2, ]
out.em[ out.em$chr == 10, ]
plot(out.em, lodcolumn=c(1), lty=1, col=c("blue"), show.marker.names=F)
plot(out.em, lodcolumn=c(1),chr=c("7"), lty=1, col=c("Blue"), show.marker.names=F)
plot(out.em, lodcolumn=c(1),chr=c("10"), lty=1, col=c("Blue"), show.marker.names=F)

#CIM
out_cim.em <- cim(mapthis, pheno.col = (42), n.marcovar=1, method=c("em"), map.function=c("kosambi") )
summary(out_cim.em, threshold=3)
max(out_cim.em)
out_cim.em[ out_cim.em$chr == 7, ]
plot(out_cim.em, show.marker.names=F, col =c("red"))
plot(out_cim.em, chr=c("7"), show.marker.names=F, col =c("red"))
attr(out_cim.em, "marker.covar")
attr(out_cim.em, "marker.covar.pos")

out_cim.em[ out_cim.em$chr == 1, ]
out_cim.em[ out_cim.em$chr == 2, ]
out_cim.em[ out_cim.em$chr == 3, ]
out_cim.em[ out_cim.em$chr == 4, ]
out_cim.em[ out_cim.em$chr == 5, ]
out_cim.em[ out_cim.em$chr == 6, ]
out_cim.em[ out_cim.em$chr == 7, ]
out_cim.em[ out_cim.em$chr == 8, ]
out_cim.em[ out_cim.em$chr == 9, ]
out_cim.em[ out_cim.em$chr == 10, ]

#combine
plot(out.mr, out.em, out_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
plot(out.mr, out.em, out_cim.em, chr="2",show.marker.names=F,col =c("Black", "blue", "red"))
plot(out.mr, out.em, out_cim.em, chr="7",show.marker.names=F,col =c("Black", "blue", "red"))
plot(out.mr, out.em, out_cim.em, chr="10",show.marker.names=F,col =c("Black", "blue", "red"))

out.mr <- scanone(mapthis, method="mr",pheno.col=42)
summary(out.mr, threshold=2)
out.em <- scanone (mapthis, pheno.col= (42),method="em")
summary(out.em, threshold=2)
out_cim.em <- cim(mapthis, pheno.col = (42), n.marcovar=3, method=c("em"), map.function=c("kosambi") )
summary(out_cim.em, threshold=2)

# combine CIM chart

#CIM all phenotypes
#SC
out1_cim.em <- cim(mapthis, pheno.col = (1), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out2_cim.em <- cim(mapthis, pheno.col = (2), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out3_cim.em <- cim(mapthis, pheno.col = (3), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out4_cim.em <- cim(mapthis, pheno.col = (4), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out5_cim.em <- cim(mapthis, pheno.col = (5), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out5_cim.em.perm <- cim(mapthis, pheno.col = (5), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 10) #replace 1000 with 10 for expediency
summary(out5_cim.em.perm)
summary(out5_cim.em, threshold=3.89)
out6_cim.em <- cim(mapthis, pheno.col = (6), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out6_cim.em.perm <- cim(mapthis, pheno.col = (6), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 1000) #replace 1000 with 10 for expediency
summary(out6_cim.em.perm)
#summary(out6_cim.em, threshold=3.89)
summary(out6_cim.em, threshold=4.07) #mcov=1
summary(out6_cim.em, threshold=4.23) #mcov=3

#DT
#out7_cim.em <- cim(mapthis, pheno.col = (7), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out8_cim.em <- cim(mapthis, pheno.col = (8), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out9_cim.em <- cim(mapthis, pheno.col = (9), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out10_cim.em <- cim(mapthis, pheno.col = (10), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out11_cim.em <- cim(mapthis, pheno.col = (11), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out12_cim.em <- cim(mapthis, pheno.col = (12), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#DPS
#out13_cim.em <- cim(mapthis, pheno.col = (13), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out14_cim.em <- cim(mapthis, pheno.col = (14), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out15_cim.em <- cim(mapthis, pheno.col = (15), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out16_cim.em <- cim(mapthis, pheno.col = (16), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out17_cim.em <- cim(mapthis, pheno.col = (17), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out17_cim.em.perm <- cim(mapthis, pheno.col = (17), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 10) #replace 1000 with 10 for expediency
summary(out17_cim.em.perm)
summary(out17_cim.em, threshold=3.89)
out18_cim.em <- cim(mapthis, pheno.col = (18), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out18_cim.em.perm <- cim(mapthis, pheno.col = (18), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 1000) #replace 1000 with 10 for expediency
summary(out18_cim.em.perm)
summary(out18_cim.em, threshold=3.89)
summary(out18_cim.em, threshold=4.07) #mcovar=1
summary(out18_cim.em, threshold=6.19) #mcovar=3

#15NT1
#out19_cim.em <- cim(mapthis, pheno.col = (19), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out20_cim.em <- cim(mapthis, pheno.col = (20), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out21_cim.em <- cim(mapthis, pheno.col = (21), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out22_cim.em <- cim(mapthis, pheno.col = (22), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out23_cim.em <- cim(mapthis, pheno.col = (23), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out24_cim.em <- cim(mapthis, pheno.col = (24), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#15NT2
#out25_cim.em <- cim(mapthis, pheno.col = (25), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out26_cim.em <- cim(mapthis, pheno.col = (26), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out27_cim.em <- cim(mapthis, pheno.col = (27), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out28_cim.em <- cim(mapthis, pheno.col = (28), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out29_cim.em <- cim(mapthis, pheno.col = (29), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out30_cim.em <- cim(mapthis, pheno.col = (30), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#15NT3
#out31_cim.em <- cim(mapthis, pheno.col = (31), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out32_cim.em <- cim(mapthis, pheno.col = (32), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out33_cim.em <- cim(mapthis, pheno.col = (33), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out34_cim.em <- cim(mapthis, pheno.col = (34), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out35_cim.em <- cim(mapthis, pheno.col = (35), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out36_cim.em <- cim(mapthis, pheno.col = (36), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#AR
#out37_cim.em <- cim(mapthis, pheno.col = (37), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out38_cim.em <- cim(mapthis, pheno.col = (38), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out39_cim.em <- cim(mapthis, pheno.col = (39), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out40_cim.em <- cim(mapthis, pheno.col = (40), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out41_cim.em <- cim(mapthis, pheno.col = (41), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out42_cim.em <- cim(mapthis, pheno.col = (42), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#PDM
#out43_cim.em <- cim(mapthis, pheno.col = (43), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out44_cim.em <- cim(mapthis, pheno.col = (44), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out45_cim.em <- cim(mapthis, pheno.col = (45), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out46_cim.em <- cim(mapthis, pheno.col = (46), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out47_cim.em <- cim(mapthis, pheno.col = (47), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out48_cim.em <- cim(mapthis, pheno.col = (48), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#PTNP
#out49_cim.em <- cim(mapthis, pheno.col = (49), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out50_cim.em <- cim(mapthis, pheno.col = (50), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out51_cim.em <- cim(mapthis, pheno.col = (51), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out52_cim.em <- cim(mapthis, pheno.col = (52), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out53_cim.em <- cim(mapthis, pheno.col = (53), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out54_cim.em <- cim(mapthis, pheno.col = (54), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#PTN
#out55_cim.em <- cim(mapthis, pheno.col = (55), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out56_cim.em <- cim(mapthis, pheno.col = (56), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out57_cim.em <- cim(mapthis, pheno.col = (57), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out58_cim.em <- cim(mapthis, pheno.col = (58), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out59_cim.em <- cim(mapthis, pheno.col = (59), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out59_cim.em.perm <- cim(mapthis, pheno.col = (59), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 10) #replace 1000 with 10 for expediency
summary(out59_cim.em.perm)
summary(out59_cim.em, threshold=3.89)
out60_cim.em <- cim(mapthis, pheno.col = (60), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out60_cim.em.perm <- cim(mapthis, pheno.col = (60), n.marcovar=1, method=c("em"), map.function=c("kosambi"), n.perm = 1000) #replace 1000 with 10 for expediency
summary(out60_cim.em.perm)
summary(out60_cim.em, threshold=3.89)
summary(out60_cim.em, threshold=3.91) #mcovar=1
summary(out60_cim.em, threshold=5.77) #mcovar=3

#GDM
#out61_cim.em <- cim(mapthis, pheno.col = (61), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out62_cim.em <- cim(mapthis, pheno.col = (62), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out63_cim.em <- cim(mapthis, pheno.col = (63), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out64_cim.em <- cim(mapthis, pheno.col = (64), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out65_cim.em <- cim(mapthis, pheno.col = (65), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out66_cim.em <- cim(mapthis, pheno.col = (66), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#GTNT
#out67_cim.em <- cim(mapthis, pheno.col = (67), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out68_cim.em <- cim(mapthis, pheno.col = (68), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out69_cim.em <- cim(mapthis, pheno.col = (69), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out70_cim.em <- cim(mapthis, pheno.col = (70), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out71_cim.em <- cim(mapthis, pheno.col = (71), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out72_cim.em <- cim(mapthis, pheno.col = (72), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#GTN
#out73_cim.em <- cim(mapthis, pheno.col = (73), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out74_cim.em <- cim(mapthis, pheno.col = (74), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out75_cim.em <- cim(mapthis, pheno.col = (75), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
#out76_cim.em <- cim(mapthis, pheno.col = (76), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out77_cim.em <- cim(mapthis, pheno.col = (77), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out78_cim.em <- cim(mapthis, pheno.col = (78), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#plot CIM stand count
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count")
plot(out5_cim.em, show.marker.names=F, col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - SC (T x B73)")
plot(out6_cim.em, show.marker.names=F, col =c("purple"))
#abline(4.07, 0, untf = FALSE) #mcov=1
abline(4.23, 0, untf = FALSE) #mcov=3
#abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - SC (T x B73)")

#plot CIM SC chrs
#plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
#abline(3.89, 0, untf = FALSE)
#title(main = "B73 - Stand Count - Chr 1")
#plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
#abline(3.89, 0, untf = FALSE)
#title(main = "B73 - Stand Count - Chr 1")
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Stand Count - Chr 2")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count - Chr 2")
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="3", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count - Chr 3")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="3", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count - Chr 3")
plot(out5_cim.em, show.marker.names=F,chr="3", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Stand Count - Chr 3")
#plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="6", col =c("Black", "blue", "red"))
#abline(3.89, 0, untf = FALSE)
#title(main = "B73 - Stand Count - Chr 6")
#plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="6", col =c("brown", "purple", "green"))
#abline(3.89, 0, untf = FALSE)
#title(main = "B73 - Stand Count - Chr 6")

summary(out1_cim.em, threshold=3.89)
summary(out2_cim.em, threshold=3.89)
summary(out3_cim.em, threshold=3.89)
summary(out4_cim.em, threshold=3.89)
summary(out5_cim.em, threshold=3.89)
summary(out6_cim.em, threshold=3.89)


#QTL effect plot for SC blupNPFC
hyper <- calc.genoprob(mapthis, step = 1)
out <- scanone(hyper, method="hk")
marker_name <- find.marker(hyper, chr=3, pos=148.6)
effectplot(hyper, mname1 = marker_name, pheno.col=6, main=paste("QTL Effect at", marker_name), ylab=("Stand Counts"))

hyper <- calc.genoprob(mapthis, step = 1)
out <- scanone(hyper, method="hk")
marker_name <- find.marker(hyper, chr=2, pos=121)
effectplot(hyper, mname1 = marker_name, pheno.col=6, main=paste("QTL Effect at", marker_name), ylab=("Stand Counts"))


#Regression SC

#blup no parents and family as covariates
summary(out6_cim.em, threshold=3.89)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_blupNPFC_chr2 <- makeqtl(mapthis, chr=c(2), pos=c(121)) 
qtl_SC_blupNPFC_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(149)) 
qtl_SC_blupNPFC <- makeqtl(mapthis, chr=c(2, 3), pos=c(121, 149))
summary(out.qtl_SC_blupNPFC_chr2 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr2, pheno.col=6, formula=y~Q1))
summary(out.qtl_SC_blupNPFC_chr3 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr3, pheno.col=6, formula=y~Q1))
summary(out.qtl_SC_blupNPFC <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC, pheno.col=6, formula=y~Q1+Q2))
plot(qtl_SC_blupNPFC_chr2)
plot(qtl_SC_blupNPFC_chr3)
plot(qtl_SC_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_SC_blupNPFC_chr2 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr2, formula=y~Q1)
summary(out.fq_SC_blupNPFC_chr2)
out.fq_SC_blupNPFC_chr3 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr3, formula=y~Q1)
summary(out.fq_SC_blupNPFC_chr3)
out.fq_SC_blupNPFC <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC, formula=y~Q1+Q2)
summary(out.fq_SC_blupNPFC)

#plot CIM DT
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence")
plot(out11_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence")

#plot CIM DT chrs
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="5", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Tassel Emergence - Chr 5")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="5", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Tassel Emergence - Chr 5")
plot(out11_cim.em, show.marker.names=F,chr="5", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Tassel Emergence - Chr 5")
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence - Chr 10")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence - Chr 10")
plot(out11_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Tassel Emergence - Chr 10")

summary(out7_cim.em, threshold=3.89)
summary(out8_cim.em, threshold=3.89)
summary(out9_cim.em, threshold=3.89)
summary(out10_cim.em, threshold=3.89)
summary(out11_cim.em, threshold=3.89)
summary(out12_cim.em, threshold=3.89)

#Regression DT

#blup no parents and family as covariate
summary(out12_cim.em, threshold=3.89)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_blupNPFC <- makeqtl(mapthis, chr=c(10), pos=c(104)) 
summary(out.qtl_DT_blupNPFC <- fitqtl(mapthis, qtl=qtl_DT_blupNPFC, pheno.col=12, formula=y~Q1))
plot(qtl_DT_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_DT_blupNPFC <- fitqtl(mapthis, qtl=qtl_DT_blupNPFC, formula=y~Q1)
summary(out.fq_DT_blupNPFC)


#plot CIM DPS
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Pollen Shed")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Pollen Shed")
plot(out17_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
#abline(4.07, 0, unft = FALSE) #mcovar=1
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - DPS (T x B73)")
plot(out18_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
#abline(4.07, 0, unft = FALSE) #mcovar=1
abline(6.19, 0, unft = FALSE) #mcovar=3
#abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - DPS (T x B73)")

#plot CIM DPS chrs
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Pollen Shed - Chr 1")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Pollen Shed - Chr 1")
plot(out17_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Days to Pollen Shed - Chr 1")
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Pollen Shed - Chr 10")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Pollen Shed - Chr 10")
plot(out17_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Days to Pollen Shed - Chr 10")

summary(out13_cim.em, threshold=3.89)
summary(out14_cim.em, threshold=3.89)
summary(out15_cim.em, threshold=3.89)
summary(out16_cim.em, threshold=3.89)
summary(out17_cim.em, threshold=3.89)
summary(out18_cim.em, threshold=3.89)

#QTL effect plot for DPS blupNPFC
hyper <- calc.genoprob(mapthis, step = 1)
out <- scanone(hyper, method="hk")
marker_name <- find.marker(hyper, chr=10, pos=104.1)
effectplot(hyper, mname1 = marker_name, pheno.col=18, main=paste("QTL Effect at", marker_name), ylab = "Days to Pollen Shed")

#Regression DPS

#blup no parents and family as covariate
summary(out18_cim.em, threshold=3.89)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_blupNPFC <- makeqtl(mapthis, chr=c(10), pos=c(104)) 
summary(out.qtl_DPS_blupNPFC <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC, pheno.col=18, formula=y~Q1))
plot(qtl_DPS_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_DPS_blupNPFC <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC, formula=y~Q1)
summary(out.fq_DPS_blupNPFC)


#plot CIM 15NT1
plot(out19_cim.em, out20_cim.em, out21_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 1")
plot(out22_cim.em, out23_cim.em, out24_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 1")
plot(out23_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 1")

#plot CIM 15NT1 chrs
plot(out19_cim.em, out20_cim.em, out21_cim.em, show.marker.names=F,chr="4", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 1 - Chr4")
plot(out22_cim.em, out23_cim.em, out24_cim.em, show.marker.names=F,chr="4", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 1 - Chr4")

summary(out19_cim.em, threshold=3.89)
summary(out20_cim.em, threshold=3.89)
summary(out21_cim.em, threshold=3.89)
summary(out22_cim.em, threshold=3.89)
summary(out23_cim.em, threshold=3.89)
summary(out24_cim.em, threshold=3.89)

#plot CIM 15NT2
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - 15N Time-point 2")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2")
plot(out29_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2")

#plot CIM 15NT2 chrs
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2 - Chr 1")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2 - Chr 1")
plot(out29_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2 - Chr 1")
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,chr="5", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - 15N Time-point 2 - Chr 5")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,chr="5", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2 - Chr 5")
plot(out29_cim.em, show.marker.names=F,chr="5", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 2 - Chr 5")

summary(out25_cim.em, threshold=3.89)
summary(out26_cim.em, threshold=3.89)
summary(out27_cim.em, threshold=3.89)
summary(out28_cim.em, threshold=3.89)
summary(out29_cim.em, threshold=3.89)
summary(out30_cim.em, threshold=3.89)


#plot CIM 15NT3
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3")
plot(out35_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3")

#plot CIM 15NT3 chrs
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,chr="2",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3 - Chr 2")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,chr="2",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3 - Chr 2")
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,chr="4",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3 - Chr 4")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,chr="4",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - 15N Time-point 3 - Chr 4")

summary(out31_cim.em, threshold=3.89)
summary(out32_cim.em, threshold=3.89)
summary(out33_cim.em, threshold=3.89)
summary(out34_cim.em, threshold=3.89)
summary(out35_cim.em, threshold=3.89)
summary(out36_cim.em, threshold=3.89)


#plot CIM aerial roots
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes")
plot(out41_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes")

#plot CIM AR chrs
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="5", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes - Chr 5")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="5", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes - Chr 5")
plot(out41_cim.em, show.marker.names=F,chr="5", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Aerial Root Nodes - Chr 5")

summary(out37_cim.em, threshold=3.89)
summary(out38_cim.em, threshold=3.89)
summary(out39_cim.em, threshold=3.89)
summary(out40_cim.em, threshold=3.89)
summary(out41_cim.em, threshold=3.89)
summary(out42_cim.em, threshold=3.89)

#plot CIM PDM
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Dry Mass")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass")
plot(out47_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass")

#plot CIM PDM chrs
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Dry Mass - Chr 2")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass - Chr 2")
plot(out47_cim.em, show.marker.names=F,chr="2", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass - Chr 2")
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="7", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass - Chr 7")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="7", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass - Chr 7")
plot(out47_cim.em, show.marker.names=F,chr="7", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Dry Mass - Chr 7")

summary(out43_cim.em, threshold=3.89)
summary(out44_cim.em, threshold=3.89)
summary(out45_cim.em, threshold=3.89)
summary(out46_cim.em, threshold=3.89)
summary(out47_cim.em, threshold=3.89)
summary(out48_cim.em, threshold=3.89)


#plot CIM PTNP
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen Percentage")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen Percentage")
plot(out53_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen Percentage")

#plot CIM PTNP chrs
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen Percentage - Chr 8")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen Percentage - Chr 8")
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen Percentage - Chr 9")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen Percentage - Chr 9")

summary(out49_cim.em, threshold=3.89)
summary(out50_cim.em, threshold=3.89)
summary(out51_cim.em, threshold=3.89)
summary(out52_cim.em, threshold=3.89)
summary(out53_cim.em, threshold=3.89)
summary(out54_cim.em, threshold=3.89)


#plot CIM PTN
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen")
plot(out59_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - PTN (T x B73)")
plot(out60_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(3.91, 0, untf = FALSE) #mcovar=1
abline(5.77, 0, untf = FALSE) #mcovar=3
#abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "QTL - PTN (T x B73)")


#plot CIM PTN chrs
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen - Chr 2")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen - Chr 2")
plot(out59_cim.em, show.marker.names=F,chr="2", col =c("purple"))
abline(3.89, 0, untf = FALSE)
abline(2.39, 0, untf = FALSE, col=c("dark gray"))
title(main = "B73 - Plant Total Nitrogen - Chr 2")
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen - Chr 8")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen - Chr 8")
plot(out59_cim.em, show.marker.names=F,chr="8", col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Plant Total Nitrogen - Chr 8")


summary(out55_cim.em, threshold=3.89)
summary(out56_cim.em, threshold=3.89)
summary(out57_cim.em, threshold=3.89)
summary(out58_cim.em, threshold=3.89)
summary(out59_cim.em, threshold=3.89)
summary(out60_cim.em, threshold=3.89)

#QTL effect plot for PTN blupNPFC
hyper <- calc.genoprob(mapthis, step = 1)
out <- scanone(hyper, method="hk")
marker_name <- find.marker(hyper, chr=2, pos=89.3)
effectplot(hyper, mname1 = marker_name, pheno.col=60, main=paste("QTL Effect at", marker_name), ylab = "Plant Total Nitrogen")

#Regression PTN

#blup no parents and family as covariate
summary(out60_cim.em, threshold=3.89)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTN_blupNPFC <- makeqtl(mapthis, chr=c(2), pos=c(89.3)) 
summary(out.qtl_PTN_blupNPFC <- fitqtl(mapthis, qtl=qtl_PTN_blupNPFC, pheno.col=60, formula=y~Q1))
plot(qtl_PTN_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_PTN_blupNPFC <- fitqtl(mapthis, qtl=qtl_PTN_blupNPFC, formula=y~Q1)
summary(out.fq_PTN_blupNPFC)


#plot CIM GDM
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass")
plot(out65_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass")

#plot CIM GDM chrs
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,chr="3",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass - Chr 3")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,chr="3",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass - Chr 3")
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,chr="5",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass - Chr 5")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,chr="5",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass - Chr 5")
plot(out65_cim.em, show.marker.names=F,chr="5",col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Dry Mass - Chr 5")

summary(out61_cim.em, threshold=3.89)
summary(out62_cim.em, threshold=3.89)
summary(out63_cim.em, threshold=3.89)
summary(out64_cim.em, threshold=3.89)
summary(out65_cim.em, threshold=3.89)
summary(out66_cim.em, threshold=3.89)

#plot CIM GTNP
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage")
plot(out71_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage")

#plot CIM GTNP chrs
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="4",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage - Chr 4")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,chr="4",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage - Chr 4")
plot(out71_cim.em, show.marker.names=F,chr="4",col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage - Chr 4")
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="7",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage - Chr 7")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,chr="7",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen Percentage - Chr 7")

summary(out67_cim.em, threshold=3.89)
summary(out68_cim.em, threshold=3.89)
summary(out69_cim.em, threshold=3.89)
summary(out70_cim.em, threshold=3.89)
summary(out71_cim.em, threshold=3.89)
summary(out72_cim.em, threshold=3.89)

#plot CIM GTN
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen")
plot(out77_cim.em, show.marker.names=F,col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen")

#plot CIM GTN chrs
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,chr="3",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen - Chr 3")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,chr="3",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen - Chr 3")
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,chr="5",col =c("Black", "blue", "red"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen - Chr 5")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,chr="5",col =c("brown", "purple", "green"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen - Chr 5")
plot(out77_cim.em, show.marker.names=F,chr="5",col =c("purple"))
abline(3.89, 0, untf = FALSE)
title(main = "B73 - Grain Total Nitrogen - Chr 5")

summary(out73_cim.em, threshold=3.89)
summary(out74_cim.em, threshold=3.89)
summary(out75_cim.em, threshold=3.89)
summary(out76_cim.em, threshold=3.89)
summary(out77_cim.em, threshold=3.89)
summary(out78_cim.em, threshold=3.89)


#lodint and bayesint for QTL SC
lodint(out1_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out1_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out4_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out4_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out5_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out5_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out6_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out6_cim.em, chr=3, expandtomarkers=TRUE)

#lodint and bayesint for QTL DT
lodint(out7_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out7_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out8_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out8_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out10_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out10_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out11_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out11_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out12_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out12_cim.em, chr=10, expandtomarkers=TRUE)

#lodint and bayesint for QTL DPS
lodint(out13_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out13_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out14_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out14_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out15_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out15_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out16_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out16_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out17_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out17_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out18_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out18_cim.em, chr=10, expandtomarkers=TRUE)

#lodint and bayesint for QTL 15NT2
lodint(out26_cim.em, chr=5, expandtomarkers=TRUE)
bayesint(out26_cim.em, chr=5, expandtomarkers=TRUE)

#lodint and bayesint for QTL 15NT3
lodint(out31_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out31_cim.em, chr=9, expandtomarkers=TRUE)

#lodint and bayesint for QTL PDM
lodint(out46_cim.em, chr=7, expandtomarkers=TRUE)
bayesint(out46_cim.em, chr=7, expandtomarkers=TRUE)

#lodint and bayesint for QTL PTNP
lodint(out50_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out50_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out51_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out51_cim.em, chr=8, expandtomarkers=TRUE)

#lodint and bayesint for QTL PTN
lodint(out55_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out55_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out57_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out57_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out58_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out58_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out59_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out59_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out60_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out60_cim.em, chr=2, expandtomarkers=TRUE)

