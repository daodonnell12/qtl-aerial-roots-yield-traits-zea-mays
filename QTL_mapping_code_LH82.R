#Input library
library(qtl)

#Input file
setwd("~/Documents/Dissertation/Dissertation_Figures/Dissertation_Figures_Ch4/Genotypes_Phenotypes/2019NewFiles")
mapthis <- read.cross("csv", "./", "genotypes_phenotypes_LH82_052320.csv", estimate.map=FALSE,genotypes=c("A","H","B"))

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
mapthis <- drop.markers(mapthis, c("S1_6850821","S1_10827206","S1_15143071","S1_15779610","S1_32577744","S1_42848278","S1_47906350","S1_47906717","S1_52161960","S1_80971388","S1_165136111","S1_179597769","S1_206305281","S1_208097466","S1_251995326","S1_262185008","S1_264959998","S1_273677663","S1_275263368","S1_276817167","S1_277349217","S1_279156239","S1_285340338","S1_285971818","S1_290856458"))
mapthis <- drop.markers(mapthis, c("S2_14805927","S2_16345553","S2_16402269","S2_26614948","S2_27645937","S2_28315702","S2_33345079","S2_39429409","S2_43709234","S2_44188941","S2_45518434","S2_46205307","S2_47111414","S2_49160585","S2_49918175","S2_54360904","S2_61245241","S2_62290816","S2_70285749","S2_83451078","S2_86076528","S2_91644011","S2_108143505","S2_111580622","S2_142632968","S2_188785420","S2_189190443","S2_198407287","S2_208377296","S2_210633818","S2_217012177","S2_224909479","S2_233598480"))
mapthis <- drop.markers(mapthis, c("S3_18077001","S3_22250501","S3_22382193","S3_23012344","S3_32250068","S3_35502433","S3_47462176","S3_48981547","S3_49148545","S3_49564943","S3_49566352","S3_50473857","S3_51901676","S3_56989982","S3_96010723","S3_102461157","S3_118968744","S3_121660277","S3_123868401","S3_124185176","S3_125548608","S3_172670386","S3_199028449","S3_199561705","S3_199562123","S3_201061921","S3_201294284","S3_219183159"))
mapthis <- drop.markers(mapthis, c("S4_120597007","S4_120887175","S4_144281175","S4_223795727","S4_233532760"))
mapthis <- drop.markers(mapthis, c("S5_3335488","S5_8000068","S5_8000418","S5_8176030","S5_18171494","S5_30218507","S5_70870679","S5_74165286","S5_76212739","S5_76216154","S5_80492787","S5_88664278","S5_111342212","S5_115585549","S5_117625205","S5_128454887","S5_111341839","S5_111341918","S5_129184341","S5_139077524","S5_144812947","S5_144815828","S5_174284687","S5_175589173","S5_193849751","S5_200292412","S5_200382626","S5_200907638","S5_210417307"))
mapthis <- drop.markers(mapthis, c("S6_31649763","S6_33638149","S6_34439965","S6_34889640","S6_59409518","S6_88983477","S6_101991870","S6_120275982","S6_120917191","S6_128735682"))
mapthis <- drop.markers(mapthis, c("S7_13408182","S7_47935582","S7_47935964","S7_50591521","S7_51800207","S7_51966842","S7_52303645","S7_80309193","S7_84701634","S7_84845749","S7_108481559","S7_108482009","S7_130238773","S7_136023171","S7_137169719","S7_149754654","S7_152064300","S7_152354373","S7_165205564","S7_167231182","S7_172158268","S7_172512783"))
mapthis <- drop.markers(mapthis, c("S8_11535401","S8_11570632","S8_15435162","S8_24523280","S8_33137479","S8_37765109","S8_52206460","S8_55996712","S8_105660390","S8_111803807","S8_121798115","S8_123530235","S8_131342237","S8_133936000","S8_133939308","S8_138939057","S8_140385582","S8_147098232","S8_148007268","S8_148008422","S8_148008651","S8_151919897","S8_156460879","S8_157097952","S8_157959908","S8_158098989","S8_170610722"))
mapthis <- drop.markers(mapthis, c("S9_20829231","S9_24597986","S9_33734463","S9_41048945","S9_53405134","S9_55347771","S9_55511805","S9_87360496","S9_87812660","S9_88649013","S9_88649541","S9_88651683","S9_91475255","S9_94038449","S9_99365183","S9_99458358","S9_100678211","S9_100684623","S9_102697880","S9_106275972","S9_124957912","S9_129518234","S9_140016203","S9_142928594","S9_143866051","S9_146458132","S9_147737403"))
mapthis <- drop.markers(mapthis, c("S10_55989558","S10_59873564","S10_60402907","S10_66014561","S10_74807759","S10_80609419","S10_83670911","S10_90244691","S10_125195158","S10_127614911","S10_133272617","S10_137830994"))

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

#Form linkage groups (checking only)
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=9)
table(lg[,2])

#Form linkage groups in actual data
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=9, reorgMarkers=TRUE)

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

mapthis <- drop.markers(mapthis, c("S6_168220601"))
mapthis <- drop.markers(mapthis, c("S4_30681","S4_2660235"))
mapthis <- drop.markers(mapthis, c("S4_240523061"))

#move markers
mn <- markernames(mapthis, chr=c(3)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(4)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 3)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 4)
mn <- markernames(mapthis, chr=c(5)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 5)
mn <- markernames(mapthis, chr=c(6)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(7)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 6)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 7)
mn <- markernames(mapthis, chr=c(8)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 8)
mn <- markernames(mapthis, chr=c(9)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 9)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
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

#Switch order of markers 
mapthis <- switch.order(mapthis, chr=1, c(2,	1,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	25,	26,	27,	28,	29,	30,	31,	32,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	44,	45,	46,	47,	48,	49,	50,	51,	52,	53,	54,	55,	56,	57,	58,	59,	60,	61,	62,	63,	64,	65,	66,	67,	68,	69,	70,	71,	72,	73,	74,	75,	76,	77,	78,	79,	80,	81,	82,	83,	84,	85,	89,	90,	88,	87,	86,	91,	92,	93,	94,	95,	96,	97,	98,	99,	100,	101,	102,	104,	103,	105,	106,	107,	108,	109,	110,	111,	112,	113,	114,	115,	116,	117,	118,	119,	120,	121,	122,	123,	124,	125,	126,	127,	128,	129,	130,	131,	132,	133,	134,	135,	136,	137,	138,	139,	140,	141,	142,	143,	144,	145,	146,	147,	148,	149,	150,	151,	152,	153,	154,	155,	156,	157,	158,	159,	160,	161,	162,	163,	164,	165,	166,	168,	167, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=1, c(168:105,103,104,102:91,87,86,88:90,85,84,80,79,81:83,2,1,3,4,78,5,77,6:10,76,75,11,12,74:71,13,70:67,14,66,15:19,65,64,20,21,63,22,62,23,61:51,24,25,50,26,49,48,27,28,47,29,30,46,31,45:43,32,33,42:38,34:37, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=1, c(167,168,166:84,80,79,81:83,78,77:67,65,66,64:3,1,2, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=2, c(2,	1,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	25,	26,	27,	28,	29,	30,	31,	32,	33,	34,	35,	36,	37,	38,	111,	110,	108,	39,	40,	105,	47,	41,	42,	48,	104,	45,	103,	101,	43,	44,	46,	49,	51,	50,	53,	112,	100,	114,	52,	54,	99,	120,	55,	98,	119,	118,	93,	56,	60,	117,	116,	65,	66,	92,	115,	113,	58,	57,	90,	86,	59,	67,	69,	61,	62,	75,	109,	107,	63,	106,	64,	68,	70,	82,	102,	80,	83,	81,	71,	72,	79,	84,	85,	97,	73,	96,	74,	95,	76,	94,	91,	89,	78,	77,	87,	88, error.prob=0.005))
#mapthis <- switch.order(mapthis, chr=2, c(1:40,42,43,46,47,45,48,41,44,49,51,50,53,52,54:56,61,60,62,57,65,66,68,69,58,59,63,70,64,71,74,75,78,79,67,80,83,84,91,72,92,94,73,76,77,98,81,82,85,99,86,102,105,87:90,108,110,112,114,93,115:117,95,96,118,97,120,119,113,100,111,101,103,104,106,107,109, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=3, c(2,1,101,3:19,100,20,99,21,22,98,97,23,96,24,25,26,95:90,27,28,89,29,30,88,31:34,87,35,36,86:84,37:42,83,43,82,44,81,80,45,79,78,46,47,77,48,76,49:51,75,74,52,73,72:70,53,54,69,55,68,56,67,57,66,65,58,59,60:64, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=4, c(46,1,45,2,3,44,4:13,43:40,14,15,39:37,16,36:34,17,18,33:31,19:30, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=5, c(100,1,101,102,99:96,2:6,95,7,94,8,9,93,10,92,11,91,12,14,13,90,15:18,89,19,88:86,20:23,85,84,24:26,83:80,27:29,79,30,78,31,77,32,76,33:35,37,36,75,74,38,73,39,72:70,40:45,69,68,46,67,47,48,66,49,65,64,50:52,63:61,53,60,59,58,55,54,56,57, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=6, c(79:65,63,64,62:1, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=7, c(1:3,81:78,4,77:70,5:9,69,10,68,67,11:13,66,14,65:62,15,61,16,60,17,59:57,18,56,19,20,55,21,54:51,22,50,23,49,48,24,47:43,25,42,26,41:38,27,37:35,28:30,34,31,33,32, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=8, c(2,1,3:38,85,84,39,40,83,42,41,43,82,81,44:46,80,78,79,47,77,48:50,76:73,51,72,52,71,70,53,54,69,55,68,67,56:62,66,65,63,64, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=9, c(75:6,4,5,3:1, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=10, c(1:21,23,22,24:54,58:55, error.prob=0.005))

#plot linkage map
plotMap(mapthis)
plotMap(mapthis, show.marker.names=F)
plotMap(mapthis, show.marker.names=F, chr=1)

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
out.mr[ out.mr$chr == 1, ]

#plot single analysis
plot(out.mr, show.marker.names=F)
plot(out.mr, show.marker.names=F, chr=c("3"))

#SIM
out.em <- scanone (mapthis, pheno.col=c(41),method="em")
out.em1 <- scantwo (mapthis, pheno.col=c(41),method="em")
summary(out.em, threshold=3.91)
max(out.em)
out.em[ out.em$chr == 1, ]
plot(out.em, lodcolumn=c(1), lty=1, col=c("blue"), show.marker.names=F)
abline(3.91, 0, untf = FALSE)
plot(out.em, lodcolumn=c(1),chr=c("1"), lty=1, col=c("Blue"), show.marker.names=F)

#CIM
out_cim.em <- cim(mapthis, pheno.col = (39), n.marcovar=1, method=c("em"), map.function=c("kosambi") )
summary(out_cim.em, threshold=3)
max(out_cim.em)
out_cim.em[ out_cim.em$chr == 8, ]
plot(out_cim.em, show.marker.names=F, col =c("red"))
plot(out_cim.em, chr=c("8"), show.marker.names=F, col =c("red"))
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
plot(out.mr, out.em, out_cim.em, chr="8",show.marker.names=F,col =c("Black", "blue", "red"))

out.mr <- scanone(mapthis, method="mr",pheno.col=40)
summary(out.mr, threshold=2)
out.em <- scanone (mapthis, pheno.col= (40),method="em")
summary(out.em, threshold=2)
out_cim.em <- cim(mapthis, pheno.col = (40), n.marcovar=3, method=c("em"), map.function=c("kosambi") )
summary(out_cim.em, threshold=2)

# combine CIM chart

#CIM all phenotypes
out1_cim.em <- cim(mapthis, pheno.col = (1), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out2_cim.em <- cim(mapthis, pheno.col = (2), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out3_cim.em <- cim(mapthis, pheno.col = (3), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out4_cim.em <- cim(mapthis, pheno.col = (4), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out5_cim.em <- cim(mapthis, pheno.col = (5), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out6_cim.em <- cim(mapthis, pheno.col = (6), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out7_cim.em <- cim(mapthis, pheno.col = (7), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out8_cim.em <- cim(mapthis, pheno.col = (8), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out9_cim.em <- cim(mapthis, pheno.col = (9), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out10_cim.em <- cim(mapthis, pheno.col = (10), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out11_cim.em <- cim(mapthis, pheno.col = (11), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out12_cim.em <- cim(mapthis, pheno.col = (12), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out13_cim.em <- cim(mapthis, pheno.col = (13), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out14_cim.em <- cim(mapthis, pheno.col = (14), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out15_cim.em <- cim(mapthis, pheno.col = (15), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out16_cim.em <- cim(mapthis, pheno.col = (16), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out17_cim.em <- cim(mapthis, pheno.col = (17), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out18_cim.em <- cim(mapthis, pheno.col = (18), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out19_cim.em <- cim(mapthis, pheno.col = (19), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out20_cim.em <- cim(mapthis, pheno.col = (20), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out21_cim.em <- cim(mapthis, pheno.col = (21), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out22_cim.em <- cim(mapthis, pheno.col = (22), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out23_cim.em <- cim(mapthis, pheno.col = (23), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out24_cim.em <- cim(mapthis, pheno.col = (24), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out25_cim.em <- cim(mapthis, pheno.col = (25), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out26_cim.em <- cim(mapthis, pheno.col = (26), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out27_cim.em <- cim(mapthis, pheno.col = (27), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out28_cim.em <- cim(mapthis, pheno.col = (28), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out29_cim.em <- cim(mapthis, pheno.col = (29), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out30_cim.em <- cim(mapthis, pheno.col = (30), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out31_cim.em <- cim(mapthis, pheno.col = (31), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out32_cim.em <- cim(mapthis, pheno.col = (32), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out33_cim.em <- cim(mapthis, pheno.col = (33), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out34_cim.em <- cim(mapthis, pheno.col = (34), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out35_cim.em <- cim(mapthis, pheno.col = (35), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out36_cim.em <- cim(mapthis, pheno.col = (36), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out37_cim.em <- cim(mapthis, pheno.col = (37), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out38_cim.em <- cim(mapthis, pheno.col = (38), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out39_cim.em <- cim(mapthis, pheno.col = (39), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out40_cim.em <- cim(mapthis, pheno.col = (40), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out41_cim.em <- cim(mapthis, pheno.col = (41), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out42_cim.em <- cim(mapthis, pheno.col = (42), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out37_sim.em <- scanone(mapthis, pheno.col = (37), method="em")
out38_sim.em <- scanone(mapthis, pheno.col = (38), method="em")
out39_sim.em <- scanone(mapthis, pheno.col = (39), method="em")
out40_sim.em <- scanone(mapthis, pheno.col = (40), method="em")
out41_sim.em <- scanone(mapthis, pheno.col = (41), method="em")
out42_sim.em <- scanone(mapthis, pheno.col = (42), method="em")
out43_cim.em <- cim(mapthis, pheno.col = (43), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out44_cim.em <- cim(mapthis, pheno.col = (44), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out45_cim.em <- cim(mapthis, pheno.col = (45), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out46_cim.em <- cim(mapthis, pheno.col = (46), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out47_cim.em <- cim(mapthis, pheno.col = (47), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out48_cim.em <- cim(mapthis, pheno.col = (48), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out49_cim.em <- cim(mapthis, pheno.col = (49), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out50_cim.em <- cim(mapthis, pheno.col = (50), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out51_cim.em <- cim(mapthis, pheno.col = (51), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out52_cim.em <- cim(mapthis, pheno.col = (52), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out53_cim.em <- cim(mapthis, pheno.col = (53), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out54_cim.em <- cim(mapthis, pheno.col = (54), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out55_cim.em <- cim(mapthis, pheno.col = (55), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out56_cim.em <- cim(mapthis, pheno.col = (56), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out57_cim.em <- cim(mapthis, pheno.col = (57), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out58_cim.em <- cim(mapthis, pheno.col = (58), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out59_cim.em <- cim(mapthis, pheno.col = (59), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out60_cim.em <- cim(mapthis, pheno.col = (60), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out61_cim.em <- cim(mapthis, pheno.col = (61), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out62_cim.em <- cim(mapthis, pheno.col = (62), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out63_cim.em <- cim(mapthis, pheno.col = (63), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out64_cim.em <- cim(mapthis, pheno.col = (64), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out65_cim.em <- cim(mapthis, pheno.col = (65), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out66_cim.em <- cim(mapthis, pheno.col = (66), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out67_cim.em <- cim(mapthis, pheno.col = (67), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out68_cim.em <- cim(mapthis, pheno.col = (68), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out69_cim.em <- cim(mapthis, pheno.col = (69), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out70_cim.em <- cim(mapthis, pheno.col = (70), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out71_cim.em <- cim(mapthis, pheno.col = (71), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out72_cim.em <- cim(mapthis, pheno.col = (72), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out73_cim.em <- cim(mapthis, pheno.col = (73), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out74_cim.em <- cim(mapthis, pheno.col = (74), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out75_cim.em <- cim(mapthis, pheno.col = (75), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out76_cim.em <- cim(mapthis, pheno.col = (76), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out77_cim.em <- cim(mapthis, pheno.col = (77), n.marcovar=1, method=c("em"), map.function=c("kosambi"))
out78_cim.em <- cim(mapthis, pheno.col = (78), n.marcovar=1, method=c("em"), map.function=c("kosambi"))

#plot CIM stand count
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#lod 2.79 based on lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
#lod 2.41 = drop 1.5 from lod significance threshold 3.89
title(main = "LH82 - Stand Count")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count")
plot(out5_cim.em, show.marker.names=F,col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count")

#plot CIM SC chrs
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count - Chr 1")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count - Chr 1")
plot(out5_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count - Chr 1")

summary(out1_cim.em, threshold=3.91)
summary(out2_cim.em, threshold=3.91)
summary(out3_cim.em, threshold=3.91)
summary(out4_cim.em, threshold=3.91)
summary(out5_cim.em, threshold=3.91)
summary(out6_cim.em, threshold=3.91)

lodint(out.qtl_SC_avg_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_avg_chr1, chr=1, expandtomarkers=TRUE)
lodint(out.qtl_SC_rep1_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_rep1_chr1, chr=1, expandtomarkers=TRUE)
lodint(out.qtl_SC_rep2_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_rep2_chr1, chr=1, expandtomarkers=TRUE)
lodint(out.qtl_SC_blup_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_blup_chr1, chr=1, expandtomarkers=TRUE)
lodint(out.qtl_SC_blupNP_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_blupNP_chr1, chr=1, expandtomarkers=TRUE)
lodint(out.qtl_SC_blupNPFC_chr1, chr=1, expandtomarkers=TRUE)
bayesint(out.qtl_SC_blupNPFC_chr1, chr=1, expandtomarkers=TRUE)

#Regression SC
#Average
summary(out1_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_avg_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(213)) 
summary(out.qtl_SC_avg_chr1 <- fitqtl(mapthis, qtl=qtl_SC_avg_chr1, pheno.col=1, formula=y~Q1))
plot(qtl_SC_avg_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_avg_chr1 <- fitqtl(mapthis, qtl=qtl_SC_avg_chr1, formula=y~Q1)
summary(out.fq_SC_avg_chr1)

#refine QTL and plot LOD profile
rqtl_SC_avg_chr1 <- refineqtl(mapthis, qtl=qtl_SC_avg_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_avg_chr1
out.fq2_SC_avg_chr1 <- fitqtl(mapthis, qtl=qtl_SC_avg_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_avg_chr1)
plotLodProfile(rqtl_SC_avg_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Average Rep1 + Rep2)")

#rep1
summary(out2_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_rep1_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(136)) 
summary(out.qtl_SC_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr1, pheno.col=2, formula=y~Q1))
plot(qtl_SC_rep1_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr1, formula=y~Q1)
summary(out.fq_SC_rep1_chr1)

#refine QTL and plot LOD profile
rqtl_SC_rep1_chr1 <- refineqtl(mapthis, qtl=qtl_SC_rep1_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_rep1_chr1
out.fq2_SC_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_rep1_chr1)
plotLodProfile(rqtl_SC_rep1_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Rep1)")

#rep2
summary(out3_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_rep2_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(216)) 
summary(out.qtl_SC_rep2_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep2_chr1, pheno.col=3, formula=y~Q1))
plot(qtl_SC_rep2_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_rep2_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep2_chr1, formula=y~Q1)
summary(out.fq_SC_rep2_chr1)

#refine QTL and plot LOD profile
rqtl_SC_rep2_chr1 <- refineqtl(mapthis, qtl=qtl_SC_rep2_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_rep2_chr1
out.fq2_SC_rep2_chr1 <- fitqtl(mapthis, qtl=qtl_SC_rep2_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_rep2_chr1)
plotLodProfile(rqtl_SC_rep2_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Rep2)")

#blup
summary(out4_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_blup_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(213)) 
summary(out.qtl_SC_blup_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blup_chr1, pheno.col=4, formula=y~Q1))
plot(qtl_SC_blup_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_blup_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blup_chr1, formula=y~Q1)
summary(out.fq_SC_blup_chr1)

#refine QTL and plot LOD profile
rqtl_SC_blup_chr1 <- refineqtl(mapthis, qtl=qtl_SC_blup_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_blup_chr1
out.fq2_SC_blup_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blup_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_blup_chr1)
plotLodProfile(rqtl_SC_blup_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out5_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_blupNP_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(213)) 
summary(out.qtl_SC_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNP_chr1, pheno.col=5, formula=y~Q1))
plot(qtl_SC_blupNP_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNP_chr1, formula=y~Q1)
summary(out.fq_SC_blupNP_chr1)

#refine QTL and plot LOD profile
rqtl_SC_blupNP_chr1 <- refineqtl(mapthis, qtl=qtl_SC_blupNP_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_blupNP_chr1
out.fq2_SC_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNP_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_blupNP_chr1)
plotLodProfile(rqtl_SC_blupNP_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariates
summary(out6_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_blupNPFC_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(213)) 
summary(out.qtl_SC_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr1, pheno.col=6, formula=y~Q1))
plot(qtl_SC_blupNPFC_chr1)

#fit the model (with QTL in fixed positions)
out.fq_SC_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr1, formula=y~Q1)
summary(out.fq_SC_blupNPFC_chr1)

#refine QTL and plot LOD profile
rqtl_SC_blupNPFC_chr1 <- refineqtl(mapthis, qtl=qtl_SC_blupNPFC_chr1, formula=y~Q1, verbose=FALSE)
rqtl_SC_blupNPFC_chr1
out.fq2_SC_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_blupNPFC_chr1)
plotLodProfile(rqtl_SC_blupNPFC_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray"))
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Stand Count (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod support interval coordinates
lodint(rqtl_SC_avg_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_rep1_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_rep2_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_blup_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_blupNP_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_blupNPFC_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_SC_avg_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_rep1_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_rep2_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_blup_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_blupNP_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_blupNPFC_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM DT
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Tassel Emergence")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Tassel Emergence")

#plot CIM DT chrs
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)

#plot CIM DPS
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.48, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from peak
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Days to Pollen Shed")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Pollen Shed")

#plot CIM DPS chrs
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.48, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from peak
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Days to Pollen Shed - Chr 10")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
summary(out13_cim.em, threshold=3.91)
summary(out14_cim.em, threshold=3.91)
summary(out15_cim.em, threshold=3.91)
summary(out16_cim.em, threshold=3.91)
summary(out17_cim.em, threshold=3.91)
summary(out18_cim.em, threshold=3.91)

#Regression DPS
#rep2
summary(out15_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_rep2 <- makeqtl(mapthis, chr=c(10), pos=c(99.8)) 
summary(out.qtl_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, pheno.col=15, formula=y~Q1))
plot(qtl_DPS_rep2)

#fit the model (with QTL in fixed positions)
out.fq_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1)
summary(out.fq_DPS_rep2)

#refine QTL and plot LOD profile
rqtl_DPS_rep2 <- refineqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1, verbose=FALSE)
rqtl_DPS_rep2
out.fq2_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_rep2)
plotLodProfile(rqtl_DPS_rep2)
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Days to Pollen Shed (Rep2)")

#obtain lod interval coordinates
lodint(rqtl_DPS_rep2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_DPS_rep2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM 15NT1
plot(out19_cim.em, out20_cim.em, out21_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 1")
plot(out22_cim.em, out23_cim.em, out24_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 1")

#plot CIM 15NT1 chrs
plot(out19_cim.em, out20_cim.em, out21_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
plot(out22_cim.em, out23_cim.em, out24_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)

#plot CIM 15NT2
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 2")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 2")

#plot CIM 15NT2 chrs

#plot CIM 15NT3
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 3")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 3")

#plot CIM 15NT3 chrs
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,chr="5", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,chr="5", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)

#plot CIM aerial roots
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes")
plot(out41_cim.em, show.marker.names=F,col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes")

#plot CIM AR chrs
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 1")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 1")
plot(out41_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 1")
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 2")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Aerial Root Nodes - Chr 2")
plot(out41_cim.em, show.marker.names=F,chr="2", col =c("purple"))
abline(3.91, 0, untf = FALSE)
#abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Aerial Root Nodes - Chr 2")
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 9")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 9")
plot(out41_cim.em, show.marker.names=F,chr="9", col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop from lod threshold
title(main = "LH82 - Aerial Root Nodes - Chr 9")

summary(out37_cim.em, threshold=3.91)
summary(out38_cim.em, threshold=3.91)
summary(out39_cim.em, threshold=3.91)
summary(out40_cim.em, threshold=3.91)
summary(out41_cim.em, threshold=3.91)
summary(out42_cim.em, threshold=3.91)

#Regression AR
#Average
summary(out37_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_avg_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(251.5))
qtl_AR_avg_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(97.1)) 
qtl_AR_avg <- makeqtl(mapthis, chr=c(1, 9), pos=c(251.1, 97.1))
summary(out.qtl_AR_avg_chr1 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr1, pheno.col=37, formula=y~Q1))
summary(out.qtl_AR_avg_chr9 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr9, pheno.col=37, formula=y~Q1))
summary(out.qtl_AR_avg <- fitqtl(mapthis, qtl=qtl_AR_avg, pheno.col=37, formula=y~Q1+Q2))
plot(qtl_AR_avg_chr1)
plot(qtl_AR_avg_chr9)
plot(qtl_AR_avg)

#fit the model (with QTL in fixed positions)
out.fq_AR_avg_chr1 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr1, formula=y~Q1)
summary(out.fq_AR_avg_chr1)
out.fq_AR_avg_chr9 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr9, formula=y~Q1)
summary(out.fq_AR_avg_chr9)
out.fq_AR_avg <- fitqtl(mapthis, qtl=qtl_AR_avg, formula=y~Q1+Q2)
summary(out.fq_AR_avg)

#refine QTL and plot LOD profile
#uses an iterative algorithm to refine the locations of QTL in a multiple-QTL model
#with the aim of obtaining the maximum likelihood estimates of the QTL positions

rqtl_AR_avg_chr1 <- refineqtl(mapthis, qtl=qtl_AR_avg_chr1, formula=y~Q1, verbose=FALSE)
rqtl_AR_avg_chr1
out.fq2_AR_avg_chr1 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_avg_chr1)
plotLodProfile(rqtl_AR_avg_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2)")

rqtl_AR_avg_chr9 <- refineqtl(mapthis, qtl=qtl_AR_avg_chr9, formula=y~Q1, verbose=FALSE)
rqtl_AR_avg_chr9
out.fq2_AR_avg_chr9 <- fitqtl(mapthis, qtl=qtl_AR_avg_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_avg_chr9)
plotLodProfile(rqtl_AR_avg_chr9)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2)")

rqtl_AR_avg <- refineqtl(mapthis, qtl=qtl_AR_avg, formula=y~Q1+Q2, verbose=FALSE)
rqtl_AR_avg
out.fq2_AR_avg <- fitqtl(mapthis, qtl=qtl_AR_avg, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_AR_avg)
plotLodProfile(rqtl_AR_avg)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2)")

#test pairwise interactions between QTL
addint(mapthis, qtl=rqtl_AR_avg, formula=y~Q1+Q2)

#rep1
summary(out38_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_rep1_chr2 <- makeqtl(mapthis, chr=c(2), pos=c(147)) 
summary(out.qtl_AR_rep1_chr2 <- fitqtl(mapthis, qtl=qtl_AR_rep1_chr2, pheno.col=38, formula=y~Q1))
plot(qtl_AR_rep1_chr2)

#fit the model (with QTL in fixed positions)
out.fq_AR_rep1_chr2 <- fitqtl(mapthis, qtl=qtl_AR_rep1_chr2, formula=y~Q1)
summary(out.fq_AR_rep1_chr2)

#refine QTL and plot LOD profile
rqtl_AR_rep1_chr2 <- refineqtl(mapthis, qtl=qtl_AR_rep1_chr2, formula=y~Q1, verbose=FALSE)
rqtl_AR_rep1_chr2
out.fq2_AR_rep1_chr2 <- fitqtl(mapthis, qtl=qtl_AR_rep1_chr2, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_rep1_chr2)
plotLodProfile(rqtl_AR_rep1_chr2)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Rep1)")

#rep2
summary(out39_cim.em, threshold=3.91)

#blup
summary(out40_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128) #QTL object containing imputed genotypes, with 128 imputations
qtl_AR_blup_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(251.5))
qtl_AR_blup_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(97.1)) 
qtl_AR_blup <- makeqtl(mapthis, chr=c(1, 9), pos=c(251.1, 97.1))
summary(out.qtl_AR_blup_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr1, pheno.col=40, formula=y~Q1))
summary(out.qtl_AR_blup_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr9, pheno.col=40, formula=y~Q1))
summary(out.qtl_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, pheno.col=40, formula=y~Q1+Q2))
plot(qtl_AR_blup_chr1)
plot(qtl_AR_blup_chr9)
plot(qtl_AR_blup)

#fit the model (with QTL in fixed positions)
out.fq_AR_blup_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr1, formula=y~Q1)
summary(out.fq_AR_blup_chr1)
out.fq_AR_blup_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr9, formula=y~Q1)
summary(out.fq_AR_blup_chr9)
out.fq_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1+Q2)
summary(out.fq_AR_blup)

#refine QTL and plot LOD profile
rqtl_AR_blup_chr1 <- refineqtl(mapthis, qtl=qtl_AR_blup_chr1, formula=y~Q1, verbose=FALSE)
rqtl_AR_blup_chr1
out.fq2_AR_blup_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blup_chr1)
plotLodProfile(rqtl_AR_blup_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP)")

rqtl_AR_blup_chr9 <- refineqtl(mapthis, qtl=qtl_AR_blup_chr9, formula=y~Q1, verbose=FALSE)
rqtl_AR_blup_chr9
out.fq2_AR_blup_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blup_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blup_chr9)
plotLodProfile(rqtl_AR_blup_chr9)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP)")

rqtl_AR_blup <- refineqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1+Q2, verbose=FALSE)
rqtl_AR_blup
out.fq2_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_AR_blup)
plotLodProfile(rqtl_AR_blup)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP)")

#test pairwise interactions between QTL
addint(mapthis, qtl=rqtl_AR_blup, formula=y~Q1+Q2)

#blup no parents
summary(out41_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_blupNP_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(251.5))
qtl_AR_blupNP_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(97.1))
qtl_AR_blupNP <- makeqtl(mapthis, chr=c(1, 9), pos=c(251.1, 97.1))
summary(out.qtl_AR_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr1, pheno.col=41, formula=y~Q1))
summary(out.qtl_AR_blupNP_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr9, pheno.col=41, formula=y~Q1))
summary(out.qtl_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, pheno.col=41, formula=y~Q1+Q2))
plot(qtl_AR_blupNP_chr1)
plot(qtl_AR_blupNP_chr9)
plot(qtl_AR_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_AR_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr1, formula=y~Q1)
summary(out.fq_AR_blupNP_chr1)
out.fq_AR_blupNP_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr9, formula=y~Q1)
summary(out.fq_AR_blupNP_chr9)
out.fq_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1+Q2)
summary(out.fq_AR_blupNP)

#refine QTL and plot LOD profile
rqtl_AR_blupNP_chr1 <- refineqtl(mapthis, qtl=qtl_AR_blupNP_chr1, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNP_chr1
out.fq2_AR_blupNP_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNP_chr1)
plotLodProfile(rqtl_AR_blupNP_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NP)")

rqtl_AR_blupNP_chr9 <- refineqtl(mapthis, qtl=qtl_AR_blupNP_chr9, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNP_chr9
out.fq2_AR_blupNP_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNP_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNP_chr9)
plotLodProfile(rqtl_AR_blupNP_chr9)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NP)")

rqtl_AR_blupNP <- refineqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1+Q2, verbose=FALSE)
rqtl_AR_blupNP
out.fq2_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_AR_blupNP)
plotLodProfile(rqtl_AR_blupNP)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NP)")

#test pairwise interactions between QTL
addint(mapthis, qtl=rqtl_AR_blupNP, formula=y~Q1+Q2)

#blup no parents and family as covariates
summary(out42_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_blupNPFC_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(251.5))
qtl_AR_blupNPFC_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(97.1))
qtl_AR_blupNPFC <- makeqtl(mapthis, chr=c(1, 9), pos=c(251.1, 97.1))
summary(out.qtl_AR_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr1, pheno.col=42, formula=y~Q1))
summary(out.qtl_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, pheno.col=42, formula=y~Q1))
summary(out.qtl_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, pheno.col=42, formula=y~Q1+Q2))
plot(qtl_AR_blupNPFC_chr1)
plot(qtl_AR_blupNPFC_chr9)
plot(qtl_AR_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_AR_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr1, formula=y~Q1)
summary(out.fq_AR_blupNPFC_chr1)
out.fq_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1)
summary(out.fq_AR_blupNPFC_chr9)
out.fq_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2)
summary(out.fq_AR_blupNPFC)

#refine QTL and plot LOD profile
rqtl_AR_blupNPFC_chr1 <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC_chr1, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNPFC_chr1
out.fq2_AR_blupNPFC_chr1 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNPFC_chr1)
plotLodProfile(rqtl_AR_blupNPFC_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NPFC)")

rqtl_AR_blupNPFC_chr9 <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNPFC_chr9
out.fq2_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNPFC_chr9)
plotLodProfile(rqtl_AR_blupNPFC_chr9)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NPFC)")

rqtl_AR_blupNPFC <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2, verbose=FALSE)
rqtl_AR_blupNPFC
out.fq2_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_AR_blupNPFC)
plotLodProfile(rqtl_AR_blupNPFC)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NPFC)")

#test pairwise interactions between QTL
addint(mapthis, qtl=rqtl_AR_blupNPFC, formula=y~Q1+Q2)

#obtain lod interval coordinates
lodint(rqtl_AR_avg_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blup_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNP_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNPFC_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_rep1_chr2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_avg_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blup_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNP_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNPFC_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_AR_avg_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blup_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNP_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNPFC_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_rep1_chr2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_avg_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blup_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNP_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNPFC_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM PDM
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Dry Mass")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.58, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Dry Mass")
plot(out47_cim.em, show.marker.names=F,col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.58, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Dry Mass")

#plot CIM PDM chrs
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Dry Mass - Chr 10")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
abline(2.58, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Dry Mass - Chr 10")
plot(out47_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.91, 0, untf = FALSE)
abline(2.58, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Dry Mass - Chr 10")
#plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
#abline(3.91, 0, untf = FALSE)
#plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
#abline(3.91, 0, untf = FALSE)
#plot(out47_cim.em, show.marker.names=F,chr="8", col =c("purple"))
#abline(3.91, 0, untf = FALSE)
#plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
#abline(3.91, 0, untf = FALSE)
#plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
#abline(3.91, 0, untf = FALSE)
#plot(out47_cim.em, show.marker.names=F,chr="9", col =c("purple"))
#abline(3.91, 0, untf = FALSE)
summary(out43_cim.em, threshold=3.91)
summary(out44_cim.em, threshold=3.91)
summary(out45_cim.em, threshold=3.91)
summary(out46_cim.em, threshold=3.91)
summary(out47_cim.em, threshold=3.91)
summary(out48_cim.em, threshold=3.91)

#Regression PDM
#blup
summary(out46_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_blup_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(185))
summary(out.qtl_PDM_blup_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blup_chr10, pheno.col=46, formula=y~Q1))
plot(qtl_PDM_blup_chr10)

#fit the model (with QTL in fixed positions)
out.fq_PDM_blup_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blup_chr10, formula=y~Q1)
summary(out.fq_PDM_blup_chr10)

#refine QTL and plot LOD profile
rqtl_PDM_blup_chr10 <- refineqtl(mapthis, qtl=qtl_PDM_blup_chr10, formula=y~Q1, verbose=FALSE)
rqtl_PDM_blup_chr10
out.fq2_PDM_blup_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blup_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_blup_chr10)
plotLodProfile(rqtl_PDM_blup_chr10)
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Dry Mass (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out47_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_blupNP_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(185))
summary(out.qtl_PDM_blupNP_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNP_chr10, pheno.col=47, formula=y~Q1))
plot(qtl_PDM_blupNP_chr10)

#fit the model (with QTL in fixed positions)
out.fq_PDM_blupNP_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNP_chr10, formula=y~Q1)
summary(out.fq_PDM_blupNP_chr10)

#refine QTL and plot LOD profile
rqtl_PDM_blupNP_chr10 <- refineqtl(mapthis, qtl=qtl_PDM_blupNP_chr10, formula=y~Q1, verbose=FALSE)
rqtl_PDM_blupNP_chr10
out.fq2_PDM_blupNP_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNP_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_blupNP_chr10)
plotLodProfile(rqtl_PDM_blupNP_chr10)
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Dry Mass (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariates
summary(out48_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_blupNPFC_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(185))
summary(out.qtl_PDM_blupNPFC_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC_chr10, pheno.col=48, formula=y~Q1))
plot(qtl_PDM_blupNPFC_chr10)

#fit the model (with QTL in fixed positions)
out.fq_PDM_blupNPFC_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC_chr10, formula=y~Q1)
summary(out.fq_PDM_blupNPFC_chr10)

#refine QTL and plot LOD profile
rqtl_PDM_blupNPFC_chr10 <- refineqtl(mapthis, qtl=qtl_PDM_blupNPFC_chr10, formula=y~Q1, verbose=FALSE)
rqtl_PDM_blupNPFC_chr10
out.fq2_PDM_blupNPFC_chr10 <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_blupNPFC_chr10)
plotLodProfile(rqtl_PDM_blupNPFC_chr10)
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Dry Mass (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_PDM_blup_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_blupNP_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_blupNPFC_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PDM_blup_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_blupNP_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_blupNPFC_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM PTNP
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Total Nitrogen Percentage")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage")

#plot CIM PTNP chrs
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 1")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 1")
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 9")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 9")
plot(out53_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 9")
summary(out49_cim.em, threshold=3.91)
summary(out50_cim.em, threshold=3.91)
summary(out51_cim.em, threshold=3.91)
summary(out52_cim.em, threshold=3.91)
summary(out53_cim.em, threshold=3.91)
summary(out54_cim.em, threshold=3.91)

#Regression PTNP
#rep1
summary(out50_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTNP_rep1_chr1 <- makeqtl(mapthis, chr=c(1), pos=c(212))
summary(out.qtl_PTNP_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_PTNP_rep1_chr1, pheno.col=50, formula=y~Q1))
plot(qtl_PTNP_rep1_chr1)

#fit the model (with QTL in fixed positions)
out.fq_PTNP_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_PTNP_rep1_chr1, formula=y~Q1)
summary(out.fq_PTNP_rep1_chr1)

#refine QTL and plot LOD profile
rqtl_PTNP_rep1_chr1 <- refineqtl(mapthis, qtl=qtl_PTNP_rep1_chr1, formula=y~Q1, verbose=FALSE)
rqtl_PTNP_rep1_chr1
out.fq2_PTNP_rep1_chr1 <- fitqtl(mapthis, qtl=qtl_PTNP_rep1_chr1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTNP_rep1_chr1)
plotLodProfile(rqtl_PTNP_rep1_chr1)
abline(3.91, 0, untf = FALSE)
abline(2.79, 0, untf = FALSE, col=c("dark gray")) #lod interval
#abline(2.41, 0, untf = FALSE, col=c("dark gray")) #1.5 drop lod threshold
title(main = "LH82 - Plant Total Nitrogen Percentage (Rep1)")

#obtain lod interval coordinates
lodint(rqtl_PTNP_rep1_chr1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PTNP_rep1_chr1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM PTN
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Total Nitrogen")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen")

#plot CIM PTN chrs
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,chr="5", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen - Chr 5")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,chr="5", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen - Chr 5")
plot(out59_cim.em, show.marker.names=F,chr="5", col =c("purple"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen - Chr 5")
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Total Nitrogen - Chr 9")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen - Chr 9")
plot(out59_cim.em, show.marker.names=F,chr="9", col =c("purple"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen - Chr 9")
summary(out55_cim.em, threshold=3.91)
summary(out56_cim.em, threshold=3.91)
summary(out57_cim.em, threshold=3.91)
summary(out58_cim.em, threshold=3.91)
summary(out59_cim.em, threshold=3.91)
summary(out60_cim.em, threshold=3.91)

#Regression PTN
#rep2
summary(out57_cim.em, threshold=3.91)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTN_rep2_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(93.3))
summary(out.qtl_PTN_rep2_chr9 <- fitqtl(mapthis, qtl=qtl_PTN_rep2_chr9, pheno.col=57, formula=y~Q1))
plot(qtl_PTN_rep2_chr9)

#fit the model (with QTL in fixed positions)
out.fq_PTN_rep2_chr9 <- fitqtl(mapthis, qtl=qtl_PTN_rep2_chr9, formula=y~Q1)
summary(out.fq_PTN_rep2_chr9)

#refine QTL and plot LOD profile
rqtl_PTN_rep2_chr9 <- refineqtl(mapthis, qtl=qtl_PTN_rep2_chr9, formula=y~Q1, verbose=FALSE)
rqtl_PTN_rep2_chr9
out.fq2_PTN_rep2_chr9 <- fitqtl(mapthis, qtl=qtl_PTN_rep2_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTN_rep2_chr9)
plotLodProfile(rqtl_PTN_rep2_chr9)
abline(3.91, 0, untf = FALSE)
abline(2.41, 0, untf = FALSE, col=c("dark gray"))
title(main = "LH82 - Plant Total Nitrogen (Rep2)")

#obtain lod interval coordinates
lodint(rqtl_PTN_rep2_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PTN_rep2_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM GDM
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass")

#plot CIM GDM chrs
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,chr="7", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass - Chr 7")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,chr="7", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass - Chr 7")

#plot CIM GTNP
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage")

#plot CIM GTNP chrs
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage - Chr 10")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage - Chr 10")

#plot CIM GTN
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen")

#plot CIM GTN chrs
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,chr="7", col =c("Black", "blue", "red"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen - Chr 7")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,chr="7", col =c("brown", "purple", "green"))
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen - Chr 7")

#plot CIM blup_noparents SC
plot(out5_cim.em, show.marker.names=F,col =c("purple")) #SC
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Stand Count")
summary(out5_cim.em, threshold=3.91)
plot(out5_cim.em, show.marker.names=F,chr="1", col =c("purple")) #SC
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Stand Count - Chr 1")

#Regression SC
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC <- makeqtl(mapthis, chr=c(1), pos=c(213)) 
summary(out.qtl_SC <- fitqtl(mapthis, qtl=qtl_SC, pheno.col=5, formula=y~Q1)) 

#plot CIM blup_noparents DT
plot(out11_cim.em, show.marker.names=F,col =c("purple")) #DT
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Tassel Emergence")
summary(out11_cim.em, threshold=3.91)
plot(out11_cim.em, show.marker.names=F,chr="8", col =c("purple")) #DT
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Tassel Emergence - Chr 8")

#Regression DT
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_DT <- makeqtl(mapthis, chr=c(5), pos=c(132))
#summary(out.qtl_DT <- fitqtl(mapthis, qtl=qtl_DT, pheno.col=11, formula=y~Q1))

#plot CIM blup_noparents DPS
plot(out17_cim.em, show.marker.names=F,col =c("purple")) #DPS
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Days to Pollen Shed")
summary(out17_cim.em, threshold=3.91)

#Regression DPS
#There were no LOD peaks above the threshold.

#plot CIM blup_noparents T1
plot(out23_cim.em, show.marker.names=F,col =c("purple")) #15NT1
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 1")
summary(out23_cim.em, threshold=3.91)
plot(out23_cim.em, show.marker.names=F,chr="2", col =c("purple")) #15NT1
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 1 - Chr 2")

#Regression T1
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_T1 <- makeqtl(mapthis, chr=c(2), pos=c(82.6)) 
#summary(out.qtl_T1 <- fitqtl(mapthis, qtl=qtl_T1, pheno.col=23, formula=y~Q1))

#plot CIM blup_noparents T2
plot(out29_cim.em, show.marker.names=F,col =c("purple")) #15NT2
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 2")
summary(out29_cim.em, threshold=3.91)

#Regression T2
#There were no LOD peaks above the threshold.

#plot CIM blup_noparents T3
plot(out35_cim.em, show.marker.names=F,col =c("purple")) #15NT3
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - 15N Time-point 3")
summary(out35_cim.em, threshold=3.91)

#Regression T3
#There were no LOD peaks above the threshold.

#plot CIM blup_noparents AR
plot(out41_cim.em, show.marker.names=F,col =c("purple")) #AR
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Aerial Root Nodes")
summary(out41_cim.em, threshold=3.91)
plot(out41_cim.em, show.marker.names=F,chr="1", col =c("purple")) #AR
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Aerial Root Nodes - Chr 1")
plot(out41_cim.em, show.marker.names=F,chr="2", col =c("purple")) #AR
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Aerial Root Nodes - Chr 2")
plot(out41_cim.em, show.marker.names=F,chr="9", col =c("purple")) #AR
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Aerial Root Nodes - Chr 9")

#Regression AR
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR <- makeqtl(mapthis, chr=c(1), pos=c(251.5))
#qtl_AR <- makeqtl(mapthis, chr=c(2), pos=c(103))
qtl_AR <- makeqtl(mapthis, chr=c(9), pos=c(97.1))
summary(out.qtl_AR <- fitqtl(mapthis, qtl=qtl_AR, pheno.col=41, formula=y~Q1)) 

#plot CIM blup_noparents PDM
plot(out47_cim.em, show.marker.names=F,col =c("purple")) #PDM
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Dry Mass")
summary(out47_cim.em, threshold=3.91)
plot(out47_cim.em, show.marker.names=F,chr="9", col =c("purple")) #PDM
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Dry Mass - Chr 9")
plot(out47_cim.em, show.marker.names=F,chr="10", col =c("purple")) #PDM
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Dry Mass - Chr 10")

#Regression PDM
mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_PDM <- makeqtl(mapthis, chr=c(9), pos=c(110)) 
qtl_PDM <- makeqtl(mapthis, chr=c(10), pos=c(185))
summary(out.qtl_PDM <- fitqtl(mapthis, qtl=qtl_PDM, pheno.col=47, formula=y~Q1)) 

#plot CIM blup_noparents PTNP
plot(out53_cim.em, show.marker.names=F,col =c("purple")) #PTNP
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage")
summary(out53_cim.em, threshold=3.91)
plot(out53_cim.em, show.marker.names=F,chr="1", col =c("purple")) #PTNP
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 1")
plot(out53_cim.em, show.marker.names=F,chr="9", col =c("purple")) #PTNP
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen Percentage - Chr 9")

#Regression PTNP
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_PTNP <- makeqtl(mapthis, chr=c(1), pos=c(212)) 
#qtl_PTNP <- makeqtl(mapthis, chr=c(9), pos=c(33.8)) 
#summary(out.qtl_PTNP <- fitqtl(mapthis, qtl=qtl_PTNP, pheno.col=53, formula=y~Q1)) 

#plot CIM blup_noparents PTN
plot(out59_cim.em, show.marker.names=F,col =c("purple")) #PTN
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Plant Total Nitrogen")
summary(out59_cim.em, threshold=3.91)

#Regression PTN
#There were no LOD peaks above the threshold.

#plot CIM blup_noparents GDM
plot(out65_cim.em, show.marker.names=F,col =c("purple")) #GDM
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass")
summary(out65_cim.em, threshold=3.91)
plot(out65_cim.em, show.marker.names=F,chr="7", col =c("purple")) #GDM
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Dry Mass - Chr 7")

#Regression GDM
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_GDM <- makeqtl(mapthis, chr=c(6), pos=c(111))
#summary(out.qtl_GDM <- fitqtl(mapthis, qtl=qtl_GDM, pheno.col=65, formula=y~Q1)) 

#plot CIM blup_noparents GTNP
plot(out71_cim.em, show.marker.names=F,col =c("purple")) #GTNP
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage")
summary(out71_cim.em, threshold=3.91)
plot(out71_cim.em, show.marker.names=F,chr="10", col =c("purple")) #GTNP
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen Percentage")

#Regression GTNP
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_GTNP <- makeqtl(mapthis, chr=c(10), pos=c(92.7)) 
#summary(out.qtl_GTNP <- fitqtl(mapthis, qtl=qtl_GTNP, pheno.col=71, formula=y~Q1))

#plot CIM blup_noparents GTN
plot(out77_cim.em, show.marker.names=F,col =c("purple")) #GTN
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen")
summary(out77_cim.em, threshold=3.91)
plot(out77_cim.em, show.marker.names=F,chr="7", col =c("purple")) #GTN
abline(3.91, 0, untf = FALSE)
title(main = "LH82 - Grain Total Nitrogen - Chr 7")

#Regression GTN
#There were no LOD peaks above the threshold.
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_GTN <- makeqtl(mapthis, chr=c(7), pos=c(110)) 
#summary(out.qtl_GTN <- fitqtl(mapthis, qtl=qtl_GTN, pheno.col=77, formula=y~Q1))

#lodint and bayesint for non-refined QTL SC
lodint(out1_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out1_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out2_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out2_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out3_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out3_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out4_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out4_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out5_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out5_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out6_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out6_cim.em, chr=1, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL DPS
lodint(out15_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out15_cim.em, chr=10, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL AR
lodint(out37_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out37_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out37_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out37_cim.em, chr=9, expandtomarkers=TRUE)
lodint(out38_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out38_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out40_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out40_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out40_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out40_cim.em, chr=9, expandtomarkers=TRUE)
lodint(out41_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out41_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out41_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out41_cim.em, chr=9, expandtomarkers=TRUE)
lodint(out42_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out42_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out42_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out42_cim.em, chr=9, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PDM
lodint(out46_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out46_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out47_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out47_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out48_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out48_cim.em, chr=10, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PTNP
lodint(out50_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out50_cim.em, chr=1, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PTN
lodint(out57_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out57_cim.em, chr=9, expandtomarkers=TRUE)

