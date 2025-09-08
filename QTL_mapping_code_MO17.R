#Input library
library(qtl)

#Input file
setwd("~/Documents/Dissertation/Dissertation_Figures/Dissertation_Figures_Ch4/Genotypes_Phenotypes/2019NewFiles/")
mapthis <- read.cross("csv", "./", "genotypes_phenotypes_MO17_052320.csv", estimate.map=FALSE,genotypes=c("A","H","B"))

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
mapthis <- drop.markers(mapthis, c("S1_15075522","S1_15779610","S1_17695130","S1_21543129","S1_23054888","S1_25376989","S1_25876187","S1_27776090","S1_31855785","S1_32088558","S1_32669088","S1_32837944","S1_38243957","S1_39855816","S1_42848244","S1_55067888","S1_58646339","S1_58739851","S1_69587711","S1_69887775","S1_70768832","S1_72022558","S1_93075334","S1_93268640","S1_93970700","S1_95112588","S1_105325600","S1_105951834","S1_107793132","S1_119551660","S1_121545002","S1_121545300","S1_148451802","S1_149321101","S1_164090748","S1_164370587","S1_167030028","S1_182197270","S1_182872518","S1_183404595","S1_183891343","S1_187588875","S1_187589178","S1_188221075","S1_189404821","S1_189588596","S1_195185504","S1_199730706","S1_198278150","S1_199339833","S1_204801636","S1_214725337","S1_216164633","S1_222205601","S1_239208621","S1_239209147","S1_240716598","S1_241001627","S1_241269788","S1_246822629","S1_248751990","S1_251995326","S1_253109100","S1_255638118","S1_257972121","S1_261412168","S1_262512982","S1_262545819","S1_262937156","S1_262998969","S1_267531250","S1_268335761","S1_279565922","S1_279566151","S1_290955698"))
mapthis <- drop.markers(mapthis, c("S2_2120469","S2_3099031","S2_6032038","S2_8287296","S2_9557024","S2_10014310","S2_10786515","S2_14267335","S2_32509770","S2_32951082","S2_33500082","S2_37575697","S2_38705600","S2_54360087","S2_57669712","S2_58365108","S2_60630391","S2_62099974","S2_67931798","S2_83555638","S2_86076528","S2_86217657","S2_122925177","S2_128586630","S2_129307459","S2_136562969","S2_144629347","S2_148879082","S2_151609073","S2_152501327","S2_152505842","S2_157236298","S2_160797696","S2_163677797","S2_164889910","S2_186642391","S2_186993990","S2_193730314","S2_193735700","S2_205944318","S2_206948327","S2_207040859","S2_216756682","S2_216989050","S2_217012177","S2_219156302","S2_219660322","S2_220400686","S2_222708077","S2_222807823","S2_222815122","S2_234229490","S2_234264560","S2_234265011"))
mapthis <- drop.markers(mapthis, c("S3_162587460","S3_191841403"))
mapthis <- drop.markers(mapthis, c("S4_178577456","S5_5378891","S5_6675642","S5_7261348","S5_7151401","S5_7588983","S5_10417884","S5_18156758","S5_21877521","S5_25490067","S5_25883072","S5_29323783","S5_30218507","S5_53005097","S5_53078943","S5_59431530","S5_61801008","S5_68512496","S5_70869281","S5_74451355","S5_76428858","S5_77994061","S5_115585549","S5_117625205","S5_143956427","S5_146885160","S5_147541972","S5_168533223","S5_170023565","S5_176876597"))
mapthis <- drop.markers(mapthis, c("S6_18157470","S6_28146415","S6_31005681","S6_31649763","S6_33638149","S6_34438516","S6_36947508","S6_37026816","S6_43101957","S6_43908386","S6_44824453","S6_45925036","S6_47017613","S6_56452361","S6_58441439","S6_62397710","S6_62449922","S6_63913485","S6_63964388","S6_66896643","S6_67259792","S6_71231565","S6_85608256","S6_87519670","S6_96922599","S6_108121295","S6_108306672","S6_110292139","S6_118827890","S6_125370528","S6_128581183","S6_135883690","S6_138501317","S6_138562829","S6_139874853","S6_146814335","S6_147061959","S6_147393371","S6_147690216","S6_147880398","S6_148593255","S6_148960226","S6_149469557","S6_149528358","S6_150939848","S6_154347667","S6_154893534","S6_155455167","S6_156056819","S6_156552473","S6_156752063","S6_156762280","S6_158190931","S6_158480726","S6_161132446","S6_164776510"))
mapthis <- drop.markers(mapthis, c("S7_41764479","S7_48162918","S7_85579298","S7_85579750","S7_101840213","S7_101842494","S7_101886029","S7_101902389","S7_105370754","S7_100524503","S7_107455959","S7_107485132","S7_108535010","S7_108535101","S7_116288791","S7_127298234","S7_127717098","S7_135204931","S7_135892383","S7_139724362","S7_141010338","S7_141011084","S7_147159640","S7_150040595","S7_150168954","S7_150176978","S7_155354972","S7_156713162","S7_158155290","S7_162309349","S7_165301200","S7_168388199","S7_170944406","S7_172718187","S7_173517584"))
mapthis <- drop.markers(mapthis, c("S8_29379714","S8_51064820","S8_60877596","S8_62836711","S8_111803807","S8_111824549","S8_116854710","S8_121798115","S8_124951681","S8_125294941","S8_128343923","S8_128428344","S8_128555667","S8_132879659","S8_133275925","S8_133948519","S8_134630475","S8_139713043","S8_143156918","S8_151919897","S8_156824249","S8_158573623","S8_160533569","S8_163214768","S8_163903758","S8_165582258","S8_166714602","S8_167603518"))
mapthis <- drop.markers(mapthis, c("S9_13821933","S9_16117941","S9_20338165","S9_20910403","S9_20906846","S9_26478204","S9_26895432","S9_27983488","S9_31283111","S9_31480908","S9_48738880","S9_55347766","S9_56143183","S9_59160364","S9_59453925","S9_59649116","S9_61814159","S9_61814464","S9_73475902","S9_77307438","S9_84224710","S9_87360496","S9_88648911","S9_88651683","S9_88734963","S9_100190820","S9_100789121","S9_101099665","S9_109668102","S9_109740653","S9_120359864","S9_123637645","S9_124149942","S9_124160253","S9_124959335","S9_127449846","S9_127646282","S9_128094690","S9_129185786","S9_129767667","S9_133546700","S9_139518470","S9_146696738","S9_147962492","S9_149382485"))
mapthis <- drop.markers(mapthis, c("S10_7156521","S10_12751979","S10_18399723","S10_31014426","S10_49982707","S10_51007457","S10_52223901","S10_66767570","S10_71471618","S10_71672350","S10_71810673","S10_72980094","S10_75809949","S10_84239903","S10_96251345","S10_111288191","S10_122252640","S10_125214335","S10_125951324","S10_130278659","S10_132702636","S10_138027125","S10_138469918","S10_139081560","S10_139877720","S10_146165279","S10_147811484","S10_147997647"))

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
lg <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=5)
table(lg[,2])

#Form linkage groups in actual data
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=5, reorgMarkers=TRUE)

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

mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 8)

mn <- markernames(mapthis, chr=c(12)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 9)

mn <- markernames(mapthis, chr=c(13)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)

mapthis <- drop.markers(mapthis, c("S4_29167364"))
mapthis <- drop.markers(mapthis, c("S3_227190097"))

#move markers
mn <- markernames(mapthis, chr=c(3)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(9)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 3)
mn <- markernames(mapthis, chr=c(4)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 9)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 4)
mn <- markernames(mapthis, chr=c(5)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)
mn <- markernames(mapthis, chr=c(6)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 5)
mn <- markernames(mapthis, chr=c(11)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 6)
mn <- markernames(mapthis, chr=c(7)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 11)
mn <- markernames(mapthis, chr=c(10)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 7)
mn <- markernames(mapthis, chr=c(8)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 10)
mn <- markernames(mapthis, chr=c(9)) 
for(marker in mn) 
  mapthis <- movemarker(mapthis, marker, 8)
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
mapthis <- switch.order(mapthis, chr=1, c(1,	2,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13,	14,	17,	15,	16,	100,	99,	101,	102,	98,	97,	96,	95,	94,	93,	92,	91,	90,	88,	89,	87,	86,	85,	84,	83,	82,	81,	80,	79,	78,	77,	76,	75,	74,	73,	72,	71,	70,	69,	66,	68,	67,	64,	65,	63,	62,	61,	60,	59,	58,	57,	55,	56,	52,	54,	53,	51,	50,	49,	48,	47,	46,	45,	44,	43,	42,	41,	40,	39,	38,	37,	36,	35,	34,	33,	18,	19,	20,	32,	31,	21,	30,	29,	28,	27,	26,	22,	23,	24,	25, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=2, c(68,	64,	67,	65,	66,	63,	62,	61,	60,	59,	58,	57,	56,	55,	54,	53,	52,	51,	50,	48,	49,	47,	46,	45,	44,	43,	42,	41,	40,	39,	38,	37,	35,	36,	34,	33,	32,	31,	30,	29,	28,	27,	26,	25,	24,	23,	22,	20,	21,	19,	18,	17,	16,	15,	14,	13,	12,	11,	10,	9,	8,	7,	6,	2,	4,	5,	3,	1, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=3, c(1,	5,	3,	6,	7,	14,	15,	16,	17,	18,	19,	20,	21,	22,	24,	23,	13,	12,	11,	10,	9,	8,	4,	2, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=4, c(1:9,11,12,10,13,14, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=5, c(1,	2,	3,	4,	5,	6,	7,	8,	9,	10,	11,	12,	13,	14,	15,	16,	17,	18,	19,	20,	21,	22,	23,	24,	25,	26,	27,	28,	29,	30,	31,	32,	34,	35,	37,	39,	40,	33,	36,	38,	46,	45,	41,	43,	42,	44, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=6, c(63,	64,	61,	62,	60,	59,	58,	57,	56,	55,	54,	53,	52,	51,	50,	49,	48,	47,	46,	45,	44,	43,	42,	41,	39,	40,	38,	37,	36,	35,	34,	33,	32,	31,	30,	29,	28,	25,	27,	26,	24,	23,	22,	21,	20,	19,	18,	17,	16,	15,	14,	13,	12,	11,	10,	9,	8,	7,	6,	5,	3,	4,	1,	2, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=7, c(48,45,44,46,47,43:4,1,3,2, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=8, c(3,2,4,1,5:47,49,48,50:52, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=9, c(42,	45,	44,	46,	43,	41,	40,	39,	38,	37,	36,	35,	34,	33,	32,	31,	30,	29,	28,	20,	21,	4,	3,	27,	26,	25,	24,	1,	2,	5,	23,	6,	7,	8,	22,	19,	9,	18,	17,	16,	15,	10,	14,	12,	11,	13, error.prob=0.005))
mapthis <- switch.order(mapthis, chr=10, c(1,	2,	3,	4,	6,	5,	7,	8,	9,	10,	11,	12,	32,	13,	14,	15,	30,	16,	29,	31,	17,	28,	27,	26,	25,	24,	23,	22,	21,	18,	20,	19,	33,	34,	35,	36,	37,	38,	39,	40,	41,	42,	43,	44,	45,	46,	47,	48,	49, error.prob=0.005))

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
out.mr <- scanone(mapthis, method="mr",pheno.col=7)
summary(out.mr, threshold=3)
max(out.mr)
out.mr[ out.mr$chr == 3, ]

#plot single analysis
plot(out.mr, show.marker.names=F)
plot(out.mr, show.marker.names=F, chr=c("3"))

#SIM
out.em <- scanone (mapthis, pheno.col=c(6),method="em")
summary(out.em, threshold=3)
max(out.em)
out.em[ out.em$chr == 3, ]
plot(out.em, lodcolumn=c(1), lty=1, col=c("blue"), show.marker.names=F)
plot(out.em, lodcolumn=c(1),chr=c("3"), lty=1, col=c("blue"), show.marker.names=F)

#CIM
out_cim.em <- cim(mapthis, pheno.col = (41), n.marcovar=1, method=c("em"), map.function=c("kosambi") )
summary(out_cim.em, threshold=3)
max(out_cim.em)
out_cim.em[ out_cim.em$chr == 5, ]
plot(out_cim.em, show.marker.names=F, col =c("red"))
plot(out_cim.em, chr=c("4"), show.marker.names=F, col =c("red"))
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
plot(out.mr, out.em, out_cim.em, chr="3",show.marker.names=F,col =c("Black", "blue", "red"))

out.mr <- scanone(mapthis, method="mr",pheno.col=6)
summary(out.mr, threshold=2)
out.em <- scanone (mapthis, pheno.col= (6),method="em")
summary(out.em, threshold=2)
out_cim.em <- cim(mapthis, pheno.col = (6), n.marcovar=3, method=c("em"), map.function=c("kosambi") )
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
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count")
plot(out5_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count")

#plot CIM SC chrs

plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="6", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 6")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="6", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 6")
plot(out5_cim.em, show.marker.names=F,chr="6", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 6")
plot(out1_cim.em, out2_cim.em, out3_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 9")
plot(out4_cim.em, out5_cim.em, out6_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 9")
plot(out5_cim.em, show.marker.names=F,chr="9", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Stand Count - Chr 9")

#SC QTL chr3

summary(out1_cim.em, threshold=3.84)
summary(out2_cim.em, threshold=3.84)
summary(out3_cim.em, threshold=3.84)
summary(out4_cim.em, threshold=3.84)
summary(out5_cim.em, threshold=3.84)
summary(out6_cim.em, threshold=3.84)

#Regression SC

#Rep 1
summary(out2_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_rep1_chr6 <- makeqtl(mapthis, chr=c(6), pos=c(26.9)) 
summary(out.qtl_SC_rep1_chr6 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr6, pheno.col=2, formula=y~Q1))
plot(qtl_SC_rep1_chr6)

#fit the model (with QTL in fixed positions)
out.fq_SC_rep1_chr6 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr6, formula=y~Q1)
summary(out.fq_SC_rep1_chr6)

#refine QTL and plot LOD profile
rqtl_SC_rep1_chr6 <- refineqtl(mapthis, qtl=qtl_SC_rep1_chr6, formula=y~Q1, verbose=FALSE)
rqtl_SC_rep1_chr6
out.fq2_SC_rep1_chr6 <- fitqtl(mapthis, qtl=qtl_SC_rep1_chr6, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_rep1_chr6)
plotLodProfile(rqtl_SC_rep1_chr6)
abline(3.84, 0, untf = FALSE)
#abline(2.79, 0, untf = FALSE, col=c("dark gray"))
abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Stand Count (Rep1)")

#blup no parents and family as covariates
summary(out6_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC_blupNPFC_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(111)) 
summary(out.qtl_SC_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr9, pheno.col=6, formula=y~Q1))
plot(qtl_SC_blupNPFC_chr9)

#fit the model (with QTL in fixed positions)
out.fq_SC_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr9, formula=y~Q1)
summary(out.fq_SC_blupNPFC_chr9)

#refine QTL and plot LOD profile
rqtl_SC_blupNPFC_chr9 <- refineqtl(mapthis, qtl=qtl_SC_blupNPFC_chr9, formula=y~Q1, verbose=FALSE)
rqtl_SC_blupNPFC_chr9
out.fq2_SC_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_SC_blupNPFC_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_SC_blupNPFC_chr9)
plotLodProfile(rqtl_SC_blupNPFC_chr9)
abline(3.84, 0, untf = FALSE)
#abline(2.79, 0, untf = FALSE, col=c("dark gray"))
abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Stand Count (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod support interval coordinates
lodint(rqtl_SC_rep1_chr6, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_SC_blupNPFC_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_SC_rep1_chr6, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_SC_blupNPFC_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)


#plot CIM DT
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence")
plot(out11_cim.em,  show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence")

#plot CIM DT chrs
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="3", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 3")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="3", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 3")
plot(out11_cim.em, show.marker.names=F,chr="3", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 3")
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 8")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 8")
plot(out11_cim.em, show.marker.names=F,chr="8", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 8")
plot(out7_cim.em, out8_cim.em, out9_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 10")
plot(out10_cim.em, out11_cim.em, out12_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 10")
plot(out11_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence - Chr 10")

summary(out7_cim.em, threshold=3.84)
summary(out8_cim.em, threshold=3.84)
summary(out9_cim.em, threshold=3.84)
summary(out10_cim.em, threshold=3.84)
summary(out11_cim.em, threshold=3.84)
summary(out12_cim.em, threshold=3.84)

#Regression DT
#Average
summary(out7_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_avg <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
summary(out.qtl_DT_avg <- fitqtl(mapthis, qtl=qtl_DT_avg, pheno.col=7, formula=y~Q1))
plot(qtl_DT_avg)

#fit the model (with QTL in fixed positions)
out.fq_DT_avg <- fitqtl(mapthis, qtl=qtl_DT_avg, formula=y~Q1)
summary(out.fq_DT_avg)

#refine QTL and plot LOD profile
rqtl_DT_avg <- refineqtl(mapthis, qtl=qtl_DT_avg, formula=y~Q1, verbose=FALSE)
rqtl_DT_avg
out.fq2_DT_avg <- fitqtl(mapthis, qtl=qtl_DT_avg, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_avg)
plotLodProfile(rqtl_DT_avg)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Average Rep1 + Rep2)")

#Rep 1
summary(out8_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_rep1 <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
summary(out.qtl_DT_rep1 <- fitqtl(mapthis, qtl=qtl_DT_rep1, pheno.col=8, formula=y~Q1))
plot(qtl_DT_rep1)

#fit the model (with QTL in fixed positions)
out.fq_DT_rep1 <- fitqtl(mapthis, qtl=qtl_DT_rep1, formula=y~Q1)
summary(out.fq_DT_rep1)

#refine QTL and plot LOD profile
rqtl_DT_rep1 <- refineqtl(mapthis, qtl=qtl_DT_rep1, formula=y~Q1, verbose=FALSE)
rqtl_DT_rep1
out.fq2_DT_rep1 <- fitqtl(mapthis, qtl=qtl_DT_rep1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_rep1)
plotLodProfile(rqtl_DT_rep1)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Rep1)")

#Rep 2
summary(out9_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_rep2_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(158)) 
qtl_DT_rep2_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(122))
qtl_DT_rep2 <- makeqtl(mapthis, chr=c(8, 10), pos=c(158, 122))
summary(out.qtl_DT_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr8, pheno.col=9, formula=y~Q1))
summary(out.qtl_DT_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr10, pheno.col=9, formula=y~Q1))
summary(out.qtl_DT_rep2 <- fitqtl(mapthis, qtl=qtl_DT_rep2, pheno.col=9, formula=y~Q1+Q2))
plot(qtl_DT_rep2_chr8)
plot(qtl_DT_rep2_chr10)
plot(qtl_DT_rep2)

#fit the model (with QTL in fixed positions)
out.fq_DT_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr8, formula=y~Q1)
summary(out.fq_DT_rep2_chr8)
out.fq_DT_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr10, formula=y~Q1)
summary(out.fq_DT_rep2_chr10)
out.fq_DT_rep2 <- fitqtl(mapthis, qtl=qtl_DT_rep2, formula=y~Q1+Q2)
summary(out.fq_DT_rep2)

#refine QTL and plot LOD profile
rqtl_DT_rep2_chr8 <- refineqtl(mapthis, qtl=qtl_DT_rep2_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DT_rep2_chr8
out.fq2_DT_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_rep2_chr8)
plotLodProfile(rqtl_DT_rep2_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Rep2)")

rqtl_DT_rep2_chr10 <- refineqtl(mapthis, qtl=qtl_DT_rep2_chr10, formula=y~Q1, verbose=FALSE)
rqtl_DT_rep2_chr10
out.fq2_DT_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DT_rep2_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_rep2_chr10)
plotLodProfile(rqtl_DT_rep2_chr10)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Rep2)")

rqtl_DT_rep2 <- refineqtl(mapthis, qtl=qtl_DT_rep2, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DT_rep2
out.fq2_DT_rep2 <- fitqtl(mapthis, qtl=qtl_DT_rep2, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DT_rep2)
plotLodProfile(rqtl_DT_rep2)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Rep2)")

#blup
summary(out10_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_blup <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
summary(out.qtl_DT_blup <- fitqtl(mapthis, qtl=qtl_DT_blup, pheno.col=10, formula=y~Q1))
plot(qtl_DT_blup)

#fit the model (with QTL in fixed positions)
out.fq_DT_blup <- fitqtl(mapthis, qtl=qtl_DT_blup, formula=y~Q1)
summary(out.fq_DT_blup)

#refine QTL and plot LOD profile
rqtl_DT_blup <- refineqtl(mapthis, qtl=qtl_DT_blup, formula=y~Q1, verbose=FALSE)
rqtl_DT_blup
out.fq2_DT_blup <- fitqtl(mapthis, qtl=qtl_DT_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_blup)
plotLodProfile(rqtl_DT_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out11_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_blupNP <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
summary(out.qtl_DT_blupNP <- fitqtl(mapthis, qtl=qtl_DT_blupNP, pheno.col=11, formula=y~Q1))
plot(qtl_DT_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_DT_blupNP <- fitqtl(mapthis, qtl=qtl_DT_blupNP, formula=y~Q1)
summary(out.fq_DT_blupNP)

#refine QTL and plot LOD profile
rqtl_DT_blupNP <- refineqtl(mapthis, qtl=qtl_DT_blupNP, formula=y~Q1, verbose=FALSE)
rqtl_DT_blupNP
out.fq2_DT_blupNP <- fitqtl(mapthis, qtl=qtl_DT_blupNP, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_blupNP)
plotLodProfile(rqtl_DT_blupNP)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariate
summary(out12_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DT_blupNPFC <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
summary(out.qtl_DT_blupNPFC <- fitqtl(mapthis, qtl=qtl_DT_blupNPFC, pheno.col=12, formula=y~Q1))
plot(qtl_DT_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_DT_blupNPFC <- fitqtl(mapthis, qtl=qtl_DT_blupNPFC, formula=y~Q1)
summary(out.fq_DT_blupNPFC)

#refine QTL and plot LOD profile
rqtl_DT_blupNPFC <- refineqtl(mapthis, qtl=qtl_DT_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_DT_blupNPFC
out.fq2_DT_blupNPFC <- fitqtl(mapthis, qtl=qtl_DT_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DT_blupNPFC)
plotLodProfile(rqtl_DT_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Tassel Emergence (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_DT_avg, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_rep1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_rep2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_rep2_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_rep2_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_blupNP, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DT_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_DT_avg, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_rep1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_rep2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_rep2_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_rep2_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_blupNP, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DT_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM DPS
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed")
plot(out17_cim.em,  show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed")

#plot CIM DPS chrs
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="3", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 3")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="3", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 3")
plot(out17_cim.em, show.marker.names=F,chr="3", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 3")
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 8")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 8")
plot(out17_cim.em, show.marker.names=F,chr="8", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 8")
plot(out13_cim.em, out14_cim.em, out15_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 10")
plot(out16_cim.em, out17_cim.em, out18_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 10")
plot(out17_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Days to Pollen Shed - Chr 10")

summary(out13_cim.em, threshold=3.84)
summary(out14_cim.em, threshold=3.84)
summary(out15_cim.em, threshold=3.84)
summary(out16_cim.em, threshold=3.84)
summary(out17_cim.em, threshold=3.84)
summary(out18_cim.em, threshold=3.84)

#Regression DPS
#Average
summary(out13_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_avg_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS_avg_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(156))
qtl_DPS_avg <- makeqtl(mapthis, chr=c(3, 8), pos=c(121, 156))
summary(out.qtl_DPS_avg_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr3, pheno.col=13, formula=y~Q1))
summary(out.qtl_DPS_avg_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr8, pheno.col=13, formula=y~Q1))
summary(out.qtl_DPS_avg <- fitqtl(mapthis, qtl=qtl_DPS_avg, pheno.col=13, formula=y~Q1+Q2))
plot(qtl_DPS_avg_chr3)
plot(qtl_DPS_avg_chr8)
plot(qtl_DPS_avg)

#fit the model (with QTL in fixed positions)
out.fq_DPS_avg_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr3, formula=y~Q1)
summary(out.fq_DPS_avg_chr3)
out.fq_DPS_avg_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr8, formula=y~Q1)
summary(out.fq_DPS_avg_chr8)
out.fq_DPS_avg <- fitqtl(mapthis, qtl=qtl_DPS_avg, formula=y~Q1+Q2)
summary(out.fq_DPS_avg)

#refine QTL and plot LOD profile
rqtl_DPS_avg_chr3 <- refineqtl(mapthis, qtl=qtl_DPS_avg_chr3, formula=y~Q1, verbose=FALSE)
rqtl_DPS_avg_chr3
out.fq2_DPS_avg_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr3, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_avg_chr3)
plotLodProfile(rqtl_DPS_avg_chr3)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2)")

rqtl_DPS_avg_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_avg_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_avg_chr8
out.fq2_DPS_avg_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_avg_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_avg_chr8)
plotLodProfile(rqtl_DPS_avg_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2)")

rqtl_DPS_avg <- refineqtl(mapthis, qtl=qtl_DPS_avg, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_avg
out.fq2_DPS_avg <- fitqtl(mapthis, qtl=qtl_DPS_avg, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_avg)
plotLodProfile(rqtl_DPS_avg)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2)")

#Rep 1
summary(out14_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_rep1_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS_rep1_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(156))
qtl_DPS_rep1 <- makeqtl(mapthis, chr=c(3, 8), pos=c(121, 156))
summary(out.qtl_DPS_rep1_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr3, pheno.col=14, formula=y~Q1))
summary(out.qtl_DPS_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr8, pheno.col=14, formula=y~Q1))
summary(out.qtl_DPS_rep1 <- fitqtl(mapthis, qtl=qtl_DPS_rep1, pheno.col=14, formula=y~Q1+Q2))
plot(qtl_DPS_rep1_chr3)
plot(qtl_DPS_rep1_chr8)
plot(qtl_DPS_rep1)

#fit the model (with QTL in fixed positions)
out.fq_DPS_rep1_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr3, formula=y~Q1)
summary(out.fq_DPS_rep1_chr3)
out.fq_DPS_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr8, formula=y~Q1)
summary(out.fq_DPS_rep1_chr8)
out.fq_DPS_rep1 <- fitqtl(mapthis, qtl=qtl_DPS_rep1, formula=y~Q1+Q2)
summary(out.fq_DPS_rep1)

#refine QTL and plot LOD profile
rqtl_DPS_rep1_chr3 <- refineqtl(mapthis, qtl=qtl_DPS_rep1_chr3, formula=y~Q1, verbose=FALSE)
rqtl_DPS_rep1_chr3
out.fq2_DPS_rep1_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr3, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_rep1_chr3)
plotLodProfile(rqtl_DPS_rep1_chr3)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep1)")

rqtl_DPS_rep1_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_rep1_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_rep1_chr8
out.fq2_DPS_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep1_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_rep1_chr8)
plotLodProfile(rqtl_DPS_rep1_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep1)")

rqtl_DPS_rep1 <- refineqtl(mapthis, qtl=qtl_DPS_rep1, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_rep1
out.fq2_DPS_rep1 <- fitqtl(mapthis, qtl=qtl_DPS_rep1, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_rep1)
plotLodProfile(rqtl_DPS_rep1)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep1)")

#Rep 2
summary(out15_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_rep2_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(158))
qtl_DPS_rep2_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(122)) 
qtl_DPS_rep2 <- makeqtl(mapthis, chr=c(8,10), pos=c(158, 122))
summary(out.qtl_DPS_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr8, pheno.col=15, formula=y~Q1))
summary(out.qtl_DPS_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr10, pheno.col=15, formula=y~Q1))
summary(out.qtl_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, pheno.col=15, formula=y~Q1+Q2))
plot(qtl_DPS_rep2_chr8)
plot(qtl_DPS_rep2_chr10)
plot(qtl_DPS_rep2)

#fit the model (with QTL in fixed positions)
out.fq_DPS_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr8, formula=y~Q1)
summary(out.fq_DPS_rep2_chr8)
out.fq_DPS_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr10, formula=y~Q1)
summary(out.fq_DPS_rep2_chr10)
out.fq_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1+Q2)
summary(out.fq_DPS_rep2)

#refine QTL and plot LOD profile
rqtl_DPS_rep2_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_rep2_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_rep2_chr8
out.fq2_DPS_rep2_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_rep2_chr8)
plotLodProfile(rqtl_DPS_rep2_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep2)")

rqtl_DPS_rep2_chr10 <- refineqtl(mapthis, qtl=qtl_DPS_rep2_chr10, formula=y~Q1, verbose=FALSE)
rqtl_DPS_rep2_chr10
out.fq2_DPS_rep2_chr10 <- fitqtl(mapthis, qtl=qtl_DPS_rep2_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_rep2_chr10)
plotLodProfile(rqtl_DPS_rep2_chr10)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep2)")

rqtl_DPS_rep2 <- refineqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_rep2
out.fq2_DPS_rep2 <- fitqtl(mapthis, qtl=qtl_DPS_rep2, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_rep2)
plotLodProfile(rqtl_DPS_rep2)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Rep2)")

#blup
summary(out16_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_blup_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS_blup_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(156))
qtl_DPS_blup <- makeqtl(mapthis, chr=c(3, 8), pos=c(121, 156))
summary(out.qtl_DPS_blup_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr3, pheno.col=16, formula=y~Q1))
summary(out.qtl_DPS_blup_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr8, pheno.col=16, formula=y~Q1))
summary(out.qtl_DPS_blup <- fitqtl(mapthis, qtl=qtl_DPS_blup, pheno.col=16, formula=y~Q1+Q2))
plot(qtl_DPS_blup_chr3)
plot(qtl_DPS_blup_chr8)
plot(qtl_DPS_blup)

#fit the model (with QTL in fixed positions)
out.fq_DPS_blup_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr3, formula=y~Q1)
summary(out.fq_DPS_blup_chr3)
out.fq_DPS_blup_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr8, formula=y~Q1)
summary(out.fq_DPS_blup_chr8)
out.fq_DPS_blup <- fitqtl(mapthis, qtl=qtl_DPS_blup, formula=y~Q1+Q2)
summary(out.fq_DPS_blup)

#refine QTL and plot LOD profile
rqtl_DPS_blup_chr3 <- refineqtl(mapthis, qtl=qtl_DPS_blup_chr3, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blup_chr3
out.fq2_DPS_blup_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr3, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blup_chr3)
plotLodProfile(rqtl_DPS_blup_chr3)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP))")

rqtl_DPS_blup_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_blup_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blup_chr8
out.fq2_DPS_blup_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blup_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blup_chr8)
plotLodProfile(rqtl_DPS_blup_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP))")

rqtl_DPS_blup <- refineqtl(mapthis, qtl=qtl_DPS_blup, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_blup
out.fq2_DPS_blup <- fitqtl(mapthis, qtl=qtl_DPS_blup, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_blup)
plotLodProfile(rqtl_DPS_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out17_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_blupNP_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS_blupNP_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(156))
qtl_DPS_blupNP <- makeqtl(mapthis, chr=c(3, 8), pos=c(121, 156))
summary(out.qtl_DPS_blupNP_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr3, pheno.col=17, formula=y~Q1))
summary(out.qtl_DPS_blupNP_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr8, pheno.col=17, formula=y~Q1))
summary(out.qtl_DPS_blupNP <- fitqtl(mapthis, qtl=qtl_DPS_blupNP, pheno.col=17, formula=y~Q1+Q2))
plot(qtl_DPS_blupNP_chr3)
plot(qtl_DPS_blupNP_chr8)
plot(qtl_DPS_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_DPS_blupNP_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr3, formula=y~Q1)
summary(out.fq_DPS_blupNP_chr3)
out.fq_DPS_blupNP_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr8, formula=y~Q1)
summary(out.fq_DPS_blupNP_chr8)
out.fq_DPS_blupNP <- fitqtl(mapthis, qtl=qtl_DPS_blupNP, formula=y~Q1+Q2)
summary(out.fq_DPS_blupNP)

#refine QTL and plot LOD profile
rqtl_DPS_blupNP_chr3 <- refineqtl(mapthis, qtl=qtl_DPS_blupNP_chr3, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blupNP_chr3
out.fq2_DPS_blupNP_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr3, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blupNP_chr3)
plotLodProfile(rqtl_DPS_blupNP_chr3)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NP))")

rqtl_DPS_blupNP_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_blupNP_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blupNP_chr8
out.fq2_DPS_blupNP_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNP_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blupNP_chr8)
plotLodProfile(rqtl_DPS_blupNP_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NP))")

rqtl_DPS_blupNP <- refineqtl(mapthis, qtl=qtl_DPS_blupNP, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_blupNP
out.fq2_DPS_blupNP <- fitqtl(mapthis, qtl=qtl_DPS_blupNP, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_blupNP)
plotLodProfile(rqtl_DPS_blupNP)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariate
summary(out18_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS_blupNPFC_chr3 <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS_blupNPFC_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(156))
qtl_DPS_blupNPFC <- makeqtl(mapthis, chr=c(3, 8), pos=c(121, 156))
summary(out.qtl_DPS_blupNPFC_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr3, pheno.col=18, formula=y~Q1))
summary(out.qtl_DPS_blupNPFC_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr8, pheno.col=18, formula=y~Q1))
summary(out.qtl_DPS_blupNPFC <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC, pheno.col=18, formula=y~Q1+Q2))
plot(qtl_DPS_blupNPFC_chr3)
plot(qtl_DPS_blupNPFC_chr8)
plot(qtl_DPS_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_DPS_blupNPFC_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr3, formula=y~Q1)
summary(out.fq_DPS_blupNPFC_chr3)
out.fq_DPS_blupNPFC_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr8, formula=y~Q1)
summary(out.fq_DPS_blupNPFC_chr8)
out.fq_DPS_blupNPFC <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC, formula=y~Q1+Q2)
summary(out.fq_DPS_blupNPFC)

#refine QTL and plot LOD profile
rqtl_DPS_blupNPFC_chr3 <- refineqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr3, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blupNPFC_chr3
out.fq2_DPS_blupNPFC_chr3 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr3, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blupNPFC_chr3)
plotLodProfile(rqtl_DPS_blupNPFC_chr3)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NPFC))")

rqtl_DPS_blupNPFC_chr8 <- refineqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr8, formula=y~Q1, verbose=FALSE)
rqtl_DPS_blupNPFC_chr8
out.fq2_DPS_blupNPFC_chr8 <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_DPS_blupNPFC_chr8)
plotLodProfile(rqtl_DPS_blupNPFC_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NPFC))")

rqtl_DPS_blupNPFC <- refineqtl(mapthis, qtl=qtl_DPS_blupNPFC, formula=y~Q1+Q2, verbose=FALSE)
rqtl_DPS_blupNPFC
out.fq2_DPS_blupNPFC <- fitqtl(mapthis, qtl=qtl_DPS_blupNPFC, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_DPS_blupNPFC)
plotLodProfile(rqtl_DPS_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Days to Pollen Shed (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_DPS_avg_chr3, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_avg_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_avg, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep1_chr3, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep1_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep2_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep2_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_rep2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blup_chr3, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blup_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNP_chr3, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNP_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNP, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNPFC_chr3, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNPFC_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_DPS_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_DPS_avg_chr3, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_avg_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_avg, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep1_chr3, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep1_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep2_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep2_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_rep2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blup_chr3, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blup_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNP_chr3, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNP_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNP, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNPFC_chr3, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNPFC_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_DPS_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM 15NT1
plot(out19_cim.em, out20_cim.em, out21_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 1")
plot(out22_cim.em, out23_cim.em, out24_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 1")
plot(out23_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 1")

#plot CIM 15NT1 chrs

summary(out19_cim.em, threshold=3.84)
summary(out20_cim.em, threshold=3.84)
summary(out21_cim.em, threshold=3.84)
summary(out22_cim.em, threshold=3.84)
summary(out23_cim.em, threshold=3.84)
summary(out24_cim.em, threshold=3.84)

#plot CIM 15NT2
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2")
plot(out29_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2")

#plot CIM 15NT2 chrs
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,chr="3", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 3")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,chr="3", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 3")
plot(out29_cim.em, show.marker.names=F,chr="3", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 3")
plot(out25_cim.em, out26_cim.em, out27_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 9")
plot(out28_cim.em, out29_cim.em, out30_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 9")
plot(out29_cim.em, show.marker.names=F,chr="9", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 2 - Chr 9")

summary(out25_cim.em, threshold=3.84)
summary(out26_cim.em, threshold=3.84)
summary(out27_cim.em, threshold=3.84)
summary(out28_cim.em, threshold=3.84)
summary(out29_cim.em, threshold=3.84)
summary(out30_cim.em, threshold=3.84)

#plot CIM 15NT3
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3")
plot(out35_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3")

#plot CIM 15NT3 chrs
plot(out31_cim.em, out32_cim.em, out33_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3 - Chr 2")
plot(out34_cim.em, out35_cim.em, out36_cim.em, show.marker.names=F,chr="2", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3 - Chr 2")
plot(out35_cim.em, show.marker.names=F,chr="2", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - 15N Time-point 3 - Chr 2")

summary(out31_cim.em, threshold=3.84)
summary(out32_cim.em, threshold=3.84)
summary(out33_cim.em, threshold=3.84)
summary(out34_cim.em, threshold=3.84)
summary(out35_cim.em, threshold=3.84)
summary(out36_cim.em, threshold=3.84)

#blup
summary(out34_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_15NT3_blup <- makeqtl(mapthis, chr=c(2), pos=c(150)) 
summary(out.qtl_15NT3_blup <- fitqtl(mapthis, qtl=qtl_15NT3_blup, pheno.col=34, formula=y~Q1))
plot(qtl_15NT3_blup)

#fit the model (with QTL in fixed positions)
out.fq_15NT3_blup <- fitqtl(mapthis, qtl=qtl_15NT3_blup, formula=y~Q1)
summary(out.fq_15NT3_blup)

#refine QTL and plot LOD profile
rqtl_15NT3_blup <- refineqtl(mapthis, qtl=qtl_15NT3_blup, formula=y~Q1, verbose=FALSE)
rqtl_15NT3_blup
out.fq2_15NT3_blup <- fitqtl(mapthis, qtl=qtl_15NT3_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_15NT3_blup)
plotLodProfile(rqtl_15NT3_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - 15N Time-point 3 (Average Rep1 + Rep2, BLUP)")

#blup NP
summary(out35_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_15NT3_blupNP <- makeqtl(mapthis, chr=c(2), pos=c(150)) 
summary(out.qtl_15NT3_blupNP <- fitqtl(mapthis, qtl=qtl_15NT3_blupNP, pheno.col=35, formula=y~Q1))
plot(qtl_15NT3_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_15NT3_blupNP <- fitqtl(mapthis, qtl=qtl_15NT3_blupNP, formula=y~Q1)
summary(out.fq_15NT3_blupNP)

#refine QTL and plot LOD profile
rqtl_15NT3_blupNP <- refineqtl(mapthis, qtl=qtl_15NT3_blupNP, formula=y~Q1, verbose=FALSE)
rqtl_15NT3_blupNP
out.fq2_15NT3_blupNP <- fitqtl(mapthis, qtl=qtl_15NT3_blupNP, formula=y~Q1, dropone=FALSE)
summary(out.fq2_15NT3_blupNP)
plotLodProfile(rqtl_15NT3_blupNP)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - 15N Time-point 3 (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariate
summary(out36_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_15NT3_blupNPFC <- makeqtl(mapthis, chr=c(2), pos=c(150)) 
summary(out.qtl_15NT3_blupNPFC <- fitqtl(mapthis, qtl=qtl_15NT3_blupNPFC, pheno.col=36, formula=y~Q1))
plot(qtl_15NT3_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_15NT3_blupNPFC <- fitqtl(mapthis, qtl=qtl_15NT3_blupNPFC, formula=y~Q1)
summary(out.fq_15NT3_blupNPFC)

#refine QTL and plot LOD profile
rqtl_15NT3_blupNPFC <- refineqtl(mapthis, qtl=qtl_15NT3_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_15NT3_blupNPFC
out.fq2_15NT3_blupNPFC <- fitqtl(mapthis, qtl=qtl_15NT3_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_15NT3_blupNPFC)
plotLodProfile(rqtl_15NT3_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - 15N Time-point 3 (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_15NT3_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_15NT3_blupNP, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_15NT3_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_15NT3_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_15NT3_blupNP, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_15NT3_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM aerial roots
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes")
plot(out41_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes")

#plot CIM AR chrs
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="7", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 7")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="7", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 7")
plot(out41_cim.em, show.marker.names=F,chr="7", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 7")
plot(out37_cim.em, out38_cim.em, out39_cim.em, show.marker.names=F,chr="9", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 9")
plot(out40_cim.em, out41_cim.em, out42_cim.em, show.marker.names=F,chr="9", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 9")
plot(out41_cim.em, show.marker.names=F,chr="9", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Aerial Root Nodes - Chr 9")

summary(out37_cim.em, threshold=3.84)
summary(out38_cim.em, threshold=3.84)
summary(out39_cim.em, threshold=3.84)
summary(out40_cim.em, threshold=3.84)
summary(out41_cim.em, threshold=3.84)
summary(out42_cim.em, threshold=3.84)

#Regression AR
#Rep 1
summary(out38_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_rep1 <- makeqtl(mapthis, chr=c(7), pos=c(191)) 
summary(out.qtl_AR_rep1 <- fitqtl(mapthis, qtl=qtl_AR_rep1, pheno.col=38, formula=y~Q1))
plot(qtl_AR_rep1)

#fit the model (with QTL in fixed positions)
out.fq_AR_rep1 <- fitqtl(mapthis, qtl=qtl_AR_rep1, formula=y~Q1)
summary(out.fq_AR_rep1)

#refine QTL and plot LOD profile
rqtl_AR_rep1 <- refineqtl(mapthis, qtl=qtl_AR_rep1, formula=y~Q1, verbose=FALSE)
rqtl_AR_rep1
out.fq2_AR_rep1 <- fitqtl(mapthis, qtl=qtl_AR_rep1, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_rep1)
plotLodProfile(rqtl_AR_rep1)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Rep1)")

#blup
summary(out40_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_blup <- makeqtl(mapthis, chr=c(9), pos=c(122)) 
summary(out.qtl_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, pheno.col=40, formula=y~Q1))
plot(qtl_AR_blup)

#fit the model (with QTL in fixed positions)
out.fq_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1)
summary(out.fq_AR_blup)

#refine QTL and plot LOD profile
rqtl_AR_blup <- refineqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1, verbose=FALSE)
rqtl_AR_blup
out.fq2_AR_blup <- fitqtl(mapthis, qtl=qtl_AR_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blup)
plotLodProfile(rqtl_AR_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out41_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_blupNP <- makeqtl(mapthis, chr=c(9), pos=c(122)) 
summary(out.qtl_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, pheno.col=41, formula=y~Q1))
plot(qtl_AR_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1)
summary(out.fq_AR_blupNP)

#refine QTL and plot LOD profile
rqtl_AR_blupNP <- refineqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNP
out.fq2_AR_blupNP <- fitqtl(mapthis, qtl=qtl_AR_blupNP, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNP)
plotLodProfile(rqtl_AR_blupNP)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariate
summary(out42_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_AR_blupNPFC_chr7 <- makeqtl(mapthis, chr=c(7), pos=c(191)) 
qtl_AR_blupNPFC_chr9 <- makeqtl(mapthis, chr=c(9), pos=c(122))
qtl_AR_blupNPFC <- makeqtl(mapthis, chr=c(7, 9), pos=c(191, 122))
summary(out.qtl_AR_blupNPFC_chr7 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr7, pheno.col=42, formula=y~Q1))
summary(out.qtl_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, pheno.col=42, formula=y~Q1))
summary(out.qtl_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, pheno.col=42, formula=y~Q1+Q2))
plot(qtl_AR_blupNPFC_chr7)
plot(qtl_AR_blupNPFC_chr9)
plot(qtl_AR_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_AR_blupNPFC_chr7 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr7, formula=y~Q1)
summary(out.fq_AR_blupNPFC_chr7)
out.fq_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1)
summary(out.fq_AR_blupNPFC_chr9)
out.fq_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2)
summary(out.fq_AR_blupNPFC)

#refine QTL and plot LOD profile
rqtl_AR_blupNPFC_chr7 <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC_chr7, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNPFC_chr7
out.fq2_AR_blupNPFC_chr7 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr7, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNPFC_chr7)
plotLodProfile(rqtl_AR_blupNPFC_chr7)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Average Rep1 + Rep2 (BLUP NPFC))")

rqtl_AR_blupNPFC_chr9 <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1, verbose=FALSE)
rqtl_AR_blupNPFC_chr9
out.fq2_AR_blupNPFC_chr9 <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC_chr9, formula=y~Q1, dropone=FALSE)
summary(out.fq2_AR_blupNPFC_chr9)
plotLodProfile(rqtl_AR_blupNPFC_chr9)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Average Rep1 + Rep2 (BLUP NPFC))")

rqtl_AR_blupNPFC <- refineqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2, verbose=FALSE)
rqtl_AR_blupNPFC
out.fq2_AR_blupNPFC <- fitqtl(mapthis, qtl=qtl_AR_blupNPFC, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_AR_blupNPFC)
plotLodProfile(rqtl_AR_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Aerial Root Nodes (Average Rep1 + Rep2 (BLUP NPFC))")

#obtain lod interval coordinates
lodint(rqtl_AR_rep1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNP, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNPFC_chr7, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNPFC_chr9, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_AR_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_AR_rep1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNP, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNPFC_chr7, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNPFC_chr9, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_AR_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM PDM
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass")
plot(out47_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass")

#plot CIM PDM chrs
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="7", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 7")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="7", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 7")
plot(out47_cim.em, show.marker.names=F,chr="7", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 7")
plot(out43_cim.em, out44_cim.em, out45_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 8")
plot(out46_cim.em, out47_cim.em, out48_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 8")
plot(out47_cim.em, show.marker.names=F,chr="8", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Dry Mass - Chr 8")

summary(out43_cim.em, threshold=3.84)
summary(out44_cim.em, threshold=3.84)
summary(out45_cim.em, threshold=3.84)
summary(out46_cim.em, threshold=3.84)
summary(out47_cim.em, threshold=3.84)
summary(out48_cim.em, threshold=3.84)

#Regression PDM
#Rep 1
summary(out44_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_rep1_chr7 <- makeqtl(mapthis, chr=c(7), pos=c(97.3)) 
qtl_PDM_rep1_chr8 <- makeqtl(mapthis, chr=c(8), pos=c(144.7))
qtl_PDM_rep1 <- makeqtl(mapthis, chr=c(7, 8), pos=c(97.3, 144.7))
summary(out.qtl_PDM_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr7, pheno.col=44, formula=y~Q1))
summary(out.qtl_PDM_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr8, pheno.col=44, formula=y~Q1))
summary(out.qtl_PDM_rep1 <- fitqtl(mapthis, qtl=qtl_PDM_rep1, pheno.col=44, formula=y~Q1+Q2))
plot(qtl_PDM_rep1_chr7)
plot(qtl_PDM_rep1_chr8)
plot(qtl_PDM_rep1)

#fit the model (with QTL in fixed positions)
out.fq_PDM_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr7, formula=y~Q1)
summary(out.fq_PDM_rep1_chr7)
out.fq_PDM_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr8, formula=y~Q1)
summary(out.fq_PDM_rep1_chr8)
out.fq_PDM_rep1 <- fitqtl(mapthis, qtl=qtl_PDM_rep1, formula=y~Q1+Q2)
summary(out.fq_PDM_rep1)

#refine QTL and plot LOD profile
rqtl_PDM_rep1_chr7 <- refineqtl(mapthis, qtl=qtl_PDM_rep1_chr7, formula=y~Q1, verbose=FALSE)
rqtl_PDM_rep1_chr7
out.fq2_PDM_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr7, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_rep1_chr7)
plotLodProfile(rqtl_PDM_rep1_chr7)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Dry Mass (Rep1)")

rqtl_PDM_rep1_chr8 <- refineqtl(mapthis, qtl=qtl_PDM_rep1_chr8, formula=y~Q1, verbose=FALSE)
rqtl_PDM_rep1_chr8
out.fq2_PDM_rep1_chr8 <- fitqtl(mapthis, qtl=qtl_PDM_rep1_chr8, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_rep1_chr8)
plotLodProfile(rqtl_PDM_rep1_chr8)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Dry Mass (Rep1)")

rqtl_PDM_rep1 <- refineqtl(mapthis, qtl=qtl_PDM_rep1, formula=y~Q1+Q2, verbose=FALSE)
rqtl_PDM_rep1
out.fq2_PDM_rep1 <- fitqtl(mapthis, qtl=qtl_PDM_rep1, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_PDM_rep1)
plotLodProfile(rqtl_PDM_rep1)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Dry Mass (Rep1)")

#blup
summary(out46_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_blup <- makeqtl(mapthis, chr=c(4), pos=c(234)) 
summary(out.qtl_PDM_blup <- fitqtl(mapthis, qtl=qtl_PDM_blup, pheno.col=46, formula=y~Q1))
plot(qtl_PDM_blup)

#fit the model (with QTL in fixed positions)
out.fq_PDM_blup <- fitqtl(mapthis, qtl=qtl_PDM_blup, formula=y~Q1)
summary(out.fq_PDM_blup)

#refine QTL and plot LOD profile
rqtl_PDM_blup <- refineqtl(mapthis, qtl=qtl_PDM_blup, formula=y~Q1, verbose=FALSE)
rqtl_PDM_blup
out.fq2_PDM_blup <- fitqtl(mapthis, qtl=qtl_PDM_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_blup)
plotLodProfile(rqtl_PDM_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Dry Mass (Average Rep1 + Rep2, BLUP)")

#blup NPFC
summary(out48_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PDM_blupNPFC <- makeqtl(mapthis, chr=c(4), pos=c(234)) 
summary(out.qtl_PDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC, pheno.col=48, formula=y~Q1))
plot(qtl_PDM_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_PDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC, formula=y~Q1)
summary(out.fq_PDM_blupNPFC)

#refine QTL and plot LOD profile
rqtl_PDM_blupNPFC <- refineqtl(mapthis, qtl=qtl_PDM_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_PDM_blupNPFC
out.fq2_PDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_PDM_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PDM_blupNPFC)
plotLodProfile(rqtl_PDM_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Dry Mass (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_PDM_rep1_chr7, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_rep1_chr8, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_rep1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PDM_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PDM_rep1_chr7, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_rep1_chr8, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_rep1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PDM_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)

#plot CIM PTNP
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage")
plot(out53_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage")

#plot CIM PTNP chrs
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 1")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 1")
plot(out53_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 1")
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 8")
plot(out52_cim.em, out53_cim.em, out54_cim.em, show.marker.names=F,chr="8", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 8")
plot(out53_cim.em, show.marker.names=F,chr="8", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 8")
plot(out49_cim.em, out50_cim.em, out51_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage - Chr 10")

summary(out49_cim.em, threshold=3.84)
summary(out50_cim.em, threshold=3.84)
summary(out51_cim.em, threshold=3.84)
summary(out52_cim.em, threshold=3.84)
summary(out53_cim.em, threshold=3.84)
summary(out54_cim.em, threshold=3.84)

#Regression PTNP
#blup
summary(out52_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTNP_blup <- makeqtl(mapthis, chr=c(1), pos=c(122)) 
summary(out.qtl_PTNP_blup <- fitqtl(mapthis, qtl=qtl_PTNP_blup, pheno.col=52, formula=y~Q1))
plot(qtl_PTNP_blup)

#fit the model (with QTL in fixed positions)
out.fq_PTNP_blup <- fitqtl(mapthis, qtl=qtl_PTNP_blup, formula=y~Q1)
summary(out.fq_PTNP_blup)

#refine QTL and plot LOD profile
rqtl_PTNP_blup <- refineqtl(mapthis, qtl=qtl_PTNP_blup, formula=y~Q1, verbose=FALSE)
rqtl_PTNP_blup
out.fq2_PTNP_blup <- fitqtl(mapthis, qtl=qtl_PTNP_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTNP_blup)
plotLodProfile(rqtl_PTNP_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage (Average Rep1 + Rep2, BLUP)")

#blup no parents
summary(out53_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTNP_blupNP <- makeqtl(mapthis, chr=c(1), pos=c(122)) 
summary(out.qtl_PTNP_blupNP <- fitqtl(mapthis, qtl=qtl_PTNP_blupNP, pheno.col=53, formula=y~Q1))
plot(qtl_PTNP_blupNP)

#fit the model (with QTL in fixed positions)
out.fq_PTNP_blupNP <- fitqtl(mapthis, qtl=qtl_PTNP_blupNP, formula=y~Q1)
summary(out.fq_PTNP_blupNP)

#refine QTL and plot LOD profile
rqtl_PTNP_blupNP <- refineqtl(mapthis, qtl=qtl_PTNP_blupNP, formula=y~Q1, verbose=FALSE)
rqtl_PTNP_blupNP
out.fq2_PTNP_blupNP <- fitqtl(mapthis, qtl=qtl_PTNP_blupNP, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTNP_blupNP)
plotLodProfile(rqtl_PTNP_blupNP)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage (Average Rep1 + Rep2, BLUP NP)")

#blup no parents and family as covariate
summary(out54_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTNP_blupNPFC <- makeqtl(mapthis, chr=c(1), pos=c(122)) 
summary(out.qtl_PTNP_blupNPFC <- fitqtl(mapthis, qtl=qtl_PTNP_blupNPFC, pheno.col=54, formula=y~Q1))
plot(qtl_PTNP_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_PTNP_blupNPFC <- fitqtl(mapthis, qtl=qtl_PTNP_blupNPFC, formula=y~Q1)
summary(out.fq_PTNP_blupNPFC)

#refine QTL and plot LOD profile
rqtl_PTNP_blupNPFC <- refineqtl(mapthis, qtl=qtl_PTNP_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_PTNP_blupNPFC
out.fq2_PTNP_blupNPFC <- fitqtl(mapthis, qtl=qtl_PTNP_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTNP_blupNPFC)
plotLodProfile(rqtl_PTNP_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen Percentage (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_PTNP_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PTNP_blupNP, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PTNP_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PTNP_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PTNP_blupNP, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PTNP_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)


#plot CIM PTN
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen")
plot(out59_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen")

#plot CIM PTN chrs
plot(out55_cim.em, out56_cim.em, out57_cim.em, show.marker.names=F,chr="4", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen - Chr 4")
plot(out58_cim.em, out59_cim.em, out60_cim.em, show.marker.names=F,chr="4", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen - Chr 4")
plot(out59_cim.em, show.marker.names=F,chr="4", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen - Chr 4")

summary(out55_cim.em, threshold=3.84)
summary(out56_cim.em, threshold=3.84)
summary(out57_cim.em, threshold=3.84)
summary(out58_cim.em, threshold=3.84)
summary(out59_cim.em, threshold=3.84)
summary(out60_cim.em, threshold=3.84)

#Regression PTN
#Rep 1
summary(out56_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTN_rep1_chr7 <- makeqtl(mapthis, chr=c(7), pos=c(111)) 
qtl_PTN_rep1_chr10 <- makeqtl(mapthis, chr=c(10), pos=c(106))
qtl_PTN_rep1 <- makeqtl(mapthis, chr=c(7, 10), pos=c(111, 106))
summary(out.qtl_PTN_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr7, pheno.col=56, formula=y~Q1))
summary(out.qtl_PTN_rep1_chr10 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr10, pheno.col=56, formula=y~Q1))
summary(out.qtl_PTN_rep1 <- fitqtl(mapthis, qtl=qtl_PTN_rep1, pheno.col=56, formula=y~Q1+Q2))
plot(qtl_PTN_rep1_chr7)
plot(qtl_PTN_rep1_chr10)
plot(qtl_PTN_rep1)

#fit the model (with QTL in fixed positions)
out.fq_PTN_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr7, formula=y~Q1)
summary(out.fq_PTN_rep1_chr7)
out.fq_PTN_rep1_chr10 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr10, formula=y~Q1)
summary(out.fq_PTN_rep1_chr10)
out.fq_PTN_rep1 <- fitqtl(mapthis, qtl=qtl_PTN_rep1, formula=y~Q1+Q2)
summary(out.fq_PTN_rep1)

#refine QTL and plot LOD profile
rqtl_PTN_rep1_chr7 <- refineqtl(mapthis, qtl=qtl_PTN_rep1_chr7, formula=y~Q1, verbose=FALSE)
rqtl_PTN_rep1_chr7
out.fq2_PTN_rep1_chr7 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr7, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTN_rep1_chr7)
plotLodProfile(rqtl_PTN_rep1_chr7)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen (Rep1)")

rqtl_PTN_rep1_chr10 <- refineqtl(mapthis, qtl=qtl_PTN_rep1_chr10, formula=y~Q1, verbose=FALSE)
rqtl_PTN_rep1_chr10
out.fq2_PTN_rep1_chr10 <- fitqtl(mapthis, qtl=qtl_PTN_rep1_chr10, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTN_rep1_chr10)
plotLodProfile(rqtl_PTN_rep1_chr10)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen (Rep1)")

rqtl_PTN_rep1 <- refineqtl(mapthis, qtl=qtl_PTN_rep1, formula=y~Q1+Q2, verbose=FALSE)
rqtl_PTN_rep1
out.fq2_PTN_rep1 <- fitqtl(mapthis, qtl=qtl_PTN_rep1, formula=y~Q1+Q2, dropone=FALSE)
summary(out.fq2_PTN_rep1)
plotLodProfile(rqtl_PTN_rep1)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen (Rep1)")

#blup
summary(out58_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_PTN_blup <- makeqtl(mapthis, chr=c(3), pos=c(79.7)) 
summary(out.qtl_PTN_blup <- fitqtl(mapthis, qtl=qtl_PTN_blup, pheno.col=58, formula=y~Q1))
plot(qtl_PTN_blup)

#fit the model (with QTL in fixed positions)
out.fq_PTN_blup <- fitqtl(mapthis, qtl=qtl_PTN_blup, formula=y~Q1)
summary(out.fq_PTN_blup)

#refine QTL and plot LOD profile
rqtl_PTN_blup <- refineqtl(mapthis, qtl=qtl_PTN_blup, formula=y~Q1, verbose=FALSE)
rqtl_PTN_blup
out.fq2_PTN_blup <- fitqtl(mapthis, qtl=qtl_PTN_blup, formula=y~Q1, dropone=FALSE)
summary(out.fq2_PTN_blup)
plotLodProfile(rqtl_PTN_blup)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Plant Total Nitrogen (Average Rep1 + Rep2, BLUP)")

#obtain lod interval coordinates
lodint(rqtl_PTN_rep1_chr7, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PTN_rep1_chr10, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PTN_rep1, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_PTN_blup, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_PTN_rep1_chr7, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PTN_rep1_chr10, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PTN_rep1, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_PTN_blup, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)


#plot CIM GDM
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass")
plot(out65_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass")

#plot CIM GDM chrs
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 1")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 1")
plot(out65_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 1")
plot(out61_cim.em, out62_cim.em, out63_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 10")
plot(out64_cim.em, out65_cim.em, out66_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 10")
plot(out65_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Dry Mass - Chr 10")

summary(out61_cim.em, threshold=3.84)
summary(out62_cim.em, threshold=3.84)
summary(out63_cim.em, threshold=3.84)
summary(out64_cim.em, threshold=3.84)
summary(out65_cim.em, threshold=3.84)
summary(out66_cim.em, threshold=3.84)

#Regression GDM
#Average
summary(out61_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_GDM_avg <- makeqtl(mapthis, chr=c(1), pos=c(26.3)) 
summary(out.qtl_GDM_avg <- fitqtl(mapthis, qtl=qtl_GDM_avg, pheno.col=61, formula=y~Q1))
plot(qtl_GDM_avg)

#fit the model (with QTL in fixed positions)
out.fq_GDM_avg <- fitqtl(mapthis, qtl=qtl_GDM_avg, formula=y~Q1)
summary(out.fq_GDM_avg)

#refine QTL and plot LOD profile
rqtl_GDM_avg <- refineqtl(mapthis, qtl=qtl_GDM_avg, formula=y~Q1, verbose=FALSE)
rqtl_GDM_avg
out.fq2_GDM_avg <- fitqtl(mapthis, qtl=qtl_GDM_avg, formula=y~Q1, dropone=FALSE)
summary(out.fq2_GDM_avg)
plotLodProfile(rqtl_GDM_avg)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Grain Dry Mass (Average Rep1 + Rep2)")

#Rep 2
summary(out63_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_GDM_rep2 <- makeqtl(mapthis, chr=c(1), pos=c(26.3)) 
summary(out.qtl_GDM_rep2 <- fitqtl(mapthis, qtl=qtl_GDM_rep2, pheno.col=63, formula=y~Q1))
plot(qtl_GDM_rep2)

#fit the model (with QTL in fixed positions)
out.fq_GDM_rep2 <- fitqtl(mapthis, qtl=qtl_GDM_rep2, formula=y~Q1)
summary(out.fq_GDM_rep2)

#refine QTL and plot LOD profile
rqtl_GDM_rep2 <- refineqtl(mapthis, qtl=qtl_GDM_rep2, formula=y~Q1, verbose=FALSE)
rqtl_GDM_rep2
out.fq2_GDM_rep2 <- fitqtl(mapthis, qtl=qtl_GDM_rep2, formula=y~Q1, dropone=FALSE)
summary(out.fq2_GDM_rep2)
plotLodProfile(rqtl_GDM_rep2)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Grain Dry Mass (Rep2)")

#blup no parents and family as covariate
summary(out66_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_GDM_blupNPFC <- makeqtl(mapthis, chr=c(1), pos=c(26.3)) 
summary(out.qtl_GDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_GDM_blupNPFC, pheno.col=66, formula=y~Q1))
plot(qtl_GDM_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_GDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_GDM_blupNPFC, formula=y~Q1)
summary(out.fq_GDM_blupNPFC)

#refine QTL and plot LOD profile
rqtl_GDM_blupNPFC <- refineqtl(mapthis, qtl=qtl_GDM_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_GDM_blupNPFC
out.fq2_GDM_blupNPFC <- fitqtl(mapthis, qtl=qtl_GDM_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_GDM_blupNPFC)
plotLodProfile(rqtl_GDM_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Grain Dry Mass (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_GDM_avg, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_GDM_rep2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_GDM_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_GDM_avg, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_GDM_rep2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_GDM_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)


#plot CIM GTNP
plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen Percentage")
plot(out70_cim.em, out71_cim.em, out72_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen Percentage")
plot(out71_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen Percentage")

#plot CIM GTNP chrs
#plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="2", col =c("Black", "blue", "red"))
#abline(3.84, 0, untf = FALSE)
#plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="3", col =c("Black", "blue", "red"))
#abline(3.84, 0, untf = FALSE)
#plot(out67_cim.em, out68_cim.em, out69_cim.em, show.marker.names=F,chr="8", col =c("Black", "blue", "red"))
#abline(3.84, 0, untf = FALSE)

summary(out67_cim.em, threshold=3.84)
summary(out68_cim.em, threshold=3.84)
summary(out69_cim.em, threshold=3.84)
summary(out70_cim.em, threshold=3.84)
summary(out71_cim.em, threshold=3.84)
summary(out72_cim.em, threshold=3.84)

#plot CIM GTN
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen")
plot(out77_cim.em, show.marker.names=F,col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen")

#plot CIM GTN chrs
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,chr="1", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 1")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,chr="1", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 1")
plot(out77_cim.em, show.marker.names=F,chr="1", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 1")
plot(out73_cim.em, out74_cim.em, out75_cim.em, show.marker.names=F,chr="10", col =c("Black", "blue", "red"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 10")
plot(out76_cim.em, out77_cim.em, out78_cim.em, show.marker.names=F,chr="10", col =c("brown", "purple", "green"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 10")
plot(out77_cim.em, show.marker.names=F,chr="10", col =c("purple"))
abline(3.84, 0, untf = FALSE)
abline(2.34, 0, untf = FALSE,col =c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen - Chr 10")

summary(out73_cim.em, threshold=3.84)
summary(out74_cim.em, threshold=3.84)
summary(out75_cim.em, threshold=3.84)
summary(out76_cim.em, threshold=3.84)
summary(out77_cim.em, threshold=3.84)
summary(out78_cim.em, threshold=3.84)

#Regression GTN
#Rep 2
summary(out75_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_GTN_rep2 <- makeqtl(mapthis, chr=c(1), pos=c(42.5)) 
summary(out.qtl_GTN_rep2 <- fitqtl(mapthis, qtl=qtl_GTN_rep2, pheno.col=75, formula=y~Q1))
plot(qtl_GTN_rep2)

#fit the model (with QTL in fixed positions)
out.fq_GTN_rep2 <- fitqtl(mapthis, qtl=qtl_GTN_rep2, formula=y~Q1)
summary(out.fq_GTN_rep2)

#refine QTL and plot LOD profile
rqtl_GTN_rep2 <- refineqtl(mapthis, qtl=qtl_GTN_rep2, formula=y~Q1, verbose=FALSE)
rqtl_GTN_rep2
out.fq2_GTN_rep2 <- fitqtl(mapthis, qtl=qtl_GTN_rep2, formula=y~Q1, dropone=FALSE)
summary(out.fq2_GTN_rep2)
plotLodProfile(rqtl_GTN_rep2)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen (Rep2)")

#blup no parents and family as covariate
summary(out78_cim.em, threshold=3.84)
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_GTN_blupNPFC <- makeqtl(mapthis, chr=c(1), pos=c(26.3)) 
summary(out.qtl_GTN_blupNPFC <- fitqtl(mapthis, qtl=qtl_GTN_blupNPFC, pheno.col=78, formula=y~Q1))
plot(qtl_GTN_blupNPFC)

#fit the model (with QTL in fixed positions)
out.fq_GTN_blupNPFC <- fitqtl(mapthis, qtl=qtl_GTN_blupNPFC, formula=y~Q1)
summary(out.fq_GTN_blupNPFC)

#refine QTL and plot LOD profile
rqtl_GTN_blupNPFC <- refineqtl(mapthis, qtl=qtl_GTN_blupNPFC, formula=y~Q1, verbose=FALSE)
rqtl_GTN_blupNPFC
out.fq2_GTN_blupNPFC <- fitqtl(mapthis, qtl=qtl_GTN_blupNPFC, formula=y~Q1, dropone=FALSE)
summary(out.fq2_GTN_blupNPFC)
plotLodProfile(rqtl_GTN_blupNPFC)
abline(3.84, 0, untf = FALSE)
#abline(2.34, 0, untf = FALSE, col=c("dark gray"))
title(main = "MO17 - Grain Total Nitrogen (Average Rep1 + Rep2, BLUP NPFC)")

#obtain lod interval coordinates
lodint(rqtl_GTN_rep2, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)
lodint(rqtl_GTN_blupNPFC, chr, qtl.index=1, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)

#obtain Bayesian credible interval
bayesint(rqtl_GTN_rep2, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)
bayesint(rqtl_GTN_blupNPFC, chr, qtl.index=1, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)


#plot CIM blup_noparents SC
plot(out5_cim.em, show.marker.names=F,col =c("purple")) #SC
abline(3.84, 0, untf = FALSE)
summary(out5_cim.em, threshold=3.84)
plot(out5_cim.em, show.marker.names=F,chr="9", col =c("purple")) #SC
abline(3.84, 0, untf = FALSE)

#Regression SC
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_SC <- makeqtl(mapthis, chr=c(9), pos=c(111)) 
summary(out.qtl_SC <- fitqtl(mapthis, qtl=qtl_SC, pheno.col=5, formula=y~Q1))

#plot CIM blup_noparents DT
plot(out11_cim.em, show.marker.names=F,col =c("purple")) #DT
abline(3.84, 0, untf = FALSE)
summary(out11_cim.em, threshold=3.84)
plot(out11_cim.em, show.marker.names=F,chr="3", col =c("purple")) #DT
abline(3.84, 0, untf = FALSE)
plot(out11_cim.em, show.marker.names=F,chr="8", col =c("purple")) #DT
abline(3.84, 0, untf = FALSE)
plot(out11_cim.em, show.marker.names=F,chr="10", col =c("purple")) #DT
abline(3.84, 0, untf = FALSE)

#Regression DT
mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_DT <- makeqtl(mapthis, chr=c(3), pos=c(156))
qtl_DT <- makeqtl(mapthis, chr=c(8), pos=c(156)) 
#qtl_DT <- makeqtl(mapthis, chr=c(10), pos=c(136)) 
summary(out.qtl_DT <- fitqtl(mapthis, qtl=qtl_DT, pheno.col=11, formula=y~Q1))

#plot CIM blup_noparents DPS
plot(out17_cim.em, show.marker.names=F,col =c("purple")) #DPS
abline(3.84, 0, untf = FALSE)
summary(out17_cim.em, threshold=3.84)
plot(out17_cim.em, show.marker.names=F,chr="3", col =c("purple")) #DPS
abline(3.84, 0, untf = FALSE)
plot(out17_cim.em, show.marker.names=F,chr="8", col =c("purple")) #DPS
abline(3.84, 0, untf = FALSE)
plot(out17_cim.em, show.marker.names=F,chr="10", col =c("purple")) #DPS
abline(3.84, 0, untf = FALSE)

#Regression DPS
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_DPS <- makeqtl(mapthis, chr=c(3), pos=c(121)) 
qtl_DPS <- makeqtl(mapthis, chr=c(8), pos=c(156))
#qtl_DPS <- makeqtl(mapthis, chr=c(10), pos=c(136))
summary(out.qtl_DPS <- fitqtl(mapthis, qtl=qtl_DPS, pheno.col=17, formula=y~Q1))

#plot CIM blup_noparents T1
plot(out23_cim.em, show.marker.names=F,col =c("purple")) #15NT1
abline(3.84, 0, untf = FALSE)
summary(out23_cim.em, threshold=3.84)
plot(out23_cim.em, show.marker.names=F,chr="9", col =c("purple")) #15NT1
abline(3.84, 0, untf = FALSE)

#Regression T1
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_T1 <- makeqtl(mapthis, chr=c(7), pos=c(208))
#summary(out.qtl_T1 <- fitqtl(mapthis, qtl=qtl_T1, pheno.col=23, formula=y~Q1))

#plot CIM blup_noparents T2
plot(out29_cim.em, show.marker.names=F,col =c("purple")) #15NT2
abline(3.84, 0, untf = FALSE)
summary(out29_cim.em, threshold=3.84)
plot(out29_cim.em, show.marker.names=F,chr="9", col =c("purple")) #15NT2
abline(3.84, 0, untf = FALSE)

#Regression T2
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_T2 <- makeqtl(mapthis, chr=c(7), pos=c(208)) 
#summary(out.qtl_T2 <- fitqtl(mapthis, qtl=qtl_T2, pheno.col=29, formula=y~Q1))

#plot CIM blup_noparents T3
plot(out35_cim.em, show.marker.names=F,col =c("purple")) #15NT3
abline(3.84, 0, untf = FALSE)
summary(out35_cim.em, threshold=3.84)
plot(out35_cim.em, show.marker.names=F,chr="2", col =c("purple")) #15NT3
abline(3.84, 0, untf = FALSE)

#Regression T3
mapthis <- sim.geno(mapthis, n.draws=128)
qtl_T3 <- makeqtl(mapthis, chr=c(2), pos=c(150))
summary(out.qtl_T3 <- fitqtl(mapthis, qtl=qtl_T3, pheno.col=35, formula=y~Q1))

#plot CIM blup_noparents AR
plot(out41_cim.em, show.marker.names=F,col =c("purple")) #AR
abline(3.84, 0, untf = FALSE)
summary(out41_cim.em, threshold=3.84)
plot(out41_cim.em, show.marker.names=F,chr="6", col =c("purple")) #AR
abline(3.84, 0, untf = FALSE)
plot(out41_cim.em, show.marker.names=F,chr="9", col =c("purple")) #AR
abline(3.84, 0, untf = FALSE)

#Regression AR
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_AR <- makeqtl(mapthis, chr=c(6), pos=c(95.3))
#qtl_AR <- makeqtl(mapthis, chr=c(9), pos=c(117.6))
#summary(out.qtl_AR <- fitqtl(mapthis, qtl=qtl_AR, pheno.col=41, formula=y~Q1))

#plot CIM blup_noparents PDM
plot(out47_cim.em, show.marker.names=F,col =c("purple")) #PDM
abline(3.84, 0, untf = FALSE)
summary(out48_cim.em, threshold=3.84)
plot(out47_cim.em, show.marker.names=F,chr="7", col =c("purple")) #PDM
abline(3.84, 0, untf = FALSE)
plot(out47_cim.em, show.marker.names=F,chr="8", col =c("purple")) #PDM
abline(3.84, 0, untf = FALSE)

#Regression PDM
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_PDM <- makeqtl(mapthis, chr=c(7), pos=c(144.7)) 
#qtl_PDM <- makeqtl(mapthis, chr=c(8), pos=c(91.3))
#summary(out.qtl_PDM <- fitqtl(mapthis, qtl=qtl_PDM, pheno.col=47, formula=y~Q1))

#plot CIM blup_noparents PTNP
plot(out55_cim.em, show.marker.names=F,col =c("purple")) #PTNP
abline(3.84, 0, untf = FALSE)
summary(out55_cim.em, threshold=3.84)
plot(out55_cim.em, show.marker.names=F,chr="7", col =c("purple")) #PTNP
abline(3.84, 0, untf = FALSE)

#Regression PTNP
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_PTNP <- makeqtl(mapthis, chr=c(7), pos=c(170)) 
#summary(out.qtl_PTNP <- fitqtl(mapthis, qtl=qtl_PTNP, pheno.col=55, formula=y~Q1))

#plot CIM blup_noparents PTN
plot(out59_cim.em, show.marker.names=F,col =c("purple")) #PTN
abline(3.84, 0, untf = FALSE)
summary(out59_cim.em, threshold=3.84)

#Regression PTN
# There were no LOD peaks above the threshold.

#plot CIM blup_noparents GDM
plot(out65_cim.em, show.marker.names=F,col =c("purple")) #GDM
abline(3.84, 0, untf = FALSE)
summary(out65_cim.em, threshold=3.84)
plot(out65_cim.em, show.marker.names=F,chr="1", col =c("purple")) #GDM
abline(3.84, 0, untf = FALSE)

#Regression GDM
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_GDM <- makeqtl(mapthis, chr=c(1), pos=c(266)) 
#summary(out.qtl_GDM <- fitqtl(mapthis, qtl=qtl_GDM, pheno.col=65, formula=y~Q1))

#plot CIM blup_noparents GTNP
plot(out71_cim.em, show.marker.names=F,col =c("purple")) #GTNP
abline(3.84, 0, untf = FALSE)
summary(out71_cim.em, threshold=3.84)

#Regression GTNP
# There were no LOD peaks above the threshold.

#plot CIM blup_noparents GTN
plot(out77_cim.em, show.marker.names=F,col =c("purple")) #GTN
abline(3.84, 0, untf = FALSE)
summary(out77_cim.em, threshold=3.84)
plot(out77_cim.em, show.marker.names=F,chr="1", col =c("purple")) #GTN
abline(3.84, 0, untf = FALSE)

#Regression GTN
#mapthis <- sim.geno(mapthis, n.draws=128)
#qtl_GTN <- makeqtl(mapthis, chr=c(1), pos=c(266)) 
#summary(out.qtl_GTN <- fitqtl(mapthis, qtl=qtl_GTN, pheno.col=77, formula=y~Q1))

#lodint and bayesint for non-refined QTL SC
lodint(out2_cim.em, chr=6, expandtomarkers=TRUE)
bayesint(out2_cim.em, chr=6, expandtomarkers=TRUE)
lodint(out6_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out6_cim.em, chr=9, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL DT
lodint(out7_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out7_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out8_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out8_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out9_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out9_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out9_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out9_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out10_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out10_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out11_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out11_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out12_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out12_cim.em, chr=8, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL DPS
lodint(out13_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out13_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out13_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out13_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out14_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out14_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out14_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out14_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out15_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out15_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out15_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out15_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out16_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out16_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out16_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out16_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out17_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out17_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out17_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out17_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out18_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out18_cim.em, chr=3, expandtomarkers=TRUE)
lodint(out18_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out18_cim.em, chr=8, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL 15NT3
lodint(out34_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out34_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out35_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out35_cim.em, chr=2, expandtomarkers=TRUE)
lodint(out36_cim.em, chr=2, expandtomarkers=TRUE)
bayesint(out36_cim.em, chr=2, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL AR
lodint(out38_cim.em, chr=7, expandtomarkers=TRUE)
bayesint(out38_cim.em, chr=7, expandtomarkers=TRUE)
lodint(out40_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out40_cim.em, chr=9, expandtomarkers=TRUE)
lodint(out41_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out41_cim.em, chr=9, expandtomarkers=TRUE)
lodint(out42_cim.em, chr=7, expandtomarkers=TRUE)
bayesint(out42_cim.em, chr=7, expandtomarkers=TRUE)
lodint(out42_cim.em, chr=9, expandtomarkers=TRUE)
bayesint(out42_cim.em, chr=9, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PDM
lodint(out44_cim.em, chr=7, expandtomarkers=TRUE)
bayesint(out44_cim.em, chr=7, expandtomarkers=TRUE)
lodint(out44_cim.em, chr=8, expandtomarkers=TRUE)
bayesint(out44_cim.em, chr=8, expandtomarkers=TRUE)
lodint(out46_cim.em, chr=4, expandtomarkers=TRUE)
bayesint(out46_cim.em, chr=4, expandtomarkers=TRUE)
lodint(out48_cim.em, chr=4, expandtomarkers=TRUE)
bayesint(out48_cim.em, chr=4, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PTNP
lodint(out52_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out52_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out53_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out53_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out54_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out54_cim.em, chr=1, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL PTN
lodint(out56_cim.em, chr=7, expandtomarkers=TRUE)
bayesint(out56_cim.em, chr=7, expandtomarkers=TRUE)
lodint(out56_cim.em, chr=10, expandtomarkers=TRUE)
bayesint(out56_cim.em, chr=10, expandtomarkers=TRUE)
lodint(out58_cim.em, chr=3, expandtomarkers=TRUE)
bayesint(out58_cim.em, chr=3, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL GDM
lodint(out61_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out61_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out63_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out63_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out66_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out66_cim.em, chr=1, expandtomarkers=TRUE)

#lodint and bayesint for non-refined QTL GTN
lodint(out75_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out75_cim.em, chr=1, expandtomarkers=TRUE)
lodint(out78_cim.em, chr=1, expandtomarkers=TRUE)
bayesint(out78_cim.em, chr=1, expandtomarkers=TRUE)