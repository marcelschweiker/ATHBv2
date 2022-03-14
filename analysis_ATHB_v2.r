
library(MASS)
library(comf)
library(ggplot2)
library(vtable)
library(effectsize)
library(ggpubr)
library(moments)
library(sjPlot)
library(psych)

################################
# load function calcPMVL from separate file
###################################


pathFiles <- "xx" # add here path to file ashrae_db2.01.csv

# ashrae_db2.01.csv is downloaded from https://www.kaggle.com/claytonmiller/ashrae-global-thermal-comfort-database-ii
dfDBII <- read.csv(paste(pathFiles, "ashrae_db2.01.csv", sep=""))

colnames(dfDBII) <- c("publication","dataContributor","Year","Season","koppenClimateClassification","Climate","City","Country","BuildingType","CoolingStrategyBuilding","CoolingStrategyOperation","HeatingStrategyBuilding","Age","Sex","ASV","Acc","Pref","AccAV","PrefAV","Comf","PMV","PPD","SET","Clo","Met","activity_10","activity_20","activity_30","activity_60","Ta","TaF","Ta_h","Ta_hF","Ta_m","Ta_mF","Ta_l","Ta_lF","Top","TopF","Tr","TrF","Tg","TgF","Tg_h","Tg_hF","Tg_m","Tg_mF","Tg_l","Tg_lF","rH","prefHum","sensHum","av","avUS", "Velocity_h","Velocity_hfpm","Velocity_m","Velocity_mfpm","Velocity_l","Velocityfpm","height","weight","Blind","Fan","Window","Door","Heater","Trm","TrmF","Database")

dfDBII$Tg <- ifelse(!is.na(dfDBII$Tg), dfDBII$Tg,
						ifelse(!is.na(dfDBII$Tg_h), dfDBII$Tg_h,
						NA))
  
dfDBIIsel <- dfDBII[,c(1,3:5,7:17,20:25,30,38,40,42,44, 50,53,61,62,68)]
dfDBIIsel$publication <- substr(dfDBIIsel$publication,1,33)


# from here transform tglobe and top to tr, to be used when tr is missing
calcTradGlobe <- function(tg, tair, va, diameter = 0.15, eg = 0.95){
	hcgn <- 1.4*(abs(tg-tair)/diameter)^0.25
	hcgf <- 6.3*((va^0.6)/(diameter^0.4))
	if(hcgn > hcgf){ # natural convection
		trcalc <- ((tg+273)^4 + (0.4*(10^8))*(abs(tg-tair)^0.25)*(tg-tair))^0.25 - 273
	} else { # forced convection
		trcalc <- ((tg+273)^4 + (2.5*(10^8))*(va^0.6)*(tg-tair))^0.25 - 273
	}
	trcalc
}

vTg <- dfDBIIsel$Tg; vTa <- dfDBIIsel$Ta; vav <- dfDBIIsel$av
trcalcTg <- NA
for(i in 1:nrow(dfDBIIsel)){
	if(!is.na(vTg[i]) & !is.na(vTa[i]) & !is.na(vav[i])){
		trcalcTg[i] <- calcTradGlobe(vTg[i], vTa[i], vav[i])
	}
}

calcTrFromTopTa <- function(ta, to, vel, met){
	varIn <- vel + 0.0052 * (met - 58)
	vara <- ifelse(varIn < 0.2, 0.5, ifelse(varIn > 0.6, 0.7, 
        0.6))
	(tr <- (to - (vara * ta)) / (1 - vara) )
}

vTop <- dfDBIIsel$Top; vMet <- dfDBIIsel$Met
trcalcTop <- NA
for(i in 1:nrow(dfDBIIsel)){
	if(!is.na(vTa[i]) & !is.na(vTop[i]) & !is.na(vav[i]) & !is.na(vMet[i])){
		trcalcTop[i] <- calcTrFromTopTa(vTa[i], vTop[i], vav[i], vMet[i])
	}
}
#summary(trcalcTop)

vTr <- dfDBIIsel$Tr
tradjoint <- NA
for(i in 1:nrow(dfDBIIsel)){
	tradjoint[i] <- ifelse(!is.na(vTr[i]), vTr[i],
						ifelse(!is.na(trcalcTg[i]), trcalcTg[i],
						ifelse(!is.na(trcalcTop[i]), trcalcTop[i],
						ifelse(!is.na(vTa[i]), vTa[i], 
						NA))))
	}
	
	
dfDBIIsel$Trad <- tradjoint
						
rm(vTg, vTa, vav, vTop, vMet, vTr)

dfDBIIselval <- dfDBIIsel[which(!is.na(dfDBIIsel$ASV)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Ta)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Trad)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$av)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$rH)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Clo)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Met)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Trm)),]
dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Season)),] #to keep sample size the same in latter analysis
nrow(dfDBIIselval)


# inspect bins and manually choose limits based on bins <49
# set up cut-off values 
breaks <- c(0.6,8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4)
# specify interval/bin labels
tags <- c("[0.6-0.8)","[0.8-1)", "[1-1.2)", "[1.2-1.4)", "[1.4-1.6)", "[1.6-1.8)","[1.8-2)", "[2-2.2)","[2.2-2.4)", "[2.4-2.6)", 
							"[2.6-2.8)","[2.8-3)", "[3-3.2)","[3.2-3.4)", "[3.4-3.6)","[3.6-3.8)", "[3.8-4)")
# bucketing values into bins
group_tags <- cut(dfDBIIselval$Met, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)

summary(group_tags)
dfDBIIselval$MetBin <- group_tags
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Met <= 2.6),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

Clo_tags <- cut(dfDBIIselval$Clo, breaks = 28)
summary(Clo_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Clo <= 2.0),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Clo >= 0.1),] # removing due to low number in bins < 49

nrow(dfDBIIselval)

Ta_tags <- cut(dfDBIIselval$Ta, breaks = 28)
summary(Ta_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Ta <= 38.5),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Ta >= 12.6),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

rH_tags <- cut(dfDBIIselval$rH, breaks = 75)
summary(rH_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$rH <= 87.7),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$rH >= 16.9),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$av <= 8.0),] # removing due to unrealistic values
av_tags <- cut(dfDBIIselval$av, breaks = 28)
summary(av_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$av <= 1.9),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Trad <= 48.0),] # removing due to unrealistic values
Trad_tags <- cut(dfDBIIselval$Trad, breaks = 28)
summary(Trad_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Trad <= 38.5),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Trad >= 12.6),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

Trad_tags <- cut((dfDBIIselval$Trad - dfDBIIselval$Ta), breaks = 28)
summary(Trad_tags)
dfDBIIselval <- dfDBIIselval[which((dfDBIIselval$Trad - dfDBIIselval$Ta) >= (-7.4)),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which((dfDBIIselval$Trad - dfDBIIselval$Ta) <= 9.2),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

Trm_tags <- cut(dfDBIIselval$Trm, breaks = 28)
summary(Trm_tags)
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Trm <= 41.3),] # removing due to low number in bins < 49
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$Trm >= (-2.7)),] # removing due to low number in bins < 49
nrow(dfDBIIselval)

dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$CoolingStrategyBuilding != "Mechanically Ventilated"),] 
dfDBIIselval <- dfDBIIselval[which(dfDBIIselval$BuildingType != "Senior center"),] #too low number for later analysis
nrow(dfDBIIselval)

dfDBIIselval <- dfDBIIselval[which(!is.na(dfDBIIselval$Season)),] #too low number for later analysis
nrow(dfDBIIselval)
# nrow with age/sex
nrow(dfDBIIselval[which(!is.na(dfDBIIselval$Age)&!is.na(dfDBIIselval$Sex)),])

summary(dfDBIIselval)

# calculate dynamic clothing
dfDBIIselval$CloO <- dfDBIIselval$Clo

ITRnude <- exp(-0.533*(dfDBIIselval$av-0.15)+0.069*(dfDBIIselval$av-0.15)^2)*0.7
corrITclothed <- exp(-0.281*(dfDBIIselval$av-0.15)+0.044*(dfDBIIselval$av-0.15)^2)
corrITclothed <- ifelse(corrITclothed>1, 1, corrITclothed)


dfDBIIselval$Clo <- ifelse(dfDBIIselval$Clo > 0.6,
				corrITclothed*dfDBIIselval$Clo,
				(((0.6-dfDBIIselval$Clo)*ITRnude+(dfDBIIselval$Clo*dfDBIIselval$Clo*corrITclothed))/0.6))
					
i <- 1
dfNew <- data.frame(pmv = 0, ppd = 0, Lraw = 0)

maxLength <- nrow(dfDBIIselval)
dfNew <- sapply(seq(maxLength), function(i) {calcPMVL(dfDBIIselval$Ta[i], dfDBIIselval$Trad[i], dfDBIIselval$av[i], dfDBIIselval$rH[i], dfDBIIselval$Clo[i], dfDBIIselval$Met[i]) } )

dfNew <- as.data.frame(t(as.data.frame(dfNew))) # change from list to data frame and transpose rows and columns
dfNew$Lraw <- as.numeric(as.character(dfNew$Lraw))
dfNew$pmv <- as.numeric(as.character(dfNew$pmv))
dfNew$ppd <- as.numeric(as.character(dfNew$ppd))
dfNewAll <- cbind(dfDBIIselval, dfNew)
#dfNewAll <- dfNewAll[-which(dfNewAll$Lraw>200),]
						
i <- 1
dfNewSET <- data.frame(pmv = 0, ppd = 0, Lraw = 0)

maxLength <- nrow(dfDBIIselval)
dfNewSET <- sapply(seq(maxLength), function(i) {calc2Node(dfDBIIselval$Ta[i], dfDBIIselval$Trad[i], dfDBIIselval$av[i], dfDBIIselval$rH[i], dfDBIIselval$Clo[i], dfDBIIselval$Met[i]) } )

dfNewSET <- as.data.frame(t(as.data.frame(dfNewSET))) # change from list to data frame and transpose rows and columns
dfNewAll$set <- as.numeric(as.character(dfNewSET$set))
dfNewAll$pts <- as.numeric(as.character(dfNewSET$pts))

dfNewAll$pmv <- ifelse(dfNewAll$pmv > 3, 3, ifelse(dfNewAll$pmv  < (-3), (-3), dfNewAll$pmv ))
dfNewAll$pts <- ifelse(dfNewAll$pts > 3, 3, ifelse(dfNewAll$pts  < (-3), (-3), dfNewAll$pts ))


nrow(dfNewAll)
summary(dfNewAll)

#https://cran.r-project.org/web/packages/vtable/vignettes/sumtable.html
st(dfNewAll[,c(3,4,7,11,12,28,29,13,20,21,22,31,26:27,30, 32:37)])

#split into training, validation and test data
cTrain <- c(seq(1, nrow(dfNewAll), 5), seq(3, nrow(dfNewAll), 5), seq(5, nrow(dfNewAll), 5))
cVal <- seq(2, nrow(dfNewAll), 5)
cTest <- seq(4, nrow(dfNewAll), 5)
#length(unique(c(cTrain, cVal, cTest))) ## check that none is double

dfTrain <- dfNewAll[cTrain,]
dfVal <- dfNewAll[cVal,]
dfTest <- dfNewAll[cTest,]
# nrow(dfTrain) + nrow(dfVal) + nrow(dfTest) ## check all rows in?


############################################
## Analysis of physiological adaptation

trm <- 0:30

HFtomet <- 0.092
metadaptPhys <- 0.193 * HFtomet
dmetadaptPhys <- ifelse(((trm - 18) * metadaptPhys) < 0, 
        0, ((trm - 18) * metadaptPhys))

# Hori adapt
dMetadaptHori <- 0.208 * (trm - 15) * 1.16222

METtmr <- 0.8 - dmetadaptPhys
BMR <- METtmr*58

BMRFroehlem <- (1436.4 +168 - 5.3 * trm) * 0.04843 / 1.7 # convert kcal/d to watt and per m²
BMRFroehlef <- (1436.4 - 5.3 * trm) * 0.04843 / 1.7 # convert kcal/d to watt and per m²

dfFroehle <- data.frame(trm = trm, BMRf = BMRFroehlef, BMRm = BMRFroehlem)

# HoriDatset was exported by https://apps.automeris.io/wpd/
dfHori <- read.csv(paste(pathFiles, "HoriDataset.csv", sep=""), sep=";", dec=",")
dfHori$BMRWatt <- dfHori$BMR*1.162

lmHoriLin <- lm(BMRWatt~trm, data=dfHori)
lmHoriPol <- lm(BMRWatt~poly(trm,3,raw=TRUE), data=dfHori)

summary(lmHoriLin)
effectsize(lmHoriLin)

summary(lmHoriPol)
effectsize(lmHoriPol)

dfTrain$MetAd <- dfTrain$Met + ((-0.23424*dfTrain$Trm)/58.2)

# check same result when converting earlier?
dfHori$MET <- dfHori$BMRWatt/58.2
lmHoriMETLin <- lm(MET~trm, data=dfHori)

summary(lmHoriMETLin)
effectsize(lmHoriMETLin)
# all fine, slope results in trm -0.0040247



############################################
## Analysis of behavioural adaptation


#summary(lm(CloO~Trm, data=dfTrain)) #R2 .13
summary(lmCloO1 <- lm(log10(CloO)~Trm, data=dfTrain)) #R2 .13
# model diagnosis á la https://www.theanalysisfactor.com/linear-models-r-diagnosing-regression-model/
par(mfrow = c(2,2)); plot(lm(CloO~Trm, data=dfTrain)) 
par(mfrow = c(1,1))
par(mfrow = c(2,2)); plot(lm(log10(CloO)~Trm, data=dfTrain)) 
par(mfrow = c(1,1))

### check for density
# skewness(dfTrain$CloO, na.rm = TRUE)
# skewness(log10(dfTrain$CloO), na.rm = TRUE)
# dfTrain$CloOLog <- log10(dfTrain$CloO)
# ggdensity(dfTrain, x = "CloO", fill = "lightgray") +
  # stat_overlay_normal_density(color = "red", linetype = "dashed")
# ggdensity(dfTrain, x = "CloOLog", fill = "lightgray") +
  # stat_overlay_normal_density(color = "red", linetype = "dashed")

# ## Box-Cox-Transformation: https://statologie.de/box-cox-transformation-r/
# #optimales Lambda für die Box-Cox-Transformation finden
# bc <- boxcox(CloO~Trm, data=dfTrain)
# (lambda <- bc$x[which.max(bc$y)])
# new_model <- lm(((CloO^lambda-1)/lambda) ~ Trm, data=dfTrain)
# par(mfrow = c(2,2)); plot(new_model); par(mfrow = c(1,1))

# bc <- boxcox(CloO~Trm*CoolingStrategyBuilding*BuildingType, data=dfTrain)
# (lambda <- bc$x[which.max(bc$y)])
# new_model <- lm(((CloO^lambda-1)/lambda) ~ Trm*CoolingStrategyBuilding*BuildingType, data=dfTrain)
# par(mfrow = c(2,2)); plot(new_model); par(mfrow = c(1,1))

# bc <- boxcox(CloO~Trm*Met*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)
# (lambda <- bc$x[which.max(bc$y)])
# new_model <- lm(((CloO^lambda-1)/lambda) ~ Trm*Met*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)
# par(mfrow = c(2,2)); plot(new_model); par(mfrow = c(1,1))

# poisson.model <- glm(CloO~Trm*Met*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, dfTrain, family = Gamma(link = "inverse"))
# summary(poisson.model)
# par(mfrow = c(2,2)); plot(poisson.model); par(mfrow = c(1,1))


# cd <- cooks.distance(new_model)
# tail(sort(cd))

# #summary(lm(CloO~Trm*CoolingStrategyBuilding*BuildingType*Season*koppenClimateClassification, data=dfTrain)) #adjR2 .52 - not justified by minor increase
# #summary(lm(log10(CloO)~Trm*CoolingStrategyBuilding*BuildingType*Season*koppenClimateClassification, data=dfTrain)) #adjR2 .56 - not justified by minor increase - potentially explaining differences in regression lines found elsewhere

# # https://www.dataquest.io/blog/tutorial-poisson-regression-in-r/
# poisson.model <- glm(CloO ~ Trm*Met*CoolingStrategyBuilding*BuildingType*Season, dfTrain, family = poisson(link = "log"))
# summary(poisson.model)
# par(mfrow = c(2,2)); plot(poisson.model); par(mfrow = c(1,1))

# poisson.model <- glm(CloO ~ Trm*CoolingStrategyBuilding*BuildingType*Season, dfTrain, family = quasipoisson(link = "log"))
# summary(poisson.model)
# par(mfrow = c(2,2)); plot(poisson.model); par(mfrow = c(1,1))


# poisson.model <- glm(CloO ~ Trm*CoolingStrategyBuilding*BuildingType*Season, dfTrain, family = Gamma(link = "inverse"))
# summary(poisson.model)
# par(mfrow = c(2,2)); plot(poisson.model); par(mfrow = c(1,1))


# poisson.model <- glm(CloO ~ Trm, dfTrain, family = poisson(link = "log"))
# summary(poisson.model)
# par(mfrow = c(2,2)); plot(poisson.model); par(mfrow = c(1,1))

# decision use log10 for comparison with Schiavon and simplicity - especially due to predictive purpose.

#summary(lm(CloO~Trm*Met, data=dfTrain)) #R2 .13
summary(lmCloO2 <- lm(log10(CloO)~Trm*Met, data=dfTrain)) #R2 .13
summary(lmCloO3 <- lm(log10(CloO)~Trm*MetAd, data=dfTrain)) #R2 .13

summary(lmCloO4 <- lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding, data=dfTrain)) #R2 .21
summary(lmCloO5 <- lm(log10(CloO)~Trm*MetAd*BuildingType, data=dfTrain)) #R2 .29
# merge certain levels
# levels(x)[c(3, 5)] <- levels(x)[c(2, 4)]

# only distinguish between office and others
dfTrain$BuildingTypeSimpl <- dfTrain$BuildingType
levels(dfTrain$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfTrain$BuildingTypeSimpl)[c(2,2,2)]
# table(dfTrain$BuildingTypeSimpl)
# table(dfTrain$BuildingType)
summary(lm(log10(CloO)~Trm*MetAd*BuildingTypeSimpl, data=dfTrain)) #R2 .28

# only distinguish between multi family and others
dfTrain$BuildingTypeSimpl <- dfTrain$BuildingType
levels(dfTrain$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfTrain$BuildingTypeSimpl)[c(3,3,2)]
# table(dfTrain$BuildingTypeSimpl)
# table(dfTrain$BuildingType)
summary(lm(log10(CloO)~Trm*MetAd*BuildingTypeSimpl, data=dfTrain)) #R2 .18 - not way forward

# combine office and classroom vs. all others
dfTrain$BuildingTypeSimpl <- dfTrain$BuildingType
levels(dfTrain$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfTrain$BuildingTypeSimpl)[c(3,2,2)]
# table(dfTrain$BuildingTypeSimpl)
# table(dfTrain$BuildingType)
summary(lm(log10(CloO)~Trm*MetAd*BuildingTypeSimpl, data=dfTrain)) #R2 .18 - not way forward

### DECISION
# only distinguish between office and others
dfTrain$BuildingTypeSimpl <- dfTrain$BuildingType
levels(dfTrain$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfTrain$BuildingTypeSimpl)[c(2,2,2)]
# table(dfTrain$BuildingTypeSimpl)
# table(dfTrain$BuildingType)
summary(lmCloO6 <- lm(log10(CloO)~Trm*MetAd*BuildingTypeSimpl, data=dfTrain)) #R2 .28

summary(lmCloO7 <- lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl, data=dfTrain)) #adjR2

summary(lmCloO8 <- lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*Season, data=dfTrain)) #adjR2 .43
par(mfrow = c(2,2)); plot(lm(log10(CloO)~Trm*CoolingStrategyBuilding*BuildingTypeSimpl*Season, data=dfTrain)) 
par(mfrow = c(1,1))

dfTrain$SeasonSimpl <- dfTrain$Season
levels(dfTrain$SeasonSimpl)[c(1,2)] <- levels(dfTrain$SeasonSimpl)[c(3,3)]
# table(dfTrain$SeasonSimpl)
# table(dfTrain$Season)
summary(lm(CloO~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)) #adjR2 .42
par(mfrow = c(2,2)); plot(lm(CloO~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)); par(mfrow = c(1,1))

summary(lmCloO9 <- lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)) #R2 .42
par(mfrow = c(2,2)); plot(lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain)); par(mfrow = c(1,1))

# --> much simpler and still good
tab_model(lmCloO1, lmCloO2, file="lmCloO12.doc")
tab_model(lmCloO3, lmCloO4, file="lmCloO34.doc")
tab_model(lmCloO5, lmCloO6, file="lmCloO56.doc")
tab_model(lmCloO7, lmCloO8, file="lmCloO78.doc")
tab_model(lmCloO9, file="lmCloO9.doc")


############################################
## Analysis of psychological adaptation PMV + SET


dfTrain$ASVpos <- dfTrain$ASV + 4

summary(lm(ASVpos~Lraw, data=dfTrain))
summary(lm(log(ASVpos)~ Lraw, data=dfTrain)) # not a better fit

summary(lm(ASVpos~Lraw+Met+Trm, data=dfTrain))

par(mfrow = c(2,2)); plot(lm(ASVpos~Lraw+Met+Trm, data=dfTrain)); par(mfrow = c(1,1))

###########################################
# LATEST HERE LOAD calcPMVL!!!
###########################################

maxLength <- nrow(dfTrain)
dfNew2 <- sapply(seq(maxLength), function(i) {calcPMVL(dfTrain$Ta[i], dfTrain$Trad[i], dfTrain$av[i], dfTrain$rH[i], dfTrain$Clo[i], dfTrain$MetAd[i]) } )

dfNew2 <- as.data.frame(t(as.data.frame(dfNew2))) # change from list to data frame and transpose rows and columns
dfNew2$Lraw <- as.numeric(as.character(dfNew2$Lraw))
dfNew2$pmv <- as.numeric(as.character(dfNew2$pmv))
dfNew2$ppd <- as.numeric(as.character(dfNew2$ppd))

dfTrain$LrawAd <- dfNew2$Lraw

summary(lm(ASVpos~Lraw, data=dfTrain))
# plot(ASVpos~Lraw, data=dfTrain)
# abline(lm(ASVpos~Lraw, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~Lraw, data=dfTrain)); par(mfrow = c(1,1))

summary(lmPsy1 <- lm(ASVpos~LrawAd, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~LrawAd, data=dfTrain)); par(mfrow = c(1,1))
plot(ASVpos~LrawAd, data=dfTrain)

summary(lmPsy2 <- lm(ASVpos~LrawAd*MetAd, data=dfTrain))
summary(lmPsy3 <- lm(ASVpos~LrawAd*Trm, data=dfTrain))

#summary(lm(ASVpos~Lraw+MetAd+Trm, data=dfTrain))
# summary(lm(ASVpos~LrawAd+MetAd+Trm, data=dfTrain))

# summary(lm(ASVpos~Lraw+MetAd+Trm, data=dfTrain))
# summary(lm(ASVpos~LrawAd+MetAd+Trm, data=dfTrain))
# par(mfrow = c(2,2)); plot(lm(ASVpos~LrawAd+MetAd+Trm, data=dfTrain)); par(mfrow = c(1,1))

# bc <- boxcox(ASVpos~LrawAd+MetAd+Trm, data=dfTrain)
# (lambda <- bc$x[which.max(bc$y)])
# new_model <- lm(((ASVpos^lambda-1)/lambda) ~ LrawAd+MetAd+Trm, data=dfTrain)
# par(mfrow = c(2,2)); plot(new_model); par(mfrow = c(1,1))

# summary(lm(log(ASVpos)~ Lraw, data=dfTrain))
# summary(lm(log(ASVpos)~ LrawAd, data=dfTrain))

summary(lmPsy4 <- lm(ASVpos~LrawAd*MetAd*Trm, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~LrawAd*MetAd*Trm, data=dfTrain)); par(mfrow = c(1,1))

summary(lm(ASVpos~LrawAd*MetAd*Trm, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~LrawAd*MetAd*Trm, data=dfTrain)); par(mfrow = c(1,1))

summary(lmPsy5 <- lm(ASVpos~LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm, data=dfTrain)); par(mfrow = c(1,1))

# end first step without building characts - use as one model for validation procedure

summary(lmPsy6 <- lm(ASVpos~LrawAd*MetAd*Trm*CoolingStrategyBuilding, data=dfTrain))
summary(lmPsy7 <- lm(ASVpos~LrawAd*MetAd*Trm*BuildingType, data=dfTrain))
summary(lmPsy8 <- lm(ASVpos~LrawAd*MetAd*Trm*BuildingTypeSimpl, data=dfTrain))

summary(lm(ASVpos~LrawAd*MetAd*Trm*CoolingStrategyBuilding*BuildingType, data=dfTrain))
summary(lmPsy9 <- lm(ASVpos~LrawAd*MetAd*Trm*CoolingStrategyBuilding*BuildingTypeSimpl, data=dfTrain))
summary(lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*CoolingStrategyBuilding, data=dfTrain))
summary(lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*BuildingTypeSimpl, data=dfTrain))
summary(lmPsy10 <- lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*CoolingStrategyBuilding*BuildingTypeSimpl, data=dfTrain))
par(mfrow = c(2,2)); plot(lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*CoolingStrategyBuilding*BuildingTypeSimpl, data=dfTrain)); par(mfrow = c(1,1))

tab_model(lmPsy1, lmPsy2, file="lmPsy12.doc")
tab_model(lmPsy3, lmPsy4, file="lmPsy34.doc")
tab_model(lmPsy5, lmPsy6, file="lmPsy56.doc")
tab_model(lmPsy7, lmPsy8, file="lmPsy78.doc")
tab_model(lmPsy9, lmPsy10, file="lmPsy10.doc")



# for discussion only
summary(lm(ASVpos~LrawAd*MetAd*Trm*Age, data=dfTrain))
summary(lm(ASVpos~LrawAd*MetAd*Trm*Sex, data=dfTrain))

summary(lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*CoolingStrategyBuilding*BuildingTypeSimpl*Age*Sex, data=dfTrain))

# end for discussion only

#### end model development
########################################################################
########################################################################


########################################################################
########################################################################
### here validation of models

dfVal$ASVkl <- cut(dfVal$ASV,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))

# classic PMV
PMVcl <- NA
for(i in 1:nrow(dfVal)){
PMVcl[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], dfVal$Clo[i], dfVal$Met[i])$pmv
}
dfVal$PMVcl <- ifelse(PMVcl < (-3), (-3), ifelse(PMVcl > 3, 3, PMVcl))
dfVal$PMVclkl <- cut(dfVal$PMVcl,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVcl),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVclkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVclkl),2)
d = dfVal$ASV-dfVal$PMVcl
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# phys Schweiker (ATHBv1)
PMVphSchw <- NA
for(i in 1:nrow(dfVal)){
metadaptPhys <- 0.017756
dmetadaptPhys <- ifelse(((dfVal$Trm[i] - 18) * metadaptPhys) < 0, 0, ((dfVal$Trm[i] - 18) * metadaptPhys))
metAd <- dfVal$Met[i] - dmetadaptPhys 
PMVphSchw[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], dfVal$Clo[i], metAd)$pmv
}
dfVal$PMVphSchw <- ifelse(PMVphSchw < (-3), (-3), ifelse(PMVphSchw > 3, 3, PMVphSchw))
dfVal$PMVphSchwkl <- cut(dfVal$PMVphSchw,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVphSchw),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVphSchwkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVphSchwkl),2)
d = dfVal$ASV-dfVal$PMVphSchw
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2


# phys HORIabs (ATHBv2)
PMVphHORI <- NA
for(i in 1:nrow(dfVal)){
metAd <- dfVal$Met[i] + ((-0.23424*dfVal$Trm[i])/58.2) 

PMVphHORI[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], dfVal$Clo[i], metAd)$pmv
}
dfVal$PMVphHORI <- ifelse(PMVphHORI < (-3), (-3), ifelse(PMVphHORI > 3, 3, PMVphHORI))
dfVal$PMVphHORIkl <- cut(dfVal$PMVphHORI,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVphHORI),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVphHORIkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVphHORIkl),2)
d = dfVal$ASV-dfVal$PMVphHORI
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2


# beh Schweiker/Wagner (ATHBv1)
dfVal$MetAd <- dfVal$Met + ((-0.23424*dfVal$Trm)/58.2)

PMVbehSchw <- NA
for(i in 1:nrow(dfVal)){
Cload <- 1.252594 + dfVal$Trm[i] * (-0.03023063) 
Cload <- ifelse(Cload > 1, 1, ifelse(Cload < 0.46, 0.46, Cload))
PMVbehSchw[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, dfVal$MetAd[i])$pmv
}
dfVal$PMVbehSchw <- ifelse(PMVbehSchw < (-3), (-3), ifelse(PMVbehSchw > 3, 3, PMVbehSchw))
dfVal$PMVbehSchwkl <- cut(dfVal$PMVbehSchw,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVbehSchw),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVbehSchwkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVbehSchwkl),2)
d = dfVal$ASV-dfVal$PMVbehSchw
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# beh trm+Co

PMVbehTrm <- NA
lmBehTrm <- lm(log10(CloO)~Trm*MetAd, data=dfTrain)
for(i in 1:nrow(dfVal)){
Cload <- coef(lmBehTrm)[1] + coef(lmBehTrm)[2]*dfVal$Trm[i] + coef(lmBehTrm)[3]*dfVal$MetAd[i] + coef(lmBehTrm)[4]*dfVal$Trm[i]*dfVal$MetAd[i]
Cload <- 10^Cload
PMVbehTrm[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, dfVal$MetAd[i])$pmv
}
dfVal$PMVbehTrm <- ifelse(PMVbehTrm < (-3), (-3), ifelse(PMVbehTrm > 3, 3, PMVbehTrm))
dfVal$PMVbehTrmkl <- cut(dfVal$PMVbehTrm,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVbehTrm),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVbehTrmkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVbehTrmkl),2)
d = dfVal$ASV-dfVal$PMVbehTrm
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# beh trm+Co+bldg

# prep variables
# only distinguish between office and others
dfVal$BuildingTypeSimpl <- dfVal$BuildingType
levels(dfVal$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfVal$BuildingTypeSimpl)[c(2,2,2)]
# table(dfVal$BuildingTypeSimpl)
# table(dfVal$BuildingType)
dfVal$SeasonSimpl <- dfVal$Season
levels(dfVal$SeasonSimpl)[c(1,2)] <- levels(dfVal$SeasonSimpl)[c(3,3)]
# table(dfVal$SeasonSimpl)
# table(dfVal$Season)

PMVbehTrmBldg <- NA
lmBehTrmBldg <- lm(log10(CloO)~Trm*MetAd*CoolingStrategyBuilding*BuildingTypeSimpl*SeasonSimpl, data=dfTrain) #R2 .42
for(i in 1:nrow(dfVal)){
Cload <- 10^predict(lmBehTrmBldg, newdata = dfVal[i,])
PMVbehTrmBldg[i] <- calcPMVPPD(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, dfVal$MetAd[i])$pmv
}
dfVal$PMVbehTrmBldg <- ifelse(PMVbehTrmBldg < (-3), (-3), ifelse(PMVbehTrmBldg > 3, 3, PMVbehTrmBldg))
dfVal$PMVbehTrmBldgkl <- cut(dfVal$PMVbehTrmBldg,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVbehTrmBldg),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVbehTrmBldgkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVbehTrmBldgkl),2)
d = dfVal$ASV-dfVal$PMVbehTrmBldg
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

#psy Trm
maxLength <- nrow(dfVal)
dfNew2 <- sapply(seq(maxLength), function(i) {calcPMVL(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], dfVal$Clo[i], dfVal$MetAd[i]) } )

dfNew2 <- as.data.frame(t(as.data.frame(dfNew2))) # change from list to data frame and transpose rows and columns
dfNew2$Lraw <- as.numeric(as.character(dfNew2$Lraw))
dfVal$LrawAd <- dfNew2$Lraw

lmPsyTrm <- lm(ASVpos~LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm, data=dfTrain)
PMVPsyTrm <- NA
for(i in 1:nrow(dfVal)){
PMVPsyTrm[i] <- predict(lmPsyTrm, newdata = dfVal[i,])-4
}
dfVal$PMVPsyTrm <- ifelse(PMVPsyTrm < (-3), (-3), ifelse(PMVPsyTrm > 3, 3, PMVPsyTrm))
dfVal$PMVPsyTrmkl <- cut(dfVal$PMVPsyTrm,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVPsyTrm),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVPsyTrmkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVPsyTrmkl),2)
d = dfVal$ASV-dfVal$PMVPsyTrm
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

#psy TrmBldg
lmPsyTrmBldg <- lm(ASVpos~(LrawAd+MetAd+Trm+LrawAd:Trm+MetAd:Trm+LrawAd:MetAd:Trm)*CoolingStrategyBuilding*BuildingTypeSimpl, data=dfTrain)
PMVPsyTrmBldg <- NA
for(i in 1:nrow(dfVal)){
PMVPsyTrmBldg[i] <- predict(lmPsyTrmBldg, newdata = dfVal[i,])-4
}
dfVal$PMVPsyTrmBldg <- ifelse(PMVPsyTrmBldg < (-3), (-3), ifelse(PMVPsyTrmBldg > 3, 3, PMVPsyTrmBldg))
dfVal$PMVPsyTrmBldgkl <- cut(dfVal$PMVPsyTrmBldg,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVPsyTrmBldg),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVPsyTrmBldgkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVPsyTrmBldgkl),2)
d = dfVal$ASV-dfVal$PMVPsyTrmBldg
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# combined
# phys + beh + psy (Trm only)

PMVcombined <- NA
for(i in 1:nrow(dfVal)){
metAd <- dfVal$Met[i] + ((-0.23424*dfVal$Trm[i])/58.2)

Cload <- coef(lmBehTrm)[1] + coef(lmBehTrm)[2]*dfVal$Trm[i] + coef(lmBehTrm)[3]*metAd + coef(lmBehTrm)[4]*dfVal$Trm[i]*metAd
Cload <- 10^Cload
Lraw <- calcPMVL(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, metAd)$Lraw
newDat <- data.frame(LrawAd = Lraw, MetAd = metAd, Trm = dfVal$Trm[i])
PMVcombined[i] <- predict(lmPsyTrm, newdata = newDat[1,])-4
}
dfVal$PMVcombined <- ifelse(PMVcombined < (-3), (-3), ifelse(PMVcombined > 3, 3, PMVcombined))
dfVal$PMVcombinedkl <- cut(dfVal$PMVcombined,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVcombined),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVcombinedkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVcombinedkl),2)
d = dfVal$ASV-dfVal$PMVcombined
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# phys + beh + psy (Trm + bldg)
PMVcombBldg <- NA
for(i in 1:nrow(dfVal)){
metAd <- dfVal$Met[i] + ((-0.23424*dfVal$Trm[i])/58.2)
Cload <- 10^predict(lmBehTrmBldg, newdata = dfVal[i,])
Lraw <- calcPMVL(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, metAd)$Lraw
newDat <- data.frame(LrawAd = Lraw, MetAd = metAd, Trm = dfVal$Trm[i], CoolingStrategyBuilding = dfVal$CoolingStrategyBuilding[i], BuildingTypeSimpl = dfVal$BuildingTypeSimpl[i])
PMVcombBldg[i] <- predict(lmPsyTrmBldg, newdata = newDat[1,])-4
}
dfVal$PMVcombBldg <- ifelse(PMVcombBldg < (-3), (-3), ifelse(PMVcombBldg > 3, 3, PMVcombBldg))
dfVal$PMVcombBldgkl <- cut(dfVal$PMVcombBldg,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVcombBldg),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVcombBldgkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVcombBldgkl),2)
d = dfVal$ASV-dfVal$PMVcombBldg
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

# phys + beh + psy (Trm only + clo schweiker)
PMVcombined <- NA
for(i in 1:nrow(dfVal)){
metAd <- dfVal$Met[i] + ((-0.23424*dfVal$Trm[i])/58.2)
Cload <- 1.252594 + dfVal$Trm[i] * (-0.03023063) 
Cload <- ifelse(Cload > 1, 1, ifelse(Cload < 0.46, 0.46, Cload))
Lraw <- calcPMVL(dfVal$Ta[i], dfVal$Trad[i], dfVal$av[i], dfVal$rH[i], Cload, metAd)$Lraw
newDat <- data.frame(LrawAd = Lraw, MetAd = metAd, Trm = dfVal$Trm[i])
PMVcombined[i] <- predict(lmPsyTrm, newdata = newDat[1,])-4
}
dfVal$PMVcombined <- ifelse(PMVcombined < (-3), (-3), ifelse(PMVcombined > 3, 3, PMVcombined))
dfVal$PMVcombinedkl <- cut(dfVal$PMVcombined,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfVal$ASV, dfVal$PMVcombined),3)
round(calcAvgAcc(dfVal$ASVkl, dfVal$PMVcombinedkl),2)
round(calcTPRTSV(dfVal$ASVkl, dfVal$PMVcombinedkl),2)
d = dfVal$ASV-dfVal$PMVcombined
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfVal$ASV-mean(dfVal$ASV))^2))
mse; mae; rmse; R2

###########################################################
## testing on test data
dfTest$ASVkl <- cut(dfTest$ASV,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))

# classic PMV
PMVcl <- NA
for(i in 1:nrow(dfTest)){
PMVcl[i] <- calcPMVPPD(dfTest$Ta[i], dfTest$Trad[i], dfTest$av[i], dfTest$rH[i], dfTest$Clo[i], dfTest$Met[i])$pmv
}
dfTest$PMVcl <- ifelse(PMVcl < (-3), (-3), ifelse(PMVcl > 3, 3, PMVcl))
dfTest$PMVclkl <- cut(dfTest$PMVcl,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfTest$ASV, dfTest$PMVcl),3)
round(calcAvgAcc(dfTest$ASVkl, dfTest$PMVclkl),2)
round(calcTPRTSV(dfTest$ASVkl, dfTest$PMVclkl),2)
d = dfTest$ASV-dfTest$PMVcl
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfTest$ASV-mean(dfTest$ASV))^2))
mse; mae; rmse; R2

# combined
# phys + beh + psy (Trm only)
# prep variables
# only distinguish between office and others
dfTest$BuildingTypeSimpl <- dfTest$BuildingType
levels(dfTest$BuildingTypeSimpl)[c(1,4,5)] <- levels(dfTest$BuildingTypeSimpl)[c(2,2,2)]
# table(dfTest$BuildingTypeSimpl)
# table(dfTest$BuildingType)
dfTest$SeasonSimpl <- dfTest$Season
levels(dfTest$SeasonSimpl)[c(1,2)] <- levels(dfTest$SeasonSimpl)[c(3,3)]
# table(dfTest$SeasonSimpl)
# table(dfTest$Season)

PMVcombined <- NA
for(i in 1:nrow(dfTest)){
metAd <- dfTest$Met[i] + ((-0.23424*dfTest$Trm[i])/58.2)
Cload <- coef(lmBehTrm)[1] + coef(lmBehTrm)[2]*dfTest$Trm[i] + coef(lmBehTrm)[3]*metAd + coef(lmBehTrm)[4]*dfTest$Trm[i]*metAd
Cload <- 10^Cload
Lraw <- calcPMVL(dfTest$Ta[i], dfTest$Trad[i], dfTest$av[i], dfTest$rH[i], Cload, metAd)$Lraw
newDat <- data.frame(LrawAd = Lraw, MetAd = metAd, Trm = dfTest$Trm[i])
PMVcombined[i] <- predict(lmPsyTrm, newdata = newDat[1,])-4
}
dfTest$PMVcombined <- ifelse(PMVcombined < (-3), (-3), ifelse(PMVcombined > 3, 3, PMVcombined))
dfTest$PMVcombinedkl <- cut(dfTest$PMVcombined,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfTest$ASV, dfTest$PMVcombined),3)
round(calcAvgAcc(dfTest$ASVkl, dfTest$PMVcombinedkl),2)
round(calcTPRTSV(dfTest$ASVkl, dfTest$PMVcombinedkl),2)
d = dfTest$ASV-dfTest$PMVcombined
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfTest$ASV-mean(dfTest$ASV))^2))
mse; mae; rmse; R2

# phys + beh + psy (Trm + bldg)

PMVcombBldg <- NA
dfTest$MetAd <- dfTest$Met + ((-0.23424*dfTest$Trm)/58.2)

for(i in 1:nrow(dfTest)){
metAd <- dfTest$Met[i] + ((-0.23424*dfTest$Trm[i])/58.2)
Cload <- 10^predict(lmBehTrmBldg, newdata = dfTest[i,])
Lraw <- calcPMVL(dfTest$Ta[i], dfTest$Trad[i], dfTest$av[i], dfTest$rH[i], Cload, metAd)$Lraw
newDat <- data.frame(LrawAd = Lraw, MetAd = metAd, Trm = dfTest$Trm[i], CoolingStrategyBuilding = dfTest$CoolingStrategyBuilding[i], BuildingTypeSimpl = dfTest$BuildingTypeSimpl[i])
PMVcombBldg[i] <- predict(lmPsyTrmBldg, newdata = newDat[1,])-4
}
dfTest$PMVcombBldg <- ifelse(PMVcombBldg < (-3), (-3), ifelse(PMVcombBldg > 3, 3, PMVcombBldg))
dfTest$PMVcombBldgkl <- cut(dfTest$PMVcombBldg,breaks=c(-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5),labels=c("-3","-2","-1","0","1","2","3"))
round(calcBias(dfTest$ASV, dfTest$PMVcombBldg),3)
round(calcAvgAcc(dfTest$ASVkl, dfTest$PMVcombBldgkl),2)
round(calcTPRTSV(dfTest$ASVkl, dfTest$PMVcombBldgkl),2)
d = dfTest$ASV-dfTest$PMVcombBldg
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)
R2 = 1-(sum((d)^2)/sum((dfTest$ASV-mean(dfTest$ASV))^2))
mse; mae; rmse; R2

