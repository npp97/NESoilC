# TODO: Add comment Author: lenovo


library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(cshapes)
library(MuMIn)
library(mgcv)
library(lattice)
require(doBy)

# ---------------------------
setwd("C:\\Users\\lenovo\\workspace\\NESoilC\\data")
load("xxal.RData")
xyplot(soil_C ~ abv_C | dominant, data = xxal)
xyplot(BAR ~ Diam_C | dominant, data = xxal)

summaryBy(soil_C + abv_C ~ dominant, data = xxal, FUN = mean)

summaryBy(soil_C ~ dominant, data = xxal, FUN = sd)
# first we fit a dummy random variable
xxal$dummy <- rep(1, nrow(xxal))

# now build in spatial autocorrelation

# change coordinates slightly since some sites have exactly the same coordinates
xxal$lat_J <- xxal$lat + (rnorm(length(xxal$lat), 0, 1e-05))
xxal$elv_J <- xxal$elv + (rnorm(length(xxal$elv), 0, 0.001))

cs1Exp <- corExp(1, form = ~lat_J + long)
cs1Exp <- Initialize(cs1Exp, xxal)
corMatrix(cs1Exp)[1:10, 1:4]

# first run some null models to check our random variable structure is appropriate
null.model2 <- lme(log(soil_C) ~ 1, data = xxal, random = ~1 | dominant)
null.model1 <- lme(log(soil_C) ~ 1, data = xxal, random = ~1 | dummy)
null.model3 <- update(null.model2, correlation = corExp(1, form = ~lat_J + long + elv))
null.model4 <- update(null.model2, correlation = corExp(1, form = ~lat_J + long))
null.model5 <- gls(log(soil_C)~elv*soil_depth*Diam_C*abv_C*lat_J*long,data=xxal,correlation = corExp(1, form = ~lat_J + long))
AICc(null.model, null.model2, null.model3, null.model4,null.model5)


# fitting a model that accounts for between study differences is better and we NEED to account for spatial autocorrelation for our results to be
# statistically valid

# now let's fit the models that Liu et al used for: precipitation - a model with a squared term for mean_precip

D1 <- lme(soil_Ccb ~ abv_C + lat + elv, data = xxal, random = ~1 | dominant)
D2 <- gls(soil_Ccb ~ abv_C + lat + elv, data = xxal)
D3 <- lme(soil_Ccb ~ abv_C + lat, data = xxal, random = ~1 | dominant)
D4 <- gls(soil_Ccb ~ abv_C + lat, data = xxal)
D5 <- glm(soil_Ccb ~ abv_C + lat, data = xxal)
D6 <- lm(soil_Ccb ~ abv_C + lat, data = xxal)  #best
D7 <- gls(soil_C ~ abv_C + I(abv_C^2) + lat, data = xxal)

AICc(D1, D2, D3, D4, D5, D6, D7)

summary(D6)
# lm(formula = soil_Ccb ~ abv_C + lat, data = xxal) Residuals: Min 1Q Median 3Q Max -2.54135 -0.76085 -0.05754 0.64515 2.48899 Coefficients: Estimate Std.
# Error t value Pr(>|t|) (Intercept) 10.42613 2.38627 4.369 1.83e-05 *** abv_C -0.03580 0.02997 -1.194 0.23347 lat -0.16149 0.04968 -3.251 0.00131 ** ---
# Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 Residual standard error: 1.012 on 250 degrees of freedom Multiple R-squared: 0.04602,
# Adjusted R-squared: 0.03839 F-statistic: 6.03 on 2 and 250 DF, p-value: 0.002769

# soil_C= (soil_Ccb*0.25+1)^4

# Diam_C_model<-lme(log(soil_C) ~ log(abv_C)+lat, data=xxal,random=~1|dominant) Diam_C_model<-gls(log(soil_C)~ abv_C+lat+elv, data=xxal) D1<-gls(soil_C ~
# abv_C+lat+elv, data=xxal, correlation = corAR1(,form = ~ 1 | dominant),weights = varPower()) D2<-gls(soil_C ~ abv_C+lat*long*elv, data=xxal, correlation
# = corAR1(form = ~ 1 | dominant)) D4<- D3<-gls(soil_C ~ abv_C, data=xxal, correlation = corAR1(form = ~ 1 | dominant)) AICc(D1,D2,D3)

# bb<-boxcox(soil_C ~ abv_C+lat+elv, data=xxal) space<-groupedData(soil_C~abv_C+lat|dominant,data=xxal) model4<-gls(soil_C~abv_C+lat,space)
# summary(model4)

# ll<-lm(residuals(Diam_C_model)~lat+long+elv,data=xxal) anova(ll)

# Diam_C_model<-lme(soil_C ~ abv_C*lat, data=xxal,random=~1|dominant)

# Diam_C_model<-lme(soil_C ~ abv_C,data=xxal,random=~1|dominant,correlation = corExp(1, form = ~ lat_J + long))

# now a global model for use in model averaging that contains all varibles that are needed I haven't included the squared and cubed terms for temp and
# precipitation becuase I don't think they make biological sense All_model<-lme(soil_C~Diam_C*lat_J*long*elv,data=xxal,random=~1|dominant,correlation =
# corExp(1, form = ~ lat_J + long)) All_model<-lme(soil_C~abv_C*lat*long*elv,data=xxal) plot(All_model) now we can dredge the model so that the value of
# each variable in predicting biomass can be assessed rather than using them in isolation
# MS1<-dredge(All_model,evaluate=T,rank=AICc,trace=T,data=xxal,REML=F) poss_mod<-get.models(MS1,subset=delta<7) modsumm<- model.sel(poss_mod, rank =
# 'AICc',fit=T) modsumm2<-subset(modsumm,modsumm$delta<7) modsumm2 averaged<-model.avg(modsumm2,fit=T) run the top model from the model averaging to get r
# squared statistic top_model<-lme(log(AGB)~Age+Mean_precip+Age*Mean_T+Age_sq,data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J +
# Long)) 
