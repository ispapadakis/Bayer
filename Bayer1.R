### R Markdown
# Intro -------------------------------------------------------------------

require(tidyverse)
require(xgboost)
require(Matrix)
require(rpart)
p_colors = c('skyblue','seagreen','tan3')

# Gather Data -------------------------------------------------------------

getwd()
subj_ = read.csv('subject.csv')
rand_ = read.csv('randomization.csv')
effi_ = read.csv('efficacy.csv')

str(subj_)
str(rand_)
str(effi_)

Bayer = subj_ %>% 
  left_join(rand_,by='subject') %>%
  left_join(effi_,by='subject')

summary(Bayer)

#Metadata
#Dependent Variable
y_var = 'nosebleeds'  # Is nosebleeds in Likert Scale?
#Exposure
exposure = 'duration'
# Possible Predictors
x_vars = c('country','eye.colour','tissue.use','previous.year','mucus.viscosity','arm','duration')


# Impute Missing
# Assign Missing to Each Own Category
Bayer$eye.colour = addNA(Bayer$eye.colour)
levels(Bayer$eye.colour)[4] = "N/A"

row_with_missing = which(apply(is.na(Bayer),1,sum)==1)
Bayer[row_with_missing,]
mv.rp = rpart(mucus.viscosity~.,Bayer[,x_vars])
#Impute Modeled to Missing mucus.viscosity
Bayer[row_with_missing,'mucus.viscosity'] = predict(mv.rp,new=Bayer[row_with_missing,x_vars])

# Data Transformation to Faciliate Reports
Bayer$arm = factor(as.character(Bayer$arm), levels=c("PLACEBO","ACTIVE"))
Bayer$tissue.use = factor(as.character(Bayer$tissue.use), levels=c("MEDIUM","HIGH"))

# EDA ---------------------------------------------------------------------

# Analysis of Dependent Variable
Bayer %>% ggplot(aes(duration,nosebleeds)) + geom_point() + geom_jitter() + 
  ggtitle('Nosebleeds by Duration (plus x-axis noise)')
Bayer %>% ggplot(aes(nosebleeds)) + geom_bar() + ggtitle('Nosebleeds Distribution')
Bayer %>% ggplot(aes(duration)) + geom_histogram(bins=10) + ggtitle('Duration Distribution')

Bayer %>% ggplot(aes(cut(duration,6),nosebleeds)) + geom_boxplot(varwidth = T)

count_mean = mean(Bayer$nosebleeds)
barplot(rbind(table(Bayer$nosebleeds),
              dpois(seq(0,5,1),count_mean)*444,
              dnbinom(seq(0,5,1),1,mu=count_mean)*444),
        col=p_colors,beside=T,legend.text=c('Observed','Poisson','Neg Binomial'),
        args.legend = list(bg="white"))
grid()

# How Many Cases with Full Year Duration
table(Bayer$duration == 365)

X = Bayer %>% subset(.,duration==365)
### country
X %>% ggplot(aes(country,nosebleeds)) + geom_boxplot(varwidth = T)
### eye.colour
X %>% ggplot(aes(eye.colour,nosebleeds,fill=eye.colour)) + geom_boxplot(varwidth = T) +
  scale_fill_manual(values=c("grey2", "lightblue2", "tan4", "white"))

Bayer %>% ggplot(aes(country,fill=eye.colour)) + geom_bar() + coord_flip() +
  scale_fill_manual(values=c("grey2", "lightblue2", "tan4", "white")) +
  ggtitle("Interaction Between Eye Color & Country")

### tissue.use
X %>% ggplot(aes(tissue.use,nosebleeds)) + geom_boxplot(varwidth = T, fill='skyblue')

### previous.year
X %>% ggplot(aes(factor(previous.year),nosebleeds)) + geom_boxplot(varwidth = T, fill='skyblue')

### mucus.viscosity
na.omit(X) %>% ggplot(aes(mucus.viscosity,nosebleeds)) + 
  geom_point() + 
  geom_smooth(method = 'loess', formula = 'y ~ x')

mv_cutpoints = c(-Inf,0.5,1,2,4,9)
X %>% ggplot(aes(nosebleeds)) + geom_histogram(bins=6,aes(fill=arm)) + facet_wrap(~cut(mucus.viscosity,mv_cutpoints))
#X %>% ggplot(aes(cut(mucus.viscosity,mv_cutpoints))) + geom_bar(aes(fill=nosebleeds))
X %>% ggplot(aes(cut(mucus.viscosity,mv_cutpoints),nosebleeds)) + geom_boxplot(varwidth = T, fill="skyblue")

# arm
X %>% ggplot(aes(arm,nosebleeds)) + geom_boxplot(fill="skyblue")


# More Placebos with < 365 duration
with(Bayer, mosaicplot(table(arm, (duration < 365)),ylab='Duration < 365', 
                       main="Full Duration by arm", col=c("skyblue","tan")))
with(Bayer, mosaicplot(table(eye.colour, (duration < 365)),ylab='Duration < 365', 
                       main='Full Duration by Eye Color', col=c("skyblue","tan")))

library(rpart)
placebo.rp = rpart(arm~.+Bayer$duration,Bayer[,x_vars])
plotcp(placebo.rp)
placebo.rp = prune(placebo.rp,.015)
plot(placebo.rp, margin=.1)
text(placebo.rp, cex=.8, all=F)
title("Predictor <arm> Not Randomly Assigned")
placebo.rp


# Model -------------------------------------------------------------------


nb.pois = glm(nosebleeds ~  country+eye.colour+tissue.use+previous.year+mucus.viscosity+arm+
               log(duration)+I(abs(mucus.viscosity-.75)<.25)
             , family=poisson, data=Bayer)
summary(nb.pois)
#plot(nb.gpois)

nb.qpois = glm(nosebleeds ~  country+eye.colour+tissue.use+previous.year+mucus.viscosity+arm+
                 log(duration)+I(abs(mucus.viscosity-.75)<.25)
               , family=quasipoisson, data=na.omit(Bayer))
summary(nb.qpois)

#library(sandwich)
#library(lmtest)
#waldtest(nb.pois)
#coef_norm = coeftest(nb.pois)
#coeftest(nb.pois, vcov = sandwich)
#coef_infl = coeftest(nb.pois, vcov=sandwich)
#plot(coef_norm,coef_infl)
#points(coef_norm,coeftest(nb.qpois),col=3)
#abline(0,c(1,1))
#grid()

library(pscl)
nb.hrd = hurdle(nosebleeds ~ country+previous.year+I(abs(mucus.viscosity-.75)<.25)|  
                  country+eye.colour+tissue.use+previous.year+mucus.viscosity+arm+
                 log(duration)+I(abs(mucus.viscosity-.75)<.25), data=Bayer)
summary(nb.hrd)
AIC(nb.hrd)

nb.zinfl = zeroinfl(nosebleeds ~  country+tissue.use+previous.year+mucus.viscosity+arm+
                  log(duration)+I(abs(mucus.viscosity-.75)<.25) |  
                    country+tissue.use+previous.year+mucus.viscosity+arm+
                    log(duration)+eye.colour, dist="negbin", data=Bayer)
summary(nb.zinfl)

# Compare Models
coef_ = merge(data.frame(zero_inflated=coef(nb.zinfl)) ,
      data.frame(hurdle=coef(nb.hrd)),by=0,all=T)
rownames(coef_) = coef_$Row.names
coef_$Row.names = NULL

pois_mod = data.frame(poisson=coef(nb.pois))
pois_mod_zero = data.frame(pois_mod,row.names=paste('zero',rownames(pois_mod),sep='_'))
pois_mod_count = data.frame(pois_mod,row.names=paste('count',rownames(pois_mod),sep='_'))
pois_mod_new = rbind(pois_mod_count,pois_mod_zero)
coef_ = merge(coef_,pois_mod_new,by=0,all=T)
rownames(coef_) = coef_$Row.names
coef_$Row.names = NULL

aic_ = data.frame(poisson=AIC(nb.pois),hurdle=AIC(nb.hrd),zero_inflated=AIC(nb.zinfl),row.names=c('AIC'))
print(rbind(coef_,aic_))
# save for presentation
write.csv(rbind(coef_,aic_),"count_models_est.csv")


# ML Model
# Not a great application of ML (only 444 cases)
numberOfClasses <- length(unique(Bayer$nosebleeds))
xgb_params <- list("objective" = "count:poisson",
                   #"eval_metric" = "mlogloss",
                   "min_child" = 5,
                   "eta"=.01)
nround    <- 450 # number of XGBoost rounds
cv.nfold  <- 5

sparse_matrix = sparse.model.matrix(nosebleeds ~ ., data = subset(Bayer,select=-subject))[,-which(names(Bayer)=='nosebleeds')]
output_vector = Bayer$nosebleeds

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
cv_model <- xgb.cv(params = xgb_params,
                   data = sparse_matrix,
                   label = output_vector,
                   nrounds = nround,
                   nfold = cv.nfold,
                   verbose = F,
                   prediction = T,
                   stratified=FALSE)
plot(cv_model$evaluation_log$test_poisson_nloglik_mean)
n_opt = which.min(cv_model$evaluation_log$test_poisson_nloglik_mean)
nb.xgb = xgboost(data = sparse_matrix, label = output_vector,
                 verbose = F, nrounds = n_opt, params = xgb_params)

importance <- xgb.importance(feature_names = colnames(sparse_matrix), model = nb.xgb)
head(importance)
xgb.plot.importance(importance_matrix = importance)


# Predictions -------------------------------------------------------------
# ML Model
nb.xgb.pred = predict(nb.xgb,newdata = sparse_matrix)

nb.xgb.pred.group = cut(nb.xgb.pred,c(.3,.38,.41,.45,.49,.57,Inf))
tapply(Bayer$nosebleeds,nb.xgb.pred.group,mean)
ggplot(Bayer,aes(nb.xgb.pred.group,nosebleeds)) + geom_boxplot(varwidth = T) +
  labs(y='Actual Nosebleed Count',
       x='Predicted Poisson Rate Group') 
cv_model$evaluation_log[n_opt,]

# Poisson Model
nosebleeds.exp = exp(predict(nb.pois))
boxplot(Bayer$nosebleeds~cut(nosebleeds.exp,c(-Inf,.2,.4,.8,1.5,4)),
        main="nosebleeds v prediction (Poisson Model)",
        ylab='nosebleeds', xlab="Binned Predictions")

# Zero-Inflated Model
nosebleeds.zpred = predict(nb.zinfl)
boxplot(Bayer$nosebleeds~cut(nosebleeds.zpred,c(-Inf,.2,.4,.8,1.5,4)),
        main="nosebleeds v prediction (Zero-Inflated Model)",
        ylab='nosebleeds', xlab="Binned Predictions")

# Zero-Inflated Model
nosebleeds.hrd_pred = predict(nb.hrd)
boxplot(Bayer$nosebleeds~cut(nosebleeds.hrd_pred,c(-Inf,.2,.4,.8,1.5,4)),
        main="nosebleeds v prediction (Hurdle Model)",
        ylab='nosebleeds', xlab="Binned Predictions")


# Specification Tests -----------------------------------------------------

nosebleeds.pois_res = resid(nb.pois)
hist(nosebleeds.pois_res)
qqnorm(nosebleeds.pois_res)
abline(0,c(1,1))
#plot(nb.pois)

nosebleeds.hrd_res = resid(nb.hrd)
hist(nosebleeds.hrd_res)
qqnorm(nosebleeds.hrd_res)
abline(0,c(1,1))

nosebleeds.zinfl_res = resid(nb.zinfl)
hist(nosebleeds.zinfl_res)
qqnorm(nosebleeds.zinfl_res)
abline(0,c(1,1))
