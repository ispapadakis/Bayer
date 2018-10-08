
# Intro -------------------------------------------------------------------

library(tidyverse)


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


# EDA ---------------------------------------------------------------------

# Dependent Variable
y_var = 'nosebleeds'

Bayer %>% ggplot(aes(duration,nosebleeds)) + geom_point() + ggtitle('Nosebleeds by Duration')
Bayer %>% ggplot(aes(nosebleeds)) + geom_bar()
Bayer %>% ggplot(aes(duration)) + geom_histogram()


# Possible Predictors
x_vars = c('country','eye.colour','tissue.use','previous.year','mucus.viscosity','arm')

# How Many Cases with Full Year Duration
table(Bayer$duration == 365)

X = Bayer %>% subset(.,duration==365)
### country
X %>% ggplot(aes(country,nosebleeds)) + geom_boxplot(varwidth = T)
### eye.colour
X %>% ggplot(aes(eye.colour,nosebleeds)) + geom_boxplot(varwidth = T)
### tissue.use
X %>% ggplot(aes(tissue.use,nosebleeds)) + geom_boxplot(varwidth = T)
### previous.year
X %>% ggplot(aes(factor(previous.year),nosebleeds)) + geom_boxplot(varwidth = T)
### mucus.viscosity
X %>% ggplot(aes(mucus.viscosity,nosebleeds)) + 
  geom_point() + 
  geom_smooth()

mv_cutpoints = c(-Inf,0.5,1,2,4,9)
X %>% ggplot(aes(nosebleeds)) + geom_histogram(bins=6,aes(fill=arm)) + facet_wrap(~cut(mucus.viscosity,mv_cutpoints))
X %>% ggplot(aes(cut(mucus.viscosity,mv_cutpoints))) + geom_bar(aes(fill=nosebleeds))
X %>% ggplot(aes(cut(mucus.viscosity,mv_cutpoints),nosebleeds)) + geom_boxplot(varwidth = T)


X %>% ggplot(aes(arm,nosebleeds)) + geom_boxplot()


# More Placebos with < 365 duration
with(Bayer, mosaicplot(table(arm, (duration < 365)),ylab='Duration < 365'))
Bayer$eye.colour = addNA(Bayer$eye.colour)
with(Bayer, mosaicplot(table(eye.colour, (duration < 365)),ylab='Duration < 365'))
library(rpart)
placebo.rp = rpart(arm~.+Bayer$duration,Bayer[,x_vars])
plotcp(placebo.rp)
prune(placebo.rp,.015)



# Model -------------------------------------------------------------------


nb.glm = glm(nosebleeds ~  country+eye.colour+tissue.use+previous.year+mucus.viscosity+arm+
               log(duration)+I(abs(mucus.viscosity-.75)<.25)
             , family=poisson, data=Bayer)
summary(nb.glm)
plot(nb.glm)

nb.glm = glm(nosebleeds ~  country+eye.colour+tissue.use+previous.year+mucus.viscosity+arm+
               log(duration)+I(abs(mucus.viscosity-.75)<.25)
             , family=quasipoisson, data=Bayer)
summary(nb.glm)

