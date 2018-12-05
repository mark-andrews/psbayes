library(readr)
library(dplyr)
library(ggplot2)
library(brms)
libray(psc)

Df <- read_csv('data/affairs.csv')

ggplot(Df,
       aes(x = yearsmarried, y = affairs, col = gender)
) + stat_smooth()

### Gender 

M <- glm(affairs ~ gender, 
         data=Df, 
         family=poisson)

M_bayes <- brm(affairs ~ gender, 
               data=Df, 
               cores = 2, 
               family = poisson(),
               prior = set_prior('normal(0, 100)'), 
               save_all_pars = T)

M_bayes_null <- brm(affairs ~ 1, 
                    data=Df, 
                    cores = 2, 
                    family = poisson(),
                    save_all_pars = T)

waic(M_bayes, M_bayes_null)
bayes_factor(M_bayes, M_bayes_null)

# Gender and age

M <- glm(affairs ~ yearsmarried + gender, 
         data=Df, 
         family=poisson)

M_bayes <- brm(affairs ~ yearsmarried + gender, 
               data=Df, 
               cores = 2, 
               family = poisson(),
               prior = set_prior('normal(0, 100)'), 
               save_all_pars = T)

M_bayes_null <- brm(affairs ~ yearsmarried,  
                    data=Df, 
                    cores = 2, 
                    family = poisson(),
                    prior = set_prior('normal(0, 100)'), 
                    save_all_pars = T)

waic(M_bayes, M_bayes_null)
bayes_factor(M_bayes_null, M_bayes)


#### All variables

M <- glm(affairs ~ gender + age + yearsmarried
         + children + religiousness + education
         + occupation + rating, 
         data=Df, 
         family=poisson)

M_bayes <- brm(affairs ~ gender + age + yearsmarried
               + children + religiousness + education
               + occupation + rating,  
               data=Df, 
               cores = 2, 
               family = poisson(),
               prior = set_prior('normal(0, 100)'), 
               save_all_pars = T)

summary(M_bayes)
marginal_effects(M_bayes)


