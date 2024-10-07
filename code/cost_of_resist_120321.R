##cost of resistance in communities?

library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)

resc <- read.csv("data/cost_res_o_050824.csv", header = T)

#malthusian growth rate

resc <- mutate(resc,
               t_1 = count * 2 * dil * dil_fact,
               m = log(t_1/t_0))

fac_labs <- c("(a) Polyculture", "(b) Monoculture")
names(fac_labs) <- c("comm", "mono")

cols <- c("white","lightgrey", "black", "white","lightgrey", "black")

a <- ggplot(resc, aes(x = treat, y = m, group = interaction(res,treat), fill = res)) +
  geom_boxplot(col = 'darkgrey') +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2) +
  theme_bw() +
  ylab(expression("Growth rate ("~italic("m")~")")) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = "bottom", title = element_text(size =11)) +
  scale_fill_manual(values = cols, labels = c('Ancestral', 'Phage\nresistant', 'Phage\nsusceptible')) +
  scale_x_discrete(labels = c('Polyculture', 'Monoculture')) +
  xlab("Treatment") +
  labs(fill = 'Strain', title = expression(paste('(a)', italic(' Ochrobactrum'), ' sp.')))

#analyses 

resc$gen2 <- as.character(seq(1:6))
resc$gen2[25:36] <- 'A'
resc$gen[is.na(resc$gen)] <- as.factor(seq(1:12))

mod1 <- lmer(m ~ res * treat + (1|gen/gen2), data = resc)
mod2 <- lmer(m ~ res + treat + (1|gen/gen2), data = resc)
anova(mod1,mod2)
mod3 <- lmer(m ~ treat + (1|gen/gen2), data = resc)
anova(mod2,mod3)
mod4 <- lmer(m ~ res + (1|gen/gen2), data = resc)
anova(mod2,mod4)
mod5 <- lmer(m ~ 1 + (1|gen/gen2), data = resc)
anova(mod3,mod5)
emmeans::emmeans(mod3, pairwise ~ treat)
summary(mod1)


### Variovorax

vario_dat <- read.csv('data/vario_monovcom_res_071024.csv', header = T)

#how density and m was calculated
#vario_dat <- mutate(vario_dat, 
                    #den = count * dil * dil_fact * froz,
                    #m = log(den/166666.7))

m1 <- lmer(m ~ rep * polm + (1|rep3), data = vario_dat)
m2 <- lmer(m ~ rep + polm + (1|rep3), data = vario_dat)
anova(m1,m2)
m3 <- lmer(m ~ polm + (1|rep3), data = vario_dat)
anova(m2,m3)
m4 <- lmer(m ~ rep + (1|rep3), data = vario_dat)
anova(m2,m4) #sig
m5 <- lmer(m ~ 1 + (1|rep3), data = vario_dat)
anova(m3,m5)
emmeans::emmeans(m3, pairwise ~ polm)

cols <- c("white","lightgrey", "black", "white","lightgrey", "black")

b <- ggplot(vario_dat, aes(x = polm, y = m, group = interaction(rep,polm), fill = rep)) +
  geom_boxplot(col = 'darkgrey') +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2) +
  theme_bw() +
  ylab(expression("Growth rate ("~italic("m")~")")) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 11, colour = 'black'), axis.title.x = element_text(size = 12, colour = 'black'), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = "bottom", axis.title.y = element_blank(), title = element_text(size =11)) +
  scale_fill_manual(values = cols, labels = c('Ancestral', 'Phage\nresistant', 'Phage\nsusceptible')) +
  scale_x_discrete(labels = c('Polyculture', 'Monoculture')) +
  xlab("Treatment") +
  labs(fill = 'Strain', title = expression(paste('(b)', italic(' Variovorax'), ' sp.')))

a + b + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
