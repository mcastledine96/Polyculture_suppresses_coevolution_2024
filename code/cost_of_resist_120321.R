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

ggplot(resc, aes(x = treat, y = m, group = interaction(res,treat), fill = res)) +
  geom_boxplot(col = 'darkgrey') +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2) +
  theme_bw() +
  ylab(expression("Growth rate ("~italic("m")~")")) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = "bottom") +
  scale_fill_manual(values = cols, labels = c('Ancestral', 'Phage\nresistant', 'Phage\nsusceptible')) +
  scale_x_discrete(labels = c('Polyculture', 'Monoculture')) +
  xlab("Treatment") +
  labs(fill = 'Strain')

#analyses 

resc$gen2 <- as.character(seq(1:6))
resc$gen2[25:36] <- 'A'
resc$gen[is.na(resc$gen)] <- as.factor(seq(1:12))

mod1 <- lmer(m ~ res * treat + (1|gen/gen2), data = resc)
mod2 <- lmer(m ~ res + treat + (1|gen/gen2), data = resc)
anova(mod1,mod2)
mod3 <- lmer(m ~ treat + (1|gen), data = resc)
anova(mod2,mod3)
mod4 <- lmer(m ~ res + (1|gen), data = resc)
anova(mod2,mod4)
mod5 <- lmer(m ~ 1 + (1|gen), data = resc)
anova(mod3,mod5)
emmeans::emmeans(mod3, pairwise ~ treat)
summary(mod1)
