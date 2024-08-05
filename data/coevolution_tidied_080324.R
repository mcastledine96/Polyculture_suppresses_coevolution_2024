#Coevolution analysis

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)

res_dat <- read.csv("data/res_tophage_280121.csv", header = T)

summary(res_dat)

#look at O

o_res <- filter(res_dat, bact == "o")

#plot

#calculate prop susc
o_res <- mutate(o_res, prop = susc / total)
o_res <- mutate(o_res, res = total - susc,
                prop2 = res/total)
o_res <- na.omit(o_res)

o_res_means <- group_by(o_res, bact, cond, time, phage)  %>%
  summarise(m_prop = mean(prop2),
            se = sd(prop2)/sqrt(length(prop2))) %>%
  ungroup()

fac_labs <- c("(a) Polyculture", "(b) Monoculture")
names(fac_labs) <- c("comm", "mono")

ggplot() +
  theme_bw() +
  geom_point(data = o_res, aes(x = time, y = prop2, group = phage, col = phage), size = 1, position = position_jitterdodge(0.8), alpha = 0.5) +
  facet_wrap(~cond, labeller = labeller(cond = fac_labs)) +
  geom_point(data = o_res_means, aes(x = time, y = m_prop, group = phage, col = phage), size = 2, position = position_dodge(0.8)) +
  geom_errorbar(data = o_res_means, aes(x = time, ymin = m_prop - se, ymax = m_prop + se, group = phage, col = phage), position = position_dodge(0.8), width = 0.5) +
  palettetown::scale_color_poke(pokemon = 'jigglypuff', spread = 4, labels = c("Anc.", "2", "4", "6")) +
  ylab("Phage resistance") +
  xlab("Bacteria time point") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  scale_x_discrete(labels = c("2", "4", "6")) +
  labs(col = "Phage time point")

#analysis

o_res$rep2 <- interaction(o_res$cond, o_res$rep)
o_res <- mutate(o_res,
                res = total - susc)

o_mod1 <- glmer(cbind(susc, res) ~ time * phage * cond + (1|rep2), family = binomial, data = o_res)
blmeco::dispersion_glmer(o_mod1) #no overdispersion
o_mod1.1 <- glmer(cbind(susc, res) ~ time + phage + cond + time:cond + time:phage + phage:cond + (1|rep2), family = binomial, data = o_res)
anova(o_mod1, o_mod1.1) #0.09331
o_mod2 <- glmer(cbind(susc, res) ~ time + phage + cond + time:phage + phage:cond + (1|rep2), family = binomial, data = o_res) #remove time:cond interaction
anova(o_mod1.1, o_mod2) #no time:cond interaction p = 0.22
o_mod3 <- glmer(cbind(susc, res) ~ time + phage + cond + time:cond + phage:cond + (1|rep2), family = binomial, data = o_res) #remove time:phage interaction
anova(o_mod1.1, o_mod3) #no time:phage interaction 0.08108
o_mod4 <- glmer(cbind(susc, res) ~ time + phage + cond + time:cond + phage:time + (1|rep2), family = binomial, data = o_res) #remove cond:phage interaction
anova(o_mod1.1, o_mod4) #no cond:phage interaction p = 0.263. Drop first
o_mod5 <- glmer(cbind(susc, res) ~ time + phage + cond + phage:time + (1|rep2), family = binomial, data = o_res) #remove time:cond interaction again
anova(o_mod4, o_mod5) #non-sig. p = 0.238
o_mod6 <- glmer(cbind(susc, res) ~ time + phage + cond + time:cond + (1|rep2), family = binomial, data = o_res) #remove phage:time
anova(o_mod4, o_mod6) #non-sig p = 0.0777. Drop time:cond first
o_mod7 <- glmer(cbind(susc, res) ~ time + phage + cond + (1|rep2), family = binomial, data = o_res)
anova(o_mod5, o_mod7) #phage:time interaction non-sig 0.0779
o_mod8 <- glmer(cbind(susc, res) ~ phage + cond + (1|rep2), family = binomial, data = o_res) #remove bact time
anova(o_mod7, o_mod8) #sig effect of bact time
o_mod9 <- glmer(cbind(susc, res) ~ time + cond + (1|rep2), family = binomial, data = o_res) #remove phage time
anova(o_mod7, o_mod9) #sig effect of phage time
o_mod10 <- glmer(cbind(susc, res) ~ time + phage + (1|rep2), family = binomial, data = o_res) #remove cond
anova(o_mod7, o_mod10) #sig effect of condition

#all effects significant as fixed non-interacting effects
summary(o_mod7)
emmeans::emmeans(o_mod7, pairwise ~ phage)
emmeans::emmeans(o_mod7, pairwise ~ cond)
emmeans::emmeans(o_mod7, pairwise ~ time)

#Analysis of resistance to ancestral phage at week two

anc_res <- filter(res_dat, phage == "anc")

anc_res_means <- group_by(anc_res, bact, time, cond) %>%
  summarise(m_prop = mean(prop),
            se = sd(prop)/sqrt(length(prop))) 

anc_res <- na.omit(anc_res)

anc_res <- mutate(anc_res,
                  res = total - susc,
                  prop = res / total)

anc_res <- filter(anc_res, time == "T2")
anc_res_means <- filter(anc_res_means, time == "T2")

ggplot() +
  theme_bw() +
  geom_point(data = anc_res, aes(x = bact, y = prop, group = cond, col = cond), position = position_dodge(0.5), alpha = 0.5) +
  geom_point(data = anc_res_means, aes(x = bact, y = m_prop, group = cond, col = cond), position = position_dodge(0.5), size = 2) +
  geom_errorbar(data = anc_res_means, aes(x = bact, ymin = m_prop - se, ymax = m_prop + se, group = cond, col = cond), position = position_dodge(0.5), width = 0.2) +
  palettetown::scale_color_poke(pokemon = 'nidoranm', spread = 2, labels = c("Polyculture", "Monoculture")) +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  labs(col = "Treatment") +
  ylab("Proportion phage resistant") +
  xlab("Species") +
  scale_x_discrete(labels = c(expression(italic(Ochrobactrum)), expression(italic(Pseudomonas)), expression(italic(Variovorax))))

m1 <- lm(asin(sqrt(prop)) ~ cond * bact, data = anc_res)
m2 <- lm(asin(sqrt(prop)) ~ cond + bact, data = anc_res)
plot(m1)
anova(m1,m2,test="F")
emmeans::emmeans(m1,pairwise~cond|bact)

