#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)
library(flextable)

phage_den <- read.csv("data/phage_densities_101220.csv", header = T)

#calculate densities
phage_den <- mutate(phage_den,
                    pfus = count * dilution * dil_fact,
                    logpfus = log10(pfus))

#sort out infinite values in logged data
phage_den$logpfus[phage_den$logpfus == -Inf] <- 0

#for reviewer 2 also consider if densities were 100 (detection limit)

phage_den$logpfus[phage_den$logpfus == 0 & phage_den$species == 'O'] <- log10(100)

#mean of phage densities

phage_den <- na.omit(phage_den)
phage_means <- group_by(phage_den, type, species, time) %>%
  summarise(mean_den = mean(logpfus),
            se = sd(logpfus)/sqrt(length(logpfus)))

#plot for individual evolution lines
phage_means$time <- substr(phage_means$time, 2, 2)
phage_den$time <- substr(phage_den$time, 2, 2)

phage_lab <- c("(a) ORM_20", "(b) CHF7MC", "(c) VAC_51")
names(phage_lab) <- c("O", "P", "V")

phage_lab <- c("(a) Ochrobactrum phage ORM_20", "(b) Pseudomonas phage CHF7MC", "(c) Variovorax phage VAC_51")

names(phage_lab) <- c("O", "P", "V")

phage_den <- mutate(phage_den,
                    species2 = forcats::fct_recode(species, `(a)~italic("Ochrobactrum")~phage~ORM_20` = 'O', `(b)~italic("Pseudomonas")~phage~CHF7MC` = 'P', `(c)~italic("Variovorax")~phage~VAC_51` = 'V'))

phage_means <- mutate(phage_means,
                      species2 = forcats::fct_recode(species, `(a)~italic("Ochrobactrum")~phage~ORM_20` = 'O', `(b)~italic("Pseudomonas")~phage~CHF7MC` = 'P', `(c)~italic("Variovorax")~phage~VAC_51` = 'V'))

week_0_density <- log10(40000/6)

phage_den$treat <- interaction(phage_den$line, phage_den$species, phage_den$type)

week_0_data <- data.frame(
  time = 0,
  logpfus = week_0_density,
  group = unique(phage_den$treat)  # Replicates and treatment average
)

week_0_data <- week_0_data %>%
  separate(col = group, into = c('line', 'species', 'type'))

week_0_data <- mutate(week_0_data,
                      species2 = forcats::fct_recode(species, `(a)~italic("Ochrobactrum")~phage~ORM_20` = 'O', `(b)~italic("Pseudomonas")~phage~CHF7MC` = 'P', `(c)~italic("Variovorax")~phage~VAC_51` = 'V'))

phage_den <- phage_den[,c(1:3,7,9,10)]

phage_den <- rbind(week_0_data, phage_den)

phage_means$treat <- interaction(phage_means$species, phage_means$type)

week_0_data <- data.frame(
  time = 0,
  mean_den = week_0_density,
  group = unique(phage_means$treat),
  se = 0# Replicates and treatment average
)

week_0_data <- week_0_data %>%
  separate(col = group, into = c('species', 'type'))

week_0_data <- mutate(week_0_data,
                      species2 = forcats::fct_recode(species, `(a)~italic("Ochrobactrum")~phage~ORM_20` = 'O', `(b)~italic("Pseudomonas")~phage~CHF7MC` = 'P', `(c)~italic("Variovorax")~phage~VAC_51` = 'V'))

phage_means <- phage_means[,-7]

phage_means <- rbind(week_0_data, phage_means)

ggplot() +
  geom_hline(yintercept = log10(100), col = 'grey') +
  facet_wrap(~species2, labeller = label_parsed)+
  geom_point(data = phage_den, aes(x = time, y = logpfus, col = type), alpha = 0.3) +
  geom_line(data = phage_den, aes(x = time, y = logpfus, linetype = type, group = interaction(type,line), col = type), alpha = 0.3) +
  theme_bw() +
  ylab("Density log10(PFU/mL)") +
  xlab("Week") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  geom_point(data = phage_means, aes(x = time, y = mean_den, group = type, col = type), size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(data = phage_means, aes(x = time, ymin = mean_den - se, ymax = mean_den + se, group = type, col = type), position = position_dodge(0.5), width = 0.3, size = 0.8) +
  geom_line(data = phage_means, aes(x = time, y = mean_den, group = type, linetype = type, col = type), position = position_dodge(0.5), size = 0.8) +
  labs(shape = "Treatment", linetype = "Treatment", col = "Treatment") +
  scale_shape_discrete(labels = c("Polyculture", "Monoculture")) +
  scale_linetype_discrete(labels = c("Polyculture", "Monoculture")) +
  palettetown::scale_color_poke(pokemon = 'nidoranm', spread = 2, labels = c("Polyculture", "Monoculture"))

phage_den <- filter(phage_den, ! time == 0)

p_phage <- filter(phage_den, species == "P")
o_phage <- filter(phage_den, species == "O")
v_phage <- filter(phage_den, species == "V")

ph_m1 <- lmer(logpfus ~ time * type + (1|line), data = p_phage)
ph_m2 <- lmer(logpfus ~ time + type + (1|line), data = p_phage)
anova(ph_m1, ph_m2) #sig interaction
emmeans::emmeans(ph_m1, pairwise ~ type | time) #sig less phage in mono lines at T2 

op_m1 <- lmer(logpfus ~ time * type + (1|line), data = o_phage)
op_m2 <- lmer(logpfus ~ time + type + (1|line), data = o_phage)
anova(op_m1, op_m2) #non-sig interaction
op_m3 <- lmer(logpfus ~ type + (1|line), data = o_phage)
anova(op_m2, op_m3) #sig effect of time
op_m4 <- lmer(logpfus ~ time + (1|line), data = o_phage)
anova(op_m2, op_m4) #sig effect of type

emmeans::emmeans(op_m2, pairwise ~ time) #phage densities tend to decline over time
emmeans::emmeans(op_m2, pairwise ~ type) #on average, phage densities in monoculture higher

vp_m1 <- lmer(logpfus ~ time * type + (1|line), data = v_phage)
vp_m2 <- lmer(logpfus ~ time + type + (1|line), data = v_phage)
anova(vp_m1, vp_m2)
vp_m3 <- lmer(logpfus ~ time + (1|line), data = v_phage)
anova(vp_m2, vp_m3)
vp_m4 <- lmer(logpfus ~ type + (1|line), data = v_phage)
anova(vp_m2, vp_m4) #sig effect of time only
vp_m5 <- lmer(logpfus ~ 1 + (1|line), data = v_phage)
anova(vp_m3, vp_m5)
emmeans::emmeans(vp_m3, pairwise ~ time)

library(flextable)

omen <- emmeans::emmeans(op_m2, pairwise ~ time)$contrasts
omen <- data.frame(omen)

omen$contrast <- gsub("time", "", omen$contrast)

omen <- omen[,c(1,2,5,6)]

omen <- mutate(omen,
               estimate = round(estimate, digits = 3),
               t.ratio = round(t.ratio, digits = 3),
               p.value = round(p.value, digits = 3))

omen$p.value[omen$p.value == 0] <- '<0.001'

emen_tab2 <- flextable(omen) %>%
  set_header_labels(contrast = "Contrast", estimate = 'Estimate', t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  align(j = c(1:4), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#Is there a correlation between density at week two and resistance?

#den_t2 <- filter(den_dat, time == 'T2')
#den_t2 <- filter(den_t2, phage == 'Y')

#res_dat <- read.csv("data/res_tophage_280121.csv", header = T)
#anc_res <- filter(res_dat, phage == "anc")

#anc_res <- mutate(anc_res,
#res = total - susc,
#prop = res / total)

#anc_res <- filter(anc_res, time == "T2")

#write.csv(anc_res, 'res_anc_t2.csv', row.names = F)
#write.csv(den_t2, 'den_t2.csv', row.names = F)

den_t2 <- read.csv('den_t2.csv', header = T)

m1 <- lm(asin(sqrt(res)) ~ logcfus, data = den_t2)
m2 <- lm(asin(sqrt(res)) ~ 1, data = den_t2)
anova(m1,m2,test='F')
