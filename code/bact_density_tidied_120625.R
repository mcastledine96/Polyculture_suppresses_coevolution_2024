### script for analysing and plotting bacterial densities

#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)
library(flextable)

#read data
den_dat <- read.csv("data/pvodensity_time_121120.csv", header = T)

#calculate population densities

den_dat <- mutate(den_dat,
                  cfus = count * dil * dil_fact * froz,
                  logcfus = log10(cfus))

#take averages

den_dat2 <- group_by(den_dat, phage, sp, time, type) %>%
  summarise(., mean_logden = mean(logcfus),
            sd = sd(logcfus),
            se = sd/sqrt(length(logcfus)))

#tidy script
den_dat <- den_dat[, c(1:5,10,12)]

#add in week 0 data / starting densities - this is for visualising the change from inoculum to week two which is displayed in the SI. Live code is for the plot in the main manuscript. 

#den_dat$treat <- interaction(den_dat$evo_line, den_dat$phage, den_dat$rep, den_dat$sp, den_dat$type)

#week_0_density <- log10(1000000/6)

#week_0_data <- data.frame(time = 0,logcfus = week_0_density,group = unique(den_dat$treat)  # Replicates and treatment average)

#week_0_data <- week_0_data %>%separate(col = group, into = c('evo_line','phage','rep', 'sp', 'type'))

#den_dat <- den_dat[,-8]

#den_dat <- rbind(week_0_data, den_dat)

#same for means
#den_dat2$treat <- interaction(den_dat2$phage, den_dat2$sp, den_dat2$type)

#week_0_density <- log10(1000000/6)

#week_0_data <- data.frame(time = 0,mean_logden = week_0_density,group = unique(den_dat2$treat),se = 0,sd = 0)

#week_0_data <- week_0_data %>%separate(col = group, into = c('phage', 'sp', 'type'))

#den_dat2 <- den_dat2[,-8]

#den_dat2 <- rbind(week_0_data, den_dat2)

#another plot of just treatments with phage
ph_bact <- filter(den_dat, phage == "Y")
ph_bact2 <- filter(den_dat2, phage == "Y")
no_phbact1 <- filter(den_dat, phage == "N")
no_phbact2 <- filter(den_dat2, phage == "N")

com_ph <- filter(den_dat, type == "comm")
com_ph2 <- filter(den_dat2, type == "comm")
mono_ph <- filter(den_dat, type == "mono")
mono_ph2 <- filter(den_dat2, type == "mono")

labels = c("0","2", "4", "6", "8")

bact_labs <- c("(i) O", "(ii) P", "(iii) V")
names(bact_labs) <- c("O", "P", "V")

bact_labs <- c("(a) Ochrobactrum sp.", "(b) Pseudomonas sp.", "(c) Variovorax sp.")

names(bact_labs) <- c("O", "P", "V")

mono_ph <- mutate(mono_ph,
                   sp2 = forcats::fct_recode(sp, `(i)~italic("Ochrobactrum")` = 'O', `(ii)~italic("Pseudomonas")` = 'P', `(iii)~italic("Variovorax")` = 'V'))

mono_ph2 <- mutate(mono_ph2,
                  sp2 = forcats::fct_recode(sp, `(i)~italic("Ochrobactrum")` = 'O', `(ii)~italic("Pseudomonas")` = 'P', `(iii)~italic("Variovorax")` = 'V'))

com_ph <- mutate(com_ph,
                  sp2 = forcats::fct_recode(sp, `(i)~italic("Ochrobactrum")` = 'O', `(ii)~italic("Pseudomonas")` = 'P', `(iii)~italic("Variovorax")` = 'V'))

com_ph2 <- mutate(com_ph2,
                 sp2 = forcats::fct_recode(sp, `(i)~italic("Ochrobactrum")` = 'O', `(ii)~italic("Pseudomonas")` = 'P', `(iii)~italic("Variovorax")` = 'V'))


###plot

bact_den1 <- ggplot() +
  geom_point(data = mono_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  facet_wrap(~sp2, labeller = label_parsed) +
  geom_line(data = mono_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  geom_point(data = mono_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 3) +
  geom_line(data = mono_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 1) +
  geom_errorbar(data = mono_ph2, aes(x = time, ymin = mean_logden - se, ymax = mean_logden+se, col = phage), width = 0.35, size = 0.8) +
  theme_bw() +
  xlab("Week") +
  ylab("Density log10(CFU/mL)") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 12)) +
  labs(title = "(a) Monoculture", col = "Phage")+
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(6.8, 8.75))

bact_den2 <- ggplot() +
  geom_point(data = com_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  facet_wrap(~sp2, labeller = label_parsed) +
  geom_line(data = com_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  geom_point(data = com_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 3) +
  geom_line(data = com_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 1) +
  geom_errorbar(data = com_ph2, aes(x = time, ymin = mean_logden - se, ymax = mean_logden+se, col = phage), width = 0.35, size = 0.8) +
  theme_bw() +
  xlab("Week") +
  ylab("Density log10(CFU/mL)") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  labs(title = "(b) Polyculture", col = "Phage")+
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4, labels = c("Absent", "Present")) +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(6.8, 8.75))

bact_den1 + bact_den2 + plot_layout(ncol = 1)

#analysis

den_dat$rando <- interaction(den_dat$evo_line, den_dat$rep)

p_den <- filter(den_dat, sp == "P")
v_den <- filter(den_dat, sp == "V")
o_den <- filter(den_dat, sp == "O")

p_mod1 <- lmer(logcfus ~ time * type * phage + (1|rando), data = p_den)
p_mod1.1 <- lmer(logcfus ~ time + type + phage + time:type + time:phage + phage:type + (1|rando), data = p_den)
anova(p_mod1, p_mod1.1) #sig 3 way interaction

#for comparisons, we are interested in comm mono / phage comparisons within a time point
emmeans::emmeans(p_mod1, pairwise ~ phage:type|time) #large differences at the start but population densities evened over time

pom <- emmeans::emmeans(p_mod1, pairwise ~ phage:type|time)$contrasts
pom <- data.frame(pom)

pom <- pom[,c(1:3, 6,7)]
pom$contrast <- gsub('comm', 'Poly.', pom$contrast)
pom$contrast <- gsub('mono', 'Mono.', pom$contrast)
pom$contrast <- gsub('N', 'Phage absent', pom$contrast)
pom$contrast <- gsub('Y', 'Phage present', pom$contrast)
pom$time <- gsub('T', '', pom$time)

pom <- mutate(pom,
              estimate = round(estimate, digits = 3),
              t.ratio = round(t.ratio, digits = 3),
              p.value = round(p.value, digits = 3))

pom$p.value[pom$p.value == 0] <- '<0.001'

pom_tab <- flextable(pom) %>%
  set_header_labels(contrast = "Contrast", time = 'Week', estimate = 'Estimate', t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  align(j = c(1:5), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(6,12,18), border = officer::fp_border(color="black")) 


#same for O
o_mod1 <- lmer(logcfus ~ time * type * phage + (1|rando), data = o_den)
o_mod1.1 <- lmer(logcfus ~ time + type + phage + time:type + time:phage + phage:type + (1|rando), data = o_den)
anova(o_mod1, o_mod1.1) #no 3 way interaction
o_mod2 <- lmer(logcfus ~ time + type + phage + time:phage + type:phage + (1|rando), data = o_den) #removing time:type interaction 
anova(o_mod1.1, o_mod2) #non-sig time:type interaction
o_mod3 <- lmer(logcfus ~ time + type + phage + time:type + type:phage + (1|rando), data = o_den) #removing time:phage interaction
anova(o_mod1.1, o_mod3) #sig time:phage interaction
o_mod4 <- lmer(logcfus ~ time + type + phage + time:phage + time:type + (1|rando), data = o_den) #remove type:phage interaction
anova(o_mod1.1, o_mod4) #sig type:phage interaction
#drop time:type interaction and continue
o_mod5 <- lmer(logcfus ~ time + type + phage + time:phage + (1|rando), data = o_den) #drop phage:type interaction
anova(o_mod2, o_mod5) #sig phage:type interaction
o_mod6 <- lmer(logcfus ~ time + type + phage + type:phage + (1|rando), data = o_den) #drop phage:time interaction
anova(o_mod2, o_mod6) #sig time:phage interaction. 
#Model cannot be simplified further

emmeans::emmeans(o_mod2, pairwise ~ phage:type|time)
emmeans::emmeans(o_mod2, pairwise ~ phage|type)

bom <- emmeans::emmeans(o_mod2, pairwise ~ phage:type|time)$contrasts
bom <- data.frame(bom)

bom <- bom[,c(1:3, 6,7)]
bom$contrast <- gsub('comm', 'Poly.', bom$contrast)
bom$contrast <- gsub('mono', 'Mono.', bom$contrast)
bom$contrast <- gsub('N', 'Phage absent', bom$contrast)
bom$contrast <- gsub('Y', 'Phage present', bom$contrast)
bom$time <- gsub('T', '', bom$time)

bom <- mutate(bom,
               estimate = round(estimate, digits = 3),
               t.ratio = round(t.ratio, digits = 3),
               p.value = round(p.value, digits = 3))

bom$p.value[bom$p.value == 0] <- '<0.001'

bom_tab <- flextable(bom) %>%
  set_header_labels(contrast = "Contrast", time = 'Week', estimate = 'Estimate', t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  align(j = c(1:5), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(6,12,18), border = officer::fp_border(color="black")) 

#same for V
v_mod1 <- lmer(logcfus ~ time * type * phage + (1|rando), data = v_den)
v_mod1.1 <- lmer(logcfus ~ time + type + phage + time:type + time:phage + phage:type + (1|rando), data = v_den)
anova(v_mod1, v_mod1.1) #non-sig 3 way interaction
v_mod2 <- lmer(logcfus ~ time + type + phage + time:phage + type:phage + (1|rando), data = v_den) #removing time:type interaction 
anova(v_mod1, v_mod2) #sig time:type interaction
v_mod3 <- lmer(logcfus ~ time + type + phage + time:type + type:phage + (1|rando), data = v_den) #removing time:phage interaction
anova(v_mod1, v_mod3) #non-sig time:phage interaction
v_mod4 <- lmer(logcfus ~ time + type + phage + time:phage + time:type + (1|rando), data = v_den) #remove type:phage interaction
anova(v_mod1, v_mod4) #non-sig type:phage interaction. This one is less sig so drop first
v_mod5 <- lmer(logcfus ~ time + type + phage + time:type + (1|rando), data = v_den) #drop time:phage
anova(v_mod4, v_mod5) #still time:phage non-sig
v_mod6 <- lmer(logcfus ~ time + type + phage + time:phage + (1|rando), data = v_den) #drop time:type
anova(v_mod4, v_mod6) #time:type significant
v_mod6.1 <- lmer(logcfus ~ time + type + phage + (1|rando), data = v_den) #drop time:type
anova(v_mod5, v_mod6.1)
v_mod7 <- lmer(logcfus ~ time + type + time:type + (1|rando), data = v_den)
anova(v_mod5, v_mod7) #phage as an independent effect significant but not as interacting

emmeans::emmeans(v_mod5, pairwise ~ type|time) #difference in poly to mono decreases over time
emmeans::emmeans(v_mod5, pairwise ~ phage) #presence of phage overall increases V density

library(flextable)
mens <- data.frame(emmeans::emmeans(v_mod5, pairwise ~ type|time)$contrasts)
mens <- mens[,c(1:3,6,7)]

mens$contrast <- "Poly. - Mono."

mens <- mutate(mens,
               estimate = round(estimate, digits = 3),
               t.ratio = round(t.ratio, digits = 3),
               p.value = round(p.value, digits = 3))

mens$p.value[mens$p.value == 0] <- "<0.001"
mens$time <- substring(mens$time, 2)

men_flex <- flextable(mens) %>%
  set_header_labels(time = "Week", contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 