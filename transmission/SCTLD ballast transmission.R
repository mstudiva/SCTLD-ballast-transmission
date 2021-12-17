#### packages ####

library(ggplot2)
library(ggpubr)
library(rcompanion)
library(MASS)
library(stringr)
library(survival)
library(survminer)


#### uv data import ####

uv <- read.csv("sctld uv transmission.csv", head=T)
str(uv)
head(uv)


#### uv data normality/transformation ####

# plots histogram and q-q plot for test of normality assumptions
# automatically exported as a pdf
pdf("sctld uv normality.pdf")
  par(mfrow=c(2,2))
  hist(uv$lesion)
  qqnorm(uv$lesion)
  qqline(uv$lesion)

  hist(uv$dose)
  qqnorm(uv$dose)
  qqline(uv$dose)
dev.off()
# inspect this plot, the first column shows the breakdown of data values and skew (histograms), the second column shows dispersion of residuals (Q-Q plots)
# the histograms should be roughly bell-shaped, and the Q-Q plots should be relatively straight lines

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(uv$lesion)
shapiro.test(uv$dose)
# both not normal

# BoxCox transformation finds best exponent value for data transformation
uv_bc_days <- boxcox(uv$lesion ~ uv$species_treatment)
uv_lambda_days <- uv_bc_days$x[which.max(uv_bc_days$y)]

uv_bc_dose <- boxcox(uv$dose ~ uv$species_treatment)
uv_lambda_dose <- uv_bc_dose$x[which.max(uv_bc_dose$y)]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test((uv$lesion^uv_lambda_days-1)/uv_lambda_days)
shapiro.test((uv$dose^uv_lambda_dose-1)/uv_lambda_dose)
# still not normal

pdf("sctld uv normality transformed.pdf")
par(mfrow=c(2,2))
hist((uv$lesion^uv_lambda_days-1)/uv_lambda_days)
qqnorm((uv$lesion^uv_lambda_days-1)/uv_lambda_days)
qqline((uv$lesion^uv_lambda_days-1)/uv_lambda_days)

hist((uv$dose^uv_lambda_dose-1)/uv_lambda_dose)
qqnorm((uv$dose^uv_lambda_dose-1)/uv_lambda_dose)
qqline((uv$dose^uv_lambda_dose-1)/uv_lambda_dose)
dev.off()
# both still not normal, but their histograms/Q-Q plots look better


#### uv statistical tests ####

# running ANOVAs on both raw and transformed data for final normality comparisons

# ANOVA on untransformed data
uv_anova_days <- aov((lesion) ~ species*treatment, data=uv)
summary(uv_anova_days)

uv_anova_dose <- aov((dose) ~ species*treatment, data=uv)
summary(uv_anova_dose)

# ANOVA on Box-Cox transformed data
uv_anova_days_bc <- aov(((lesion^uv_lambda_days-1)/uv_lambda_days) ~ species*treatment, data=uv)
summary(uv_anova_days_bc)

uv_anova_dose_bc <- aov(((dose^uv_lambda_dose-1)/uv_lambda_dose) ~ species*treatment, data=uv)
summary(uv_anova_dose_bc)

# Q-Q plots for raw and transformed data
pdf("sctld uv normality ANOVA.pdf")
par(mfrow=c(2,2))
qqnorm(uv_anova_days$residuals, main="days"); qqline(uv_anova_days$residuals);
qqnorm(uv_anova_dose$residuals, main="dose"); qqline(uv_anova_dose$residuals)
qqnorm(uv_anova_days_bc$residuals, main="days transformed"); qqline(uv_anova_days_bc$residuals)
qqnorm(uv_anova_dose_bc$residuals, main="dose transformed"); qqline(uv_anova_dose_bc$residuals)
dev.off()
# despite deviations from normality in Box-Cox transformed data, Q-Q plot of ANOVA residuals appears more normal than with raw data
# proceeding with ANOVA of transformed data, so copy the outputs from the second set of tests above

# Tukey post hoc tests
uv_tukey_days <- TukeyHSD(uv_anova_days_bc)
uv_tukey_days

uv_tukey_dose <- TukeyHSD(uv_anova_dose_bc)
uv_tukey_dose

# creating dataframes of the pairwise comparisons needed for plots and doing a bit of table reformatting
uv_letters_days <- data.frame(uv_tukey_days$`species:treatment`)
uv_letters_days$Var <- rownames(uv_letters_days)
names(uv_letters_days)[5] <- "comparison"
uv_letters_days$comparison = str_replace_all(uv_letters_days$comparison,":","_")
uv_letters_days$p.adj[is.na(uv_letters_days$p.adj)] <- 1
uv_letters_days

uv_letters_dose <- data.frame(uv_tukey_dose$`species:treatment`)
uv_letters_dose$Var <- rownames(uv_letters_dose)
names(uv_letters_dose)[5] <- "comparison"
uv_letters_dose$comparison = str_replace_all(uv_letters_dose$comparison,":","_")
uv_letters_dose$p.adj[is.na(uv_letters_dose$p.adj)] <- 1
uv_letters_dose

# creates compact letter display of significant pairwise differences for figure
uv_cld_days <- cldList(p.adj ~ comparison, data = uv_letters_days, threshold = 0.05)
uv_cld_days

uv_cld_dose <- cldList(p.adj ~ comparison, data = uv_letters_dose, threshold = 0.05)
uv_cld_dose


#### uv time to transmission figures ####

# keeps treatment order as imported
uv$species_treatment=factor(uv$species_treatment, levels=unique(uv$species_treatment)) 

# boxplots comparing time to transmission among species/treatments
uv_transmission_days <-
  ggboxplot(
    uv,
    x = "species_treatment",
    y = "lesion",
    color = "grey30",
    palette = c("#80cdc1","#018571","#a6611a","#dfc27d", "#018571","#a6611a","#dfc27d"),
    fill = "species_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Treatment",
           y = "Time to Transmission (d)",
           fill = 'Treatment') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=uv_cld_days, aes(x = Group, y=0, label=Letter)) 
uv_transmission_days

ggsave("sctld uv transmission days.pdf", plot= uv_transmission_days, width=6, height=4, units="in", dpi=300)

uv_transmission_dose <-
  ggboxplot(
    uv,
    x = "species_treatment",
    y = "dose",
    color = "grey30",
    palette = c("#80cdc1","#018571","#a6611a","#dfc27d","#80cdc1","#018571","#a6611a","#dfc27d"),
    fill = "species_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Treatment",
           y = "Estimated Dose (L)",
           fill = 'Treatment') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=uv_cld_dose, aes(x = Group, y=0, label=Letter)) 
uv_transmission_dose

ggsave("sctld uv transmission dose.pdf", plot= uv_transmission_dose, width=6, height=4, units="in", dpi=300)


#### uv transmission rate figure ####

uv_rate <- read.csv("sctld uv rate.csv", head = T)
uv_rate
uv_rate$treatment=factor(uv_rate$treatment, levels=unique(uv_rate$treatment)) 
uv_transrate <- ggplot(uv_rate, aes(fill=forcats::fct_rev(condition), y=rate, x=treatment)) + 
  scale_fill_manual(values=c("#80cdc1", "#80cdc1")) +
  geom_col(width = 0.5) +
  theme_bw() 
uv_transrate

ggsave("sctld uv rate.pdf", plot= uv_transrate, width=7.25, height=1.5, units="in", dpi=300)


#### uv survivorship ####

# subset dataframe by species
uv_ofav <- subset(uv, species=="Of")
uv_ofav$treatment=factor(uv_ofav$treatment, levels=c("HW","DW","UV","DC")) 

# subsetting O. faveolata dataframe to remove healthy water for comparisons between disease water and healthy water
uv_ofav2 <- subset(uv_ofav, treatment!="HW")
uv_ofav2$treatment=factor(uv_ofav2$treatment, levels=c("DW","UV","DC")) 
uv_ofav2

uv_pstr <- subset(uv, species=="Ps")
# subsetting P. strigosa dataframe to remove healthy water since no disease observed
uv_pstr <- subset(uv_pstr, treatment!="HW")
uv_pstr$treatment=factor(uv_pstr$treatment, levels=c("DW","UV","DC")) 
uv_pstr

# create survival objects for each species (using the Kaplan-Meier method)
uv_survOf <- Surv(time = uv_ofav$days, event = uv_ofav$status)
uv_survOf

uv_survOf2 <- Surv(time = uv_ofav2$days, event = uv_ofav2$status)
uv_survOf2

uv_survPs <- Surv(time = uv_pstr$days, event = uv_pstr$status)
uv_survPs

# run survival model for each species
uv_fitOf <- survfit(uv_survOf ~ treatment, data = uv_ofav)
summary(uv_fitOf)

uv_fitOf2 <- survfit(uv_survOf2 ~ treatment, data = uv_ofav2)
summary(uv_fitOf2)

uv_fitPs <- survfit(uv_survPs ~ treatment, data = uv_pstr)
summary(uv_fitPs)

# Kaplan-Meier plots for each species
fill.color_ofav<-c("#80cdc1","#a6611a","#dfc27d","#018571")

uv_survival_Ofav<-ggsurvplot(uv_fitOf, data = uv_ofav, pval = TRUE, xlab="Days", ylab="Health probability",
                          conf.int = T, risk.table=T, palette=fill.color_ofav,
                          break.time.by=5, xlim=c(0,42), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
uv_survival_Ofav

fill.color_ofav2<-c("#a6611a","#dfc27d","#018571")

uv_survival_Ofav2<-ggsurvplot(uv_fitOf2, data = uv_ofav2, pval = TRUE, xlab="Days", ylab="Health probability",
                             conf.int = T, risk.table=T, palette=fill.color_ofav2,
                             break.time.by=5, xlim=c(0,42), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
uv_survival_Ofav2

fill.color_Pstr<-c("#a6611a","#dfc27d","#018571")

uv_survival_Pstr<-ggsurvplot(uv_fitPs, data = uv_pstr, pval = TRUE, xlab="Days", ylab="Health probability",
                          conf.int = T, risk.table=T, palette=fill.color_Pstr,
                          break.time.by=5, xlim=c(0,42), risk.table.y.text = FALSE) + ggtitle("Pseudodiplora strigosa") 
uv_survival_Pstr

# hazard ratio by disease treatments
uv_hazOf <- coxph(uv_survOf ~ treatment, data = uv_ofav)
summary(uv_hazOf)

uv_hazOf2 <- coxph(uv_survOf2 ~ treatment, data = uv_ofav2)
summary(uv_hazOf2)

uv_hazPs <- coxph(uv_survPs ~ treatment, data = uv_pstr)
summary(uv_hazPs)

# plots
uv_hazard_Ofav <- ggforest(uv_hazOf, data = uv_ofav)

uv_hazard_Ofav2 <- ggforest(uv_hazOf2, data = uv_ofav2)

uv_hazard_Pstr <- ggforest(uv_hazPs, data = uv_pstr)

#### uv survivorship figures ####

uv_survival_multiplot<-ggarrange(uv_survival_Ofav$plot,
                              uv_survival_Pstr$plot, 
                              uv_survival_Ofav$table, 
                              uv_survival_Pstr$table,
                              uv_hazard_Ofav,
                              uv_hazard_Pstr,
                              heights = c(2, 0.5, 0.75),
                              ncol = 2, nrow = 3)
uv_survival_multiplot

ggsave("sctld uv survivorship.pdf", uv_survival_multiplot, width=10, height=10,dpi = 300)

uv_survival_multiplot2<-ggarrange(uv_survival_Ofav2$plot,
                                 uv_survival_Pstr$plot, 
                                 uv_survival_Ofav2$table, 
                                 uv_survival_Pstr$table,
                                 uv_hazard_Ofav2,
                                 uv_hazard_Pstr,
                                 heights = c(2, 0.5, 0.75),
                                 ncol = 2, nrow = 3)
uv_survival_multiplot2

ggsave("sctld uv survivorship_noHW.pdf", uv_survival_multiplot2, width=10, height=10,dpi = 300)


#### ballast data import ####

ballast <- read.csv("sctld ballast transmission.csv", head=T)
str(ballast)
head(ballast)


#### ballast data normality/transformation ####

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(ballast$lesion)
# not normal

# BoxCox transformation finds best exponent value for data transformation
ballast_bc_days <- boxcox(ballast$lesion ~ ballast$species_treatment)
ballast_lambda_days <- ballast_bc_days$x[which.max(ballast_bc_days$y)]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test((ballast$lesion^ballast_lambda_days-1)/ballast_lambda_days)
# still not normal

pdf("sctld ballast normality.pdf")
par(mfrow=c(2,2))
hist(ballast$lesion)
hist((ballast$lesion^ballast_lambda_days-1)/ballast_lambda_days)
qqnorm(ballast$lesion)
qqline(ballast$lesion)
qqnorm((ballast$lesion^ballast_lambda_days-1)/ballast_lambda_days)
qqline((ballast$lesion^ballast_lambda_days-1)/ballast_lambda_days)
dev.off()
# still not normal


#### ballast statistical tests ####

# running ANOVAs on both raw and transformed data for final normality comparisons

# ANOVA on untransformed data
ballast_anova_days <- aov((lesion) ~ species*treatment, data=ballast)
summary(ballast_anova_days)

# ANOVA on Box-Cox transformed data
ballast_anova_days_bc <- aov(((lesion^ballast_lambda_days-1)/ballast_lambda_days) ~ species*treatment, data=ballast)
summary(ballast_anova_days_bc)

# Q-Q plots for raw and transformed data
pdf("sctld ballast normality ANOVA.pdf")
par(mfrow=c(2,2))
qqnorm(ballast_anova_days$residuals, main="days"); qqline(ballast_anova_days$residuals);
qqnorm(ballast_anova_days_bc$residuals, main="days transformed"); qqline(ballast_anova_days_bc$residuals)
dev.off()
# not really any difference between raw and Box-Cox transformed, so proceeding with transformed analysis for consistency with UV experiment

# Tukey post hoc tests
ballast_tukey_days <- TukeyHSD(ballast_anova_days_bc)
ballast_tukey_days

# creating dataframes of the pairwise comparisons needed for plots and doing a bit of table reformatting
ballast_letters_days <- data.frame(ballast_tukey_days$`species:treatment`)
ballast_letters_days$Var <- rownames(ballast_letters_days)
names(ballast_letters_days)[5] <- "comparison"
ballast_letters_days$comparison = str_replace_all(ballast_letters_days$comparison,":","_")
ballast_letters_days$p.adj[is.na(ballast_letters_days$p.adj)] <- 1
ballast_letters_days

# creates compact letter display of significant pairwise differences for figure
ballast_cld_days <- cldList(p.adj ~ comparison, data = ballast_letters_days, threshold = 0.05)
ballast_cld_days
ballast_cld_days$Group <- c("Ps_DC","Of_DW120","Ps_DW120","Of_DW24","Ps_DW24","Of_HW120","Ps_HW120","Of_DC")

#### ballast time to transmission figures ####

# keeps treatment order as imported
ballast$species_treatment=factor(ballast$species_treatment, levels=unique(ballast$species_treatment)) 

# boxplots comparing time to transmission among species/treatments
ballast_transmission_days <-
  ggboxplot(
    ballast,
    x = "species_treatment",
    y = "lesion",
    color = "grey30",
    palette = c("#80cdc1","#018571","#dfc27d","#a6611a", "#018571","#a6611a","#dfc27d"),
    fill = "species_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.5),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Treatment",
           y = "Time to Transmission (d)",
           fill = 'Treatment') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=ballast_cld_days, aes(x = Group, y=0, label=Letter)) 
ballast_transmission_days

ggsave("sctld ballast transmission days.pdf", plot= ballast_transmission_days, width=6, height=4, units="in", dpi=300)


#### ballast transmission rate figure ####

ballast_rate <- read.csv("sctld ballast rate.csv", head = T)
ballast_rate
ballast_rate$treatment=factor(ballast_rate$treatment, levels=unique(ballast_rate$treatment)) 
ballast_transrate <- ggplot(ballast_rate, aes(fill=forcats::fct_rev(condition), y=rate, x=treatment)) + 
  scale_fill_manual(values=c("#80cdc1", "#80cdc1")) +
  geom_col(width = 0.5) +
  theme_bw() 
ballast_transrate

ggsave("sctld ballast rate.pdf", plot= ballast_transrate, width=7.25, height=1.5, units="in", dpi=300)


#### ballast survivorship ####

# subset dataframe by species
ballast_ofav <- subset(ballast, species=="Of")
ballast_ofav$treatment=factor(ballast_ofav$treatment, levels=c("HW120","DW24","DW120","DC")) 

# subsetting O. faveolata dataframe to remove healthy water for comparisons between ballast water treatments
ballast_ofav2 <- subset(ballast_ofav, treatment!="HW120")
ballast_ofav2$treatment=factor(ballast_ofav2$treatment, levels=c("DW24","DW120","DC")) 
ballast_ofav2

ballast_pstr <- subset(ballast, species=="Ps")
# subsetting P. strigosa dataframe to remove healthy water and disease water 24hr since no disease observed
ballast_pstr <- subset(ballast_pstr, treatment!="HW120")
ballast_pstr <- subset(ballast_pstr, treatment!="DW24")
ballast_pstr$treatment=factor(ballast_pstr$treatment, levels=c("DW120","DC")) 
ballast_pstr

# create survival objects for each species (using the Kaplan-Meier method)
ballast_survOf <- Surv(time = ballast_ofav$days, event = ballast_ofav$status)
ballast_survOf

ballast_survOf2 <- Surv(time = ballast_ofav2$days, event = ballast_ofav2$status)
ballast_survOf2

ballast_survPs <- Surv(time = ballast_pstr$days, event = ballast_pstr$status)
ballast_survPs

# run survival model for each species
ballast_fitOf <- survfit(ballast_survOf ~ treatment, data = ballast_ofav)
summary(ballast_fitOf)

ballast_fitOf2 <- survfit(ballast_survOf2 ~ treatment, data = ballast_ofav2)
summary(ballast_fitOf2)

ballast_fitPs <- survfit(ballast_survPs ~ treatment, data = ballast_pstr)
summary(ballast_fitPs)

# Kaplan-Meier plots for each species
fill.color_ofav<-c("#80cdc1","#dfc27d","#a6611a","#018571")

ballast_survival_Ofav<-ggsurvplot(ballast_fitOf, data = ballast_ofav, pval = TRUE, xlab="Days", ylab="Health probability",
                             conf.int = T, risk.table=T, palette=fill.color_ofav,
                             break.time.by=5, xlim=c(0,28), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
ballast_survival_Ofav

fill.color_ofav2<-c("#dfc27d","#a6611a","#018571")

ballast_survival_Ofav2<-ggsurvplot(ballast_fitOf2, data = ballast_ofav2, pval = TRUE, xlab="Days", ylab="Health probability",
                                  conf.int = T, risk.table=T, palette=fill.color_ofav2,
                                  break.time.by=5, xlim=c(0,28), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
ballast_survival_Ofav2

fill.color_Pstr<-c("#a6611a","#018571")

ballast_survival_Pstr<-ggsurvplot(ballast_fitPs, data = ballast_pstr, pval = TRUE, xlab="Days", ylab="Health probability",
                             conf.int = T, risk.table=T, palette=fill.color_Pstr,
                             break.time.by=5, xlim=c(0,28), risk.table.y.text = FALSE) + ggtitle("Pseudodiplora strigosa") 
ballast_survival_Pstr

# hazard ratio by disease treatments
ballast_hazOf <- coxph(ballast_survOf ~ treatment, data = ballast_ofav)
summary(ballast_hazOf)

ballast_hazOf2 <- coxph(ballast_survOf2 ~ treatment, data = ballast_ofav2)
summary(ballast_hazOf2)

ballast_hazPs <- coxph(ballast_survPs ~ treatment, data = ballast_pstr)
summary(ballast_hazPs)

# plots
ballast_hazard_Ofav <- ggforest(ballast_hazOf, data = ballast_ofav)

ballast_hazard_Ofav2 <- ggforest(ballast_hazOf2, data = ballast_ofav2)

ballast_hazard_Pstr <- ggforest(ballast_hazPs, data = ballast_pstr)

#### ballast survivorship figures ####

ballast_survival_multiplot<-ggarrange(ballast_survival_Ofav$plot,
                                 ballast_survival_Pstr$plot, 
                                 ballast_survival_Ofav$table, 
                                 ballast_survival_Pstr$table,
                                 ballast_hazard_Ofav,
                                 ballast_hazard_Pstr,
                                 heights = c(2, 0.5, 0.75),
                                 ncol = 2, nrow = 3)
ballast_survival_multiplot

ggsave("sctld ballast survivorship.pdf", ballast_survival_multiplot, width=10, height=10,dpi = 300)

ballast_survival_multiplot2<-ggarrange(ballast_survival_Ofav2$plot,
                                      ballast_survival_Pstr$plot, 
                                      ballast_survival_Ofav2$table, 
                                      ballast_survival_Pstr$table,
                                      ballast_hazard_Ofav2,
                                      ballast_hazard_Pstr,
                                      heights = c(2, 0.5, 0.75),
                                      ncol = 2, nrow = 3)
ballast_survival_multiplot2

ggsave("sctld ballast survivorship noHW.pdf", ballast_survival_multiplot2, width=10, height=10,dpi = 300)
