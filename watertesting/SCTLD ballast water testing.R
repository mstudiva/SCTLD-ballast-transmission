#### packages ####

library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(stringr)
library(MASS)
library(ape)
library(reshape)
library(ggplot2)
library(Rmisc)


#### collimated beam data import ####

beam <- read.csv(file = "sctld collimated beam.csv", head=T)
str(beam)
head(beam)

beam_mean <- summarySE(beam, measurevar="cfu", groupvars=c("treatment","dose"))
beam_mean

beam_plot <-  ggplot(beam_mean, aes(x=dose, y=cfu, color=treatment)) +
  geom_errorbar(aes(ymin=cfu-se, ymax=cfu+se), width=10, position=position_dodge(0.9), size = 0.75) +
  geom_line(position=position_dodge(0.9), size = 0.5) +
  geom_point(position=position_dodge(0.9), size = 3) +
  scale_color_manual(values=c("black","blue")) +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  scale_x_continuous(breaks = seq(0, 300, by = 50)) +
  theme_classic() +
  xlab("Estimated dose (mWs-1 cm-2)") +
  ylab("CFU mL-1")
beam_plot

ggsave("sctld collimatedbeam plot.pdf", plot= beam_plot, width=6, height=4, units="in", dpi=300)


#### UV cell counts data import ####

uv_counts <- read.csv(file = "sctld uv cellcounts.csv", head = T)
str(uv_counts)
head(uv_counts)

# histograms
pdf(file = "sctld uv cellcounts hist.pdf", width = 6, height = 10)
par(mfrow=c(4,2))
hist(uv_counts$live, breaks = 50)
hist(uv_counts$flag, breaks = 50)
hist(uv_counts$dino, breaks = 50)
hist(uv_counts$cili, breaks = 50)
hist(uv_counts$anne, breaks = 50)
hist(uv_counts$pdia, breaks = 50)
hist(uv_counts$cdia, breaks = 50)
hist(uv_counts$other, breaks = 50)
dev.off()
# mostly right-skewed, conducting square-root transformation

uv_counts$live_sqrt <- sqrt(uv_counts$live)
uv_counts$flag_sqrt <- sqrt(uv_counts$flag)
uv_counts$dino_sqrt <- sqrt(uv_counts$dino)
uv_counts$cili_sqrt <- sqrt(uv_counts$cili)
uv_counts$anne_sqrt <- sqrt(uv_counts$anne)
uv_counts$pdia_sqrt <- sqrt(uv_counts$pdia)
uv_counts$cdia_sqrt <- sqrt(uv_counts$cdia)
uv_counts$other_sqrt <- sqrt(uv_counts$other)

# histograms of transformed variables
pdf(file = "sctld uv cellcounts hist transformed.pdf", width = 6, height = 10)
par(mfrow=c(4,2))
hist(uv_counts$live_sqrt, breaks = 50)
hist(uv_counts$flag_sqrt, breaks = 50)
hist(uv_counts$dino_sqrt, breaks = 50)
hist(uv_counts$cili_sqrt, breaks = 50)
hist(uv_counts$anne_sqrt, breaks = 50)
hist(uv_counts$pdia_sqrt, breaks = 50)
hist(uv_counts$cdia_sqrt, breaks = 50)
hist(uv_counts$other_sqrt, breaks = 50)
dev.off()


#### UV cell counts PERMANOVA ####

# Setting seed allows randomized processes to be repeated later
set.seed(1337)

# creates a dissimilarity matrix of square-root transformed variables, based on Euclidean distance
uv_counts_dist <- vegdist(uv_counts[c(14:20)], method="bray", na.rm= TRUE)

# testing for homogeneity of variance among locations
uv_counts_disp <- betadisper(uv_counts_dist, group=uv_counts$treatment)
permutest(uv_counts_disp, bias.adjust = TRUE, perm = 9999)
# insignificant test value indicates homogeneous variance

# running the PERMANOVA
# time is accounted as a random factor using 'strata'
uv_counts_perm <- adonis(uv_counts_dist ~ time*treatment, data=uv_counts, permutations = 9999, parallel = getOption("mc.cores"))
uv_counts_perm
write.csv(uv_counts_perm$aov.tab, file = "sctld uv cellcounts PERMANOVA.csv")
# Significant time effect, insignificant treatment effect

# pairwise PERMANOVA
uv_counts_perm_time <- pairwise.adonis2(uv_counts_dist ~ time, data = uv_counts)
uv_counts_perm_time

# exporting pairwise results
uv_counts_perm_out <- bind_rows(uv_counts_perm_time$t0_vs_t1, uv_counts_perm_time$t0_vs_t2, uv_counts_perm_time$t0_vs_t3, uv_counts_perm_time$t1_vs_t2, uv_counts_perm_time$t1_vs_t3, uv_counts_perm_time$t2_vs_t3,  .id = "Comparison")
uv_counts_perm_out
write.csv(uv_counts_perm_out, file = "sctld uv cellcounts PERMANOVA pair.csv")


#### UV total cell counts ####

shapiro.test(uv_counts$live)
shapiro.test(uv_counts$live_sqrt)

uv_counts_anova <- aov(live_sqrt ~ time*treatment, data=uv_counts)
summary(uv_counts_anova)

uv_counts_tukey <- TukeyHSD(uv_counts_anova, which = "time")
uv_counts_tukey
write.csv(uv_counts_tukey$time, file = "sctld uv cellcounts Tukey.csv")


#### UV PAM data import ####

uv_pam <- read.csv(file = "sctld uv PAM.csv", head = T)
str(uv_pam)
head(uv_pam)

# histograms
pdf(file = "sctld uv PAM hist.pdf", width = 6, height = 6)
par(mfrow=c(2,2))
hist(uv_pam$Fo, breaks = 50)
hist(uv_pam$Fv, breaks = 50)
qqnorm(uv_pam$Fo)
qqline(uv_pam$Fo)
qqnorm(uv_pam$Fv)
qqline(uv_pam$Fv)
dev.off()


#### UV PAM ANOVAs ####

shapiro.test(uv_pam$Fo)
shapiro.test(uv_pam$Fv)

# BoxCox transformation finds best exponent value for data transformation
uv_pam_Fo_bc <- boxcox(uv_pam$Fo ~ uv_pam$time)
uv_pam_Fo_lambda <- uv_pam_Fo_bc$x[which.max(uv_pam_Fo_bc$y)]

uv_pam_Fv_bc <- boxcox(uv_pam$Fv+1 ~ uv_pam$time)
uv_pam_Fv_lambda <- uv_pam_Fv_bc$x[which.max(uv_pam_Fv_bc$y)]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test((uv_pam$Fo^uv_pam_Fo_lambda-1)/uv_pam_Fo_lambda)
shapiro.test((uv_pam$Fv^uv_pam_Fv_lambda-1)/uv_pam_Fv_lambda)

pdf("sctld uv PAM hist transformed.pdf")
par(mfrow=c(2,2))
hist((uv_pam$Fo^uv_pam_Fo_lambda-1)/uv_pam_Fo_lambda)
qqnorm((uv_pam$Fo^uv_pam_Fo_lambda-1)/uv_pam_Fo_lambda)
qqline((uv_pam$Fo^uv_pam_Fo_lambda-1)/uv_pam_Fo_lambda)

hist((uv_pam$Fv^uv_pam_Fv_lambda-1)/uv_pam_Fv_lambda)
qqnorm((uv_pam$Fv^uv_pam_Fv_lambda-1)/uv_pam_Fv_lambda)
qqline((uv_pam$Fv^uv_pam_Fv_lambda-1)/uv_pam_Fv_lambda)
dev.off()

# ANOVA on untransformed data
uv_pam_anova_Fo <- aov(Fo ~ time*treatment, data=uv_pam)
summary(uv_pam_anova_Fo)

uv_pam_anova_Fv <- aov(Fv ~ time*treatment, data=uv_pam)
summary(uv_pam_anova_Fv)

# ANOVA on Box-Cox transformed data
uv_pam_anova_Fo_bc <- aov(((Fo^uv_pam_Fo_lambda-1)/uv_pam_Fo_lambda) ~ time*treatment, data=uv_pam)
summary(uv_pam_anova_Fo_bc)

uv_pam_anova_Fv_bc <- aov(((Fv^uv_pam_Fv_lambda-1)/uv_pam_Fv_lambda) ~ time*treatment, data=uv_pam)
summary(uv_pam_anova_Fv_bc)

# Q-Q plots for raw and transformed data
pdf("sctld uv PAM normality ANOVA.pdf")
par(mfrow=c(2,2))
qqnorm(uv_pam_anova_Fo$residuals, main="Fo"); qqline(uv_pam_anova_Fo$residuals);
qqnorm(uv_pam_anova_Fv$residuals, main="Fv"); qqline(uv_pam_anova_Fv$residuals)
qqnorm(uv_pam_anova_Fo_bc$residuals, main="Fo transformed"); qqline(uv_pam_anova_Fo_bc$residuals)
qqnorm(uv_pam_anova_Fv_bc$residuals, main="Fv transformed"); qqline(uv_pam_anova_Fv_bc$residuals)
dev.off()

uv_pam_tukey_Fo <- TukeyHSD(uv_pam_anova_Fo_bc)
uv_pam_tukey_Fo
write.csv(uv_pam_tukey_Fo$`time:treatment`, file = "sctld uv PAM Fo Tukey.csv")

uv_pam_tukey_Fv <- TukeyHSD(uv_pam_anova_Fv_bc)
uv_pam_tukey_Fv
write.csv(uv_pam_tukey_Fv$`time:treatment`, file = "sctld uv PAM Fv Tukey.csv")


#### UV figures ####

library(rcompanion)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(gridExtra)

# modifying datasets with blank data rows of time points so cell count and PAM dataframes line up in multiplot
uv_counts <- uv_counts %>% add_row(time = "t4", treatment = "DW")
uv_counts <- uv_counts %>% add_row(time = "t4", treatment = "UV")
uv_counts <- uv_counts %>% add_row(time = "t4", treatment = "HW")
uv_counts

uv_pam <- uv_pam %>% add_row(time = "t0", treatment = "HW",  .before = T)
uv_pam <- uv_pam %>% add_row(time = "t0", treatment = "UV", .before = T)
uv_pam <- uv_pam %>% add_row(time = "t0", treatment = "DW", .before = T)
uv_pam

# creates a new column of combined time/treatment
uv_counts$time_treatment <- paste(uv_counts$time,uv_counts$treatment, sep = "_")
uv_pam$time_treatment <- paste(uv_pam$time,uv_pam$treatment, sep = "_")

# PCoA of cell counts by taxon 
uv_counts_pcoa <- pcoa(uv_counts_dist)
scores <- uv_counts_pcoa$values
vectors <- uv_counts_pcoa$vectors

pdf(file="sctld uv cellcounts PCoA.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(vectors[,1], vectors[,2],col=c("#d8b365","#f6e8c3","#5ab4ac")[as.numeric(as.factor(uv_counts$treatment))],pch=c(15,1,19,17)[as.numeric((as.factor(uv_counts$time)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
# ordispider(vectors, uv_counts$treatment, label=F, col=c("#d8b365","#f6e8c3","#5ab4ac"))
legend("topright", legend=c("DW", "UV","HW"), fill = c("#d8b365","#f6e8c3","#5ab4ac"), bty="n")
legend("bottomright", legend=c("t0","t1","t2","t3"), pch=c(15,1,19,17), bty="n")
plot(vectors[,1], vectors[,2],col=c("black","grey75","grey50","grey25")[as.numeric(as.factor(uv_counts$time))],pch=c(15,18,25,3)[as.numeric(as.factor(uv_counts$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
# ordispider(vectors, uv_counts$condition, label=F, col=c("black","grey75","grey50","grey25"))
legend("topright", legend=c("t0","t1","t2","t3"), fill = c("black","grey75","grey50","grey25"), bty="n")
legend("bottomright", legend=c("DW", "UV","HW"), pch=c(15,18,25,3), bty="n")
dev.off()

# stacked column chart

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
uv_counts_perm_letters <- data.frame(cbind(uv_counts_perm_out$Comparison,uv_counts_perm_out$'Pr(>F)'))
uv_counts_perm_letters <- na.omit(uv_counts_perm_letters)
uv_counts_perm_letters <- dplyr::rename(uv_counts_perm_letters, comparison = X1, p.adj = X2)
uv_counts_perm_letters$comparison = c("t1-t0","t2-t0","t3-t0","t2-t1","t3-t1","t3-t2")
uv_counts_perm_letters$p.adj <- as.numeric(paste(uv_counts_perm_letters$p.adj))
uv_counts_perm_letters

# creates compact letter display of significant pairwise differences for figure
uv_counts_perm_cld <- cldList(p.adj ~ comparison, data = uv_counts_perm_letters, threshold = 0.05, remove.zero = FALSE)
uv_counts_perm_cld=uv_counts_perm_cld[order(uv_counts_perm_cld$Group),] 
uv_counts_perm_cld$time <- uv_counts_perm_cld$Group
uv_counts_perm_cld <- uv_counts_perm_cld %>% add_row(Group = "t4", time = "t4", Letter = "no data")
uv_counts_perm_cld <- uv_counts_perm_cld %>% add_column(taxa="other")
uv_counts_perm_cld

# transposing and reformatting dataframe to make abundance a single column
uv_counts_percent <- dplyr::select(uv_counts, 2:3,6:12,21)
uv_counts_percent <- melt(uv_counts_percent, id = c("time","treatment","time_treatment"))
uv_counts_percent$time_treatment=factor(uv_counts_percent$time_treatment, levels=unique(uv_counts_percent$time_treatment)) 
uv_counts_percent <- cast(uv_counts_percent, time_treatment~variable, mean)
uv_counts_percent <- melt(uv_counts_percent, id = time_treatment)
uv_counts_percent <- dplyr::rename(uv_counts_percent, taxa = variable, abundance = value)
uv_counts_percent$time_treatment2 <- uv_counts_percent$time_treatment
uv_counts_percent <- uv_counts_percent %>% separate(time_treatment2, c('time', 'treatment'))
uv_counts_percent

# creates dummy dataframes for custom text labels of stats outputs
uv_counts_perm_stats <- data.frame(time_treatment = "t4_UV",abundance = 0,lab = "time: F3,35=28.091, p<0.001",time = "t4",taxa="other")
uv_counts_perm_stats2 <- data.frame(time_treatment = "t4_UV",abundance = 0,lab = "treat: F2,35=0.573, p=0.71",time = "t4",taxa="other")

uv_counts_percent$taxa=factor(uv_counts_percent$taxa, levels=unique(uv_counts_percent$taxa)) 

# percent stacked barplot
uv_count_abundance <- ggplot(uv_counts_percent, aes(fill=taxa, y=abundance, x=time_treatment)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Time",
       y = "Percent abundance",
       fill = 'Taxon') + 
  scale_fill_manual(values=c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#56B4E9", "#CC79A7","#999999")) +
  facet_grid(. ~ time, scales="free", switch = c("both")) +
  geom_text(data = uv_counts_perm_stats,label = uv_counts_perm_stats$lab, aes(x=3, y=0.95)) +
  geom_text(data = uv_counts_perm_stats2,label = uv_counts_perm_stats2$lab, aes(x=2.85, y=0.85)) +
  geom_text(data = uv_counts_perm_cld, aes(x=2, y=-0.025, label=Letter)) +
  theme_classic() 
uv_count_abundance 

ggsave("sctld uv counts abundance plot.pdf", plot= uv_count_abundance, width=10, height=3, units="in", dpi=300)

# cell counts
# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
uv_counts_letters <- data.frame(uv_counts_tukey$time)
uv_counts_letters$comparison <- rownames(uv_counts_letters)
uv_counts_letters$p.adj[is.na(uv_counts_letters$p.adj)] <- 1
uv_counts_letters

# creates compact letter display of significant pairwise differences for figure
uv_counts_cld <- cldList(p.adj ~ comparison, data = uv_counts_letters, threshold = 0.05, remove.zero = F)
uv_counts_cld=uv_counts_cld[order(uv_counts_cld$Group),] 
uv_counts_cld$time <- uv_counts_cld$Group
uv_counts_cld <- uv_counts_cld %>% add_row(Group = "t4", time = "t4", Letter = "no data")
uv_counts_cld

# creates dummy dataframes for custom text labels of stats outputs
uv_counts_stats <- data.frame(time_treatment = "t4_UV",lab = "time: F3,35=69.614, p<0.001",time = "t4")
uv_counts_stats2 <- data.frame(time_treatment = "t4_UV",lab = "treat: F2,35=0.164, p=0.85",time = "t4")

# keeps order as imported
uv_counts$time_treatment=factor(uv_counts$time_treatment, levels=unique(uv_counts$time_treatment)) 

# boxplot of live cell counts among time/treatments
uv_counts_plot <-
  ggboxplot(
    uv_counts,
    x = "time_treatment",
    y = "live",
    color = "grey30",
    palette = c("#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Time",
           y = "Live cell counts (L-1)",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right") +
  geom_text(data = uv_counts_cld, aes(x=2, y=600, label=Letter)) +
  geom_text(data = uv_counts_stats,label = uv_counts_stats$lab, aes(x=3, y=3900)) +
  geom_text(data = uv_counts_stats2,label = uv_counts_stats2$lab, aes(x=2.85, y=3500)) 
uv_counts_plot

ggsave("sctld uv counts plot.pdf", plot= uv_counts_plot, width=10, height=3, units="in", dpi=300)

# PAM Fo
uv_pam_letters_Fo <- data.frame(uv_pam_tukey_Fo$'time:treatment')
uv_pam_letters_Fo$comparison <- rownames(uv_pam_letters_Fo)
uv_pam_letters_Fo$comparison = str_replace_all(uv_pam_letters_Fo$comparison,":","_")
uv_pam_letters_Fo$p.adj[is.na(uv_pam_letters_Fo$p.adj)] <- 1
uv_pam_letters_Fo

# creates compact letter display of significant pairwise differences for figure
uv_pam_cld_Fo <- cldList(p.adj ~ comparison, data = uv_pam_letters_Fo, threshold = 0.05)
uv_pam_cld_Fo$Group  <- factor(uv_pam_cld_Fo$Group , levels = c("t1_DW", "t1_UV","t1_HW","t2_DW", "t2_UV","t2_HW","t3_DW", "t3_UV","t3_HW","t4_DW", "t4_UV","t4_HW"))
uv_pam_cld_Fo=uv_pam_cld_Fo[order(uv_pam_cld_Fo$Group),]
uv_pam_cld_Fo$time_treatment <- uv_pam_cld_Fo$Group
uv_pam_cld_Fo$time <- c("t1","t1","t1","t2","t2","t2","t3","t3","t3","t4","t4","t4")
uv_pam_cld_Fo <- uv_pam_cld_Fo %>% add_row(Group = "t0_HW", time_treatment = "t0_HW", time = "t0", .before = T)
uv_pam_cld_Fo <- uv_pam_cld_Fo %>% add_row(Group = "t0_UV", time_treatment = "t0_UV", time = "t0", Letter = "no data", .before = T)
uv_pam_cld_Fo <- uv_pam_cld_Fo %>% add_row(Group = "t0_DW", time_treatment = "t0_DW", time = "t0", .before = T)
uv_pam_cld_Fo

# creates dummy dataframes for custom text labels of stats outputs
uv_pam_Fo_stats <- data.frame(time_treatment = "t0_UV",lab = "time: F3,179=207.5, p<0.001",time = "t0")
uv_pam_Fo_stats2 <- data.frame(time_treatment = "t0_UV",lab = "treat: F2,179=200.9, p<0.001",time = "t0")
uv_pam_Fo_stats3 <- data.frame(time_treatment = "t0_UV",lab = "int: F6,179=261.4, p<0.001",time = "t0")

uv_pam$time_treatment=factor(uv_pam$time_treatment, levels=unique(uv_pam$time_treatment)) 

# boxplot of Fo among time/treatments
uv_pam_Fo_plot <-
  ggboxplot(
    uv_pam,
    x = "time_treatment",
    y = "Fo",
    color = "grey30",
    palette = c("#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Time",
           y = "Fo",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right") +
  geom_text(data = uv_pam_cld_Fo, aes(x=time_treatment, y=50, label=Letter)) +
  geom_text(data = uv_pam_Fo_stats,label = uv_pam_Fo_stats$lab, aes(x=3, y=400)) +
  geom_text(data = uv_pam_Fo_stats2,label = uv_pam_Fo_stats2$lab, aes(x=3, y=370)) +
  geom_text(data = uv_pam_Fo_stats3,label = uv_pam_Fo_stats3$lab, aes(x=2.85, y=340)) 
uv_pam_Fo_plot

ggsave("sctld uv PAM Fo plot.pdf", plot= uv_pam_Fo_plot, width=10, height=3, units="in", dpi=300)

# PAM Fv
uv_pam_letters_Fv <- data.frame(uv_pam_tukey_Fv$'time:treatment')
uv_pam_letters_Fv$comparison <- rownames(uv_pam_letters_Fv)
uv_pam_letters_Fv$comparison = str_replace_all(uv_pam_letters_Fv$comparison,":","_")
uv_pam_letters_Fv$p.adj[is.na(uv_pam_letters_Fv$p.adj)] <- 1
uv_pam_letters_Fv

# creates compact letter display of significant pairwise differences for figure
uv_pam_cld_Fv <- cldList(p.adj ~ comparison, data = uv_pam_letters_Fv, threshold = 0.05)
uv_pam_cld_Fv$Group  <- factor(uv_pam_cld_Fv$Group , levels = c("t1_DW", "t1_UV","t1_HW","t2_DW", "t2_UV","t2_HW","t3_DW", "t3_UV","t3_HW","t4_DW", "t4_UV","t4_HW"))
uv_pam_cld_Fv=uv_pam_cld_Fv[order(uv_pam_cld_Fv$Group),]
uv_pam_cld_Fv$time_treatment <- uv_pam_cld_Fv$Group
uv_pam_cld_Fv$time <- c("t1","t1","t1","t2","t2","t2","t3","t3","t3","t4","t4","t4")
uv_pam_cld_Fv <- uv_pam_cld_Fv %>% add_row(Group = "t0_HW", time_treatment = "t0_HW", time = "t0", .before = T)
uv_pam_cld_Fv <- uv_pam_cld_Fv %>% add_row(Group = "t0_UV", time_treatment = "t0_UV", time = "t0", Letter = "no data", .before = T)
uv_pam_cld_Fv <- uv_pam_cld_Fv %>% add_row(Group = "t0_DW", time_treatment = "t0_DW", time = "t0", .before = T)
uv_pam_cld_Fv

# creates dummy dataframes for custom text labels of stats outputs
uv_pam_Fv_stats <- data.frame(time_treatment = "t0_UV",lab = "time: F3,179=618.17, p<0.001",time = "t0")
uv_pam_Fv_stats2 <- data.frame(time_treatment = "t0_UV",lab = "treat: F2,179=298.2, p<0.001",time = "t0")
uv_pam_Fv_stats3 <- data.frame(time_treatment = "t0_UV",lab = "int: F6,179=40.79, p<0.001",time = "t0")

# boxplot of Fv among time/treatments
uv_pam_Fv_plot <-
  ggboxplot(
    uv_pam,
    x = "time_treatment",
    y = "Fv",
    color = "grey30",
    palette = c("#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac","#d8b365","#f6e8c3", "#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Time",
           y = "Fv",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right") +
  geom_text(data = uv_pam_cld_Fv, aes(x=time_treatment, y=-20, label=Letter)) +
  geom_hline(yintercept=400, linetype="dashed", color = "grey") +
  geom_text(data = uv_pam_Fv_stats,label = uv_pam_Fv_stats$lab, aes(x=3, y=690)) +
  geom_text(data = uv_pam_Fv_stats2,label = uv_pam_Fv_stats2$lab, aes(x=2.95, y=630)) +
  geom_text(data = uv_pam_Fv_stats3,label = uv_pam_Fv_stats3$lab, aes(x=2.8, y=570)) 
uv_pam_Fv_plot

ggsave("sctld uv PAM Fv plot.pdf", plot= uv_pam_Fv_plot, width=10, height=3, units="in", dpi=300)

# cell counts/PAM multiplot
multiplot=grid.arrange(uv_count_abundance, uv_counts_plot, uv_pam_Fo_plot, uv_pam_Fv_plot, ncol=1, nrow=4, widths=c(10), heights=c(3,3,3,3))

ggsave("sctld uv watertesting panel.pdf", plot= multiplot, width=10, height=12, units="in", dpi=300)


#### ballast cell counts data import ####

ballast_counts <- read.csv(file = "sctld ballast cellcounts.csv", head = T)
str(ballast_counts)
head(ballast_counts)

# histograms
pdf(file = "sctld ballast cellcounts hist.pdf", width = 6, height = 10)
par(mfrow=c(4,2))
hist(ballast_counts$live, breaks = 50)
hist(ballast_counts$flag, breaks = 50)
hist(ballast_counts$dino, breaks = 50)
hist(ballast_counts$cili, breaks = 50)
hist(ballast_counts$anne, breaks = 50)
hist(ballast_counts$pdia, breaks = 50)
hist(ballast_counts$cdia, breaks = 50)
hist(ballast_counts$other, breaks = 50)
dev.off()

# mostly right-skewed, conducting square-root transformation
ballast_counts$live_sqrt <- sqrt(ballast_counts$live)
ballast_counts$flag_sqrt <- sqrt(ballast_counts$flag)
ballast_counts$dino_sqrt <- sqrt(ballast_counts$dino)
ballast_counts$cili_sqrt <- sqrt(ballast_counts$cili)
ballast_counts$anne_sqrt <- sqrt(ballast_counts$anne)
ballast_counts$pdia_sqrt <- sqrt(ballast_counts$pdia)
ballast_counts$cdia_sqrt <- sqrt(ballast_counts$cdia)
ballast_counts$other_sqrt <- sqrt(ballast_counts$other)

# histograms of transformed variables
pdf(file = "sctld ballast cellcounts hist transformed.pdf", width = 6, height = 10)
par(mfrow=c(4,2))
hist(ballast_counts$live_sqrt, breaks = 50)
hist(ballast_counts$flag_sqrt, breaks = 50)
hist(ballast_counts$dino_sqrt, breaks = 50)
hist(ballast_counts$cili_sqrt, breaks = 50)
hist(ballast_counts$anne_sqrt, breaks = 50)
hist(ballast_counts$pdia_sqrt, breaks = 50)
hist(ballast_counts$cdia_sqrt, breaks = 50)
hist(ballast_counts$other_sqrt, breaks = 50)
dev.off()


#### ballast cell counts PERMANOVA ####

# Setting seed allows randomized processes to be repeated later
set.seed(1337)

# creates a dissimilarity matrix of square-root transformed variables, based on Euclidean distance
ballast_counts_dist <- vegdist(ballast_counts[c(14:20)], method="bray", na.rm= TRUE)

# testing for homogeneity of variance among locations
ballast_counts_disp <- betadisper(ballast_counts_dist, group=ballast_counts$treatment)
permutest(ballast_counts_disp, bias.adjust = TRUE, perm = 9999)
# insignificant test value indicates homogeneous variance

# running the PERMANOVA
# time is accounted as a random factor using 'strata'
ballast_counts_perm <- adonis(ballast_counts_dist ~ time*treatment, data=ballast_counts, permutations = 9999, parallel = getOption("mc.cores"))
ballast_counts_perm

write.csv(ballast_counts_perm$aov.tab, file = "sctld ballast cellcounts PERMANOVA.csv")
# Significant time effect, insignificant treatment effect

# pairwise PERMANOVA
ballast_counts$time_treatment <- paste(ballast_counts$time,ballast_counts$treatment, sep = "_")

ballast_counts_perm_int <- pairwise.adonis2(ballast_counts_dist ~ time_treatment, data = ballast_counts)
ballast_counts_perm_int

# exporting pairwise results
ballast_counts_perm_out <- bind_rows(ballast_counts_perm_int$t0_HW120_vs_t0_DW24, ballast_counts_perm_int$t0_HW120_vs_t0_DW120, ballast_counts_perm_int$t0_HW120_vs_t1_HW120, ballast_counts_perm_int$t0_HW120_vs_t1_DW24, ballast_counts_perm_int$t0_HW120_vs_t1_DW120, ballast_counts_perm_int$t0_DW24_vs_t0_DW120, ballast_counts_perm_int$t0_DW24_vs_t1_HW120, ballast_counts_perm_int$t0_DW24_vs_t1_DW24, ballast_counts_perm_int$t0_DW24_vs_t1_DW120, ballast_counts_perm_int$t0_DW120_vs_t1_HW120, ballast_counts_perm_int$t0_DW120_vs_t1_DW24, ballast_counts_perm_int$t0_DW120_vs_t1_DW120, ballast_counts_perm_int$t1_HW120_vs_t1_DW24, ballast_counts_perm_int$t1_HW120_vs_t1_DW120, ballast_counts_perm_int$t1_DW24_vs_t1_DW120,  .id = "Comparison")
ballast_counts_perm_out
write.csv(ballast_counts_perm_out, file = "sctld ballast cellcounts PERMANOVA pair.csv")


#### ballast total cell counts ####

shapiro.test(ballast_counts$live)
shapiro.test(ballast_counts$live_sqrt)

ballast_counts_anova <- aov(live_sqrt ~ time*treatment, data=ballast_counts)
summary(ballast_counts_anova)

ballast_counts_tukey <- TukeyHSD(ballast_counts_anova)
ballast_counts_tukey
write.csv(ballast_counts_tukey$`time:treatment`, file = "sctld ballast cellcounts Tukey.csv")


#### ballast PAM data import ####

ballast_pam <- read.csv(file = "sctld ballast PAM.csv", head = T)
str(ballast_pam)
head(ballast_pam)

# histograms
pdf(file = "sctld ballast PAM hist.pdf", width = 6, height = 6)
par(mfrow=c(2,2))
hist(ballast_pam$Fo, breaks = 50)
hist(ballast_pam$Fv, breaks = 50)
qqnorm(ballast_pam$Fo)
qqline(ballast_pam$Fo)
qqnorm(ballast_pam$Fv)
qqline(ballast_pam$Fv)
dev.off()


#### ballast PAM ANOVAs ####

shapiro.test(ballast_pam$Fo)
shapiro.test(ballast_pam$Fv)

# BoxCox transformation finds best exponent value for data transformation
ballast_pam_Fo_bc <- boxcox(ballast_pam$Fo ~ ballast_pam$time)
ballast_pam_Fo_lambda <- ballast_pam_Fo_bc$x[which.max(ballast_pam_Fo_bc$y)]

ballast_pam_Fv_bc <- boxcox(ballast_pam$Fv+1 ~ ballast_pam$time)
ballast_pam_Fv_lambda <- ballast_pam_Fv_bc$x[which.max(ballast_pam_Fv_bc$y)]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test((ballast_pam$Fo^ballast_pam_Fo_lambda-1)/ballast_pam_Fo_lambda)
shapiro.test((ballast_pam$Fv^ballast_pam_Fv_lambda-1)/ballast_pam_Fv_lambda)

pdf("sctld ballast PAM hist transformed.pdf")
par(mfrow=c(2,2))
hist((ballast_pam$Fo^ballast_pam_Fo_lambda-1)/ballast_pam_Fo_lambda)
qqnorm((ballast_pam$Fo^ballast_pam_Fo_lambda-1)/ballast_pam_Fo_lambda)
qqline((ballast_pam$Fo^ballast_pam_Fo_lambda-1)/ballast_pam_Fo_lambda)

hist((ballast_pam$Fv^ballast_pam_Fv_lambda-1)/ballast_pam_Fv_lambda)
qqnorm((ballast_pam$Fv^ballast_pam_Fv_lambda-1)/ballast_pam_Fv_lambda)
qqline((ballast_pam$Fv^ballast_pam_Fv_lambda-1)/ballast_pam_Fv_lambda)
dev.off()

# ANOVA on untransformed data
ballast_pam_anova_Fo <- aov(Fo ~ time*treatment, data=ballast_pam)
summary(ballast_pam_anova_Fo)

ballast_pam_anova_Fv <- aov(Fv ~ time*treatment, data=ballast_pam)
summary(ballast_pam_anova_Fv)

# ANOVA on Box-Cox transformed data
ballast_pam_anova_Fo_bc <- aov(((Fo^ballast_pam_Fo_lambda-1)/ballast_pam_Fo_lambda) ~ time*treatment, data=ballast_pam)
summary(ballast_pam_anova_Fo_bc)

ballast_pam_anova_Fv_bc <- aov(((Fv^ballast_pam_Fv_lambda-1)/ballast_pam_Fv_lambda) ~ time*treatment, data=ballast_pam)
summary(ballast_pam_anova_Fv_bc)

# Q-Q plots for raw and transformed data
pdf("sctld ballast PAM normality ANOVA.pdf")
par(mfrow=c(2,2))
qqnorm(ballast_pam_anova_Fo$residuals, main="Fo"); qqline(ballast_pam_anova_Fo$residuals);
qqnorm(ballast_pam_anova_Fv$residuals, main="Fv"); qqline(ballast_pam_anova_Fv$residuals)
qqnorm(ballast_pam_anova_Fo_bc$residuals, main="Fo transformed"); qqline(ballast_pam_anova_Fo_bc$residuals)
qqnorm(ballast_pam_anova_Fv_bc$residuals, main="Fv transformed"); qqline(ballast_pam_anova_Fv_bc$residuals)
dev.off()

ballast_pam_tukey_Fo <- TukeyHSD(ballast_pam_anova_Fo_bc)
ballast_pam_tukey_Fo
write.csv(ballast_pam_tukey_Fo$`time:treatment`, file = "sctld ballast PAM Fo Tukey.csv")

ballast_pam_tukey_Fv <- TukeyHSD(ballast_pam_anova_Fv_bc)
ballast_pam_tukey_Fv
write.csv(ballast_pam_tukey_Fv$`time:treatment`, file = "sctld ballast PAM Fv Tukey.csv")


#### ballast figures ####

library(rcompanion)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(gridExtra)

# PCoA of cell counts by taxon 
ballast_counts_pcoa <- pcoa(ballast_counts_dist)
scores <- ballast_counts_pcoa$values
vectors <- ballast_counts_pcoa$vectors

pdf(file="sctld ballast cellcounts PCoA.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(vectors[,1], vectors[,2],col=c("#f6e8c3","#d8b365","#5ab4ac")[as.numeric(as.factor(ballast_counts$treatment))],pch=c(15,1)[as.numeric((as.factor(ballast_counts$time)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Treatment")
# ordispider(vectors, ballast_counts$treatment, label=F, col=c("#f6e8c3","#d8b365","#5ab4ac"))
legend("topright", legend=c("DW24","DW120","HW120"), fill = c("#f6e8c3","#d8b365","#5ab4ac"), bty="n")
legend("bottomright", legend=c("t0","t1"), pch=c(15,1), bty="n")
plot(vectors[,1], vectors[,2],col=c("black","grey75")[as.numeric(as.factor(ballast_counts$time))],pch=c(15,18,25)[as.numeric(as.factor(ballast_counts$treatment))], xlab="Coordinate 1", ylab="Coordinate 2", main="Time")
# ordispider(vectors, ballast_counts$time, label=F, col=c("black","grey75"))
legend("topright", legend=c("t0","t1"), fill = c("black","grey75"), bty="n")
legend("bottomright", legend=c("DW24","DW120","HW120"), pch=c(15,18,25), bty="n")
dev.off()

# stacked column chart

# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
ballast_counts_perm_letters <- data.frame(cbind(ballast_counts_perm_out$Comparison,ballast_counts_perm_out$'Pr(>F)'))
ballast_counts_perm_letters <- na.omit(ballast_counts_perm_letters)
ballast_counts_perm_letters <- dplyr::rename(ballast_counts_perm_letters, comparison = X1, p.adj = X2)
ballast_counts_perm_letters$comparison = c("t0_HW120-t0_DW24","t0_HW120-t0_DW120","t0_HW120-t1_HW120", "t0_HW120-t1_DW24", "t0_HW120-t1_DW120", "t0_DW24-t0_DW120", "t0_DW24-t1_HW120", "t0_DW24-t1_DW24", "t0_DW24-t1_DW120", "t0_DW120-t1_HW120", "t0_DW120-t1_DW24", "t0_DW120-t1_DW120", "t1_HW120-t1_DW24", "t1_HW120-t1_DW120", "t1_DW24-t1_DW120")
ballast_counts_perm_letters$p.adj <- as.numeric(paste(ballast_counts_perm_letters$p.adj))
ballast_counts_perm_letters

# creates compact letter display of significant pairwise differences for figure
ballast_counts_perm_cld <- cldList(p.adj ~ comparison, data = ballast_counts_perm_letters, threshold = 0.05, remove.zero = FALSE)
ballast_counts_perm_cld=ballast_counts_perm_cld[order(ballast_counts_perm_cld$Group),] 
ballast_counts_perm_cld$time <- ballast_counts_perm_cld$Group
ballast_counts_perm_cld <- ballast_counts_perm_cld %>% add_column(taxa="other")
ballast_counts_perm_cld

# transposing and reformatting dataframe to make abundance a single column
ballast_counts_percent <- dplyr::select(ballast_counts, 2:3,6:12,21)
ballast_counts_percent <- melt(ballast_counts_percent, id = c("time","treatment","time_treatment"))
ballast_counts_percent$time_treatment=factor(ballast_counts_percent$time_treatment, levels=unique(ballast_counts_percent$time_treatment)) 
ballast_counts_percent <- cast(ballast_counts_percent, time_treatment~variable, mean)
ballast_counts_percent <- melt(ballast_counts_percent, id = time_treatment)
ballast_counts_percent <- dplyr::rename(ballast_counts_percent, taxa = variable, abundance = value)
ballast_counts_percent$time_treatment2 <- ballast_counts_percent$time_treatment
ballast_counts_percent <- ballast_counts_percent %>% separate(time_treatment2, c('time', 'treatment'))
ballast_counts_percent

# creates dummy dataframes for custom text labels of stats outputs
ballast_counts_perm_stats <- data.frame(time_treatment = "t1_DW24",abundance = 0,lab = "time: F1,17=2.976, p=0.021",time = "t1",taxa="other")
ballast_counts_perm_stats2 <- data.frame(time_treatment = "t1_DW24",abundance = 0,lab = "treat: F2,17=4.071, p=0.001",time = "t1",taxa="other")
ballast_counts_perm_stats3 <- data.frame(time_treatment = "t1_DW24",abundance = 0,lab = "int: F2,17=2.729, p=0.012",time = "t1",taxa="other")

ballast_counts_percent$taxa=factor(ballast_counts_percent$taxa, levels=unique(ballast_counts_percent$taxa)) 

# percent stacked barplot
ballast_count_abundance <- ggplot(ballast_counts_percent, aes(fill=taxa, y=abundance, x=time_treatment)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x = "Time",
       y = "Percent abundance",
       fill = 'Taxon') + 
  scale_fill_manual(values=c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#56B4E9", "#CC79A7","#999999")) +
  facet_grid(. ~ time, scales="free", switch = c("both")) +
  geom_text(data = ballast_counts_perm_stats,label = ballast_counts_perm_stats$lab, aes(x=3, y=0.95)) +
  geom_text(data = ballast_counts_perm_stats2,label = ballast_counts_perm_stats2$lab, aes(x=2.85, y=0.85)) +
  geom_text(data = ballast_counts_perm_stats3,label = ballast_counts_perm_stats3$lab, aes(x=2.85, y=0.75)) +
  # geom_text(data = ballast_counts_perm_cld, aes(x=2, y=-0.025, label=Letter)) +
  theme_classic() 
ballast_count_abundance 

ggsave("sctld ballast counts abundance plot.pdf", plot= ballast_count_abundance, width=6, height=3, units="in", dpi=300)

# cell counts
# creating dataframe of the pairwise comparisons needed for plots and doing a bit of table reformatting
ballast_counts_letters <- data.frame(ballast_counts_tukey$`time:treatment`)
ballast_counts_letters$comparison <- rownames(ballast_counts_letters)
ballast_counts_letters$comparison = str_replace_all(ballast_counts_letters$comparison,":","_")
ballast_counts_letters$p.adj[is.na(ballast_counts_letters$p.adj)] <- 1
ballast_counts_letters

# creates compact letter display of significant pairwise differences for figure
ballast_counts_cld <- cldList(p.adj ~ comparison, data = ballast_counts_letters, threshold = 0.05, remove.zero = F)
ballast_counts_cld$Group  <- factor(ballast_counts_cld$Group , levels = c("t0_DW24", "t0_DW120","t0_HW120","t1_DW24", "t1_DW120","t1_HW120"))
ballast_counts_cld=ballast_counts_cld[order(ballast_counts_cld$Group),]
ballast_counts_cld$time_treatment <- ballast_counts_cld$Group
ballast_counts_cld$time <- c("t0","t0","t0","t1","t1","t1")
ballast_counts_cld

# creates dummy dataframes for custom text labels of stats outputs
ballast_counts_stats <- data.frame(time_treatment = "t1_DW24",lab = "time: F1,17=4.895, p=0.047",time = "t1")
ballast_counts_stats2 <- data.frame(time_treatment = "t1_DW24",lab = "treat: F2,17=12.54, p=0.001",time = "t1")
ballast_counts_stats3 <- data.frame(time_treatment = "t1_DW24",lab = "int: F2,17=6.433, p=0.013",time = "t1")

# keeps order as imported
ballast_counts$time_treatment=factor(ballast_counts$time_treatment, levels=unique(ballast_counts$time_treatment)) 

# boxplot of live cell counts among time/treatments
ballast_counts_plot <-
  ggboxplot(
    ballast_counts,
    x = "time_treatment",
    y = "live",
    color = "grey30",
    palette = c("#f6e8c3","#d8b365","#5ab4ac","#f6e8c3","#d8b365","#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Time",
           y = "Live cell counts (L-1)",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", legend.text = element_blank()) +
  geom_text(data = ballast_counts_cld, aes(x=time_treatment, y=0, label=Letter)) +
  geom_text(data = ballast_counts_stats,label = ballast_counts_stats$lab, aes(x=2.5, y=1800)) +
  geom_text(data = ballast_counts_stats2,label = ballast_counts_stats2$lab, aes(x=2.55, y=1600)) +
  geom_text(data = ballast_counts_stats3,label = ballast_counts_stats3$lab, aes(x=2.45, y=1400))
ballast_counts_plot

ggsave("sctld ballast counts plot.pdf", plot= ballast_counts_plot, width=6, height=3, units="in", dpi=300)

# PAM Fo
ballast_pam_letters_Fo <- data.frame(ballast_pam_tukey_Fo$'time:treatment')
ballast_pam_letters_Fo$comparison <- rownames(ballast_pam_letters_Fo)
ballast_pam_letters_Fo$comparison = str_replace_all(ballast_pam_letters_Fo$comparison,":","_")
ballast_pam_letters_Fo$p.adj[is.na(ballast_pam_letters_Fo$p.adj)] <- 1
ballast_pam_letters_Fo

# creates compact letter display of significant pairwise differences for figure
ballast_pam_cld_Fo <- cldList(p.adj ~ comparison, data = ballast_pam_letters_Fo, threshold = 0.05, remove.zero = F)
ballast_pam_cld_Fo$Group  <- factor(ballast_pam_cld_Fo$Group , levels = c("t0_DW24", "t0_DW120","t0_HW120","t1_DW24", "t1_DW120","t1_HW120"))
ballast_pam_cld_Fo=ballast_pam_cld_Fo[order(ballast_pam_cld_Fo$Group),]
ballast_pam_cld_Fo$time_treatment <- ballast_pam_cld_Fo$Group
ballast_pam_cld_Fo$time <- c("t0","t0","t0","t1","t1","t1")
ballast_pam_cld_Fo

# creates dummy dataframes for custom text labels of stats outputs
ballast_pam_Fo_stats <- data.frame(time_treatment = "t1_DW24",lab = "time: F1,89=19.29, p<0.001",time = "t1")
ballast_pam_Fo_stats2 <- data.frame(time_treatment = "t1_DW24",lab = "treat: F2,89=119.65, p<0.001",time = "t1")
ballast_pam_Fo_stats3 <- data.frame(time_treatment = "t1_DW24",lab = "int: F2,89=97.17, p<0.001",time = "t1")

ballast_pam$time_treatment <- paste(ballast_pam$time,ballast_pam$treatment, sep = "_")

ballast_pam$time_treatment=factor(ballast_pam$time_treatment, levels=unique(ballast_pam$time_treatment)) 

# boxplot of Fo among time/treatments
ballast_pam_Fo_plot <-
  ggboxplot(
    ballast_pam,
    x = "time_treatment",
    y = "Fo",
    color = "grey30",
    palette = c("#f6e8c3","#d8b365","#5ab4ac","#f6e8c3","#d8b365","#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Time",
           y = "Fo",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", legend.text = element_blank()) +
  geom_text(data = ballast_pam_cld_Fo, aes(x=time_treatment, y=15, label=Letter)) +
  geom_text(data = ballast_pam_Fo_stats,label = ballast_pam_Fo_stats$lab, aes(x=2.5, y=85)) +
  geom_text(data = ballast_pam_Fo_stats2,label = ballast_pam_Fo_stats2$lab, aes(x=2.55, y=80)) +
  geom_text(data = ballast_pam_Fo_stats3,label = ballast_pam_Fo_stats3$lab, aes(x=2.43, y=75))
ballast_pam_Fo_plot

ggsave("sctld ballast PAM Fo plot.pdf", plot= ballast_pam_Fo_plot, width=6, height=3, units="in", dpi=300)

# PAM Fv
ballast_pam_letters_Fv <- data.frame(ballast_pam_tukey_Fv$'time:treatment')
ballast_pam_letters_Fv$comparison <- rownames(ballast_pam_letters_Fv)
ballast_pam_letters_Fv$comparison = str_replace_all(ballast_pam_letters_Fv$comparison,":","_")
ballast_pam_letters_Fv$p.adj[is.na(ballast_pam_letters_Fv$p.adj)] <- 1
ballast_pam_letters_Fv

# creates compact letter display of significant pairwise differences for figure
ballast_pam_cld_Fv <- cldList(p.adj ~ comparison, data = ballast_pam_letters_Fv, threshold = 0.05, remove.zero = F)
ballast_pam_cld_Fv$Group  <- factor(ballast_pam_cld_Fv$Group , levels = c("t0_DW24", "t0_DW120","t0_HW120","t1_DW24","t1_DW120","t1_HW120"))
ballast_pam_cld_Fv=ballast_pam_cld_Fv[order(ballast_pam_cld_Fv$Group),]
ballast_pam_cld_Fv$time_treatment <- ballast_pam_cld_Fv$Group
ballast_pam_cld_Fv$time <- c("t0","t0","t0","t1","t1","t1")
ballast_pam_cld_Fv

# creates dummy dataframes for custom text labels of stats outputs
ballast_pam_Fv_stats <- data.frame(time_treatment = "t1_DW24",lab = "time: F1,89=174.7, p<0.001",time = "t1")
ballast_pam_Fv_stats2 <- data.frame(time_treatment = "t1_DW24",lab = "treat: F2,899=262.3, p<0.001",time = "t1")
ballast_pam_Fv_stats3 <- data.frame(time_treatment = "t1_DW24",lab = "int: F2,89=61.2, p<0.001",time = "t1")

# boxplot of Fv among time/treatments
ballast_pam_Fv_plot <-
  ggboxplot(
    ballast_pam,
    x = "time_treatment",
    y = "Fv",
    color = "grey30",
    palette = c("#f6e8c3","#d8b365","#5ab4ac","#f6e8c3","#d8b365","#5ab4ac"),
    fill = "time_treatment",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + 
  labs(x = "Time",
           y = "Fv",
           fill = 'Treatment') + 
  facet_grid(. ~ time, scales="free", space="free", switch = c("both")) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right", legend.text = element_blank()) +
  geom_text(data = ballast_pam_cld_Fv, aes(x=time_treatment, y=-20, label=Letter)) +
  geom_hline(yintercept=400, linetype="dashed", color = "grey") +
  geom_text(data = ballast_pam_Fv_stats,label = ballast_pam_Fv_stats$lab, aes(x=2.5, y=690)) +
  geom_text(data = ballast_pam_Fv_stats2,label = ballast_pam_Fv_stats2$lab, aes(x=2.55, y=630)) +
  geom_text(data = ballast_pam_Fv_stats3,label = ballast_pam_Fv_stats3$lab, aes(x=2.4, y=570)) 
ballast_pam_Fv_plot

ggsave("sctld ballast PAM Fv plot.pdf", plot= ballast_pam_Fv_plot, width=6, height=3, units="in", dpi=300)

# cell counts/PAM multiplot
multiplot2=grid.arrange(ballast_count_abundance, ballast_counts_plot, ballast_pam_Fo_plot, ballast_pam_Fv_plot, ncol=1, nrow=4, widths=c(10), heights=c(3,3,3,3))

ggsave("sctld ballast watertesting panel.pdf", plot= multiplot2, width=6, height=12, units="in", dpi=300)
