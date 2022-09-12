#Analysis of images collected 8/11/22 of B73 and BTX623 GPS2 and GPS1-NR treated and untreated with GA3
# One of the samples in the data set was exposed to microscope conditions for 24h with images taken every 4 hours, 
# for a total of 7 images, while the other unexposed sample was only imaged at hour 0 after GA treatment, and 24h aftter GA treatment after the end of the exposed sample timecourse
library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(agricolae)
library(emmeans)
library(ggbeeswarm)
library(readr)
library(stringr)


Exposed_Sample_Before_GA3_Results <- read_csv("/Volumes/Seagate/B73_BTX623_before_and_afterexposure_noexposure/Exposed_Sample_Before_GA3/data/Exposed_Sample_Before_GA3_Results.csv")
Exposed_0h_thru_24h_Results <- read_csv("/Volumes/Seagate/B73_BTX623_before_and_afterexposure_noexposure/Exposed_0h_thru_24h_AfterGA3/data/Exposed_0h_thru_24h_Results.csv")

Exposed_0h_thru_24h_Results <- subset(Exposed_0h_thru_24h_Results, Area <= 1000)
Exposed_Sample_Before_GA3_Results <- subset(Exposed_Sample_Before_GA3_Results, Area <= 1000)

Exposed_sample <- rbind(Exposed_Sample_Before_GA3_Results, Exposed_0h_thru_24h_Results)

# Now I want to do some plotting and statistics on the exposed sample as a rep for GPS2 and GPS1-NR treatments

Exposed_sample$TIME <- "value"
Exposed_sample <- Exposed_sample %>% mutate(TIME = case_when(
  grepl(pattern = "Result of Before", x = Label) ~ "-10min",
  grepl(pattern = "slice0001", x = Label) ~ "0h",
  grepl(pattern = "slice0002", x = Label) ~ "4h",
  grepl(pattern = "slice0003", x = Label) ~ "8h",
  grepl(pattern = "slice0004", x = Label) ~ "12h",
  grepl(pattern = "slice0005", x = Label) ~ "16h",
  grepl(pattern = "slice0006", x = Label) ~ "20h",
  grepl(pattern = "slice0007", x = Label) ~ "24h",
  grepl(pattern = "t:1/7", x = Label) ~ "0h",
  grepl(pattern = "t:2/7", x = Label) ~ "4h",
  grepl(pattern = "t:3/7", x = Label) ~ "8h",
  grepl(pattern = "t:4/7", x = Label) ~ "12h",
  grepl(pattern = "t:5/7", x = Label) ~ "16h",
  grepl(pattern = "t:6/7", x = Label) ~ "20h",
  grepl(pattern = "t:7/7", x = Label) ~ "24h"
))

Exposed_sample$TREATMENT <- "value"  
Exposed_sample <- Exposed_sample %>% mutate(TREATMENT = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A2 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A3 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A4 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "B1 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B2 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B3 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B4 Region", x = Label) ~ "0uM GA3",
))

Exposed_sample$construct <- "value"  
Exposed_sample <- Exposed_sample %>% mutate(construct = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "GPS2",
  grepl(pattern = "A2 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "A3 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "A4 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B1 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "B2 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B3 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B4 Region", x = Label) ~ "GPS1-NR",
))

Exposed_sample$Species <- "value"  
Exposed_sample <- Exposed_sample %>% mutate(Species = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "B73",
  grepl(pattern = "A2 Region", x = Label) ~ "B73",
  grepl(pattern = "A3 Region", x = Label) ~ "BTX623",
  grepl(pattern = "A4 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B1 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B2 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B3 Region", x = Label) ~ "B73",
  grepl(pattern = "B4 Region", x = Label) ~ "B73",
))

Exposed_sample$trconspec <- "value"  
Exposed_sample <- Exposed_sample %>% mutate(trconspec = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "B73_GPS2_100uMGA3",
  grepl(pattern = "A2 Region", x = Label) ~ "B73_GPS1-NR_100uMGA3",
  grepl(pattern = "A3 Region", x = Label) ~ "BTX623_GPS1-NR_100uMGA3",
  grepl(pattern = "A4 Region", x = Label) ~ "BTX623_GPS2_100uMGA3",
  grepl(pattern = "B1 Region", x = Label) ~ "BTX623_GPS1-NR_0uMGA3",
  grepl(pattern = "B2 Region", x = Label) ~ "BTX623_GPS2_0uMGA3",
  grepl(pattern = "B3 Region", x = Label) ~ "B73_GPS2_0uMGA3",
  grepl(pattern = "B4 Region", x = Label) ~ "B73_GPS1-NR_0uMGA3",
))

colnames(Exposed_sample)[4] <- "Ratio"
Exposed_sample <- within(Exposed_sample, {
  TIME <- factor(TIME, levels = c("-10min", "0h", "2h", "4h", "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h", "22h", "24h"))
  TREATMENT <- factor(TREATMENT, levels = c("0uM GA3", "100uM GA3"))})

plottitle <- "Maize B73 and Sorghum BTX623 GPS2 and GPS1-NR FRET response to 100GA3"

Exposed_sample2 <- Exposed_sample %>% group_by(trconspec, TIME) %>% summarise(n = n(), mean = mean(Ratio), median = median(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem) 

# This is a decent plot, but now I want to clean it up a bit
ggplot(Exposed_sample2, aes(x = TIME, y = mean, color = trconspec, group = trconspec)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, size = 1, position = position_dodge(width = 0.3)) + 
  geom_point(lwd = 2, position = position_dodge(width = 0.3)) + labs(x = "TIME(minutes)", y = "mean(FRET/CFP)", title = plottitle)


Exposed_sample2$Species_construct <- "value"  
Exposed_sample2 <- Exposed_sample2 %>% mutate(Species_construct = case_when(
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) ~ "B73_GPS2",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) ~ "B73_GPS1-NR",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) ~ "BTX623_GPS1-NR",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) ~ "BTX623_GPS2",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) ~ "BTX623_GPS1-NR",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) ~ "BTX623_GPS2",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) ~ "B73_GPS2",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) ~ "B73_GPS1-NR",
))

Exposed_sample2$Species_TREATMENT <- "value"  
Exposed_sample2 <- Exposed_sample2 %>% mutate(TREATMENT = case_when(
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) ~ "0uMGA3",
))

# This plot is a little better, but I would like to try a violin plot too
ggplot(Exposed_sample2, aes(x = TIME, y = mean, color = Species_construct, group = trconspec, shape = TREATMENT)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, size = 1, position = position_dodge(width = 0.3)) + 
  geom_point(lwd = 5, position = position_dodge(width = 0.3)) + labs(x = "TIME", y = "FRET Ratio", title = plottitle)

library(ggpubr)
BTX623_Exposed <- subset(Exposed_sample, Species == "BTX623")
BTX623_Exposed <- subset(BTX623_Exposed, Ratio > 0)
ggviolin(BTX623_Exposed, x = BTX623_Exposed$TIME, y = BTX623_Exposed$Ratio, fill = BTX623_Exposed$trconspec, palette = c("chartreuse4", "cornflowerblue", "darkmagenta", "darkorange1"))
# write out statistical comparissons you want to make in violin plot
mean_compare <- compare_means(Ratio~trconspec, BTX623_Exposed, group.by = "TIME", method = "wilcox.test", p.adjust.method = "fdr")
mycomparisons <- list(c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"), c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3"))
BTX_violin <- ggviolin(BTX623_Exposed, x = "TIME", y = "Ratio", fill = "trconspec", title = "BTX623 GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
BTX_violn_2 <- ggpar(BTX_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16)

#Comparing means for BTX constructs
BTX623_Exposed_GPS2 <- subset(BTX623_Exposed, construct == "GPS2")
BTX623_GPS2_means_compared <- compare_means(Ratio ~ TREATMENT, BTX623_Exposed_GPS2, group.by = "TIME", method = "wilcox.test", p.adjust.method = "fdr")
BTX623_GPS2_means_summary <- BTX623_Exposed_GPS2 %>% group_by(TREATMENT, TIME) %>% summarise(n = n(), mean = mean(Ratio), median = median(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem, CI = qnorm(0.975) * sem) 
BTX_GPS2_100uM_treatment <- subset(BTX623_Exposed_GPS2, TREATMENT == "100uM GA3")
BTX_GPS2_0uM_treatment <- subset(BTX623_Exposed_GPS2, TREATMENT == "0uM GA3")
BTX_GPS2_100uM_compare <- compare_means(Ratio ~ TIME, BTX_GPS2_100uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")
BTX_GPS2_0uM_compare <- compare_means(Ratio ~ TIME, BTX_GPS2_0uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")


  
BTX623_Exposed_GPS1NR <- subset(BTX623_Exposed, construct == "GPS1-NR")
BTX623_GPS1NR_means_compared <- compare_means(Ratio ~ TREATMENT, BTX623_Exposed_GPS1NR, group.by = "TIME", method = "wilcox.test", p.adjust.method = "fdr")
BTX623_GPS1NR_means_summary <- BTX623_Exposed_GPS1NR %>% group_by(TREATMENT, TIME) %>% summarise(n = n(), mean = mean(Ratio), median = median(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem, CI = qnorm(0.975) * sem) 
BTX623_GPS1NR_100uM_treatment <- subset(BTX623_Exposed_GPS1NR, TREATMENT == "100uM GA3")
BTX_GPS1NR_0uM_treatment <- subset(BTX623_Exposed_GPS1NR, TREATMENT == "0uM GA3")
BTX_GPS1NR_100uM_compare <- compare_means(Ratio ~ TIME, BTX623_GPS1NR_100uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")
BTX_GPS1NR_0uM_compare <- compare_means(Ratio ~ TIME, BTX_GPS1NR_0uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")



# Comparing means for B73 constructs
B73_Exposed_GPS2 <- subset(B73_Exposed, construct == "GPS2")
B73_GPS2_means_compared <- compare_means(Ratio ~ TREATMENT, B73_Exposed_GPS2, group.by = "TIME", method = "wilcox.test", p.adjust.method = "fdr")
B73_GPS2_means_summary <- B73_Exposed_GPS2 %>% group_by(TREATMENT, TIME) %>% summarise(n = n(), mean = mean(Ratio), median = median(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem, CI = qnorm(0.975) * sem) 
B73_GPS2_100uM_treatment <- subset(B73_Exposed_GPS2, TREATMENT == "100uM GA3")
B73_GPS2_0uM_treatment <- subset(B73_Exposed_GPS2, TREATMENT == "0uM GA3")
B73_GPS2_100uM_compare <- compare_means(Ratio ~ TIME, B73_GPS2_100uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")
B73_GPS2_0uM_compare <- compare_means(Ratio ~ TIME, B73_GPS2_0uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")



B73_Exposed_GPS1NR <- subset(B73_Exposed, construct == "GPS1-NR")
B73_GPS1NR_means_compared <- compare_means(Ratio ~ TREATMENT, B73_Exposed_GPS1NR, group.by = "TIME", method = "wilcox.test", p.adjust.method = "fdr")
B73_GPS1NR_means_summary <- B73_Exposed_GPS1NR %>% group_by(TREATMENT, TIME) %>% summarise(n = n(), mean = mean(Ratio), median = median(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem, CI = qnorm(0.975) * sem) 
B73_GPS1NR_100uM_treatment <- subset(B73_Exposed_GPS1NR, TREATMENT == "100uM GA3")
B73_GPS1NR_0uM_treatment <- subset(B73_Exposed_GPS1NR, TREATMENT == "0uM GA3")
B73_GPS1NR_100uM_compare <- compare_means(Ratio ~ TIME, B73_GPS1NR_100uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")
B73_GPS1NR_0uM_compare <- compare_means(Ratio ~ TIME, B73_GPS1NR_0uM_treatment, method = "wilcox.test", p.adjust.method = "fdr")





# Trying to get this stupid code to work was taking forever so I am going to give up for now and just fix it in word
stat_compare_means(mapping = aes(x = TIME, y = Ratio, group = TREATMENT), comparisons = mycomparisons, method = "t.test")
BTX_violn_2 + facet_wrap(~TIME) + stat_compare_means(comparisons = list(c("BTX623_GPS1-NR_100uMGA3", "BTX623_GPS1-NR_0uMGA3"), c("BTX623_GPS2_100uMGA3", "BTX623_GPS2_0uMGA3")))

# check the results of the compare_means code statistical significance
aov_BTXex <- aov(Ratio~trconspec*TIME, BTX623_Exposed)
summary(aov_BTXex)
HSD.test(aov_BTXex, "TIME", group = T, console = T)
BTXex_mod <- lm(Ratio~trconspec*TIME, BTX623_Exposed)                                               
BTXex_emm <- emmeans(BTXex_mod, ~ trconspec | TIME)
pairs(BTXex_emm)
# comparing the p-values shows both emmeans() and compare_means() methods came to the same conclusion about significance of p-values


# Now do all the above to the exposed B73 data!
B73_Exposed <- subset(Exposed_sample, Species == "B73")
B73_Exposed <- subset(B73_Exposed, Ratio > 0)
mean_compare_B73 <- compare_means(Ratio~trconspec, B73_Exposed, group.by = "TIME", method = "wilcox.test")
B73_violin <- ggviolin(B73_Exposed, x = "TIME", y = "Ratio", fill = "trconspec", palette = c("darkorange1", "chartreuse", "cornflowerblue", "darkmagenta"), title = "B73 GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
B73_violin
B73_violn_2 <- ggpar(B73_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16)
B73_violn_2
aov_B73ex <- aov(Ratio~trconspec*TIME, B73_Exposed)
summary(aov_B73ex)
HSD.test(aov_B73ex, "TIME", group = T, console = T)
B73ex_mod <- lm(Ratio~trconspec*TIME, B73_Exposed)                                               
B73ex_emm <- emmeans(B73ex_mod, ~ trconspec | TIME)
pairs(B73ex_emm)
#comparing the 2 methods for getting p-vals still results in no significance at all timepoints of B73

# Now read in the non-exposed datasets and combine them
nonexposed_0h_ResultsA1_A4 <- read_csv("/Volumes/Seagate/B73_BTX623_before_and_afterexposure_noexposure/Unexposed_Sample_0h/data/nonexposed_0h_ResultsA1_A4.csv")
nonexposed_0h_ResultsB1_B4 <- read_csv("/Volumes/Seagate/B73_BTX623_before_and_afterexposure_noexposure/Unexposed_Sample_0h/data/nonexposed_0h_ResultsB1_B4.csv")
Unexposed_24h_Results <- read_csv("/Volumes/Seagate/B73_BTX623_before_and_afterexposure_noexposure/Unexposed_Sample_24h/data/Unexposed_24h_Results.csv")

nonexposed <- rbind(nonexposed_0h_ResultsA1_A4, nonexposed_0h_ResultsB1_B4, Unexposed_24h_Results)

nonexposed <- subset(nonexposed, Area <= 1000)

nonexposed$TIME <- "value"
nonexposed <- nonexposed %>% mutate(TIME = case_when(
  grepl(pattern = "Result of 0h", x = Label) ~ "0h",
  grepl(pattern = "Result of unexposed_24h", x = Label) ~ "24h"
))

nonexposed$TREATMENT <- "value"  
nonexposed <- nonexposed %>% mutate(TREATMENT = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A2 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A3 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "A4 Region", x = Label) ~ "100uM GA3",
  grepl(pattern = "B1 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B2 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B3 Region", x = Label) ~ "0uM GA3",
  grepl(pattern = "B4 Region", x = Label) ~ "0uM GA3",
))

nonexposed$construct <- "value"  
nonexposed <- nonexposed %>% mutate(construct = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "GPS2",
  grepl(pattern = "A2 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "A3 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "A4 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B1 Region", x = Label) ~ "GPS1-NR",
  grepl(pattern = "B2 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B3 Region", x = Label) ~ "GPS2",
  grepl(pattern = "B4 Region", x = Label) ~ "GPS1-NR",
))

nonexposed$Species <- "value"  
nonexposed <- nonexposed %>% mutate(Species = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "B73",
  grepl(pattern = "A2 Region", x = Label) ~ "B73",
  grepl(pattern = "A3 Region", x = Label) ~ "BTX623",
  grepl(pattern = "A4 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B1 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B2 Region", x = Label) ~ "BTX623",
  grepl(pattern = "B3 Region", x = Label) ~ "B73",
  grepl(pattern = "B4 Region", x = Label) ~ "B73",
))

nonexposed$trconspec <- "value"  
nonexposed <- nonexposed %>% mutate(trconspec = case_when(
  grepl(pattern = "A1 Region", x = Label) ~ "B73_GPS2_100uMGA3",
  grepl(pattern = "A2 Region", x = Label) ~ "B73_GPS1-NR_100uMGA3",
  grepl(pattern = "A3 Region", x = Label) ~ "BTX623_GPS1-NR_100uMGA3",
  grepl(pattern = "A4 Region", x = Label) ~ "BTX623_GPS2_100uMGA3",
  grepl(pattern = "B1 Region", x = Label) ~ "BTX623_GPS1-NR_0uMGA3",
  grepl(pattern = "B2 Region", x = Label) ~ "BTX623_GPS2_0uMGA3",
  grepl(pattern = "B3 Region", x = Label) ~ "B73_GPS2_0uMGA3",
  grepl(pattern = "B4 Region", x = Label) ~ "B73_GPS1-NR_0uMGA3",
))

colnames(nonexposed)[4] <- "Ratio"
nonexposed <- within(nonexposed, {
  TIME <- factor(TIME, levels = c("0h", "24h"))
  TREATMENT <- factor(TREATMENT)
  trconspec <- factor(trconspec, levels = c("B73_GPS1-NR_0uMGA3", "B73_GPS1-NR_100uMGA3", "B73_GPS2_0uMGA3", "B73_GPS2_100uMGA3", "BTX623_GPS1-NR_0uMGA3", "BTX623_GPS1-NR_100uMGA3", "BTX623_GPS2_0uMGA3", "BTX623_GPS2_100uMGA3"))})

nonexposed$Ratio <- as.numeric(nonexposed$Ratio)

nonexposed_naomit <- na.omit(nonexposed)
#if you have problems getting the summary line of code to work, you have to use the na.omit, because you probably have nas in the dataset

nonexposed_sample2 <- nonexposed_naomit %>% group_by(trconspec, TIME) %>% summarise(n = n(), mean = mean(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem) 

ggplot(nonexposed_sample2, aes(x = TIME, y = mean, color = trconspec, group = trconspec)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, size = 1, position = position_dodge(width = 0.3)) + 
  geom_point(lwd = 2, position = position_dodge(width = 0.3)) + labs(x = "TIME(minutes)", y = "mean(FRET/CFP)", title = "Nonexposed sample")

nonexposed_sample2$Species_construct <- "value"  
nonexposed_sample2 <- nonexposed_sample2 %>% mutate(Species_construct = case_when(
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) ~ "B73_GPS2",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) ~ "B73_GPS1-NR",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) ~ "BTX623_GPS1-NR",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) ~ "BTX623_GPS2",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) ~ "BTX623_GPS1-NR",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) ~ "BTX623_GPS2",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) ~ "B73_GPS2",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) ~ "B73_GPS1-NR",
))

nonexposed_sample2$Species_TREATMENT <- "value"  
nonexposed_sample2 <- nonexposed_sample2 %>% mutate(TREATMENT = case_when(
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) ~ "100uMGA3",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) ~ "0uMGA3",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) ~ "0uMGA3",
))

library(ggpubr)
#BTX623 nonexposed
BTX623_nonexposed <- subset(nonexposed, Species == "BTX623")
BTX623_nonexposed <- subset(BTX623_nonexposed, Ratio > 0)
ggviolin(BTX623_nonexposed, x = BTX623_nonexposed$TIME, y = BTX623_nonexposed$Ratio, fill = BTX623_nonexposed$trconspec, palette = c("chartreuse4", "cornflowerblue", "darkmagenta", "darkorange1"))
# write out statistical comparissons you want to make in violin plot
mean_compare_BTX623_nonexposed <- compare_means(Ratio~trconspec, BTX623_nonexposed, group.by = "TIME", method = "wilcox.test")
BTX_nonexposed_violin <- ggviolin(BTX623_nonexposed, x = "TIME", y = "Ratio", fill = "trconspec", title = "BTX623 Nonexposed GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
BTX_nonexposed_violn_2 <- ggpar(BTX_nonexposed_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16)
BTX_nonexposed_violn_2
aov_BTXunex <- aov(Ratio~trconspec*TIME, BTX623_nonexposed)
summary(aov_BTXunex)
HSD.test(aov_BTXunex, "TIME", group = T, console = T)
BTXunex_mod <- lm(Ratio~trconspec*TIME, BTX623_nonexposed)                                               
BTXunex_emm <- emmeans(BTXunex_mod, ~ trconspec | TIME)
pairs(BTXunex_emm)

#B73 nonexposed
B73_nonexposed <- subset(nonexposed, Species == "B73")
B73_nonexposed <- subset(B73_nonexposed, Ratio > 0)
ggviolin(B73_nonexposed, x = B73_nonexposed$TIME, y = B73_nonexposed$Ratio, fill = B73_nonexposed$trconspec, palette = c("chartreuse4", "cornflowerblue", "darkmagenta", "darkorange1"))
# write out statistical comparissons you want to make in violin plot
mean_compare_B73_nonexposed <- compare_means(Ratio~trconspec, B73_nonexposed, group.by = "TIME", method = "wilcox.test")
B73_nonexposed_violin <- ggviolin(B73_nonexposed, x = "TIME", y = "Ratio", fill = "trconspec", title = "B73 Nonexposed GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
B73_nonexposed_violn_2 <- ggpar(B73_nonexposed_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16)
B73_nonexposed_violn_2
aov_B73unex <- aov(Ratio~trconspec*TIME, B73_nonexposed)
summary(aov_B73unex)
HSD.test(aov_B73unex, "TIME", group = T, console = T)
B73unex_mod <- lm(Ratio~trconspec*TIME, B73_nonexposed)                                               
BTXunex_emm <- emmeans(BTXunex_mod, ~ trconspec | TIME)
pairs(BTXunex_emm)

#Now I want to have the nonexposed and exposed 0h and 24h plotted together to compare
Exposed_sample_0h_24h <- subset(Exposed_sample, TIME == "0h" | TIME == "24h")
bothexp_unexp <- rbind(Exposed_sample_0h_24h, nonexposed)

bothexp_unexp$Exposure <- "value"  
bothexp_unexp <- bothexp_unexp %>% mutate(Exposure = case_when(
  grepl(pattern = "nonlightexposed", x = Label) ~ "Unexposed",
  grepl(pattern = "Result of 0h_to_24h_timecourse", x = Label) ~ "Exposed",
  grepl(pattern = "unexposed", x = Label) ~ "Unexposed",
  grepl(pattern = "0h_to_24h_timecourse", x = Label) ~ "Exposed"
))


bothexp_unexp$trconexposure <- "value"  
bothexp_unexp <- bothexp_unexp %>% mutate(trconexposure = case_when(
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "B73_GPS2_100uMGA3_Exposed",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "B73_GPS2_0uMGA3_Exposed",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "B73_GPS1-NR_100uMGA3_Exposed",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "B73_GPS1-NR_0uMGA3_Exposed",
  grepl(pattern = "B73_GPS2_100uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "B73_GPS2_100uMGA3_Unexposed",
  grepl(pattern = "B73_GPS2_0uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "B73_GPS2_0uMGA3_Unexposed",
  grepl(pattern = "B73_GPS1-NR_100uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "B73_GPS1-NR_100uMGA3_Unexposed",
  grepl(pattern = "B73_GPS1-NR_0uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "B73_GPS1-NR_0uMGA3_Unexposed",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "BTX623_GPS2_100uMGA3_Exposed",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "BTX623_GPS2_0uMGA3_Exposed",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "BTX623_GPS1-NR_100uMGA3_Exposed",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) & grepl(pattern = "Exposed", x = Exposure) ~ "BTX623_GPS1-NR_0uMGA3_Exposed",
  grepl(pattern = "BTX623_GPS2_100uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "BTX623_GPS2_100uMGA3_Unexposed",
  grepl(pattern = "BTX623_GPS2_0uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "BTX623_GPS2_0uMGA3_Unexposed",
  grepl(pattern = "BTX623_GPS1-NR_100uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "BTX623_GPS1-NR_100uMGA3_Unexposed",
  grepl(pattern = "BTX623_GPS1-NR_0uMGA3", x = trconspec) & grepl(pattern = "Unexposed", x = Exposure) ~ "BTX623_GPS1-NR_0uMGA3_Unexposed",
))

#B73 exposed and unexposed compared
bothexposures_B73_0h <- subset(bothexp_unexp, Species == "B73" & TIME == "0h") 
B73_0h_exposed_and_unexposed <- ggviolin(bothexposures_B73_0h, x = "Exposure", y = "Ratio", fill = "trconspec", palette = c("firebrick2", "green", "dodgerblue1", "mediumorchid1"), add = "boxplot")
ggpar(B73_0h_exposed_and_unexposed, legend.title = "Construct_Treatment", font.legend = 18, legend = "right", font.x = 16, font.y = 16)
B73_ex_unex_stats_0h <- compare_means(Ratio~trconspec, bothexposures_B73_0h, group.by = "Exposure", method = "wilcox.test")

B73_ex_unex_stats_0h_trconexp <- compare_means(Ratio ~ trconexposure, bothexposures_B73_0h, group.by = "trconspec", method = "wilcox.test")

bothexposures_B73_24h <- subset(bothexp_unexp, Species == "B73" & TIME == "24h")
B73_ex_unex_stats_24h <- compare_means(Ratio~trconexposure, bothexposures_B73_24h, group.by = "trconspec", method = "wilcox.test")

aov_B73_0h_exposure <- aov(Ratio~trconexposure, bothexposures_B73_0h)
summary(aov_B73_0h_exposure)
HSD.test(aov_B73_0h_exposure, "trconexposure", group = T, console = T)
B730hunexex_mod <- lm(Ratio~trconexposure, bothexposures_B73_0h)                                               
B730hunexex_emm <- emmeans(B730hunexex_mod, ~ trconexposure)
pairs(B730hunexex_emm)

# Try facet wrapping, this is good!
facet_wrap_data_B73 <- subset(bothexp_unexp, Species == "B73")
B73_exposed_unexposed_facet_violin <- ggviolin(facet_wrap_data_B73, x = "Exposure", y = "Ratio", fill = "trconspec", palette = c("firebrick2", "green", "dodgerblue1", "mediumorchid1"), add = "boxplot") + facet_wrap(~TIME) 
ggpar(B73_exposed_unexposed_facet_violin, legend.title = "Construct_Treatment", font.legend = 18, legend = "right", font.x = 16, font.y = 16) + theme(strip.text.x = element_text(size = 20))

B73_ex_unex_stats <- compare_means(Ratio ~ TIME, facet_wrap_data_B73, group.by = "trconexposure", method = "wilcox.test")
bothexposures_B73_0h <- subset(bothexp_unexp, Species == "B73" & TIME == "0h")


#BTX623 exposed and unexposed compared
bothexposures_BTX623_0h <- subset(bothexp_unexp, Species == "BTX623" & TIME == "0h") 
BTX623_0h_exposed_and_unexposed <- ggviolin(bothexposures_BTX623_0h, x = "Exposure", y = "Ratio", fill = "trconspec", palette = c("firebrick4", "darkolivegreen", "darkturquoise", "darkviolet"), add = "boxplot")
ggpar(BTX623_0h_exposed_and_unexposed, legend.title = "Construct_Treatment", font.legend = 18, legend = "right", font.x = 16, font.y = 16)
BTX623_ex_unex_stats_0h <- compare_means(Ratio~trconspec, bothexposures_BTX623_0h, group.by = "Exposure", method = "wilcox.test")

BTX623_ex_unex_stats_0h_trconexpo <- compare_means(Ratio ~ trconexposure, bothexposures_BTX623_0h, group.by = "trconspec", method = "wilcox.test")

bothexposures_BTX623_24h <- subset(bothexp_unexp, Species == "BTX623" & TIME == "24h")
BTX623_ex_unex_stats_24h <- compare_means(Ratio~trconexposure, bothexposures_BTX623_24h, group.by = "trconspec", method = "wilcox.test")

aov_BTX623_0h_exposure <- aov(Ratio~trconexposure, bothexposures_BTX623_0h)
summary(aov_BTX623_0h_exposure)
HSD.test(aov_BTX623_0h_exposure, "trconexposure", group = T, console = T)
BTX6230hunexex_mod <- lm(Ratio~trconexposure, bothexposures_BTX623_0h)                                               
BTX6230hunexex_emm <- emmeans(BTX6230hunexex_mod, ~ trconexposure)
pairs(BTX6230hunexex_emm)

# Try facet wrapping, this is good!
facet_wrap_data_BTX623 <- subset(bothexp_unexp, Species == "BTX623")
BTX623_exposed_unexposed_facet_violin <- ggviolin(facet_wrap_data_BTX623, x = "Exposure", y = "Ratio", fill = "trconspec", palette = c("firebrick4", "darkolivegreen", "darkturquoise", "darkviolet"), add = "boxplot") + facet_wrap(~TIME) 
ggpar(BTX623_exposed_unexposed_facet_violin, legend.title = "Construct_Treatment", font.legend = 18, legend = "right", font.x = 16, font.y = 16) + theme(strip.text.x = element_text(size = 20))

BTX623_ex_unex_stats <- compare_means(Ratio ~ TIME, facet_wrap_data_BTX623, group.by = "trconexposure", method = "wilcox.test")

#For my presentation I want just the exposed -10min and 24h data to compare
BTX_10min_24h_exposed <- subset(BTX623_Exposed, TIME == "-10min" | TIME == "24h")
BTX_10min_24h_violin <- ggviolin(BTX_10min_24h_exposed, x = "TIME", y = "Ratio", fill = "trconspec", title = "BTX623 GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
BTX_10min_24h_violin_2 <- ggpar(BTX_10min_24h_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16) + guides(fill = guide_legend(nrow = 2)) 
BTX_10min_24h_violin_2

B73_10min_24h_exposed <- subset(B73_Exposed, TIME == "-10min" | TIME == "24h")
B73_10min_24h_violin <- ggviolin(B73_10min_24h_exposed, x = "TIME", y = "Ratio", fill = "trconspec", palette = c("firebrick2", "green", "dodgerblue1", "mediumorchid1"), title = "B73 GPS2 and GPS1-NR FRET Ratio with Exogenous GA3", add = "boxplot") 
B73_10min_24h_violin_2 <- ggpar(B73_10min_24h_violin, legend.title = "Construct_Treatment", font.legend = 18, font.title = 22, legend = "bottom", font.x = 16, font.y = 16) + guides(fill = guide_legend(nrow = 2)) 
B73_10min_24h_violin_2





