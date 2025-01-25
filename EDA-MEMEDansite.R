#### Data Manipulations ####

library(ggplot2)
library(dplyr)
library(tidyr)
# install.packages("gplots")
# devtools::install_github("gmonette/spida2")
library(spida2)
library(gplots)
library(devtools)
library(car)
library(DescTools)
library(tidyverse)
library(readxl)

# Reading the data

MemeDansite <- read_excel("Downloads/MemeDansite veriler SON-1.xlsx")

# 0 means control group

MemeDansite$AltTip

MemeDansite$Tanı[is.na(MemeDansite$Tanı)] <- 0

MemeDansite[MemeDansite$Dansite1 == "b","Dansite1"] = "B"

MemeDansite[MemeDansite$Dansite2 == "b","Dansite2"] = "B"

MemeDansite[MemeDansite$Dansite3 == "b","Dansite3"] = "B"

MemeDansite[MemeDansite$Dansite4 == "b","Dansite4"] = "B"

MemeDansite[MemeDansite$Dansite5 == "b","Dansite5"] = "B"


# Assigning random numbers as name 


MemeDansite$İsim = 1453:(1452+nrow(MemeDansite)) 

# Pivot Tarih columns to long format

long_Meme_Dans <- MemeDansite %>%
  pivot_longer(
    cols = starts_with("Tarih"),
    names_to = "Tarih_Kontrol",
    values_to = "Tarih",
    names_pattern = "Tarih(\\d+)Ay"
  )


# Pivot Dansite columns to long format

long_dansite <- MemeDansite %>%
  pivot_longer(
    cols = starts_with("Dansite"),
    names_to = "Dansite_Kontrol",
    values_to = "Dansite",
    names_pattern = "Dansite(\\d+)"
  )


# Add the long form Dansite Variable 

long_Meme_Dans$Dansite = long_dansite$Dansite

long_Meme_Dans[long_Meme_Dans$Dansite == "b","Dansite"] = "B"

long_Meme_Dans$Dansite = as.factor(long_Meme_Dans$Dansite)

# Turning Diagnosis into a Factor Variable

long_Meme_Dans[491:nrow(long_Meme_Dans),"AltTip"] = "0"

long_Meme_Dans$AltTip = as.factor(long_Meme_Dans$AltTip)

long_Meme_Dans$Grup = as.factor(long_Meme_Dans$Grup)

long_Meme_Dans$Grup = ifelse(long_Meme_Dans$Grup==1,"Control","Patient")



# Spaghetti Plots # 

# Based on Control and Patient Groups 

plot_C_P <- ggplot(long_Meme_Dans, aes(x = Tarih, y = Dansite,group = İsim)) +
  geom_line(alpha=1,aes(colour = Grup),linewidth=0.9,orientation = T) +
  labs(title = "Mammographic Density vs Time Before Diagnosis for Control and Patient Subject Groups", 
       x = "Time before Diagnosis", y = "Mammographic Density") +
  theme_minimal() +
  facet_wrap(~ Grup, scales = "fixed") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 3.4), linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = 6.5, y = 3.6, label = "Diagnosis", vjust = 2, hjust = 0.1, size = 5, color = "red") + 
  scale_x_reverse() +
  scale_color_manual(values = c("Control" = "gray", "Patient" = "navy"), name = "Group") +  # Line colors
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"),  # Increase title size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16),  # Increase y-axis title size
        axis.text.x = element_text(size = 16),  # Increase x-axis text size
        axis.text.y = element_text(size = 16,face = "bold"),
        strip.text = element_text(size = 16,face = "bold"))

# Alternative Line Plot 


# Data Preparation for LongCat

MemeDansite[99:nrow(MemeDansite),"AltTip"] = "0"
Na_MemeDansite = MemeDansite[!is.na(MemeDansite$AltTip),]
Na_MemeDansite$AltTip = as.factor(Na_MemeDansite$AltTip)

y = MemeDansite[,c("Dansite1","Dansite2","Dansite3","Dansite4","Dansite5")]
times = MemeDansite[,c("Tarih1Ay","Tarih2Ay","Tarih3Ay","Tarih4Ay","Tarih5Ay")]
times = cbind(times,0)
labels = levels(long_Meme_Dans$Dansite)
groupLabels = c("","","")

cols = c("purple4","#377EB8","#4DAF4A")
cols = c("purple4","#377EB8","darkgreen")
cols = brewer.pal(4,"Spectral")
cols = c(cols[1:2],"yellow3")
cols = c("tan2","tomato2","brown3")
"Violin"

# Alt-tip 

lc = longCat(y = (y),times = -(times),id = MemeDansite$İsim,Labels = labels)
y = as.matrix(y)
y[y=="B"] = 1
y[y=="C"] = 2
y[y=="D"] = 3
y  = matrix(as.numeric(y),nrow = nrow(y),ncol = ncol(y))
lc$y = y

par(bg='transparent', mar=c(4, 4, 4, 4), xpd=TRUE,cex=1.5)
lc.g <- sorter(lc, group=Na_MemeDansite$AltTip, groupLabels=groupLabels)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', lwd=8, legendBuffer=0,
            llwd = 8,ylab = "Mammographic Density",main="Mammographic Density vs Time Before Diagnosis for Sub-type Cancers \n N = 242") 

text(x = 3, y = 265,cex=1.2 ,labels = "3", srt = 270, adj = 1, col = "black")
text(x = 3, y = 240,cex=1.2 ,labels = "2", srt = 270, adj = 1, col = "black")
text(x = 3, y = 190,cex=1.2 ,labels = "1", srt = 270, adj = 1, col = "black")
text(x = 3, y = 70,cex=1.2 ,labels = "0", srt = 270, adj = 1, col = "black")

legend(-84, 320, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = T,xjust = 0)

#legend(-115, 260, legend=labels, col=cols, lty=1, lwd=8,
#       bty = "n", 
#       inset = c(0, 0),x.intersp = 0.2,
#       y.intersp = 0.5,horiz = F,xjust = 0)

#legend(-72, 295, legend=labels, col=cols, lty=1, lwd=8,
#       bty = "n", 
#       inset = c(0, 0),x.intersp = 0.2,
#       y.intersp = 0.5,horiz = T,xjust = 0)

# Data Preparation for longCat 

y = MemeDansite[,c("Dansite1","Dansite2","Dansite3","Dansite4","Dansite5")]
times = MemeDansite[,c("Tarih1Ay","Tarih2Ay","Tarih3Ay","Tarih4Ay","Tarih5Ay")]
times = cbind(times,0)
labels = levels(long_Meme_Dans$Dansite)
groupLabels = levels(as.factor(long_Meme_Dans$Grup))

#cols <- brewer.pal(3, "Set1")

#cols = c("purple4","#377EB8","#4DAF4A")

#cols = c("purple4","#377EB8","darkgreen")

y$Dansite1 = as.factor(y$Dansite1)
y$Dansite2 = as.factor(y$Dansite2)
y$Dansite3 = as.factor(y$Dansite3)
y$Dansite4 = as.factor(y$Dansite4)
y$Dansite5 = as.factor(y$Dansite5)

lc = longCat(y = y,times = -(times),id = MemeDansite$İsim,Labels = labels)
y = as.matrix(y)
y[y=="B"] = 3
y[y=="C"] = 2
y[y=="D"] = 1
y  = matrix(as.numeric(y),nrow = nrow(y),ncol = ncol(y))
lc$y = y

# Overall 

par(bg='transparent', mar = c(5, 4, 4, 2), xpd=T,cex=1.5)
longCatPlot(lc, cols=cols, xlab='Time Before Diagnosis', lwd=10, 
            legendBuffer=0,lcex = 1,ylab = "Mammographic Density",
            llwd = 8,main="Mammographic Density vs Time Before Diagnosis \n N = 242")


legend(-115, 260, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = F,xjust = 0)

# legend(105, 200, legend=labels, col=cols, lty=1, lwd=4)

# Control vs Patient 

grDevices::pdf("Heatmap_Mammographic_Density.pdf",width = 14,height = 10)
grouplabels = c("","")
par(bg='transparent', mar = c(4, 4, 4, 4), xpd=T,cex=1.5)
lc.g <- sorter(lc, group=MemeDansite$Grup, groupLabels=grouplabels)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', 
            lwd=8, legendBuffer=0,groupBuffer = 0,gcex = 1.2,
            llwd = 8,ylab = "Mammographic Density",main="Mammographic Density vs Time Before Diagnosis for Control and Patient Subject Groups \n N = 643")

text(x = 5, y = 410,cex=1.2 ,labels = "Control", srt = 270, adj = 1, col = "black")
text(x = 5, y = 60,cex=1.2 ,labels = "Patient", srt = 270, adj = 1, col = "black")

legend(-71, 750, legend=c("D","C","B"), col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.5,
       y.intersp = 0.5,horiz = T,xjust = 0,seg.len = 4)

dev.off()

# Menapoz 

grDevices::pdf("Heatmap_Mammographic_Density_Menopause.pdf",width = 14,height = 10)
grouplabels = c("","")
par(bg='transparent', mar = c(4, 4, 4, 4), xpd=T,cex=1.5)
lc.g <- sorter(lc, group=MemeDansite$Menapoz, groupLabels=grouplabels)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', 
            lwd=8, legendBuffer=0,groupBuffer = 0,gcex = 1.2,
            llwd = 8,ylab = "Mammographic Density",main="Mammographic Density vs Time Before Diagnosis for Control and Patient Subject Groups \n N = 643")

text(x = 5, y = 400,cex=1.2 ,labels = "Control", srt = 270, adj = 1, col = "black")
text(x = 5, y = 60,cex=1.2 ,labels = "Patient", srt = 270, adj = 1, col = "black")

legend(-65, 750, legend=c("D","C","B"), col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.5,
       y.intersp = 0.5,horiz = T,xjust = 0,seg.len = 4)

dev.off()

# Aile Öyküsü 

grouplabels = c("","")
par(bg='transparent', mar = c(4, 4, 4, 4), xpd=T,cex=1.5)
lc.g <- sorter(lc, group=MemeDansite$AileÖyküsü, groupLabels=grouplabels)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', 
            lwd=8, legendBuffer=0,groupBuffer = 0,gcex = 1.2,
            llwd = 8,ylab = "Mammographic Density",main="Mammographic Density vs Time Before Diagnosis for Aile Öyküsü \n N = 242")

text(x = 5, y = 230,cex=1.2 ,labels = "1", srt = 270, adj = 1, col = "black")
text(x = 5, y = 90,cex=1.2 ,labels = "0", srt = 270, adj = 1, col = "black")


legend(-115, 260, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = F,xjust = 0)

legend(-84, 295, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = T,xjust = 0)


# DoğumSayısı 

grouplabels = c("","","","","","")
par(bg='transparent', mar = c(4, 4, 4, 4), xpd=T,cex=1.5)
lc.g <- sorter(lc, group=MemeDansite$DoğumSayısı, groupLabels=grouplabels)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', 
            lwd=8, legendBuffer=0,groupBuffer = 0,gcex = 1.2,
            llwd = 8,ylab = "Mammographic Density",main="Mammographic Density vs Time Before Diagnosis for Doğum Sayısı \n N = 242")

text(x = 4, y = 301,cex=1.2 ,labels = "5", srt = 270, adj = 1, col = "black")
text(x = 4, y = 273,cex=1.2 ,labels = "4", srt = 270, adj = 1, col = "black")
text(x = 4, y = 219,cex=1.2 ,labels = "3", srt = 270, adj = 1, col = "black")
text(x = 4, y = 130,cex=1.2 ,labels = "2", srt = 270, adj = 1, col = "black")
text(x = 4, y = 65,cex=1.2 ,labels = "1", srt = 270, adj = 1, col = "black")
text(x = 4, y = 18,cex=1.2 ,labels = "0", srt = 270, adj = 1, col = "black")


legend(-115, 260, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = F,xjust = 0)


# For Each Cancer Type 

na_long_meme = long_Meme_Dans[!is.na(long_Meme_Dans$AltTip),]
custom_colors <- colorRampPalette(c("gray", "navy"))(length(levels(na_long_meme$AltTip)))
custom_colors = c("gray","navy","navy","navy")

plot_Tani <- ggplot(na_long_meme, aes(x = Tarih, y = Dansite,group = İsim)) +
  geom_line(alpha=1,aes(colour = AltTip),linewidth=0.9) +
  labs(title = "Mammographic Density vs Time Before Diagnosis for Cancer Sub-Types and Control Group", 
       x = "Time before Diagnosis", y = "Mammographic Density") +
  theme_minimal() +
  facet_wrap(~ AltTip, scales = "fixed") +
  scale_x_reverse() +
  scale_color_manual(values = custom_colors) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 3.3), linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = 6.5, y = 3.6, label = "Diagnosis", vjust = 2, hjust = 0.1, size = 5, color = "red") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"),  # Increase title size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16),  # Increase y-axis title size
        axis.text.x = element_text(size = 16),  # Increase x-axis text size
        axis.text.y = element_text(size = 16,face = "bold"),
        strip.text = element_text(size = 16,face = "bold"))

plot_Tani <- ggplot(na_long_meme, aes(x = Tarih, y = Dansite,group = İsim,colour = as.factor(İsim))) +
  geom_line(alpha=1,aes(colour = as.factor(İsim)),linewidth=0.9) +
  labs(title = "Mammographic Density vs Time Before Diagnosis for Cancer Sub-Types and Control Group", 
       x = "Time before Diagnosis", y = "Mammographic Density") +
  theme_minimal() +
  facet_wrap(~ AltTip, scales = "fixed") +
  scale_x_reverse() +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 3.3), linetype = "dashed", color = "red", size = 0.8) +
  annotate("text", x = 6.5, y = 3.6, label = "Diagnosis", vjust = 2, hjust = 0.1, size = 5, color = "red") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"),  # Increase title size
        axis.title.x = element_text(size = 16),  # Increase x-axis title size
        axis.title.y = element_text(size = 16),  # Increase y-axis title size
        axis.text.x = element_text(size = 16),  # Increase x-axis text size
        axis.text.y = element_text(size = 16,face = "bold"),
        strip.text = element_text(size = 16,face = "bold"))




# Frequency Plot for Cancer Type

p1 = ggplot(na_long_meme, aes(x = AltTip)) + geom_bar(fill="navyblue") + 
  theme_minimal() + ggtitle("Bar Plot of Sub-types") + 
  labs(x="Cancer Type",y="Frequency")

# Frequency Plot for Dansite Type

p2 = ggplot(long_Meme_Dans, aes(x = Dansite)) + geom_bar(fill="skyblue") + 
  theme_minimal() + ggtitle("Bar Plot of Mammographic Density") + 
  labs(x="Mammographic Density",y="Frequency")

grid.arrange(p1,p2)

# Spaghetti Plots # 

# Aile Öyküsü

plot_AileO <- ggplot(long_Meme_Dens, aes(x = Tarih, y = Densite,group = İsim)) +
  geom_line(alpha=0.2,aes(colour = AileÖyküsü)) +
  labs(title = "Mammographic Density vs Time Before Diagnosis among AileÖyküsü Type", 
       x = "Time before Diagnosis (0 = Diagnosis Happened)", y = "Mammographic Density") +
  theme_minimal() +
  facet_wrap(~ AileÖyküsü, scales = "fixed") +
  scale_color_manual(values = c("0"="gray","1"="navy")) +
  theme(legend.position = "none")

# Doğum Sayısı 

custom_colors <- colorRampPalette(c("gray", "navy"))(length(levels(as.factor(long_Meme_Dans$DoğumSayısı))))

plot_Dogum <- ggplot(long_Meme_Dens, aes(x = Tarih, y = Densite,group = İsim)) +
  geom_line(alpha=0.2,aes(colour = as.factor(DoğumSayısı))) +
  labs(title = "Mammographic Density vs Time Before Diagnosis for DoğumSayısı", 
       x = "Time before Diagnosis (0 = Diagnosis Happened)", y = "Mammographic Density") +
  theme_minimal() +
  facet_wrap(~DoğumSayısı, scales = "fixed") +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "none")

# Doğum Sayısı for Each Diagnosis Type 

p41 = ggplot(na_long_meme, aes(x = AltTip, y = DoğumSayısı,fill=AltTip)) + 
  geom_boxplot() + theme_minimal() + ggtitle("Box Plot of Doğum Sayısı by Sub-Types")

# Doğum Sayısı for Each Dansite Type

p42 = ggplot(long_Meme_Dans, aes(x = Dansite, y = DoğumSayısı,fill=Dansite)) + 
  geom_boxplot() + theme_minimal() + ggtitle("Box Plot of Doğum Sayısı by Mammographic Density")

grid.arrange(p41,p42,ncol=2)

### Statistical Tests ###

# Kruskal Wallis test for Dansite Type among Diagnosis Type

kruskal_test_result <- kruskal.test(Dansite ~ AltTip, data = long_Meme_Dans)
print(kruskal_test_result)

# Perform Dunn's test for post-hoc analysis with p-value adjustment (Benjamini-Hochberg method)

library(FSA)
dunn_test_result <- dunnTest(as.numeric(Dansite) ~ AltTip, data = long_Meme_Dans, method = "bh")
print(dunn_test_result)


# Create a contingency table for Dansite Type vs Diagnosis Type 

contingency_table <- table(long_Meme_Dans$Dansite, long_Meme_Dans$AltTip)

# Print the contingency table

print(contingency_table)

# Plotting the Contingency Table with Barplot

custom_colors <- colorRampPalette(c("skyblue", "blue3"))(length(levels(long_Meme_Dans$Dansite)))

ggplot(na_long_meme, aes(x = AltTip, fill = Dansite)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution of Dansite Type by Mammographic Density", x = "Sub-Types", y = "Frequency") +
  theme_minimal()

# Perform the Chi-square test for Dansite Type vs Diagnosis Type 

chi_square_result <- chisq.test(contingency_table)

# Print the Chi-square test result

print(chi_square_result)

### Last and First Observation Difference for Dansite Type ###

# Convert the ordinal variable to numeric

long_Meme_Dans <- long_Meme_Dans %>%
  mutate(Dansite_numeric = as.numeric(Dansite))


# Calculate the first and last observations

long_meme_diff <- long_Meme_Dans %>%
  group_by(İsim) %>%
  summarize(
    FirstObservation = first(Dansite),
    LastObservation = last(Dansite),
    AltTip = first(AltTip)
  )


# Calculate the difference between the last and first observations

long_meme_diff_numerical <- long_Meme_Dans %>%
  group_by(İsim) %>%
  summarize(
    FirstObservation = first(Dansite_numeric),
    LastObservation = last(Dansite_numeric),
    Difference = LastObservation - FirstObservation,
    AltTip = first(AltTip)
  )

# Considering a numerical difference between the Dansite Type 

# ANOVA Test for the Difference of the first and last Dansite Observation among Diagnosis Type

aov_test_result <- aov(Difference ~ AltTip, data = long_meme_diff_numerical)
summary(aov_test_result)

# Fisher LSD-Test 

PostHocTest(aov_test_result,method = "lsd")

# Considering a categorical difference between the Dansite Type

long_diff <- long_meme_diff %>%
  mutate(DifferenceCategory = map2_chr(FirstObservation, LastObservation, map_diff_to_label))

long_diff$DifferenceCategory = as.factor(long_diff$DifferenceCategory)


# Kruskal Wallis Test for the Difference of the first and last Dansite Observation among Diagnosis Type

kruskal_test_result <- kruskal.test(DifferenceCategory ~ AltTip, data = long_diff)

kruskal_test_result

# Dunn Test for pair-wise difference 

dunn_test_result <- dunnTest(as.numeric(DifferenceCategory) ~ AltTip, data = long_diff, method = "bh")

dunn_test_result


# Create a contingency table - Aile Oykusu vs Dansite

contingency_table <- table(long_Meme_Dans$Dansite, long_Meme_Dans$AileÖyküsü)

# Print the contingency table

print(contingency_table)

# Plotting the contingency table for AileÖyküsü vs Dansite 

custom_colors <- colorRampPalette(c("skyblue", "blue3"))(length(levels(long_Meme_Dans$Dansite)))

ggplot(long_Meme_Dans, aes(x = AileÖyküsü, fill = Dansite)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution of Mammographic Density by AileÖyküsü", x = "AileÖyküsü", y = "Frequency") +
  theme_minimal()

# Perform the Chi-square test for AileÖyküsü vs Dansite

chi_square_result <- chisq.test(contingency_table)

# Print the Chi-square test result

print(chi_square_result)

# Create a contingency table - Aile Oykusu vs Diagnosis

contingency_table <- table(long_Meme_Dans$AileÖyküsü, long_Meme_Dans$AltTip)

# Print the contingency table

print(contingency_table)

# Plotting the contingency table for AileÖyküsü vs Diagnosis 

custom_colors <- colorRampPalette(c("skyblue", "blue3"))(length(levels(long_Meme_Dans$AltTip)))

ggplot(na_long_meme, aes(x = AileÖyküsü, fill = AltTip)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution of Mammographic Density by AileÖyküsü", x = "AileÖyküsü", y = "Frequency") +
  theme_minimal()

# Perform the Chi-square test

chi_square_result <- chisq.test(contingency_table)

# Print the Chi-square test result

print(chi_square_result)


# ANOVA Test For Doğumsayısı vs Diagnosis

aov_test_result <- aov(DoğumSayısı ~ AltTip, data = long_Meme_Dans)
summary(aov_test_result)
PostHocTest(aov_test_result,method = "lsd")

# Anova Test for Doğumsayısı vs Dansite 

aov_test_result <- aov(DoğumSayısı ~ Dansite, data = long_Meme_Dans)
summary(aov_test_result)
PostHocTest(aov_test_result,method = "lsd")

### Age Group Analysis - Before 50 - After 50 ###

# Subsetting Age Groups 

before_50 <- subset(na_long_meme,Yaş< 50)
after_50s <- subset(na_long_meme,Yaş>= 50)

# Extracting the Proportion of the Age Groups 

before_50s_prop <- nrow(before_50)/nrow(long_Meme_Dans)
after_50s_prop <- nrow(after_50s)/nrow(long_Meme_Dans)

Prop_age <- data.frame("Before 50" = before_50s_prop,"After 50" = after_50s_prop)

rownames(Prop_age) <- "Proportion"

colnames(Prop_age) <- c("Before 50","After 50")

Prop_age <- Prop_age *100

# Proportion of Age for Each Group

Prop_age

# Range of Age

range_of_age <- range(long_Meme_Dans$Yaş)
range_of_age


# Line plots of Dansite Type over time for different patients

# Before 50

p11 = ggplot(before_50, aes(x = Tarih, y = Dansite, group = İsim)) + 
  geom_line(alpha = 0.2,col="#4682b4") + theme_minimal() + 
  ggtitle("Dansite Type Over Time - Before 50") + 
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0, 100))

# After 50

p12 = ggplot(after_50s, aes(x = Tarih, y = Dansite, group = İsim)) + 
  geom_line(alpha = 0.2,col="#00008b") + theme_minimal() + 
  ggtitle("Dansite Type Over Time - After 50") + 
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(0, 100))+
  scale_y_discrete(limits = c("b", "B", "C", "D"))


# Plotting together

grid.arrange(p11,p12,ncol=2)


# Frequency of Breast Cancer Type for Each Age Group

p33 = ggplot(before_50, aes(x = AltTip)) + geom_bar(fill="#4682b4") + theme_minimal() + ggtitle("Frequency of Sub-types")
p34 = ggplot(after_50s, aes(x = AltTip)) + geom_bar(fill="#00008b") + theme_minimal() + ggtitle("Frequency of Sub-types")

grid.arrange(p33,p34)

# Frequency of Dansite Type for Each Age Group

p31 = ggplot(before_50, aes(x = Dansite)) + geom_bar(fill="#4682b4") + theme_minimal() + ggtitle("Frequency of Mammographic Density")
p32 = ggplot(after_50s, aes(x = Dansite)) + geom_bar(fill="#00008b") + theme_minimal() + ggtitle("Frequency of Mammographic Density")

grid.arrange(p31,p32)

# Adding a column to categorize age

long_Meme_Dans <- long_Meme_Dans %>%
  mutate(age_group = case_when(
    Yaş <= 50 ~ "Before 50",
    Yaş > 50  ~ "After 50"
  ))

long_Meme_Dans$age_group = as.factor(long_Meme_Dans$age_group)

# Kruska Wallis Test for each Age Group vs Dansite

kruskal_test_result <- kruskal.test(Dansite ~ age_group, data = long_Meme_Dans)
print(kruskal_test_result)

# Kruska Wallis Test for each Age Group vs Diagnosis Type

kruskal_test_result <- kruskal.test(AltTip ~ age_group, data = long_Meme_Dans)
print(kruskal_test_result)


lc = longCat(y = y,times = -(times),id = MemeDansite$İsim,Labels = labels)

# Age - Groups - Heat-mapperoni

grouplabels = c("","")
MemeDansite$Age = ifelse(MemeDansite$Yaş>50, "After 50","Before 50")
lc.g <- sorter(lc, group=MemeDansite$Age, groupLabels=grouplabels)
par(bg='transparent', mar = c(4, 4, 4, 2), xpd=T,cex=1.5)
longCatPlot(lc.g, cols=cols, xlab='Time Before Diagnosis', lwd=10, 
            legendBuffer=0,lcex = 1,ylab = "Mammographic Density",
            llwd = 8,main="Mammographic Density vs Time Before Diagnosis \n N = 242")

text(x = 5, y = 180,cex=1.2 ,labels = "Before 50", srt = 270, adj = 1, col = "black")
text(x = 5, y = 60,cex=1.2 ,labels = "After 50", srt = 270, adj = 1, col = "black")

legend(-115, 260, legend=labels, col=cols, lty=1, lwd=8,
       bty = "n", 
       inset = c(0, 0),x.intersp = 0.2,
       y.intersp = 0.5,horiz = F,xjust = 0)



library(longpower)
library(simr)

# Example of setting up a generalized linear mixed model
# Define a simple GLMM model

long_Meme_Dans$Dansite = as.factor(long_Meme_Dans$Dansite)
model <- glmer(Dansite ~ 0 + (1 + Grup*Tarih || İsim), data = long_Meme_Dans, family = binomial)
summary(model)

a = summary(model)
var(a$residuals)
a$residuals
t = seq(1,5,1)

table(MemeDansite$Grup)

# Based on Alzheimer's Disease Example 

edland.linear.power(delta=1.5, t=t, sig.level=0.05,
                    power = NULL,alternative = "two.sided",n = 98,lambda = 98/(242-98),
                    sig2.int = 12,sig2.s = 24,sig2.e = 10)


# Simulation Study 
library(lme4)

model <- glmer(Dansite ~ Grup + (1 + Tarih || İsim), data = long_Meme_Dans, family = binomial)
summary(model)
model2 <- extend(model, along="Grup", n=1200)
pc2 <- powerCurve(model2)
print(pc2)

powerSim(model)

# 3)	Hypothetical clinical trial – Treatment Effect – longpower package

Ra <- matrix(0.25, nrow = 5, ncol = 5)
diag(Ra) <- 1 #exchangeable correlation matrix for group A
ra <- c(1,2,3,4,5) #retention in group A
sigmaa <- 10  #standard deviation for group A

power.mmrm(Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 1.5, power = NULL,N = 242,sig.level = 0.05
           ,lambda = 98/(242-98))




