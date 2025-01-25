# Some of the below libraries are required 

library(GLMMadaptive)
library(ordinal)
library(questionr)
library(dplyr)
library(spida2)
library(nlme)
library(tibble)
library(multipol)
library(pracma)
library(matconv)
library(matlab)
library(Matrix)
library(AlgDesign)
library(MASS)
library(lme4)
library(readxl)
library(MCMCglmm)
library(brms)
library(tidyr)

### Two Phase Model ### 

two_phase_model = function(data, ydep, wdep, zvar, dvar, id,Random_Intercept = FALSE){
  
  stop <- 0
  
  if (length(data) == 0) {
    stop <- 1
    stop("ERROR: DATA is not provided.")
  }
  
  if (length(ydep) == 0) {
    stop <- 1
    stop("ERROR: Dependent variable YDEP is not provided.")
  }
  
  if (length(wdep) == 0) {
    stop <- 1
    stop("ERROR: Dependent variable WDEP is not provided.")
  }
  
  if (length(zvar) == 0) {
    stop <- 1
    stop("ERROR: Independent variable ZVAR is not provided.")
  }
  
  if (length(dvar) == 0) {
    stop <- 1
    stop("ERROR: Independent variable DVAR is not provided.")
  }
  
  if (length(id) == 0) {
    stop <- 1
    stop("ERROR: ID is not provided.")
  }
  
  if (stop > 0) {
    cat("NOTE: Macro stopped because of errors!")
    return("")
  }
  
  # NAIVE ESTIMATES
  
  p <- length(zvar)
  q <- length(dvar)
  
  towork <- data
  
  # Assigning the name for the id variable
  
  colnames(towork)[colnames(towork)==id] = "id"
  
  # Arranging with id order
  
  towork <- towork %>% arrange(id)  # Might Adjust later 
  
  # First Phase Model #
  
  formula_r1 <- paste(dvar,collapse = " + ")
    
  formula_f <- bf(paste(wdep," ~ (1 +", formula_r1,"||","id",")"))
    
  # Dansite ~ (1 + Tarih || id)
    
  model_brm_m1 <- brm(formula_f,chain=1, data = towork,family = "cumulative",iter = 20000) # Sample Size 600
    
  model_brm_m2 <- brm(formula_f,chain=1, data = towork,family = "cumulative") # Lower Iteration 
    
  model_brm_mcs <- brm(formula_f,chain=1, data = towork,family = "cumulative",iter = 20000) # Match Case Study
  
  model_brm_menapoz <- brm(formula_f,chain=1, data = towork,family = "cumulative",iter = 20000) # Menapoz = 1
  
  model_brm_menapoz2 <- brm(formula_f,chain=1, data = towork,family = "cumulative",iter = 20000) # Menapoz = 0
    
  # Random Effects 
    
  random_df_brm = brms::ranef(model_brm_menapoz)
  random_df_brm = random_df_brm$id
  random_df = random_df_brm[,1,]
  random_df = as.data.frame(random_df)
    
  # Covariance Structure 
    
  G_matrix = diag(brms::VarCorr(model_brm_menapoz)$id$sd[,1])^2
    
  # Name Structure 
    
  rownames(random_df) = unique(towork[,"id"])
    
  colnames(random_df) = c("Mean_Density_Category","Density_Change_Over_Time")
  
  # First Phase Random Information Matrix 
  
  I_M_1 = solve(G_matrix)
  
  # Second Phase Model Preparation 
  
  # Extracting the Random Parameters 
  
  mx <- rownames_to_column(as.data.frame(random_df),"id")
  
  # Combining Random Parameters to a data frame of zvar, ydep
  
  zy <- towork %>%
    group_by(id) %>%
    filter(row_number() == 1) %>%
    dplyr::select(id, zvar, ydep) %>%
    ungroup()
  
  zymx <- merge(zy, mx, by = "id")
  
  # FOR MULTICATEGORY RESPONSE, I.E. FOR MULTINOMIAL RESPONSE
  
  zymx[,ydep] = as.factor(zymx[,ydep])
  
  # Second Phase Model Formula Preperation 
  
  formula1 <- paste(zvar, collapse = " + ")
  formula2 <- paste(colnames(random_df), collapse = " + ")
  formula3 <- paste(formula1,"+",formula2)
  formula4 <- as.formula(paste(ydep, "~",formula3))
  
  # Second Phase Model | Proportional Odds Assumption is made
  
  b_or_m = nlevels(zymx$ydep)
  
  logistic_model = glm(formula = formula4, data = zymx, family = "binomial")
  
  log_model_summary = summary(logistic_model)
  
  Output = log_model_summary$coefficients
  
  Output[colnames(random_df),2] = sqrt( Output[colnames(random_df),2]^2 + diag(solve(I_M_1)) )
  
  Output[,3] = Output[,1,drop=F] / Output[,2,drop=F] 
  
  p_values <- round(2 * (1 - pnorm(abs(Output[,3,drop=F]))),6)
  
  colnames(p_values) = "p-value"
  
  # Merging the p-value
  
  Output = cbind(Output[,1:3],p_values)
  
  # Output 
  
  return(list("Table" = Output,"AIC" = AIC(logistic_model),"Residual Deviance" = logistic_model$deviance))
  
}

# Importing the dataset

MemeDansite <- read_excel("Downloads/MemeDansite veriler SON-1.xlsx")

# Data Preparation

# 0 means control group

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

long_Meme_Dans$Dansite = as.factor(long_Meme_Dans$Dansite)

long_Meme_Dans$Dansite = factor(long_Meme_Dans$Dansite,ordered = T)


# Running the Function

df = as.data.frame(long_Meme_Dans[,c("İsim","Yaş","Dansite","Tarih","Grup","AileÖyküsü","DoğumSayısı","Menapoz","VKI","AltTip","MenapozYaş")])
df$Menapoz = as.factor(df$Menapoz)
df$Grup = as.factor(df$Grup)
df$Dansite = factor(df$Dansite,ordered = T,levels = c("D","C","B"))

head(df)
df_menapoz = subset(df,Menapoz == 1)

df_menapoz

# Random Intercept Needs to be specified if needed

# With Random Intercept 

Result_TPM = two_phase_model(df,"Grup","Dansite",zvar = c("Yaş","Menapoz","VKI"),"Tarih","İsim",Random_Intercept = TRUE)

# Without Random Intercept 

two_phase_model(df,"Grup","Dansite",zvar = c("Yaş","Menapoz","VKI"),"Tarih","İsim")

# To play with manually

ydep = "Grup"
wdep = "Dansite"
dvar = "Tarih"
df$AileÖyküsü = as.factor(df$AileÖyküsü)
zvar = c("Yaş","Menapoz","VKI")
id = "İsim"
data = df
data = data[!is.na(data$VKI),]
data$Grup = as.factor(as.numeric(data$Grup) -1)

data$"1" = rep(1,nrow(data))

data = mcs_df

data = df_menapoz

model_brm_trial <- brm(Dansite ~ Grup + Yaş + Menapoz + VKI + (1 + Tarih || id),
                       chain=1, 
                       data = towork,
                       family = "cumulative",iter = 20000)


# Permutation-Based P-value 

B = 1000
tpm = list()
z_tpm = list()
bs_df = df

for(i in 1:B){
  
  bs_df$Grup = sample(bs_df$Grup)
  tpm[[i]] = two_phase_model(bs_df,"Grup","Dansite","Yaş","Tarih","İsim",Random_Intercept = TRUE)
  z_tpm[[i]] = tpm[[i]]$Table[,3]
  
}

Output = Result_TPM$Table

Int_z_variables <- lapply(z_tpm, function(x) x[[1]])
Int_z_variables <- unlist(Int_z_variables)
exact_p_int = length(Int_z_variables[z_variables > Output[1,3]]) / 94

Age_z_variables <- lapply(z_tpm, function(x) x[[2]])
Age_z_variables <- unlist(Age_z_variables)
exact_p_age = length(Age_z_variables[Age_z_variables > Output[2,3]]) / 94

MDC_z_variables <- lapply(z_tpm, function(x) x[[3]])
MDC_z_variables <- unlist(MDC_z_variables)
exact_p_ri = length(MDC_z_variables[MDC_z_variables > Output[3,3]]) / 94

DCOT_z_variables <- lapply(z_tpm, function(x) x[[4]])
DCOT_z_variables <- unlist(DCOT_z_variables)
exact_p_rs = length(DCOT_z_variables[DCOT_z_variables > Output[4,3]]) / 94

exact_p = rbind(exact_p_int, exact_p_age, exact_p_ri, exact_p_rs)

# Calculate Wald confidence intervals (95%)

z_value <- 1.96  # For a 95% confidence level

wald_conf_int <- cbind(
  round(Result_TPM$Table[,1] - z_value * Result_TPM$Table[,2],2),
  round(Result_TPM$Table[,1] + z_value * Result_TPM$Table[,2],2)
)

# Print Wald confidence intervals

colnames(wald_conf_int) <- c("2.5 %", "97.5 %")
print(wald_conf_int)

variabla_names_eng = c("Age","Menopause","BMI","Mean Density Category","Density Change Over Time","Age:Menopause")
variabla_names_tr = c("Yaş","Menapoz","VKI","Ortalama Dansite Kategorisi","Dansitenin Zamanla Değişimi","Yaş:Menapoz")

Col_name_eng = c("Risk factor", "Estimate", "95% CI", "P value")
Col_name_tr = c("İstatiksel Tahmin", "95% Güven Aralığı", "P değeri")

CI_results = paste0("(",wald_conf_int[,1], " to ", wald_conf_int[,2], ")")
Table_2 = data.frame("0" = round(Result_TPM$Table[,1],3),"1" = CI_results,"2" = round(Result_TPM$Table[,4],3))
Table_2 = Table_2[2:nrow(Table_2),]
rownames(Table_2) = variabla_names_tr
colnames(Table_2) = Col_name_tr
Table_2



z_value <- 1.96  # For a 95% confidence level

wald_conf_int <- cbind(
  round(Result_TPM$Table[,1] - z_value * Result_TPM$Table[,2],2),
  round(Result_TPM$Table[,1] + z_value * Result_TPM$Table[,2],2)
)

# Print Wald confidence intervals

colnames(wald_conf_int) <- c("2.5 %", "97.5 %")
print(wald_conf_int)

variabla_names_eng = c("Age","BMI","Mean Density Category","Density Change Over Time")
variabla_names_tr = c("Yaş","VKI","Ortalama Dansite Kategorisi","Dansitenin Zamanla Değişimi")

Col_name_eng = c("Risk factor", "Estimate", "95% CI", "P value")
Col_name_tr = c("İstatiksel Tahmin", "95% Güven Aralığı", "P değeri")

CI_results = paste0("(",wald_conf_int[,1], " to ", wald_conf_int[,2], ")")
Table_2 = data.frame("0" = round(Result_TPM$Table[,1],3),"1" = CI_results,"2" = round(Result_TPM$Table[,4],3))
Table_2 = Table_2[2:nrow(Table_2),]
rownames(Table_2) = variabla_names_tr
colnames(Table_2) = Col_name_tr
Table_2


hist(z_variables_vector,breaks = 30,probability = T)
f0 = function(x) dnorm(x)
plot(f0,xlim=c(-2,2),col="red",add=T)
shapiro.test(z_variables_vector)
qqnorm(z_variables_vector)

# Match-case Study

mcs = as.data.frame(subset(MemeDansite,Grup == 1))

n1 = as.data.frame(subset(mcs,Menapoz == 1 & Yaş < 53))
n2 = as.data.frame(subset(mcs,Menapoz == 0 & Yaş < 53))
n3 = as.data.frame(subset(mcs,Menapoz == 0 & Yaş >= 53))
n4 = as.data.frame(subset(mcs,Menapoz == 1 & Yaş >= 53))

nrow(n1)
nrow(n2)
nrow(n3)
nrow(n4)

case_n = nrow(n1) + nrow(n2) + nrow(n3) + nrow(n4)

mcs2 = as.data.frame(subset(MemeDansite,Grup == 2))

m1 = as.data.frame(subset(mcs2,Menapoz == 1 & Yaş < 53))
m2 = as.data.frame(subset(mcs2,Menapoz == 0 & Yaş < 53))
m3 = as.data.frame(subset(mcs2,Menapoz == 0 & Yaş >= 53))
m4 = as.data.frame(subset(mcs2,Menapoz == 1 & Yaş >= 53))

nrow(m1)
nrow(m2)
nrow(m3)
nrow(m4)

control_n = nrow(m1) + nrow(m2) + nrow(m3) + nrow(m4)

mcs_df = rbind(n1,n2,n3,n4,m1, m2[sample(nrow(m2), size = nrow(n2)),],
               m3[sample(nrow(m3), size = nrow(n3)),],m4[sample(nrow(m4), size = nrow(n4)),])

nrow(mcs_df)

# Pivot Tarih columns to long format

mcs_long_Meme_Dans <- mcs_df %>%
  pivot_longer(
    cols = starts_with("Tarih"),
    names_to = "Tarih_Kontrol",
    values_to = "Tarih",
    names_pattern = "Tarih(\\d+)Ay"
  )


# Pivot Dansite columns to long format

mcs_long_dansite <- mcs_df %>%
  pivot_longer(
    cols = starts_with("Dansite"),
    names_to = "Dansite_Kontrol",
    values_to = "Dansite",
    names_pattern = "Dansite(\\d+)"
  )

# Add the long form Dansite Variable 

mcs_long_Meme_Dans$Dansite = mcs_long_dansite$Dansite

mcs_long_Meme_Dans$Dansite = as.factor(mcs_long_Meme_Dans$Dansite)

mcs_long_Meme_Dans$Dansite = factor(mcs_long_Meme_Dans$Dansite,ordered = T)


# Running the Function

mcs_df = as.data.frame(mcs_long_Meme_Dans[,c("İsim","Yaş","Dansite","Tarih","Grup","AileÖyküsü","DoğumSayısı","Menapoz","VKI","AltTip")])
mcs_df$Menapoz = as.factor(mcs_df$Menapoz)
mcs_df$Grup = as.factor(mcs_df$Grup)
mcs_df$Dansite = factor(mcs_df$Dansite,ordered = T,levels = c("D","C","B"))


# Stratified Sampling 

hist(MemeDansite$Yaş)
MemeDansite$Menapoz = as.factor(MemeDansite$Menapoz)

n1 = as.data.frame(subset(MemeDansite,Menapoz == 1 & Yaş < 50))
n2 = as.data.frame(subset(MemeDansite,Menapoz == 0 & Yaş < 50))
n3 = as.data.frame(subset(MemeDansite,Menapoz == 1 & Yaş >= 50))
n4 = as.data.frame(subset(MemeDansite,Menapoz == 0 & Yaş >= 50))

n11 = floor(250 * nrow(n1) / nrow(MemeDansite))
n12 = floor(250 * nrow(n2) / nrow(MemeDansite))
n21 = floor(250 * nrow(n3) / nrow(MemeDansite))
n22 = floor(250 * nrow(n4) / nrow(MemeDansite))

library(dplyr)

N11 = sample_n(tbl = n1,size = n11)
N12 = sample_n(tbl = n2,size = n12)
N21 = sample_n(tbl = n3,size = n21)
N22 = sample_n(tbl = n4,size = n22)

ss_df = rbind(N11,N12,N21,N22)

# Pivot Tarih columns to long format

ss_long_Meme_Dans <- ss_df %>%
  pivot_longer(
    cols = starts_with("Tarih"),
    names_to = "Tarih_Kontrol",
    values_to = "Tarih",
    names_pattern = "Tarih(\\d+)Ay"
  )


# Pivot Dansite columns to long format

ss_long_dansite <- ss_df %>%
  pivot_longer(
    cols = starts_with("Dansite"),
    names_to = "Dansite_Kontrol",
    values_to = "Dansite",
    names_pattern = "Dansite(\\d+)"
  )

# Add the long form Dansite Variable 

ss_long_Meme_Dans$Dansite = ss_long_dansite$Dansite

ss_long_Meme_Dans$Dansite = as.factor(ss_long_Meme_Dans$Dansite)

ss_long_Meme_Dans$Dansite = factor(ss_long_Meme_Dans$Dansite,ordered = T)


# Running the Function

ss = as.data.frame(ss_long_Meme_Dans[,c("İsim","Yaş","Dansite","Tarih","Grup","AileÖyküsü","DoğumSayısı","Menapoz","VKI","AltTip")])
ss$Menapoz = as.factor(ss$Menapoz)
ss$Grup = as.factor(ss$Grup)
ss$Dansite = factor(ss$Dansite,ordered = T,levels = c("D","C","B"))
ss$AltTip[491:nrow(ss)] = 0 
ss$AltTip = as.factor(ss$AltTip)

# Random Intercept Needs to be specified if needed

# With Random Intercept 

SS_Result_TPM = two_phase_model(ss,"Grup","Dansite",zvar = c("Yaş","Menapoz","VKI"),"Tarih","İsim",Random_Intercept = TRUE)

# Creating Tables 

MemeDansite$Grup = as.factor(MemeDansite$Grup)
MemeDansite$AileÖyküsü = as.factor(MemeDansite$AileÖyküsü)
MemeDansite$Dansite1 = as.factor(MemeDansite$Dansite1)
MemeDansite$Dansite2 = as.factor(MemeDansite$Dansite2)
MemeDansite$Dansite3 = as.factor(MemeDansite$Dansite3)
MemeDansite$Dansite4 = as.factor(MemeDansite$Dansite4)
MemeDansite$Dansite5 = as.factor(MemeDansite$Dansite5)

# Creating Subsets 

Case = subset(MemeDansite,Grup == 1)
Control = subset(MemeDansite,Grup == 2)

# Age - Mean & SD 

Mean_Age_Case = round(mean(Case$Yaş),2)
Mean_Age_Control = round(mean(Control$Yaş),2)

t.test(x = Case$Yaş,y = Control$Yaş)

Sd_Age_Case = round(sd(Case$Yaş),2)
Sd_Age_Control = round(sd(Control$Yaş),2)

age_case = paste0(Mean_Age_Case, " (", Sd_Age_Case,")")
age_control = paste0(Mean_Age_Control, " (", Sd_Age_Control,")")

# BMI - Mean & SD

Mean_BMI_Case = round(mean(Case$VKI),2)
control_vki = na.omit(Control$VKI)
Mean_BMI_Control = round(mean(control_vki),2)

t.test(x = Case$VKI,y = Control$VKI)

Sd_BMI_Case = round(sd(Case$VKI),2)
  Sd_BMI_Control = round(sd(control_vki),2)

BMI_case = paste0(Mean_BMI_Case, " (", Sd_BMI_Case,")")
BMI_control = paste0(Mean_BMI_Control, " (", Sd_BMI_Control,")")

# No. of longitudinal mammograms, mean (SD)c

Case_NM = paste0(5," (1)")
Control_NM = paste0(5," (1)")

# No. of months between mammograms

# Case 

D1_Case = Case$Tarih1Ay - Case$Tarih2Ay
D2_Case = Case$Tarih2Ay - Case$Tarih3Ay
D3_Case = Case$Tarih3Ay - Case$Tarih4Ay
D4_Case = Case$Tarih4Ay - Case$Tarih5Ay

D_Case = c(D1_Case, D2_Case, D3_Case, D4_Case)

D_Case_Mean = round(mean(D_Case),2) # Mean
D_Case_Sd = round(sd(D_Case),2) # SD

D_Case_v = paste0(D_Case_Mean," (",D_Case_Sd,")") 

# Control 

D1_Control = Control$Tarih1Ay - Control$Tarih2Ay
D2_Control = Control$Tarih2Ay - Control$Tarih3Ay
D3_Control = Control$Tarih3Ay - Control$Tarih4Ay
D4_Control = Control$Tarih4Ay - Control$Tarih5Ay

D_Control = c(D1_Control, D2_Control, D3_Control, D4_Control)

D_Control_Mean = round(mean(D_Control),2) # Mean
D_Control_Sd = round(sd(D_Control),2) # SD

D_Control_v = paste0(D_Control_Mean," (",D_Control_Sd,")") 

t.test(x = D_Case,y = D_Control)

# Time between last mammogram and diagnosis date

Last_MM_Diagnosis_Mean = round(mean(Case$Tarih5Ay),2) 
Last_MM_Diagnosis_Sd = round(sd(Case$Tarih5Ay),2)

Last_MM_Diagnosis = paste0(Last_MM_Diagnosis_Mean," (",Last_MM_Diagnosis_Sd,")")
Last_MM_Diagnosis_Control = NA

# Mammogram Types At First Visit

First_Dansite_Case = Case$Dansite1
table(First_Dansite_Case) 
table(First_Dansite_Case) / nrow(Case) * 100

First_Dansite_Control = Control$Dansite1
table(First_Dansite_Control) 
table(First_Dansite_Control) / nrow(Control)

# Mammogram Types At Last Visit

Last_Dansite_Case = Case$Dansite5
table(Last_Dansite_Case) 
round(table(Last_Dansite_Case) / nrow(Case) * 100,4)

Last_Dansite_Control = Control$Dansite5
table(Last_Dansite_Control) 
round(table(Last_Dansite_Control) / nrow(Control)*100,2)

# Menopause Ratio 

Case_Menopause = subset(Case, Menapoz == 1)
nrow(Case_Menopause) / nrow(Case)

n_case_menopause_1 = nrow(Case_Menopause)
n_case_menopause_0 = nrow(Case) - nrow(Case_Menopause)

Control_Menopause = subset(Control, Menapoz == 1)
nrow(Control_Menopause) / nrow(Control)

n_control_menopause_1 = nrow(Control_Menopause)
n_control_menopause_0 = nrow(Control) - nrow(Control_Menopause)

menapoz_chisq = data.frame("Menapoz=0"=c(n_case_menopause_0,n_control_menopause_0),
                           "Menapoz=1"=c(n_case_menopause_1,n_control_menopause_1))

rownames(menapoz_chisq) = c("Case","Control")

menapoz_chisq_result <- chisq.test(menapoz_chisq)

menapoz_chisq_result

# Aile Öyküsü Ratio 

Case_AÖ = subset(Case, AileÖyküsü == 1)
nrow(Case_AÖ) / nrow(Case)

n_case_fs_1 = nrow(Case_AÖ)
n_case_fs_0 = nrow(Case) - nrow(Case_AÖ)

Control_AÖ = subset(Control, AileÖyküsü == 1)
nrow(Control_AÖ) / nrow(Control)

n_control_fs_1 = nrow(Control_AÖ)
n_control_fs_0 = nrow(Control) - nrow(Control_AÖ)

fs_chisq = data.frame("FS=0"=c(n_case_fs_0,n_control_fs_0),
                           "FS=1"=c(n_case_fs_1,n_control_fs_1))

rownames(fs_chisq) = c("Case","Control")

fs_chisq_result <- chisq.test(fs_chisq)

fs_chisq_result

# Doğum Sayısı Ratio 

Mean_DS_Case = round(mean(Case$DoğumSayısı),2) 
Mean_DS_Control = round(mean(Control$DoğumSayısı),2)

Sd_DS_Case = round(sd(Case$DoğumSayısı),2)
Sd_DS_Control = round(sd(Control$DoğumSayısı),2)

DS_Case = paste0(Mean_DS_Case," (",Sd_DS_Case,")")
DS_Control = paste0(Mean_DS_Control," (",Sd_DS_Control,")")

t.test(Case$DoğumSayısı,Control$DoğumSayısı)

r_names_eng = c("Age, mean (SD)", 
                "BMI, mean (SD)",
                "No. of longitudinal mammograms, mean (SD)", 
                "No. of months between mammograms", 
                "Time between last mammogram and diagnosis date, mean (SD)", 
                "Mammogram types at first visit",
                "Mammogram types at last visit", "Menopause status", 
                "Family story status",
                "Birth rate, mean (SD)")

c_names_eng = c("Cases (n = 98)", "Controls (n = 144)")

r_names_tr = c("Yaş, ortalama (SS)", 
               "VKI, ortalama (SS)",
               "Uzunlamasına mamogram kontrol sayısı, ortalama (SS)",
               "Mamogram kontrolleri arası geçen ay sayısı, ortalama (SS)",
               "Son mammogram kontrolü ile tanı arasında geçen süre, ortalama (SS)",
               "İlk kontroldeki mamogram türleri",
               "B","C","D",
               "Son kontroldeki mamogram türleri",
               "B1","C1","D1",
               "Menapoz Statüsü",
               "Aile Öyküsü Statüsü",
               "Doğum Sayısı, ortalama (SS)")

c_names_tr = c("Vaka (n = 98)", "Kontrol (n = 144)")

Case_variables = c(age_case,BMI_case,Case_NM,D_Case_v,Last_MM_Diagnosis,"",0,39,59,"",54,38,6,58,22,DS_Case)
Control_variables = c(age_control,BMI_control,Control_NM,D_Control_v,"NA","",0,66,78,"",114,28,2,59,33,DS_Control)
Table_1 = data.frame("1" = Case_variables,"2" = Control_variables)
colnames(Table_1) = c_names_tr
rownames(Table_1) = r_names_tr

library(xtable)

xtable(Table_1,caption = "Table 1. Risk Factors at Mammography by Case-Control Status")
xtable(Table_2,caption = "Risk Factors at Mammography by Case-Control Status")


# Density Change over case and control groups

Case = subset(MemeDansite,Grup == 1)
Control = subset(MemeDansite,Grup == 2)

# Case

Case_d_d = subset(Case,Dansite1 == "D" & Dansite5 == "D")
Case_d_c = subset(Case,Dansite1 == "D" & Dansite5 == "C")
Case_d_b = subset(Case,Dansite1 == "D" & Dansite5 == "B")

Case_c_c = subset(Case,Dansite1 == "C" & Dansite5 == "C")
Case_c_b = subset(Case,Dansite1 == "C" & Dansite5 == "B")


# Control

Control_d_d = subset(Control,Dansite1 == "D" & Dansite5 == "D")
Control_d_c = subset(Control,Dansite1 == "D" & Dansite5 == "C")
Control_d_b = subset(Control,Dansite1 == "D" & Dansite5 == "B")

Control_c_c = subset(Control,Dansite1 == "C" & Dansite5 == "C")
Control_c_b = subset(Control,Dansite1 == "C" & Dansite5 == "B")


# Number of Subjects (%)

nrow(Case_d_d) / nrow(Case) * 100
nrow(Case_d_c) / nrow(Case) * 100
nrow(Case_d_b) / nrow(Case) * 100
nrow(Case_c_c) / nrow(Case) * 100
nrow(Case_c_b) / nrow(Case) * 100


nrow(Control_d_d) / nrow(Control) * 100
nrow(Control_d_c) / nrow(Control) * 100
nrow(Control_d_b) / nrow(Control) * 100
nrow(Control_c_c) / nrow(Control) * 100
nrow(Control_c_b) / nrow(Control) * 100

# Contingency Table - Overall

cont_table_density = data.frame("Case" = c(nrow(Case_d_d),nrow(Case_d_c),nrow(Case_d_b), nrow(Case_c_c), nrow(Case_c_b)), 
                                "Control" = c(nrow(Control_d_d),nrow(Control_d_c),nrow(Control_d_b),nrow(Control_c_c), nrow(Control_c_b)))

rownames(cont_table_density) = c("D-D","D-C","D-B","C-C","C-B")

chi_square_result_overall <- chisq.test(cont_table_density)
chi_square_result_overall

# Before Menopause

Before_Menopause_df = subset(MemeDansite,Menapoz == 0)

Case_BM = subset(Before_Menopause_df,Grup == 1)
Control_BM = subset(Before_Menopause_df,Grup == 2)

# Case

Case_BM_d_d = subset(Case_BM,Dansite1 == "D" & Dansite5 == "D")
Case_BM_d_c = subset(Case_BM,Dansite1 == "D" & Dansite5 == "C")
Case_BM_d_b = subset(Case_BM,Dansite1 == "D" & Dansite5 == "B")

Case_BM_c_c = subset(Case_BM,Dansite1 == "C" & Dansite5 == "C")
Case_BM_c_b = subset(Case_BM,Dansite1 == "C" & Dansite5 == "B")


# Control

Control_BM_d_d = subset(Control_BM,Dansite1 == "D" & Dansite5 == "D")
Control_BM_d_c = subset(Control_BM,Dansite1 == "D" & Dansite5 == "C")
Control_BM_d_b = subset(Control_BM,Dansite1 == "D" & Dansite5 == "B")

Control_BM_c_c = subset(Control_BM,Dansite1 == "C" & Dansite5 == "C")
Control_BM_c_b = subset(Control_BM,Dansite1 == "C" & Dansite5 == "B")


# Number of Subjects (%)

nrow(Case_BM_d_d) / nrow(Case_BM) * 100
nrow(Case_BM_d_c) / nrow(Case_BM) * 100
nrow(Case_BM_d_b) / nrow(Case_BM) * 100
nrow(Case_BM_c_c) / nrow(Case_BM) * 100
nrow(Case_BM_c_b) / nrow(Case_BM) * 100


nrow(Control_BM_d_d) / nrow(Control_BM) * 100
nrow(Control_BM_d_c) / nrow(Control_BM) * 100
nrow(Control_BM_d_b) / nrow(Control_BM) * 100
nrow(Control_BM_c_c) / nrow(Control_BM) * 100
nrow(Control_BM_c_b) / nrow(Control_BM) * 100

# Contingency Table 

cont_table_density_BM = data.frame("Case" = c(nrow(Case_BM_d_d),nrow(Case_BM_d_c),nrow(Case_BM_d_b), nrow(Case_BM_c_c), nrow(Case_BM_c_b)), 
                                "Control" = c(nrow(Control_BM_d_d),nrow(Control_BM_d_c),nrow(Control_BM_d_b),nrow(Control_BM_c_c), nrow(Control_BM_c_b)))

rownames(cont_table_density_BM) = c("D-D","D-C","D-B","C-C","C-B")

chi_square_result_BM <- chisq.test(cont_table_density_BM[-4,])
chi_square_result_BM


# After Menopause

After_Menopause_df = subset(MemeDansite,Menapoz == 1)

Case_AM = subset(After_Menopause_df,Grup == 1)
Control_AM = subset(After_Menopause_df,Grup == 2)

# Case

Case_AM_d_d = subset(Case_AM,Dansite1 == "D" & Dansite5 == "D")
Case_AM_d_c = subset(Case_AM,Dansite1 == "D" & Dansite5 == "C")
Case_AM_d_b = subset(Case_AM,Dansite1 == "D" & Dansite5 == "B")

Case_AM_c_c = subset(Case_AM,Dansite1 == "C" & Dansite5 == "C")
Case_AM_c_b = subset(Case_AM,Dansite1 == "C" & Dansite5 == "B")


# Control

Control_AM_d_d = subset(Control_AM,Dansite1 == "D" & Dansite5 == "D")
Control_AM_d_c = subset(Control_AM,Dansite1 == "D" & Dansite5 == "C")
Control_AM_d_b = subset(Control_AM,Dansite1 == "D" & Dansite5 == "B")

Control_AM_c_c = subset(Control_AM,Dansite1 == "C" & Dansite5 == "C")
Control_AM_c_b = subset(Control_AM,Dansite1 == "C" & Dansite5 == "B")


# Number of Subjects (%)

nrow(Case_AM_d_d) / nrow(Case_AM) * 100
nrow(Case_AM_d_c) / nrow(Case_AM) * 100
nrow(Case_AM_d_b) / nrow(Case_AM) * 100
nrow(Case_AM_c_c) / nrow(Case_AM) * 100
nrow(Case_AM_c_b) / nrow(Case_AM) * 100


nrow(Control_AM_d_d) / nrow(Control_AM) * 100
nrow(Control_AM_d_c) / nrow(Control_AM) * 100
nrow(Control_AM_d_b) / nrow(Control_AM) * 100
nrow(Control_AM_c_c) / nrow(Control_AM) * 100
nrow(Control_AM_c_b) / nrow(Control_AM) * 100

# Contingency Table - Overall

cont_table_density_AM = data.frame("Case" = c(nrow(Case_AM_d_d),nrow(Case_AM_d_c),nrow(Case_AM_d_b), nrow(Case_AM_c_c), nrow(Case_AM_c_b)), 
                                   "Control" = c(nrow(Control_AM_d_d),nrow(Control_AM_d_c),nrow(Control_AM_d_b),nrow(Control_AM_c_c), nrow(Control_AM_c_b)))

rownames(cont_table_density_AM) = c("D-D","D-C","D-B","C-C","C-B")

chi_square_result_AM <- chisq.test(cont_table_density_AM)
chi_square_result_AM


# All Results
chi_square_result_overall
chi_square_result_BM
chi_square_result_AM


cont_table_density_MD1 = data.frame("Case" = c(88,165), 
                                   "Control" = c(182,208))

cont_table_density_MD2 = data.frame("Case" = c(133,108,12), 
                                   "Control" = c(310,78,2))

rownames(cont_table_density_MD1) = c("C","D")
rownames(cont_table_density_MD2) = c("B","C","D")

chisq.test(cont_table_density_MD1)
chisq.test(cont_table_density_MD2)
