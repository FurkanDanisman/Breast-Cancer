# Required Libraries 

library(SimCorMultRes)
library(brms)
library(multgee)
library(MCMCglmm)
library(ggplot2)

# Number of Simulation 

B = 1000

# Original Test 

p_olg = rep(0,B)
p_nlg = rep(0,B)
power_brm = rep(0,B)

# Effect Size Variation

effect_size = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
exp(effect_size)

effect_size1 = log(seq(1.1,2,0.1))
effect_size2 = log(seq(2.2,2.8,0.2))
es = log(2.1)


p_olg1 = rep(0,B)
p_olg2 = rep(0,B)
p_olg3 = rep(0,B)
p_olg4 = rep(0,B)
p_olg5 = rep(0,B)
p_olg6 = rep(0,B)
p_olg7 = rep(0,B)
p_olg8 = rep(0,B)
p_olg9 = rep(0,B)
p_olg10 = rep(0,B)
p_olg11 = rep(0,B)
p_olg12 = rep(0,B)
p_olg13 = rep(0,B)
p_olg14 = rep(0,B)
p_olg15 = rep(0,B)


p_nlg1 = rep(0,B)
p_nlg2 = rep(0,B)
p_nlg3 = rep(0,B)
p_nlg4 = rep(0,B)
p_nlg5 = rep(0,B)
p_nlg6 = rep(0,B)
p_nlg7 = rep(0,B)
p_nlg8 = rep(0,B)
p_nlg9 = rep(0,B)
p_nlg10 = rep(0,B)
p_nlg11 = rep(0,B)
p_nlg12 = rep(0,B)
p_nlg13 = rep(0,B)
p_nlg14 = rep(0,B)
p_nlg15 = rep(0,B)


# Sample Size Variation 

sample_sizes = c(75,100,150,200,250,300,350,400,450,500)

# First Plot 

p_olg_n11 = rep(0,B)
p_olg_n21 = rep(0,B)
p_olg_n31 = rep(0,B)
p_olg_n41 = rep(0,B)
p_olg_n51 = rep(0,B)
p_olg_n61 = rep(0,B)
p_olg_n71 = rep(0,B)
p_olg_n81 = rep(0,B)
p_olg_n91 = rep(0,B)
p_olg_n101 = rep(0,B)

# Second Plot 

p_olg_n12 = rep(0,B)
p_olg_n22 = rep(0,B)
p_olg_n32 = rep(0,B)
p_olg_n42 = rep(0,B)
p_olg_n52 = rep(0,B)
p_olg_n62 = rep(0,B)
p_olg_n72 = rep(0,B)
p_olg_n82 = rep(0,B)
p_olg_n92 = rep(0,B)
p_olg_n102 = rep(0,B)

# Third Plot 

p_olg_n13 = rep(0,B)
p_olg_n23 = rep(0,B)
p_olg_n33 = rep(0,B)
p_olg_n43 = rep(0,B)
p_olg_n53 = rep(0,B)
p_olg_n63 = rep(0,B)
p_olg_n73 = rep(0,B)
p_olg_n83 = rep(0,B)
p_olg_n93 = rep(0,B)
p_olg_n103 = rep(0,B)

# Fourth Plot 

p_olg_n14 = rep(0,B)
p_olg_n24 = rep(0,B)
p_olg_n34 = rep(0,B)
p_olg_n44 = rep(0,B)
p_olg_n54 = rep(0,B)
p_olg_n64 = rep(0,B)
p_olg_n74 = rep(0,B)
p_olg_n84 = rep(0,B)
p_olg_n94 = rep(0,B)
p_olg_n104 = rep(0,B)

# Fifth Plot 

p_olg_n15 = rep(0,B)
p_olg_n25 = rep(0,B)
p_olg_n35 = rep(0,B)
p_olg_n45 = rep(0,B)
p_olg_n55 = rep(0,B)
p_olg_n65 = rep(0,B)
p_olg_n75 = rep(0,B)
p_olg_n85 = rep(0,B)
p_olg_n95 = rep(0,B)
p_olg_n105 = rep(0,B)

# Sixth Plot 

p_olg_n16 = rep(0,B)
p_olg_n26 = rep(0,B)
p_olg_n36 = rep(0,B)
p_olg_n46 = rep(0,B)
p_olg_n56 = rep(0,B)
p_olg_n66 = rep(0,B)
p_olg_n76 = rep(0,B)
p_olg_n86 = rep(0,B)
p_olg_n96 = rep(0,B)
p_olg_n106 = rep(0,B)

# First Plot 

p_nlg_n11 = rep(0,B)
p_nlg_n21 = rep(0,B)
p_nlg_n31 = rep(0,B)
p_nlg_n41 = rep(0,B)
p_nlg_n51 = rep(0,B)
p_nlg_n61 = rep(0,B)
p_nlg_n71 = rep(0,B)
p_nlg_n81 = rep(0,B)
p_nlg_n91 = rep(0,B)
p_nlg_n101 = rep(0,B)


# Second Plot 

p_nlg_n12 = rep(0,B)
p_nlg_n22 = rep(0,B)
p_nlg_n32 = rep(0,B)
p_nlg_n42 = rep(0,B)
p_nlg_n52 = rep(0,B)
p_nlg_n62 = rep(0,B)
p_nlg_n72 = rep(0,B)
p_nlg_n82 = rep(0,B)
p_nlg_n92 = rep(0,B)
p_nlg_n102 = rep(0,B)

# Third Plot 

p_nlg_n13 = rep(0,B)
p_nlg_n23 = rep(0,B)
p_nlg_n33 = rep(0,B)
p_nlg_n43 = rep(0,B)
p_nlg_n53 = rep(0,B)
p_nlg_n63 = rep(0,B)
p_nlg_n73 = rep(0,B)
p_nlg_n83 = rep(0,B)
p_nlg_n93 = rep(0,B)
p_nlg_n103 = rep(0,B)

# Fourth Plot 

p_nlg_n14 = rep(0,B)
p_nlg_n24 = rep(0,B)
p_nlg_n34 = rep(0,B)
p_nlg_n44 = rep(0,B)
p_nlg_n54 = rep(0,B)
p_nlg_n64 = rep(0,B)
p_nlg_n74 = rep(0,B)
p_nlg_n84 = rep(0,B)
p_nlg_n94 = rep(0,B)
p_nlg_n104 = rep(0,B)

# Fifth Plot 

p_nlg_n15 = rep(0,B)
p_nlg_n25 = rep(0,B)
p_nlg_n35 = rep(0,B)
p_nlg_n45 = rep(0,B)
p_nlg_n55 = rep(0,B)
p_nlg_n65 = rep(0,B)
p_nlg_n75 = rep(0,B)
p_nlg_n85 = rep(0,B)
p_nlg_n95 = rep(0,B)
p_nlg_n105 = rep(0,B)

# Sixth Plot 

p_nlg_n16 = rep(0,B)
p_nlg_n26 = rep(0,B)
p_nlg_n36 = rep(0,B)
p_nlg_n46 = rep(0,B)
p_nlg_n56 = rep(0,B)
p_nlg_n66 = rep(0,B)
p_nlg_n76 = rep(0,B)
p_nlg_n86 = rep(0,B)
p_nlg_n96 = rep(0,B)
p_nlg_n106 = rep(0,B)

# 6 Plot Together 

p_olg_es1 = list()
p_olg_es2 = list()
p_olg_es3 = list()
p_olg_es4 = list()
p_olg_es5 = list()
p_olg_es6 = list()

p_nlg_es1 = list()
p_nlg_es2 = list()
p_nlg_es3 = list()
p_nlg_es4 = list()
p_nlg_es5 = list()
p_nlg_es6 = list()

# Simulation Process
  
  for (i in 402:B) { 
    
    # Sample size
    
    sample_size <- 615
    
    # Cluster size
    
    cluster_size <- 5
    
    # Category-specific intercepts - Thresholds
    
    beta_intercepts <- c(0.5, 3)
    
    # Time-varying regression parameters associated with covariates
    
    b1 = rep(0.05,5) # Covariate 1 - Time 
    b2 = rep(es,5) # Covariate 2 - Group
    
    beta_coefficients <- matrix(c(b1,b2),5,2)
    
    # Time-stationary covariate
    
    t = list()
    
    for (k in 1:sample_size) {
      
      t5 = rt(1,5,4)*4+40
      t4 = t5 - abs(-rt(1,4,0.5)*2 - 11)
      t3 = t4 - abs(-rt(1,4,0.5)*3 - 11)
      t3[t3<21] = 21
      t2 = t3 - abs(-rt(1,4,1)*3 - 10)
      t2[t2<3] = 3
      t1 = t2 - abs(-rt(1,4,1)*1.7 - 11)
      t1[t1<0] = 0
      
      t[[k]] = round(c(t5,t4,t3,t2,t1))
      
    }
    
    t = unlist(t)
    t = as.vector(t)
    
    # Plot Comparison for each iteration 
    
    # par(mfrow=c(1,2))
    #hist(t,breaks = 80)
    #hist(df$Tarih,breaks = 80)
    
    # Indicating which iteration it is on
    
    cat("Iteration",i,"\n")
    
    # Case vs Control Probability and Simulation
    
    p1 = 0.4
    p2 = 1 - p1
    
    case_control <- rbinom(sample_size, 1, p1)
    case_control = rep(case_control, each = cluster_size)
    
    # Latent correlation matrix for the NORTA method
    
    latent_correlation_matrix <- toeplitz(c(1, 0.75, 0.75, 0.75, 0.75))
    
    # Simulation of ordinal responses
    
    simulated_ordinal_dataset <- rmult.clm(clsize = cluster_size, # Number of observation
                                           intercepts = beta_intercepts, # Cut-points for Dansite
                                           betas = beta_coefficients, # Time-varying regression parameters 
                                           xformula = ~t + case_control, # Formula
                                           cor.matrix = latent_correlation_matrix, # Latent Correlation Matrix
                                           link = "logit", # Preffered Link
                                           xdata = cbind(as.data.frame(t),as.data.frame(case_control))) # Data set
    
    # Simulated dataframe
    
    sim_df = simulated_ordinal_dataset$simdata
    
    # Multgee Package - NomLORgee Model
    
    fit_nlg = nomLORgee(factor(y) ~ t + case_control,id = id,data = sim_df)
    
    # P-value for NomLORgee Model
    
    sum_nlg = wald(fit_nlg,"case")
    p_nlg15[i] = sum_nlg$case$anova$`p-value`
    
    # Multgee Package - OrdLORgee Model
    
    sim_df$y = factor(sim_df$y,ordered = T)
    fit_olg = ordLORgee(y ~ t + case_control,id = id,data = sim_df)
    
    # P-value for OrdLORgee 
    
    sum_olg = summary(fit_olg)
    p_olg15[i] = sum_olg$coefficients[4,4]
    
    }


# Power for OrdLORgee Model 

mean(p_olg <= alpha)

# Power for NomLORgee Model

mean(p_nlg <= alpha)

# Effect Size 0.75

# Ordinal Model 

# Rho = 0.75 -> Power = 0.95
# Rho = 0.85 -> Power = 0.93
# Rho = 0.95 -> Power = 0.91

# Nominal Model - Warnings

# Rho = 0.75 -> Power = 0.54
# Rho = 0.85 -> Power = 0.53
# Rho = 0.95 -> Power = 0.50

# Effect Size 0.5

# Ordinal Model 

# Rho = 0.75 -> Power = 0.68
# Rho = 0.85 -> Power = 0.60
# Rho = 0.95 -> Power = 0.53

# Nominal Model - Warnings 

# Rho = 0.75 -> Power = 0.57
# Rho = 0.85 -> Power = 0.46
# Rho = 0.95 -> Power = 0.43

# BRM MODEL

model_brm <- brm(y ~ (0 + t||id) + case_control,chain=1, data = sim_df,family = "cumulative")

# P-value for brm

sum_brm = summary(model_brm)

power_brm[i] = !(sum_brm$fixed[3,4]>0 & sum_brm$fixed[3,3] < 0)

# Power for BRM Model

mean(power_brm)

# Alpha Value

alpha = 0.05

# Original Data set

df$Dansite = factor(df$Dansite,ordered = T)

model_olg = ordLORgee(Dansite ~ Tarih + Grup,id = İsim,data = df)
model_olg = ordLORgee(Dansite ~ Tarih + Grup,id = İsim,data = mcs_df)
summary(model_olg)

# Power Curve for Effect Size - Ordinal 

PC = c(mean(p_olg1 <= alpha), mean(p_olg2 <= alpha), mean(p_olg3 <= alpha), mean(p_olg4 <= alpha), 
       mean(p_olg5 <= alpha), mean(p_olg6 <= alpha), mean(p_olg7 <= alpha), mean(p_olg8 <= alpha), 
       mean(p_olg9 <= alpha), mean(p_olg10 <= alpha),mean(p_olg15 <= alpha),mean(p_olg11 <= alpha),
       mean(p_olg12 <= alpha),mean(p_olg13 <= alpha),mean(p_olg14 <= alpha)) 

PC_Order = c(mean(p_olg1 <= alpha), mean(p_olg2 <= alpha), mean(p_olg3 <= alpha), mean(p_olg4 <= alpha), 
       mean(p_olg5 <= alpha), mean(p_olg6 <= alpha), mean(p_olg7 <= alpha), mean(p_olg8 <= alpha), 
       mean(p_olg9 <= alpha), mean(p_olg10 <= alpha),mean(p_olg15 <= alpha),
       mean(p_olg12 <= alpha),mean(p_olg13 <= alpha)) 

# Create a data frame

data_ef <- data.frame(Odds = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3), Power = PC_Order)

# Plot the power curve

library(Unicode)

pc_olg_ef = ggplot(data_ef, aes(x = Odds, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Power Curve for Ordinal Longitudinal Data: Case vs Control Test with n = 643 & n_case = 253", x = "Odds Ratio", y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = data_ef$Odds) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_ef_tr = ggplot(data_ef, aes(x = Odds, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Uzunlamasına Sıralı Veri Seti için Vaka vs Kontrol Testi: n = 643 & n_vaka = 253", x = "Odds Oranı", y = "İstatiksel Kuvvet") +
  theme_minimal() +
  scale_x_continuous(breaks = data_ef$Odds) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )


# Power Curve for Effect Size - Nominal 

PC_nlg = c(mean(p_nlg1 <= alpha), mean(p_nlg2 <= alpha), mean(p_nlg3 <= alpha), mean(p_nlg4 <= alpha), 
           mean(p_nlg5 <= alpha), mean(p_nlg6 <= alpha), mean(p_nlg7 <= alpha), mean(p_nlg8 <= alpha), 
           mean(p_nlg9 <= alpha), mean(p_nlg10 <= alpha),mean(p_nlg15 <= alpha),mean(p_nlg11 <= alpha),mean(p_nlg12 <= alpha))

# Create a data frame

data_ef_nlg <- data.frame(Odds = seq(1.1,2.3,0.1), Power = PC_nlg)

# Plot the power curve

library(Unicode)

Pc_ef_nlg = ggplot(data_ef_nlg, aes(x = Odds, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Power Curve for Nominal Longitudinal Data: Case vs Control Test with n = 643 n_case = 253", x = "Odds Ratio", y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = data_ef_nlg$Odds) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

Pc_ef_nlg_tr = ggplot(data_ef_nlg, aes(x = Odds, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Uzunlamasına Nominal Veri Seti için Vaka vs Kontrol Testi n = 643 n_vaka = 253", x = "Odds Ratio", y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = data_ef_nlg$Odds) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

library(gridExtra)

grid.arrange(pc_olg_ef,Pc_ef_nlg,pc_olg_ef_tr,Pc_ef_nlg_tr)


grDevices::pdf("PC_Ordinal_Eng.pdf",width = 12)
pc_olg_ef
dev.off()

grDevices::pdf("PC_Nominal_Eng.pdf",width = 12)
Pc_ef_nlg
dev.off()

library(Cairo)

CairoPDF("PC_Ordinal_TR.pdf",width = 12)
pc_olg_ef_tr
dev.off()

CairoPDF("PC_Nominal_TR.pdf",width = 12)
Pc_ef_nlg_tr
dev.off()


# Power Curve for Sample Size - Ordinal

PC_n1 = c(mean(p_olg_es1[[1]] <= alpha), mean(p_olg_es1[[2]] <= alpha), mean(p_olg_es1[[3]] <= alpha), mean(p_olg_es1[[4]] <= alpha), 
          mean(p_olg_es1[[5]] <= alpha), mean(p_olg_es1[[6]] <= alpha),    mean(p_olg_es1[[7]] <= alpha), mean(p_olg_es1[[8]] <= alpha),
          mean(p_olg_es1[[9]] <= alpha),mean(p_olg_es1[[10]] <= alpha)) 

PC_n2 = c(mean(p_olg_es2[[1]] <= alpha), mean(p_olg_es2[[2]] <= alpha), mean(p_olg_es2[[3]] <= alpha), mean(p_olg_es2[[4]] <= alpha), 
          mean(p_olg_es2[[5]] <= alpha), mean(p_olg_es2[[6]] <= alpha),    mean(p_olg_es2[[7]] <= alpha), mean(p_olg_es2[[8]] <= alpha),
          mean(p_olg_es2[[9]] <= alpha),mean(p_olg_es2[[10]] <= alpha)) 

PC_n3 = c(mean(p_olg_es3[[1]] <= alpha), mean(p_olg_es3[[2]] <= alpha), mean(p_olg_es3[[3]] <= alpha), mean(p_olg_es3[[4]] <= alpha), 
          mean(p_olg_es3[[5]] <= alpha), mean(p_olg_es3[[6]] <= alpha),    mean(p_olg_es3[[7]] <= alpha), mean(p_olg_es3[[8]] <= alpha),
          mean(p_olg_es3[[9]] <= alpha),mean(p_olg_es3[[10]] <= alpha)) 

PC_n4 = c(mean(p_olg_es4[[1]] <= alpha), mean(p_olg_es4[[2]] <= alpha), mean(p_olg_es4[[3]] <= alpha), mean(p_olg_es4[[4]] <= alpha), 
          mean(p_olg_es4[[5]] <= alpha), mean(p_olg_es4[[6]] <= alpha),    mean(p_olg_es4[[7]] <= alpha), mean(p_olg_es4[[8]] <= alpha),
          mean(p_olg_es4[[9]] <= alpha),mean(p_olg_es4[[10]] <= alpha)) 

PC_n5 = c(mean(p_olg_es5[[1]] <= alpha), mean(p_olg_es5[[2]] <= alpha), mean(p_olg_es5[[3]] <= alpha), mean(p_olg_es5[[4]] <= alpha), 
          mean(p_olg_es5[[5]] <= alpha), mean(p_olg_es5[[6]] <= alpha),    mean(p_olg_es5[[7]] <= alpha), mean(p_olg_es5[[8]] <= alpha),
          mean(p_olg_es5[[9]] <= alpha),mean(p_olg_es5[[10]] <= alpha)) 

PC_n6 = c(mean(p_olg_es6[[1]] <= alpha), mean(p_olg_es6[[2]] <= alpha), mean(p_olg_es6[[3]] <= alpha), mean(p_olg_es6[[4]] <= alpha), 
          mean(p_olg_es6[[5]] <= alpha), mean(p_olg_es6[[6]] <= alpha),    mean(p_olg_es6[[7]] <= alpha), mean(p_olg_es6[[8]] <= alpha),
          mean(p_olg_es6[[9]] <= alpha),mean(p_olg_es6[[10]] <= alpha)) 

# Create a data frame

data_n1 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n1)
data_n2 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n2)
data_n3 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n3)
data_n4 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n4)
data_n5 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n5)
data_n6 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n6)

# Plot the power curve
round(exp(c(0.1,0.25,0.5,0.7,0.75,1)),2)

pc_olg_n1 = ggplot(data_n1, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "a) Odds Ratio = 1.11", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_n2 = ggplot(data_n2, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "b) Odds Ratio = 1.28", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_n3 = ggplot(data_n3, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "c) Odds Ratio = 1.65", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_n4 = ggplot(data_n4, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "d) Odds Ratio = 2.00", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_n5 = ggplot(data_n5, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "e) Odds Ratio = 2.12", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_olg_n6 = ggplot(data_n6, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "f) Odds Ratio = 2.72", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

library(patchwork)
combined_plot <- (pc_olg_n1 + pc_olg_n2) / (pc_olg_n3 + pc_olg_n4) / (pc_olg_n5 + pc_olg_n6) +
  plot_annotation(title = "Power Curve for Ordinal Longitudinal Data\n Case vs Control Test \n",
                  theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5)))

combined_plot <- (pc_olg_n1 + pc_olg_n2) / (pc_olg_n3 + pc_olg_n4) / (pc_olg_n5 + pc_olg_n6) 

# Print the combined plot
print(combined_plot)

# Power Curve for Sample Size - Nominal

PC_n1 = c(mean(p_nlg_es1[[1]] <= alpha), mean(p_nlg_es1[[2]] <= alpha), mean(p_nlg_es1[[3]] <= alpha), mean(p_nlg_es1[[4]] <= alpha), 
          mean(p_nlg_es1[[5]] <= alpha), mean(p_nlg_es1[[6]] <= alpha),    mean(p_nlg_es1[[7]] <= alpha), mean(p_nlg_es1[[8]] <= alpha),
          mean(p_nlg_es1[[9]] <= alpha),mean(p_nlg_es1[[10]] <= alpha)) 

PC_n2 = c(mean(p_nlg_es2[[1]] <= alpha), mean(p_nlg_es2[[2]] <= alpha), mean(p_nlg_es2[[3]] <= alpha), mean(p_nlg_es2[[4]] <= alpha), 
          mean(p_nlg_es2[[5]] <= alpha), mean(p_nlg_es2[[6]] <= alpha),    mean(p_nlg_es2[[7]] <= alpha), mean(p_nlg_es2[[8]] <= alpha),
          mean(p_nlg_es2[[9]] <= alpha),mean(p_nlg_es2[[10]] <= alpha)) 


PC_n3 = c(mean(p_nlg_es3[[1]] <= alpha), mean(p_nlg_es3[[2]] <= alpha), mean(p_nlg_es3[[3]] <= alpha), mean(p_nlg_es3[[4]] <= alpha), 
          mean(p_nlg_es3[[5]] <= alpha), mean(p_nlg_es3[[6]] <= alpha),    mean(p_nlg_es3[[7]] <= alpha), mean(p_nlg_es3[[8]] <= alpha),
          mean(p_nlg_es3[[9]] <= alpha),mean(p_nlg_es3[[10]] <= alpha)) 

PC_n4 = c(mean(p_nlg_es4[[1]] <= alpha), mean(p_nlg_es4[[2]] <= alpha), mean(p_nlg_es4[[3]] <= alpha), mean(p_nlg_es4[[4]] <= alpha), 
          mean(p_nlg_es4[[5]] <= alpha), mean(p_nlg_es4[[6]] <= alpha),    mean(p_nlg_es4[[7]] <= alpha), mean(p_nlg_es4[[8]] <= alpha),
          mean(p_nlg_es4[[9]] <= alpha),mean(p_nlg_es4[[10]] <= alpha)) 

PC_n5 = c(mean(p_nlg_es5[[1]] <= alpha), mean(p_nlg_es5[[2]] <= alpha), mean(p_nlg_es5[[3]] <= alpha), mean(p_nlg_es5[[4]] <= alpha), 
          mean(p_nlg_es5[[5]] <= alpha), mean(p_nlg_es5[[6]] <= alpha),    mean(p_nlg_es5[[7]] <= alpha), mean(p_nlg_es5[[8]] <= alpha),
          mean(p_nlg_es5[[9]] <= alpha),mean(p_nlg_es5[[10]] <= alpha)) 

PC_n6 = c(mean(p_nlg_es6[[1]] <= alpha), mean(p_nlg_es6[[2]] <= alpha), mean(p_nlg_es6[[3]] <= alpha), mean(p_nlg_es6[[4]] <= alpha), 
          mean(p_nlg_es6[[5]] <= alpha), mean(p_nlg_es6[[6]] <= alpha),    mean(p_nlg_es6[[7]] <= alpha), mean(p_nlg_es6[[8]] <= alpha),
          mean(p_nlg_es6[[9]] <= alpha),mean(p_nlg_es6[[10]] <= alpha)) 

# Create a data frame

data_n1 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n1)
data_n2 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n2)
data_n3 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n3)
data_n4 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n4)
data_n5 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n5)
data_n6 <- data.frame(Sample_sizes = sample_sizes, Power = PC_n6)

# Plot the power curve

pc_nlg_n1 = ggplot(data_n1, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "a) Odds Ratio = 1.11", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_nlg_n2 = ggplot(data_n2, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "b) Odds Ratio = 1.28", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_nlg_n3 = ggplot(data_n3, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "c) Odds Ratio = 1.65", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_nlg_n4 = ggplot(data_n4, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "d) Odds Ratio = 2.00", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_nlg_n5 = ggplot(data_n5, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "e) Odds Ratio = 2.12", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

pc_nlg_n6 = ggplot(data_n6, aes(x = Sample_sizes, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "f) Odds Ratio = 2.72", x = "Sample Size", y = "Power") +
  theme_minimal() + 
  scale_x_continuous(breaks = sample_sizes) +
  ylim(c(0,1)) +
  theme(
    plot.title = element_text(size = 18),  # Center and bold title
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text = element_text(size = 14),      # Increase axis text size
    axis.text.x = element_text(size = 14),    # Increase x-axis text size
    axis.text.y = element_text(size = 14)     # Increase y-axis text size
  )

combined_plot_nlg <- (pc_nlg_n1 + pc_nlg_n2) / (pc_nlg_n3 + pc_nlg_n4) / (pc_nlg_n5 + pc_nlg_n6)

# Print the combined plot

print(combined_plot_nlg)
print(combined_plot)

# Creating t variable

par(mfrow=c(1,2))
hist(rt(1000,5,4)*4+40) 
hist(MemeDansite$Tarih1Ay)

t5 = rt(1,5,4)*4+40
t4 = t5 - abs(-rt(1,4,0.5)*2 - 11)
t3 = t4 - abs(-rt(1,4,0.5)*3 - 11)
t3[t3<21] = 21
t2 = t3 - abs(-rt(1,4,1)*3 - 10)
t2[t2<3] = 3
t1 = t2 - abs(-rt(1,4,1)*1.7 - 11)
t1[t1<0] = 0

round(c(t5,t4,t3,t2,t1))

par(mfrow=c(2,2))
hist(MemeDansite$Tarih5Ay - MemeDansite$Tarih4Ay)
hist(MemeDansite$Tarih4Ay - MemeDansite$Tarih3Ay)
hist(MemeDansite$Tarih3Ay - MemeDansite$Tarih2Ay)
hist(MemeDansite$Tarih2Ay - MemeDansite$Tarih1Ay)

par(mfrow=c(1,2))

hist(-abs(-rt(1000,4,0.5)*2 - 11))
hist(MemeDansite$Tarih2Ay - MemeDansite$Tarih1Ay)

hist(-abs(-rt(1000,4,0.5)*3 - 11))
hist(MemeDansite$Tarih3Ay - MemeDansite$Tarih2Ay)

hist(-abs(-rt(1000,4,1)*3 - 11))
hist(MemeDansite$Tarih4Ay - MemeDansite$Tarih3Ay)

hist(-abs(-rt(1000,4,1)*1.7 - 11))
hist(MemeDansite$Tarih5Ay - MemeDansite$Tarih4Ay)


library(tidyr)

df_wide <- pivot_wider(sim_df, names_from = time, values_from = c(y, t))
head(df_wide)

time = as.data.frame(df_wide[,8:12])


# Creating t variable - 2

par(mfrow=c(1,2))
hist(rt(1000,5,4)*4+40) 
hist(MemeDansite$Tarih1Ay)

t5 = rt(1,5,4)*4+40
t4 = t5 - abs(-rt(1,4,0.5)*2 - 11)
t3 = t4 - abs(-rt(1,4,0.5)*3 - 11)
t2 = t3 - abs(-rt(1,4,1)*3 - 10)
t1 = t2 - abs(-rt(1,4,1)*1.7 - 11)
t1[t1<0] = 0

round(c(t5,t4,t3,t2,t1))

# OR 

t1 = rpois(1,1)*2
t2 = t1 + abs(-rt(1,4,1)*1.7 - 11)
t3 = t2 + abs(-rt(1,4,1)*3 - 10)
t4 = t3 + abs(-rt(1,4,0.5)*3 - 11)
t5 = t4 + abs(-rt(1,4,0.5)*2 - 11)

round(c(t5,t4,t3,t2,t1))

par(mfrow=c(2,2))
hist(MemeDansite$Tarih5Ay - MemeDansite$Tarih4Ay)
hist(MemeDansite$Tarih4Ay - MemeDansite$Tarih3Ay)
hist(MemeDansite$Tarih3Ay - MemeDansite$Tarih2Ay)
hist(MemeDansite$Tarih2Ay - MemeDansite$Tarih1Ay)

par(mfrow=c(1,2))

hist(-abs(-rt(1000,4,0.5)*2 - 11))
hist(MemeDansite$Tarih2Ay - MemeDansite$Tarih1Ay)

hist(-abs(-rt(1000,4,0.5)*3 - 11))
hist(MemeDansite$Tarih3Ay - MemeDansite$Tarih2Ay)

hist(-abs(-rt(1000,4,1)*3 - 11))
hist(MemeDansite$Tarih4Ay - MemeDansite$Tarih3Ay)

hist(-abs(-rt(1000,4,1)*1.7 - 11))
hist(MemeDansite$Tarih5Ay - MemeDansite$Tarih4Ay)


library(tidyr)

df_wide <- pivot_wider(sim_df, names_from = time, values_from = c(y, t))
head(df_wide)

time = as.data.frame(df_wide[,8:12])


syntetic_PC = c(0.08,0.14,0.315,0.425,0.645,0.82,0.92,0.985,0.99,0.995)
syntetic_PC_n = c(0.125/2,0.12,0.27,0.33,0.52,0.74,0.83,0.94,0.97,0.995)
es = seq(0.1,1,0.1)
or = round(or,2)
es_or = c("1.11","1.22","1.35","1.49","1.65","1.82","2.00","2.23",
          "2.46","2.71")

data_ef = data.frame("Effect_Size" = es,"Power" = syntetic_PC)

# Ordinal

pc_olg_ef = ggplot(data_ef, aes(x = Effect_Size, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Ordinal Longitudinal Data Set for Case vs Control Test: n = 242 & n_case = 98", x = "Odds Ratio", y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.1,1,0.1),labels = es_or) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 24),  # Center and bold title
    axis.title.x = element_text(size = 22),  # Increase x-axis label size
    axis.title.y = element_text(size = 22),  # Increase y-axis label size
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 20),    # Increase legend text size
    axis.text = element_text(size = 20),      # Increase axis text size
    axis.text.x = element_text(size = 20),    # Increase x-axis text size
    axis.text.y = element_text(size = 20)     # Increase y-axis text size
  )

data_ef = data.frame("Effect_Size" = es,"Power" = syntetic_PC_n)

# Nominal 

pc_nlg_ef = ggplot(data_ef, aes(x = Effect_Size, y = Power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue") +
  ggplot2::labs(title = "Nominal Longitudinal Data Set for Case vs Control Test: n = 242 & n_case = 98", x = "Odds Ratio", y = "Power") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0.1,1,0.1),labels = es_or) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",size = 24),  # Center and bold title
    axis.title.x = element_text(size = 22),  # Increase x-axis label size
    axis.title.y = element_text(size = 22),  # Increase y-axis label size
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 20),    # Increase legend text size
    axis.text = element_text(size = 20),      # Increase axis text size
    axis.text.x = element_text(size = 20),    # Increase x-axis text size
    axis.text.y = element_text(size = 20)     # Increase y-axis text size
  )

library(extrafont)
cairo_pdf("Ordinal_PC_Eng.pdf",height = 10,width = 14)
pc_olg_ef
dev.off()

cairo_pdf("Nominal_PC_Eng.pdf",height = 10,width = 14)
pc_nlg_ef
dev.off()

cairo_pdf("Ordinal_PCs_Eng.pdf",height = 10,width = 12)
combined_plot
dev.off()

cairo_pdf("Nominal_PCs_Eng.pdf",height = 10,width = 12)
combined_plot_nlg
dev.off()

getMaxArrayCorrelation = function(a,b){
  
  n = length(a)
  
  if( n < 1 | n > 2*10^5){
    stop("n is out of the acceptable range")
  }
  
  for(i in 1:length(b)){
    if( a[i] < 1 | a[i] > 2*10^5 | b[i] < 1 | b[i] > 2*10^5){
      stop("input value is out of the acceptable range")
    }
  }
  
  possiblities = rep(0,1000)
  
  for(i in 1:1000){
    sample_attempt = sample(b,size = length(b))
    s_a = sample_attempt - a
    possiblities[i] = sum(s_a[s_a > 0])
  }
  
  return(max(possiblities))
  
}

getMaxArrayCorrelation(a = c(1,2,3,4,5),b = c(2,3,4,5,6)) 
