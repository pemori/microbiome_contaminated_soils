#Scripts for the article "Microbial refuge under shrubs limits soil Cu contamination effects"
#Submitted to Soil Biology and Biochemistry
#
#December 2nd, 2025

# Inquiries to Dr. Pedro Mondaca > pedro.mondaca@pucv.cl

setwd("C:\\Users\\pedro\\Research\\microbiome_contaminated_soils\\Github")

library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(multcomp)
library(car)
library(lsmeans)
library(openxlsx)
library(emmeans)
library(userfriendlyscience)

#Load data ####
#Mapping
sample<-factor(c("Bulk Soil","Bulk Soil","Bulk Soil","Soil surrounding root","Soil surrounding root","Soil surrounding root"))
contamination<-factor(c("High Cont","Mid Cont","Uncont","High Cont","Mid Cont","Uncont"))
sampxcont<-factor(c("BS HC","BS MC","BS UC","SSR HC","SSR MC","SSR UC"))
mapping<-data.frame(sample,contamination,sampxcont)
mapping <-mapping %>% slice(rep(1:n(),each=10))
summary(mapping)

#BGQ
BGQ <- read_excel("BGQ3.xlsx", 
                  sheet = "All")
names <- c('pH','CE','SOM','N','P_disp','K_disp','Cu_sol','pCu','N_disp')
BGQ[,names] <- lapply(BGQ[,names] , as.numeric)
BGQ1<-BGQ[,-c(1:4)]
BGQ_map <- cbind(BGQ1,mapping)
rownames(BGQ_map)<-BGQ$ID
str(BGQ_map)

#Summary table
# Group by the factor and summarize mean and standard error for each numeric variable
#fx for standard error
calculate_SE <- function(x) {
  n <- sum(!is.na(x))
  if (n > 1) {
    return(sd(x, na.rm = TRUE) / sqrt(n))
  } else {
    return(NA)
  }
}
BGQ_map
summary_table <- BGQ_map %>%
  group_by(sampxcont) %>%
  summarise_at(
    vars(2:11),
    list(mean = ~ mean(., na.rm = TRUE), SE = ~ calculate_SE(.))
  )
summary_table<-as.data.frame(summary_table)
View(summary_table)

# Adjust summarise for all environmental variables
BGQ_map %>%
  group_by(contamination) %>%
  summarise(
    n = n(),
    missing_pCu = sum(is.na(pCu)),
    mean_pCu = mean(pCu, na.rm = TRUE),
    
    missing_pH = sum(is.na(pH)),
    mean_pH = mean(pH, na.rm = TRUE),
    
    missing_CE = sum(is.na(CE)),
    mean_CE = mean(CE, na.rm = TRUE),
    
    missing_SOM = sum(is.na(SOM)),
    mean_SOM = mean(SOM, na.rm = TRUE),
    
    missing_N = sum(is.na(N)),
    mean_N = mean(N, na.rm = TRUE),
    
    missing_P_disp = sum(is.na(P_disp)),
    mean_P_disp = mean(P_disp, na.rm = TRUE),
    
    missing_K_disp = sum(is.na(K_disp)),
    mean_K_disp = mean(K_disp, na.rm = TRUE),
    
    missing_Cu = sum(is.na(Cu)),
    mean_Cu = mean(Cu, na.rm = TRUE),
    
    missing_Cu_sol = sum(is.na(Cu_sol)),
    mean_Cu_sol = mean(Cu_sol, na.rm = TRUE)
  )


# Create a filtered dataset excluding missing values
BGQ_map_nomiss <- BGQ_map %>%
  filter(
    !is.na(pCu),
    !is.na(pH),
    !is.na(CE),
    !is.na(SOM),
    !is.na(N),
    !is.na(P_disp),
    !is.na(K_disp),
    !is.na(Cu),
    !is.na(Cu_sol)
  )

BGQ_map_nomiss

# Subset directo por contamination
BGQ_UC <- subset(BGQ_map_nomiss, contamination == "Uncont")
BGQ_MC <- subset(BGQ_map_nomiss, contamination == "Mid Cont")
BGQ_HC <- subset(BGQ_map_nomiss, contamination == "High Cont")

# Subset directo por sample
BGQ_BS <- subset(BGQ_map_nomiss, sample == "Bulk Soil")
BGQ_SSR <- subset(BGQ_map_nomiss, sample == "Soil surrounding root")


# ANOVA-like test (Welch ANOVA)
oneway.test(pH ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$pH, BGQ_UC$sample, method="games-howell")

oneway.test(pH ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$pH, BGQ_MC$sample, method="games-howell")

oneway.test(pH ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$pH, BGQ_HC$sample, method="games-howell")

oneway.test(pH ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$pH, BGQ_BS$contamination, method="games-howell")

oneway.test(pH ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$pH, BGQ_SSR$contamination, method="games-howell")

#CE
# ANOVA-like test (Welch ANOVA)
oneway.test(CE ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$CE, BGQ_UC$sample, method="games-howell")

oneway.test(CE ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$CE, BGQ_MC$sample, method="games-howell")

oneway.test(CE ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$CE, BGQ_HC$sample, method="games-howell")

oneway.test(CE ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$CE, BGQ_BS$contamination, method="games-howell")

oneway.test(CE ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$CE, BGQ_SSR$contamination, method="games-howell")


#SOM
# ANOVA-like test (Welch ANOVA)
oneway.test(SOM ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$SOM, BGQ_UC$sample, method="games-howell")

oneway.test(SOM ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$SOM, BGQ_MC$sample, method="games-howell")

oneway.test(SOM ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$SOM, BGQ_HC$sample, method="games-howell")

oneway.test(SOM ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$SOM, BGQ_BS$contamination, method="games-howell")

oneway.test(SOM ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$SOM, BGQ_SSR$contamination, method="games-howell")

#N_disp
# ANOVA-like test (Welch ANOVA)
oneway.test(N_disp ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$N_disp, BGQ_UC$sample, method="games-howell")

oneway.test(N_disp ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$N_disp, BGQ_MC$sample, method="games-howell")

oneway.test(N_disp ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$N_disp, BGQ_HC$sample, method="games-howell")

oneway.test(N_disp ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$N_disp, BGQ_BS$contamination, method="games-howell")

oneway.test(N_disp ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$N_disp, BGQ_SSR$contamination, method="games-howell")

#P_disp
# ANOVA-like test (Welch ANOVA)
oneway.test(P_disp ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$P_disp, BGQ_UC$sample, method="games-howell")

oneway.test(P_disp ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$P_disp, BGQ_MC$sample, method="games-howell")

oneway.test(P_disp ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$P_disp, BGQ_HC$sample, method="games-howell")

oneway.test(P_disp ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$P_disp, BGQ_BS$contamination, method="games-howell")

oneway.test(P_disp ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$P_disp, BGQ_SSR$contamination, method="games-howell")

#K_disp
# ANOVA-like test (Welch ANOVA)
oneway.test(K_disp ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$K_disp, BGQ_UC$sample, method="games-howell")

oneway.test(K_disp ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$K_disp, BGQ_MC$sample, method="games-howell")

oneway.test(K_disp ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$K_disp, BGQ_HC$sample, method="games-howell")

oneway.test(K_disp ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$K_disp, BGQ_BS$contamination, method="games-howell")

oneway.test(K_disp ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$K_disp, BGQ_SSR$contamination, method="games-howell")

#Cu
# ANOVA-like test (Welch ANOVA)
oneway.test(Cu ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$Cu, BGQ_UC$sample, method="games-howell")

oneway.test(Cu ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$Cu, BGQ_MC$sample, method="games-howell")

oneway.test(Cu ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$Cu, BGQ_HC$sample, method="games-howell")

oneway.test(Cu ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$Cu, BGQ_BS$contamination, method="games-howell")

oneway.test(Cu ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$Cu, BGQ_SSR$contamination, method="games-howell")

#Cu_sol
# ANOVA-like test (Welch ANOVA)
oneway.test(Cu_sol ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$Cu_sol, BGQ_UC$sample, method="games-howell")

oneway.test(Cu_sol ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$Cu_sol, BGQ_MC$sample, method="games-howell")

oneway.test(Cu_sol ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$Cu_sol, BGQ_HC$sample, method="games-howell")

oneway.test(Cu_sol ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$Cu_sol, BGQ_BS$contamination, method="games-howell")

oneway.test(Cu_sol ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$Cu_sol, BGQ_SSR$contamination, method="games-howell")

# ANOVA-like test (Welch ANOVA)
oneway.test(pCu ~ sample, data = BGQ_UC, var.equal = FALSE)
# Games-Howell post hoc test
posthocTGH(BGQ_UC$pCu, BGQ_UC$sample, method="games-howell")

oneway.test(pCu ~ sample, data = BGQ_MC, var.equal = FALSE)
posthocTGH(BGQ_MC$pCu, BGQ_MC$sample, method="games-howell")

oneway.test(pCu ~ sample, data = BGQ_HC, var.equal = FALSE)
posthocTGH(BGQ_HC$pCu, BGQ_HC$sample, method="games-howell")

oneway.test(pCu ~ contamination, data = BGQ_BS, var.equal = FALSE)
posthocTGH(BGQ_BS$pCu, BGQ_BS$contamination, method="games-howell")

oneway.test(pCu ~ contamination, data = BGQ_SSR, var.equal = FALSE)
posthocTGH(BGQ_SSR$pCu, BGQ_SSR$contamination, method="games-howell")


