# 
# ---- Environment ----
setwd("/home/jana/Documents/PhD/Projects/inProgress/AlphaPart/Revision3/SimulationExample")

library(package = "AlphaSimR")
library(package = "tidyverse")
library(package = "AlphaPart")
library(package = "tidyr")


# ---- General parameters ----
# Set the number of males and females selected in each population
# This simulation includes only selection on males
nMales_pop1   =  5
nFemales_pop1 = 500
nMales_pop2   =  10
nFemales_pop2 = 1000
nMales_pop3   =  10
nFemales_pop3 = 2000
# Set the number of burn-in and evaluation generations
nGenerationBurn = 10
nGenerationEval = 10

# Set the traits to record - we are simulating three correlated traits
GenMeanCols = c("GenMeanT1","GenMeanT2", "GenMeanT3")
GenVarCols  = c("GenVarT1", "GenVarT2", "GenVarT3")

# Simulate the founder population
founderPop = runMacs(nInd = (nMales_pop1+nMales_pop2+nMales_pop3+nFemales_pop1+nFemales_pop2+nFemales_pop3),
                     nChr = 10,
                     segSites = 1000,
                     nThreads = 4,
                     species = "GENERIC")
                     #species = "CATTLE")

save.image(file = "~/Documents/PhD/Projects/inProgress/AlphaPart/Revision3/SimulationExample/FounderPop_generic.R")
load("~/Documents/PhD/Projects/inProgress/AlphaPart/Revision3/SimulationExample/FounderPop_generic.R")
###################################################################################
###################################################################################
# ---- Simulation/Base population parameters ----
SP = SimParam$new(founderPop)
# Set the genetic variance with a high correlation between traits
VarA = matrix(data = c(1.0, 0.9, 0.8,
                       0.9, 1.0, 0.8,
                       0.8, 0.8, 1.0), nrow = 3); cov2cor(VarA)
# Set the environmental variances to achieve 0.7 and 0.9 h2 for the traits
VarE = matrix(data = c(1.0/0.7-1, 0.0, 0.0, 
                       0.0, 1.0/0.9-1, 0.0,
                       0.0, 0.0, 1.0/0.8-1), nrow = 3); cov2cor(VarE)
# Set the phenotypic variance and check the h2
VarP = VarA + VarE; diag(VarA) / diag(VarP)
# Add an additive trait to the population with the var-cov structure
SP$addTraitA(nQtlPerChr = 1000, mean = c(0, 0, 0), var = diag(VarA), cor = cov2cor(VarA))
# Set gender as systematic (half F - half M)
SP$setGender(gender = "yes_sys")

# ---- Base GN population ----
Base = newPop(founderPop)[1:(nMales_pop1+nMales_pop2+nMales_pop3+nFemales_pop1+nFemales_pop2+nFemales_pop3)]
Base@gender <- sample(c(rep("F", nFemales_pop1+nFemales_pop2+nFemales_pop3), rep("M", nMales_pop1+nMales_pop2+nMales_pop3)))
# Here select the base females and males for both populations - pop1 and pop2
BaseMales_pop1   = Base[Base@gender == "M"][1:nMales_pop1]
BaseFemales_pop1 = Base[Base@gender == "F"][1:nFemales_pop1]
BaseMales_pop2   = Base[Base@gender == "M"][(nMales_pop1+1):(nMales_pop1+nMales_pop2)]
BaseFemales_pop2 = Base[Base@gender == "F"][(nFemales_pop1+1):(nFemales_pop1 + nFemales_pop2)]
BaseMales_pop3   = Base[Base@gender == "M"][(nMales_pop1+nMales_pop2+1):(nMales_pop1 + nMales_pop2 + nMales_pop3)]
BaseFemales_pop3 = Base[Base@gender == "F"][(nFemales_pop1+nFemales_pop2+1):(nFemales_pop1 + nFemales_pop2 + nFemales_pop3)]
# Check the number of individuals in each population
BaseFemales_pop1@nInd
BaseFemales_pop2@nInd
BaseMales_pop1@nInd
BaseMales_pop2@nInd
sum(BaseFemales_pop1@id %in% BaseFemales_pop2@id)
sum(BaseMales_pop1@id %in% BaseMales_pop2@id)
sum(BaseMales_pop2@id %in% BaseMales_pop3@id)
rm(Base)
  
###################################################################################
###################################################################################

# ---- Burn-in ----
# In the burn-in you create two populations - pop1 and pop2
# Select only the males - pop1 on trait 1 and pop2 on trait 2 (r=0.9)
DataBurn = data.frame(Generation = rep(1:nGenerationBurn, each=6), 
                      Gender = rep(c("F", "M"), nGenerationBurn*3),
                      Population = rep(c("Pop1", "Pop1", "Pop2", "Pop2", "Pop3", "Pop3"), nGenerationBurn),
                      GenMeanT1 = NA, GenMeanT2 = NA, GenMeanT3 = NA, 
                      GenVarT1  = NA, GenVarT2  = NA, GenVarT3  = NA)
PedBurnIn <- tibble()

for (Generation in 1:nGenerationBurn) {
  for (pop in c("Pop1", "Pop2", "Pop3")) {
  
    # Set breeding individuals according to the population
    if (pop == "Pop1") {
      BaseFemales = BaseFemales_pop1
      BaseMales = BaseMales_pop1
    } else if (pop == "Pop2") {
      BaseFemales = BaseFemales_pop2
      BaseMales = BaseMales_pop2
    } else if (pop == "Pop3") {
      BaseFemales = BaseFemales_pop3
      BaseMales = BaseMales_pop3
    }
    
    # Mate the individuals
    SelCand = randCross2(females = BaseFemales, males = BaseMales,
                         nCrosses = BaseFemales@nInd, nProgeny = 2)
    
    # Save metrics
    for (gender in c("F", "M")) {
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender) & (DataBurn$Pop == pop), GenMeanCols] =
        colMeans(SelCand@gv[SelCand@gender == gender,])
      DataBurn[(DataBurn$Generation == Generation) & (DataBurn$Gender == gender) & (DataBurn$Pop == pop), GenVarCols] =
        diag(var(SelCand@gv[SelCand@gender == gender,]))
    }
    
    
    # Phenotype
    SelCand = setPheno(pop = SelCand, varE = VarE)
    # if (Generation == 1) {
    #   VarA <- varG(SelCand)
    #   diag(VarE) <- diag(varP(SelCand)) - diag(varG(SelCand))
    # }
    # 
    # Track the pedigree and related info
    PedBurnIn = rbind(PedBurnIn,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "BurnIn",
                           Population = pop,
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           PhenoT3    = SelCand@pheno[,3],
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2],
                           TbvT3      = SelCand@gv[, 3]))
    
    # Select new parents according to the population
    if (pop == "Pop1") {
      BaseMales_pop1   = selectInd(pop = SelCand, nInd = nMales_pop1,   gender = "M",
                              use = "pheno", trait = 1)
      BaseFemales_pop1 = selectInd(pop = SelCand, nInd = nFemales_pop1, gender = "F")
    } else if (pop == "Pop2") {
      BaseMales_pop2   = selectInd(pop = SelCand, nInd = nMales_pop2,   gender = "M",
                                   use = "pheno", trait = 2)
      BaseFemales_pop2 = selectInd(pop = SelCand, nInd = nFemales_pop2, gender = "F")
    } else if (pop == "Pop3") {
      BaseMales_pop3   = selectInd(pop = SelCand, nInd = nMales_pop3,   gender = "M",
                                   use = "pheno", trait = 3)
      BaseFemales_pop3 = selectInd(pop = SelCand, nInd = nFemales_pop3, gender = "F")
    }

  }
}

table(PedBurnIn$Generation, PedBurnIn$Population)
table(PedBurnIn$Gender, PedBurnIn$Population)
# Plot genetic means
DataBurn %>%
gather(key = "Metric", value = "Value", GenMeanCols) %>%
ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(rows = vars(Population), cols = vars(Gender)) + 
  ylab(label = "Genetic mean")

# Plot only the one correlated trait
OneTrait <- DataBurn %>% gather(key = "Metric", value = "Value", GenMeanCols)
rbind(OneTrait[OneTrait$Population == "Pop1" & OneTrait$Metric == "GenMeanT1",], 
      OneTrait[OneTrait$Population == "Pop2" & OneTrait$Metric == "GenMeanT2",],
      OneTrait[OneTrait$Population == "Pop3" & OneTrait$Metric == "GenMeanT3",]) %>% 
ggplot(., aes(Generation, Value, color = Metric)) +
geom_line() + facet_grid(cols = vars(Gender)) + 
ylab(label = "Genetic mean")


# Plot genetic variances
DataBurn %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
ggplot(., aes(Generation, Value, color = Metric)) +
geom_line() +  facet_grid(rows = vars(Population), cols = vars(Gender)) + 
ylab(label = "Genetic variance")


# ---- Import ----

PedEval <- PedBurnIn
DataEval = data.frame(Generation = rep((nGenerationBurn+1):(nGenerationBurn + nGenerationEval), each=6), 
                      Gender = rep(c("F", "M"), nGenerationEval*3),
                      Population = rep(c("Pop1", "Pop1", "Pop2", "Pop2", "Pop3", "Pop3"), nGenerationBurn),
                      GenMeanT1 = NA, GenMeanT2 = NA, GenMeanT3 = NA, 
                      GenVarT1  = NA, GenVarT2  = NA, GenVarT3  = NA)

accuracies <- data.frame(Population = NA, Generation = NA, Trait = NA, Cor = NA)
##############################################################################3
##############################################################################3
import = 0.2
for (Generation in (1 + nGenerationBurn):(nGenerationEval + nGenerationBurn)) {
  for (pop in c("Pop1", "Pop2", "Pop3")) {

    # If this if the first generation of evaluation, select from burn-in animals
    if (Generation == (1 + nGenerationBurn)) {
      # If this is population 1, create a mixture of pop1, pop2 and pop3 fathers
      if (pop == "Pop1") {
        Females = BaseFemales_pop1
        Males = c(BaseMales_pop1, BaseMales_pop2, BaseMales_pop3)
        
        # If pop == Pop1, mate with a mating plan - the males contains half the import % of males from Pop2 and half the import % of males from Pop3
        # If pop == Pop2, randomly mate
        matingPlan = cbind(rep(Females@id, 2), 
                           c(sample(BaseMales_pop1@id, size = Females@nInd*2*(1-import), replace=T), 
                             sample(BaseMales_pop2@id, size = Females@nInd*2*(import/2), replace=T),
                             sample(BaseMales_pop3@id, size = Females@nInd*2*(import/2), replace=T)
                           ))
        SelCand = makeCross2(females = Females, males = Males, crossPlan = matingPlan, 
                             nProgeny = 2)
      } else if (pop == "Pop2") {
        #If this is Pop2, just select the parents
        Females = BaseFemales_pop2
        Males = BaseMales_pop2
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      } else if (pop == "Pop3") {
        #If this is Pop2, just select the parents
        Females = BaseFemales_pop3
        Males = BaseMales_pop3
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      }
    }
    
    #I this is not the first generation of evaluation, select from previous round selected parens
    if (Generation > (1 + nGenerationBurn)) {
      #If this is pop 1, create a mixture of pop1 and pop2 fathers
      if (pop == "Pop1") {
        Females = Females_pop1
        Males = c(Males_pop1, Males_pop2, Males_pop3)
        
        # If pop == Pop1, mate with a mating plan - the males contain the import % of males from Pop2
        # If pop == Pop2, randomly mate
        matingPlan = cbind(rep(Females@id, 2), 
                           c(sample(Males_pop1@id, size = Females@nInd*2*(1-import), replace=T), 
                             sample(c(Males_pop2@id, Males_pop3@id), size = Females@nInd*2*import, replace=T)))
        SelCand = makeCross2(females = Females, males = Males, crossPlan = matingPlan, 
                             nProgeny = 2)
        
      } else if (pop == "Pop2") {
        Females = Females_pop2
        Males = Males_pop2
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      } else if (pop == "Pop3") {
        Females = Females_pop3
        Males = Males_pop3
        SelCand = randCross2(females = Females, males = Males, nCrosses = Females@nInd,
                             nProgeny = 2)
      }
    }

    # meanG(SelCand)
    
    # Save metrics
    for (gender in c("F", "M")) {
      DataEval[(DataEval$Generation == Generation) & (DataEval$Gender == gender) & (DataEval$Pop == pop), GenMeanCols] =
        colMeans(SelCand@gv[SelCand@gender == gender,])
      DataEval[(DataEval$Generation == Generation) & (DataEval$Gender == gender) & (DataEval$Pop == pop), GenVarCols] =
        diag(var(SelCand@gv[SelCand@gender == gender,]))
    }
    
    # Phenotype
    SelCand = setPheno(pop = SelCand, varE = VarE)

    # Track pedigree
    PedEval = rbind(PedEval,
                    tibble(Generation = Generation,
                           IId        = SelCand@id,
                           FId        = SelCand@father,
                           MId        = SelCand@mother,
                           Gender     = SelCand@gender,
                           Program    = "Eval",
                           Population = pop,
                           PhenoT1    = SelCand@pheno[,1],
                           PhenoT2    = SelCand@pheno[,2],
                           PhenoT3    = SelCand@pheno[,3],
                           TbvT1      = SelCand@gv[, 1],
                           TbvT2      = SelCand@gv[, 2],
                           TbvT3      = SelCand@gv[, 3]))
   
    
    # Compute Accuracies (pheno-TGV)
    accuracies <- rbind(accuracies, c(pop, Generation, 1, cor(SelCand@pheno[,1], SelCand@gv[,1]) ))
    accuracies <- rbind(accuracies, c(pop, Generation, 2, cor(SelCand@pheno[,2], SelCand@gv[,2]) ))
    accuracies <- rbind(accuracies, c(pop, Generation, 3, cor(SelCand@pheno[,3], SelCand@gv[,3]) ))
    
    # Select
    if (pop == "Pop1") {
      Males_pop1   = selectInd(pop = SelCand, nInd = nMales_pop1,   gender = "M",
                                   use = "pheno", trait = 1)
      Females_pop1 = selectInd(pop = SelCand, nInd = nFemales_pop1, gender = "F")
    } else if (pop == "Pop2") {
      Males_pop2   = selectInd(pop = SelCand, nInd = nMales_pop2,   gender = "M",
                                   use = "pheno", trait = 2)
      Females_pop2 = selectInd(pop = SelCand, nInd = nFemales_pop2, gender = "F")
    } else if (pop == "Pop3") {
      Males_pop3   = selectInd(pop = SelCand, nInd = nMales_pop3,   gender = "M",
                                   use = "pheno", trait = 3)
      Females_pop3 = selectInd(pop = SelCand, nInd = nFemales_pop3, gender = "F")
    }
    # Clean
    rm(SelCand)
  }
}

    
# Plot genetic means
rbind(DataBurn, DataEval) %>%
  gather(key = "Metric", value = "Value", GenMeanCols) %>%
  ggplot(., aes(Generation, Value, color = Population)) +
  geom_line() + facet_grid(rows = vars(Metric), cols = vars(Gender)) + 
  ylab(label = "Genetic mean")


# Plot only the one correlated trait
OneTrait <- rbind(DataBurn, DataEval) %>% gather(key = "Metric", value = "Value", GenMeanCols)
rbind(OneTrait[OneTrait$Population == "Pop1" & OneTrait$Metric == "GenMeanT1",], 
      OneTrait[OneTrait$Population == "Pop2" & OneTrait$Metric == "GenMeanT2",],
      OneTrait[OneTrait$Population == "Pop3" & OneTrait$Metric == "GenMeanT3",]) %>% 
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(cols = vars(Gender)) + 
  ylab(label = "Genetic mean")

# Plot genetic variances
DataEval %>% gather(key = "Metric", value= "Value", GenVarCols) %>%
  ggplot(., aes(Generation, Value, color = Metric)) +
  geom_line() + facet_grid(rows = vars(Population), cols = vars(Gender)) + 
ylab(label = "Genetic variance")  


#PedEval <- pedSetBase(PedEval, keep = PedEval$Generation > 9, colId = "IId", colFid = "FId", colMid = "MId")
PedEval$Tbv <- ifelse(PedEval$Population == "Pop1", PedEval$TbvT1, ifelse(PedEval$Population == "Pop2", PedEval$TbvT2, PedEval$TbvT3))
#PedEval$PopGender <- paste(PedEval$Pop, PedEval$Gender, sep="_")
library(gridExtra)
library(grid)
png("InputDF.png", res=600, width=85, height=40, units="mm")
p<-gridExtra::tableGrob(head(as.data.frame(PedEval[, c("Generation", "IId", "FId", "MId", "Population", "Tbv")])), 
                        theme=ttheme_minimal(base_size = 8),
                        rows=NULL
                        )
separators <- replicate(ncol(p)-1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
## add vertical lines on the left side of columns (after 2nd)
p <- gtable::gtable_add_grob(p, grobs = separators,
                             t = 2, b = nrow(p), l = seq_len(ncol(p)-1)+1)
gridExtra::grid.arrange(p)
dev.off()


Part = AlphaPart(x = as.data.frame(PedEval[, c("Generation", "IId", "FId", "MId", "Population", "Tbv")]),
                  colId = "IId", colFid = "FId", colMid = "MId",
                  colPath = "Population", colBV = "Tbv")
#############
#figures for paper
png("PartDF.png", res=600, width=170, height=40, units="mm")
p<-gridExtra::tableGrob(head(Part$Tbv), 
                        theme=ttheme_minimal(base_size = 8),
                        rows=NULL
)
separators <- replicate(ncol(p)-1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
## add vertical lines on the left side of columns (after 2nd)
p <- gtable::gtable_add_grob(p, grobs = separators,
                             t = 2, b = nrow(p), l = seq_len(ncol(p)-1)+1)
gridExtra::grid.arrange(p)
dev.off()
####################
sumPart <- summary(Part, by="Generation")
#plot(sumPart)

Part1 = Part
Part1$value <- Part$value[Part$value$Pop == "Pop1",]
#Part1$TbvT1 <- Part$TbvT1[Part$TbvT1$Pop == "Pop1",]
#Part1$TbvT2 <- Part$TbvT2[Part$TbvT2$Pop == "Pop1",]
#Part1$TbvT3 <- Part$TbvT3[Part$TbvT3$Pop == "Pop1",]
Part2 = Part
Part2$value <- Part$value[Part$value$Pop == "Pop2",]
# Part2$TbvT1 <- Part$TbvT1[Part$TbvT1$Pop == "Pop2",]
# Part2$TbvT2 <- Part$TbvT2[Part$TbvT2$Pop == "Pop2",]
# Part2$TbvT3 <- Part$TbvT3[Part$TbvT3$Pop == "Pop2",]
Part3 = Part
Part3$value <- Part$value[Part$value$Pop == "Pop3",]
# Part3$TbvT1 <- Part$TbvT1[Part$TbvT1$Pop == "Pop3",]
# Part3$TbvT2 <- Part$TbvT2[Part$TbvT2$Pop == "Pop3",]
# Part3$TbvT3 <- Part$TbvT3[Part$TbvT3$Pop == "Pop3",]

Part1Summary = summary(object = Part1, by = "Generation")
#############
#figures for paper
png("ParSumtDF.png", res=600, width=85, height=40, units="mm")
p<-gridExtra::tableGrob(head(round(Part1Summary$Tbv,5)), 
                        theme=ttheme_minimal(base_size = 8),
                        rows=NULL
)
separators <- replicate(ncol(p)-1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
## add vertical lines on the left side of columns (after 2nd)
p <- gtable::gtable_add_grob(p, grobs = separators,
                             t = 2, b = nrow(p), l = seq_len(ncol(p)-1)+1)
gridExtra::grid.arrange(p)
dev.off()
####################
Part2Summary = summary(object = Part2, by = "Generation")
Part3Summary = summary(object = Part3, by = "Generation")


plot(summary(Part, by="Generation"))
plot(Part1Summary)
plot(Part2Summary)
plot(Part3Summary)


PartFix = pedFixBirthYear(x = as.data.frame(PedEval[, c("Generation", "IId", "FId", "MId", "Population", "Tbv")]),
                 interval=1,
                 colId = "IId", colFid = "FId", colMid = "MId")

PedEval <- as.data.frame(PedEval[, c("Generation", "IId", "FId", "MId", "Population", "Tbv")])
PedRebase = pedSetBase(x = as.data.frame(PedEval[, c("Generation", "IId", "FId", "MId", "Population", "Tbv")]),
                          keep = PedEval$Generation > 10,
                          colId = "IId", colFid = "FId", colMid = "MId")


PedEvalSubset <- AlphaPartSubset(Part1Summary, paths = c("Pop1", "Pop2"))
PedEvalSum <- AlphaPartSum(Part1Summary, map=list(c("Domestic", "Pop1"), c("Import", "Pop2", "Pop3")))
plot(PedEvalSum)
