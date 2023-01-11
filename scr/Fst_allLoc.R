# LOAD LIBRARIES & FUNCTIONS ----

source("scr/libraries.R")
source("scr/xtrafunctions.R")
source("scr/ggplot.R")
source("scr/genind.R")

library(assigner)
library(radiator)

# set up to run in parallel ====
library(parallel)
library(foreach)
library(doMC)

# set number of cores to run in parallel
registerDoMC(20)

# set how numbers are printed
options(scipen=999)


# LOAD & FORMAT DATA SET ----

# import genotypes as genind ====

gen <- read.genepop(file = "data/POPGEN/BMA_by_pop_genepop.gen",
                    ncode = 3L, quiet = FALSE)

Inds <- as.data.frame(indNames(gen)) %>%
  rename(LIB_ID = `indNames(gen)`) %>%
  separate(LIB_ID, into = c("SP", "PLATE", "WELL", "SAMPLE_ID"), sep = "-", remove = FALSE, extra = "merge") %>%
  select(-SP)

# Load as strata ====
SampleInfo <- read_delim("data/POPGEN/SampleInfo.txt", delim = "\t")

strata <- left_join(Inds, SampleInfo) %>%
  distinct() %>%
  mutate(POP = ordered(POP, levels = pops)) %>%
  mutate(REGION = case_when(POP %in% c("FLA") ~ "SWATL",
                            POP %in% c("CAMP") ~ "SGULF",
                            POP %in% c("FLGS", "FLGN") ~ "EGULF",
                            POP %in% c("MB", "MISS", "CS", "LA") ~ "CGULF",
                            POP %in% c("CC") ~ "WGULF",),
         REGION = ordered(REGION, levels = reg),
         OCEAN = ifelse(POP == "FLA", "ATL", "GULF"))

strata(gen) <- strata

# define populations using defined stratification
setPop(gen) <- ~POP


# create neutral and outlier data sets ====
outl <- list()

outl[[1]] <- read_delim("results/arlequin01.outlier", delim = "\t")
outl[[2]] <- read_delim("results/bayescan.pop.outlier", delim = "\t")

outl <- ldply(outl) %>%
  mutate(LOCUS = str_trim(LOCUS, side = "both")) %>%
  distinct()

# neutral data set
removeloc <- outl$LOCUS

gen_n <- genind.rem.loci(gen, removeloc)

# outlier data set
keeploc <- outl$LOCUS
gen_o <- gen[loc = keeploc]

rm(SampleInfo)


# CALCULATE GLOBAL FST ---- 

# group genetic data for comparison
setPop(gen) <- ~POP

pop <- popNames(gen)

# calculate global Fst and CIs ====
tidy <- tidy_genomic_data(data = gen, filename = NULL, parallel = 65)

fst <- assigner::fst_WC84(data = tidy, 
                          pop.levels = pop,
                          holdout.samples = NULL,
                          pairwise = TRUE,
                          ci = TRUE, 
                          iteration.ci = 1000,
                          quantiles.ci = c(0.025, 0.975), 
                          digits = 9,
                          verbose = TRUE)

# write results
df <- fst$fst.overall

write_delim(df, "results/estuaries_allLoc.globalfst.ci", delim = "\t")

# global fst and permuted p-values ==== 

dat <- genind2hierfstat(gen)
dat_rand <- dat

# number of permutations
nperm <- 1000

# create vector with locus names
loc <- locNames(gen)

# calculate F-statistics for genotypes permuted betw pop ====

fst.glob.sim <- foreach(1:nperm) %dopar%  {
  
  dat_rand$pop <- sample(dat_rand$pop, replace = FALSE)
  
  wc(dat_rand)
  
}

# calculate F-statistics for empirical data ====
fst.glob.obs <- wc(dat)

# parse permuted Fst (all loci) ====
sim_fst <- list()

for(i in 1:nperm){
  
  sim_fst[[i]] <- as.data.frame(fst.glob.sim[[i]]$FST) %>%
    rename(SIM_FST = `fst.glob.sim[[i]]$FST`)
  
}

# calculate p-values ====
sim_fst <- ldply(sim_fst, data.frame)

obs_fst <- fst.glob.obs$FST

larger <- filter(sim_fst, SIM_FST > obs_fst)
larger <- nrow(larger)

pval <- larger/nperm

# write to file ====
STAT <- c("OBS_FST", "PVAL", "NPERM")
VALUE <- c(obs_fst, pval, nperm)

results <- data.frame(STAT, VALUE)

write_delim(results, "results/estuary_allLoc.globalfst.pval", delim = "\t")


# PAIRWISE FST OCEANS ----

# set groups to compare ====
setPop(gen) <- ~OCEAN

pop <- popNames(gen)

# number of groups being compared
n <- length(pop)

# calculate pairwise Fst & CIs ====
tidy <- tidy_genomic_data(data = gen, filename = NULL, parallel = 10)

fst <- assigner::fst_WC84(data = tidy, 
                          pop.levels = pop,
                          holdout.samples = NULL,
                          pairwise = TRUE,
                          ci = TRUE, 
                          iteration.ci = 10000,
                          quantiles.ci = c(0.025, 0.975), 
                          digits = 9,
                          verbose = TRUE,
                          parallel.core = 55)

# write results
df <- fst$pairwise.fst

write_delim(df, "results/ocean_allLoc.fst.ci", delim = "\t")


# compute Fst matrix ====

dat <- genind2hierfstat(gen)
mat.obs <- pairwise.WCfst(dat)

temp <- as.data.frame(mat.obs) %>%
  rownames_to_column("GRP1") %>%
  gather(key = "GRP2", value = "obsFst", 2:(n+1)) %>%
  filter(GRP1 != GRP2)

write_delim(temp, "results/ocean_allLoc.fst.WC84", delim = "\t")

# calculate pairwise Fst for individuals permuted between groups ====

# create list with NPERM matrices of permuted Fst values
NBPERM <- 1000

# permute individuals between groups for each pairwise comparison
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.WCfst(mutate(dat, pop = sample(pop, replace = FALSE))), mc.cores = 55)

# create data frame with permuted values
fst_perm <- list()

for(i in 1:length(mat.perm)){
  
  fst <- as.data.frame(mat.perm[[i]]) %>%
    rownames_to_column("GRP1") %>%
    gather(key = "GRP2", value = "ppFST", 2:(n+1)) %>%
    filter(GRP1 != GRP2)
  
  fst_perm[[i]] <- fst
  
}

fst_perm <- ldply(fst_perm, data.frame)

write_delim(fst_perm, "results/ocean_allLoc.fst.perm", delim = "\t")

# get p-values for each pairwise comparison ====

# use randtest to determine p-value (i.e. is observed value different from permuted values)
# p.globs.p<-sum(gglobs.p>=gglobs.p[nperm+1])/(nperm+1)  p-val is sum(times observed value is > permuted value / total permutations)

ppfst_pval <- list()

for(i in 1:(nrow(mat.obs)-1)){
  
  for(j in 2:nrow(mat.obs)){
    
    ppfst_pval[[paste(rownames(mat.obs)[i], rownames(mat.obs)[j], sep = "-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter = "greater")
    
  }
  
}

# create data frame with p-values
COMP <- names(ppfst_pval)

PVAL <- rep(NA, length(ppfst_pval))

for (i in 1:length(PVAL)) {
  
  PVAL[i] <- ppfst_pval[[i]]$pvalue
  
}

pval <- data.frame(COMP, PVAL)

write_delim(pval, "results/ocean_allLoc.fst.pval", delim = "\t")


# PAIRWISE FST REGIONS ----

# set groups to compare ====
setPop(gen) <- ~REGION

pop <- popNames(gen)

# number of groups to compare
n <- length(pop)


# pairwise Fst and CIs ====
tidy <- tidy_genomic_data(data = gen, filename = NULL, parallel = 10)

fst <- assigner::fst_WC84(data = tidy, 
                          pop.levels = pop,
                          holdout.samples = NULL,
                          pairwise = TRUE,
                          ci = TRUE, 
                          iteration.ci = 10000,
                          quantiles.ci = c(0.025, 0.975), 
                          digits = 9,
                          verbose = TRUE,
                          parallel.core = 55)

df <- fst$pairwise.fst

write_delim(df, "results/region_allLoc.fst.ci", delim = "\t")


# compute Fst matrix ====
dat <- genind2hierfstat(gen)
mat.obs <- pairwise.WCfst(dat)

temp <- as.data.frame(mat.obs) %>%
  rownames_to_column("GRP1") %>%
  gather(key = "GRP2", value = "obsFst", 2:(n+1)) %>%
  filter(GRP1 != GRP2)

write_delim(temp, "results/region_allLoc.fst", delim = "\t")

# calculate pairwise Fst for individuals permuted between groups ====

# create list with NPERM matrices of permuted Fst values
NBPERM <- 1000

# permute individuals between groups for each pairwise comparison
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.WCfst(mutate(dat, pop = sample(pop, replace = FALSE))), mc.cores = 55)

# create data frame with permuted values
fst_perm <- list()

for(i in 1:length(mat.perm)){
  
  fst <- as.data.frame(mat.perm[[i]]) %>%
    rownames_to_column("GRP1") %>%
    gather(key = "GRP2", value = "ppFST", 2:(n+1)) %>%
    filter(GRP1 != GRP2)
  
  fst_perm[[i]] <- fst
  
}

fst_perm <- ldply(fst_perm, data.frame)

write_delim(fst_perm, "results/region_allLoc.fst.perm", delim = "\t")

# get p-values for each pairwise comparison ====

# use randtest to determine p-value (i.e. is observed value different from permuted values)
# p.globs.p<-sum(gglobs.p>=gglobs.p[nperm+1])/(nperm+1)  p-val is sum(times observed value is > permuted value / total permutations)

ppfst_pval <- list()

for(i in 1:(nrow(mat.obs)-1)){
  
  for(j in 2:nrow(mat.obs)){
    
    ppfst_pval[[paste(rownames(mat.obs)[i], rownames(mat.obs)[j], sep = "-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter = "greater")
    
  }
  
}

# create data frame with p-values
COMP <- names(ppfst_pval)

PVAL <- rep(NA, length(ppfst_pval))

for (i in 1:length(PVAL)) {
  
  PVAL[i] <- ppfst_pval[[i]]$pvalue
  
}

pval <- data.frame(COMP, PVAL)

write_delim(pval, "results/region_allLoc.fst.pval", delim = "\t")



# PAIRWISE FST ESTUARIES ----

# set groups to compare ====
setPop(gen) <- ~POP

pop <- popNames(gen)

# number of groups to compare
n <- length(pop)


# pairwise Fst and CIs ====
tidy <- tidy_genomic_data(data = gen, filename = NULL, parallel = 20)

fst <- assigner::fst_WC84(data = tidy, 
                          pop.levels = pop,
                          holdout.samples = NULL,
                          pairwise = TRUE,
                          ci = TRUE, 
                          iteration.ci = 10000,
                          quantiles.ci = c(0.025, 0.975), 
                          digits = 9,
                          verbose = TRUE,
                          parallel.core = 50)
# pairwise fst
df <- fst$pairwise.fst

write_delim(df, "results/estuaries_allLoc.fst.ci", delim = "\t")


# compute Fst matrix ====
dat <- genind2hierfstat(gen)
mat.obs <- pairwise.WCfst(dat)

gen <- as.data.frame(mat.obs) %>%
  rownames_to_column("GRP1") %>%
  gather(key = "GRP2", value = "obsFst", 2:(n+1)) %>%
  filter(GRP1 != GRP2)

write_delim(gen, "results/estuary_allLoc.fst.WC84", delim = "\t")


# calculate pairwise Fst for individuals permuted between groups ====

# create list with NPERM matrices of permuted Fst values
NBPERM <- 1000

# permute individuals between groups for each pairwise comparison
mat.perm <- mclapply(1:NBPERM, function(i) pairwise.WCfst(mutate(dat, pop = sample(pop, replace = FALSE))), mc.cores = 55)

# create data frame with permuted values
fst_perm <- list()

for(i in 1:length(mat.perm)){
  
  fst <- as.data.frame(mat.perm[[i]]) %>%
    rownames_to_column("GRP1") %>%
    gather(key = "GRP2", value = "ppFST", 2:(n+1)) %>%
    filter(GRP1 != GRP2)
  
  fst_perm[[i]] <- fst
  
}

fst_perm <- ldply(fst_perm, data.frame)

write_delim(fst_perm, "results/estuary_allLoc.fst.perm", delim = "\t")

# get p-values for each pairwise comparison ====

# use randtest to determine p-value (i.e. is observed value different from permuted values)
# p-val is sum(times observed value is > permuted value / total permutations)

ppfst_pval <- list()

for(i in 1:(nrow(mat.obs)-1)){
  
  for(j in 2:nrow(mat.obs)){
    
    ppfst_pval[[paste(rownames(mat.obs)[i], rownames(mat.obs)[j], sep = "-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter = "greater")
    
  }
  
}

# create data frame with p-values
COMP <- names(ppfst_pval)

PVAL <- rep(NA, length(ppfst_pval))

for (i in 1:length(PVAL)) {
  
  PVAL[i] <- ppfst_pval[[i]]$pvalue
  
}

pval <- data.frame(COMP, PVAL)

write_delim(pval, "results/estuary_allLoc.fst.pval", delim = "\t")


