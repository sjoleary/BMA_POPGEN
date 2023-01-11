## Converting a genind object from R to a DIYABC object with the SNPs intact##
#gen.net is my neutral data genind object
#haplotyping/codes.SNP.TRS.F10.gen is my data from the last haplotyping of my data before I read it into R for my analysis

# Read in data set

gen.net <- read.genepop(file = "data/POPGEN/BMA_by_pop_thin_neutral_snps_genepop.gen",
                    ncode = 3L, quiet = FALSE)

Inds <- as.data.frame(indNames(gen)) %>%
  rename(LIB_ID = `indNames(gen)`) %>%
  separate(LIB_ID, into = c("SP", "PLATE", "WELL", "SAMPLE_ID"), sep = "-", remove = FALSE, extra = "merge") %>%
  select(-SP)

# Import sample meta-data
SampleInfo <- read_delim("data/POPGEN/SampleInfo.txt", delim = "\t")

strata <- left_join(Inds, SampleInfo) %>%
  distinct() %>%
  mutate(POP = ordered(POP, levels = pops)) %>%
  mutate(ESTUARY = ordered(ESTUARY, levels = est_levels),
         REGION = case_when(POP %in% c("FLA") ~ "SWATL",
                            POP %in% c("CAMP") ~ "SGULF",
                            POP %in% c("FLGS", "FLGN") ~ "EGULF",
                            POP %in% c("MB", "MISS", "CS", "LA") ~ "CGULF",
                            POP %in% c("CC") ~ "WGULF",),
         REGION = ordered(REGION, levels = reg),
         OCEAN = ifelse(POP == "FLA", "ATL", "GULF"))

strata(gen.net) <- strata


# Read in the sequence data
seqs <- read.table("data/Haplotyping/codes.BMA.gen", head=F)

keep <- as.data.frame(locNames(gen.net)) %>%
  separate(`locNames(gen.net)`, into = c("dDocent", "Contig", "Locus", "SNP"), sep = "_") %>%
  unite(LOCUS, dDocent, Contig, Locus, sep = "_")

seqs <- seqs[seqs$V1 %in% keep$LOCUS,]

# Define variables
indv <- indNames(gen.net)
loci <- as.character(seqs[,1])

# Make a list as long as the number of loci
seq.list <- vector("list",nrow(seqs))

# Pull the allele number and variable DNA sequences for each locus
for(i in 1:nrow(seqs)){
  
  seqs2 <- data.frame(matrix(unlist(strsplit(as.character(gsub(":", ",", seqs[i,2])),",")),ncol=2, byrow=T))
  colnames(seqs2) <- c("Seq", "Allele")
  seq.list[[i]] <- seqs2
  
}

seq.list <- setNames(seq.list, as.character(seqs[,1]))

# Make the final data table
final <- matrix(nrow=length(indv), ncol=2*length(loci))
rownames(final)<-indNames(gen.net)
colnames(final)<-rep(locNames(gen.net), each=2)

#Fill in the final table
rm.loci <- vector()

for(j in loci){
  
  geno <- gen.net@tab[ , grep(paste(j,".",sep=""),colnames(gen.net@tab), fixed=T)]
  
  if(class(geno)=="integer"){final <- final[,-(grep(j,colnames(final)))]; print(paste("Locus", j, "removed because it is monomorphic")); rm.loci <- c(rm.loci,j); next}
  colnames(geno) <- gsub(paste(j,".",sep=""), "", colnames(geno))
  
  print(paste("processing locus", j))
  
  for(i in indv){
    
    if(sum(is.na(geno[rownames(geno)==i,]))>0){next}
    
    tmp <- geno[rownames(geno)==i,]
    tmp <- tmp[tmp > 0]
    
    if(length(tmp) == 1){final[rownames(final)==i,colnames(final)==j] <- rep(as.character(seq.list[[j]][seq.list[[j]][,2]==names(tmp),1]),2)
    
    } else if(length(tmp) == 2){a <- as.character(seq.list[[j]][seq.list[[j]][,2]==names(tmp)[1],1])
    
    b<-as.character(seq.list[[j]][seq.list[[j]][,2]==names(tmp)[2],1])
    
    final[rownames(final)==i,colnames(final)==j] <- c(a,b)}
  }}

# Redefining loci
loci.end <- unique(colnames(final))

#Define populations

#If your strata file has the population names, then the top part of this code is not necessary
setPop(gen.net) <- ~POP

# table(pop(gen.net), gen.net@strata$POP)
# EAST <- indNames(gen.net)[pop(gen.net)==1]
# ATL <- indNames(gen.net)[pop(gen.net)==2]
# WEST <- indNames(gen.net)[pop(gen.net)==3]

#This part is necessary
#This could probably be written as a loop to run over all the haplotypes, but I only had 3 groups so hard coding it made sense

for (p in strata$POP){
  
  
  
  
  write.table(final[rownames(final) %in% p,], glue("data/POPGEN/{p}_haplotypes.txt"), 
              sep="\t", quote=F, col.names=F, row.names=T)
  
}


```{bash}

#Title the DIYABC file
#This includes the sex ratio
echo "Angel Shark DIYABC by KMEANS < NM = 1.0NF >" > DIYABC.file

#Add the loci to the DIYABC file
awk 'BEGIN{OFS="\t"}{print $0, "<A>"}' loci.txt >> DIYABC.file

#Add the 1st population to the DIYABC file
#It might be possible to loop this across each population, but I didn't take the time to try
echo "Pop" >> DIYABC.file
awk '{for (i = 2; i <= NF; i=i+2) print "<["$i"]["$(i+1)"]> ";}' WEST_haplotypes.txt > WEST.tmp

for i in $(seq 1 $(wc -l WEST_haplotypes.txt| cut -d " " -f 1)); do
echo $i
LOCI=$(wc -l loci.txt | cut -f1 -d" ")
TAIL=$(echo "$i*$LOCI" | bc)
head -n $TAIL WEST.tmp | tail -n $LOCI | sed -z 's:\n:\t:g' > tmp.geno
awk -v var=$i 'NR==var{print $1, ","}' WEST_haplotypes.txt > tmp.indv
cat tmp.indv tmp.geno | sed -z 's:\n:\t:g' >> DIYABC.file
done

#Add the 2nd population to the DIYABC file
echo -e "\nPop" >> DIYABC.file
awk '{for (i = 2; i <= NF; i=i+2) print "<["$i"]["$(i+1)"]> ";}' EAST_haplotypes.txt > EAST.tmp

for i in $(seq 1 $(wc -l EAST_haplotypes.txt| cut -d " " -f 1)); do
echo $i
LOCI=$(wc -l loci.txt | cut -f1 -d" ")
TAIL=$(echo "$i*$LOCI" | bc)
head -n $TAIL EAST.tmp | tail -n $LOCI | sed -z 's:\n:\t:g' > tmp.geno
awk -v var=$i 'NR==var{print $1, ","}' EAST_haplotypes.txt > tmp.indv
cat tmp.indv tmp.geno | sed -z 's:\n:\t:g' >> DIYABC.file
done

#Add the 3rd population to the DIYABC file
echo -e "\nPop" >> DIYABC.file
awk '{for (i = 2; i <= NF; i=i+2) print "<["$i"]["$(i+1)"]> ";}' ATL_haplotypes.txt > ATL.tmp

for i in $(seq 1 $(wc -l ATL_haplotypes.txt| cut -d " " -f 1)); do
echo $i
LOCI=$(wc -l loci.txt | cut -f1 -d" ")
TAIL=$(echo "$i*$LOCI" | bc)
head -n $TAIL ATL.tmp | tail -n $LOCI | sed -z 's:\n:\t:g' > tmp.geno
awk -v var=$i 'NR==var{print $1, ","}' ATL_haplotypes.txt > tmp.indv
cat tmp.indv tmp.geno | sed -z 's:\n:\t:g' >> DIYABC.file
done

#Put each sample on its own line
#Code below should correct this without doing it by hand
sed -i -e 's/> \tWEST/>\nWEST/g' DIYABC.file
sed -i -e 's/> \tEAST/>\nEAST/g' DIYABC.file
sed -i -e 's/> \tAtl/>\nAtl/g' DIYABC.file

#Remove all NA values
sed -i 's/NA//g' DIYABC.file

# DIYABC.file is complete
