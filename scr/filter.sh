# filter LQ SNP calls
vcftools --vcf data/VCF/temp/SFL.F0.recode.vcf --out data/VCF/temp/SFL.F1 --minQ 20 --minGQ 20 --minDP 3 --mac 3 --recode --recode-INFO-all

# query stats
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --depth
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --site-mean-depth
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --missing-indv
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --missing-site
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --het
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --geno-depth
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --singletons
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --indv-freq-burden
vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/SFL.F1 --freq 

vcftools --vcf data/VCF/temp/SFL.F1.recode.vcf --out data/VCF/temp/SFL.F2a --max-missing 0.5 --min-meanDP 10 --recode --recode-INFO-all

# library SFL-4
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/SFL4 --keep data/VCF/SFL4.ind  --recode --recode-INFO-all

vcftools --vcf data/VCF/temp/SFL4.recode.vcf --out data/VCF/SFL4 --missing-site

# library SFL-5
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/SFL5 --keep data/VCF/SFL5.ind  --recode --recode-INFO-all

vcftools --vcf data/VCF/temp/SFL5.recode.vcf --out data/VCF/SFL5 --missing-site

# library SFL-6
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/SFL6 --keep data/VCF/SFL6.ind  --recode --recode-INFO-all

vcftools --vcf data/VCF/temp/SFL6.recode.vcf --out data/VCF/SFL6 --missing-site

# Alabama
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/AL --keep data/VCF/AL.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/AL.recode.vcf --out data/VCF/AL --missing-site

# Florida Atlantic
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/FLA --keep data/VCF/FLA.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/FLA.recode.vcf --out data/VCF/FLA --missing-site

# Florida Gulf
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/FLG --keep data/VCF/FLG.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/FLG.recode.vcf --out data/VCF/FLG --missing-site

# Louisiana
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/LA --keep data/VCF/LA.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/LA.recode.vcf --out data/VCF/LA --missing-site

# Mississippi
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/MS --keep data/VCF/MS.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/MS.recode.vcf --out data/VCF/MS --missing-site

# South Carolina
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/SC --keep data/VCF/SC.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/SC.recode.vcf --out data/VCF/SC --missing-site

# Texas
vcftools --vcf data/VCF/temp/SFL.F2a.recode.vcf --out data/VCF/temp/TX --keep data/VCF/TX.ind  --recode --recode-INFO-all
vcftools --vcf data/VCF/temp/TX.recode.vcf --out data/VCF/TX --missing-site
