# Script for running GWAS on PPMI and Columbia, meta-analyzing, and running COJO 
## Remove outliers, based on the R script outlier files
awk '{print $2}' NY_outliers.txt > temp.txt
sed s/'"'/''/g temp.txt | sed '1d' > temp2.txt
rev temp2.txt | cut -c 1-6 | rev > temp.txt
grep -f temp.txt ny_gcase.fam | awk '{print $1,$2}' > remove_outliers.txt
plink --bfile ny_gcase --remove remove_outliers.txt --make-bed --out ny_gcase_noOut

awk '{print $2}' PPMI_outliers.txt > temp.txt
sed s/'"'/''/g temp.txt | sed '1d' > temp2.txt
rev temp2.txt | cut -c 1-9 | rev > temp.txt
grep -f temp.txt ppmi_gcase.fam | awk '{print $1,$2}' > remove_outliers.txt
plink --bfile ppmi_gcase --remove remove_outliers.txt --make-bed --out ppmi_gcase_noOut

## General regression example - PPMI
plink --bfile ppmi_no_gcase_outl --linear --maf 0.01 --ci .95 \
 --covar ppmi_covar_new.txt \
 --covar-name Age,Sex,Pheno,E326K,N370S,T369M,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out results/ppmi_age_sex_pheno_GBA

awk '{if ($12!="NA") print $0}' results/ppmi_age_sex_pheno_GBA.assoc.linear | grep 'ADD' | sort -gk12 > results/p.ppmi_age_sex_pheno_GBA.assoc.linear

## General regression example - Columbia
plink --bfile ny_gcase_noOut --linear --maf 0.01 --ci .95 \
 --covar covar_ny.txt \
 --covar-name Sex,Age,Pheno,N370S,E326K,T369M,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out results/ny_age_sex_pheno_GBA

 awk '{if ($12!="NA") print $0}' Columbia/results/ny_age_sex_pheno_GBA.assoc.linear | grep 'ADD' | sort -gk12 > Columbia/results/p.ny_age_sex_pheno_GBA.assoc.linear

## Generate sumstats for each analysis 
# Run script to get full, metal, and cojo summary stats 
# assoc.list is just a list of the analysis names you want to meta analyze
# Run in directory above cohorts 

cat ~/runs/emsom/projects/GCase_GWAS/assoc.list | while read line 
do
    echo $line > marker.txt
    cp ~/runs/emsom/projects/GCase_GWAS/Columbia/results/p.ny_$line.assoc.linear assoc
    cp ~/runs/emsom/projects/GCase_GWAS/Columbia/allChrs_NY_addsam.info infos
    R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
    paste lambda.txt >> Columbia/results/plots/NY_lambdas.tab
    mv QQ.tiff Columbia/results/plots/QQ_NY_$line.tiff    
    mv Metal.tab Columbia/results/sumstats/Metal_NY_$line.tab
    mv COJO.tab Columbia/results/sumstats/COJO_NY_$line.tab
    mv fullSTATS.tab Columbia/results/sumstats/fullSTATS_NY_$line.tab
    mv ManH.tiff Columbia/results/plots/ManH_NY_$line.tiff
done

cat ~/runs/emsom/projects/GCase_GWAS/assoc.list | while read line 
do
    echo $line > marker.txt
    cp ~/runs/emsom/projects/GCase_GWAS/PPMI/results/p.ppmi_$line.assoc.linear assoc
    cp ~/runs/emsom/projects/GCase_GWAS/PPMI/allChrsPPMI.info infos
    R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
    paste lambda.txt >> PPMI/results/plots/PPMI_lambdas.tab
    mv QQ.tiff PPMI/results/plots/QQ_PPMI_$line.tiff    
    mv Metal.tab PPMI/results/sumstats/Metal_PPMI_$line.tab
    mv COJO.tab PPMI/results/sumstats/COJO_PPMI_$line.tab
    mv fullSTATS.tab PPMI/results/sumstats/fullSTATS_PPMI_$line.tab
    mv ManH.tiff PPMI/results/plots/ManH_PPMI_$line.tiff
done

rm QQ.tiff
rm Metal.tab 
rm COJO.tab
rm fullSTATS.tab
rm ManH.tiff

# Run script to generate annotated significant SNP files 
# keep.PDsnps.txt file contains significant SNPs from Nalls et al. 2019 PD GWAS 
mkdir Columbia/results/sigSNPs
cat ~/runs/emsom/projects/GCase_GWAS/assoc.list | while read line 
do
    head -n 1 Columbia/results/sumstats/fullSTATS_NY_$line.tab > header.txt
    awk '{if ($10<1e-4) print}' Columbia/results/sumstats/fullSTATS_NY_$line.tab > sigSNPs.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs.tab > input.annovar.txt

    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish

    cat header.txt sigSNPs.txt > temp.txt
    paste temp.txt annotatedtable.hg19_multianno.txt > temp2.txt
    sort -gk10 temp2.txt > Columbia/results/sigSNPs/Annotated_NY_$line.tab

    grep -f keep.PDsnps.txt Columbia/results/sumstats/fullSTATS_NY_$line.tab > sigSNPs_PD.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs_PD.txt > input.annovar_PD.txt

    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar_PD.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable_PD -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish

    cat header.txt sigSNPs_PD.txt > temp.txt 
    paste temp.txt annotatedtable_PD.hg19_multianno.txt > temp2.txt
    sort -gk10 temp2.txt > Columbia/results/sigSNPs/Annotated_NY_PD-snps_$line.tab
done

mkdir PPMI/results/sigSNPs
cat ~/runs/emsom/projects/GCase_GWAS/assoc.list | while read line 
do
    head -n 1 PPMI/results/sumstats/fullSTATS_PPMI_$line.tab > header.txt
    awk '{if ($10<1e-4) print}' PPMI/results/sumstats/fullSTATS_PPMI_$line.tab > sigSNPs.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs.tab > input.annovar.txt

    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish

    cat header.txt sigSNPs.txt > temp.txt
    paste temp.txt annotatedtable.hg19_multianno.txt > temp2.txt
    sort -gk10 temp2.txt > PPMI/results/sigSNPs/Annotated_PPMI_$line.tab

    grep -f keep.PDsnps.txt PPMI/results/sumstats/fullSTATS_PPMI_$line.tab > sigSNPs_PD.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs_PD.txt > input.annovar_PD.txt

    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar_PD.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable_PD -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish

    cat header.txt sigSNPs_PD.txt > temp.txt 
    paste temp.txt annotatedtable_PD.hg19_multianno.txt > temp2.txt
    sort -gk10 temp2.txt > PPMI/results/sigSNPs/Annotated_PPMI_PD-snps_$line.tab
done

rm temp.txt
rm temp2.txt
rm sigSNPs.tab
rm input.annovar.txt 
rm annotatedtable.*
rm header.txt 

## META
module load metal 
cat ~/runs/emsom/projects/GCase_GWAS/assoc.list | while read line 
do
    cp PPMI/results/sumstats/Metal_PPMI_$line.tab meta1.tab
    cp Columbia/results/sumstats/Metal_NY_$line.tab meta2.tab
    metal ~/runs/emsom/scripts/metal.txt 
    awk 'BEGIN{FS=OFS="\t"}{split($1,snp,":"); print snp[1],snp[2]}' metal1.tbl > chrbp.txt
    paste chrbp.txt metal1.tbl > meta_analysis/meta_GCase_$line.tab 

    cp meta_analysis/meta_GCase_$line.tab meta.tab
    R < ~/runs/emsom/scripts/sumstats_from_meta.R --no-save
    mv fullSTATS.meta.tab meta_analysis/fullSTATS_meta_$line.tab
    mv QQ.meta.tiff meta_analysis/QQ_meta_$line.tiff
    mv ManH.meta.tiff meta_analysis/ManH_meta_$line.tiff

    head -n 1 meta_analysis/fullSTATS_meta_$line.tab > header.txt
    awk '{if ($11<1e-5) print}' meta_analysis/fullSTATS_meta_$line.tab > sigSNPs_meta.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs_meta.txt > input.annovar.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt  ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish
    cat header.txt sigSNPs_meta.txt > temp.txt 
    paste temp.txt annotatedtable.hg19_multianno.txt | sort -gk 11 > meta_analysis/sigSNPs_meta_$line.tab

    grep -f keep.PDsnps.txt meta_analysis/fullSTATS_meta_$line.tab > sigSNPs_PD.txt
    awk '{print $1,$2,$2,$5,$4}' sigSNPs_PD.txt > input.annovar_PD.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar_PD.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg19 -out annotatedtable_PD -remove -protocol refGene,avsnp147,clinvar_20190305,gnomad_genome \
    -operation g,f,f,f -nastring . -polish
    head -n 1 meta_analysis/fullSTATS_meta_$line.tab > header.txt 
    cat header.txt sigSNPs_PD.txt > temp.txt 
    paste temp.txt annotatedtable_PD.hg19_multianno.txt > temp2.txt
    sort -gk10 temp2.txt > meta_analysis/Annotated_meta_PD-snps_$line.tab
done

## Conditional and joint analysis in Columbia 
gcta64 --cojo-file Columbia/results/sumstats/COJO_NY_$line.tab --bfile ny_gcase_noOut --cojo-slct --cojo-p 5e-8 --out NY_gcase_COJO_$line
### only hit is 1:155952732:G:A 
plink --bfile --bfile ny_gcase_noOut --linear  --geno 0.05 --maf 0.01 --ci 0.95 --freq --condition 1:155952732:G:A --hide-covar --covar  --covar covar_ny.txt --covar-name Sex,Age,Pheno,N370S,E326K,T369M,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --out NY_gcase_COJO_$line 
awk '{if ($12!="NA") print $0}' NY_gcase_COJO_$line.assoc.linear | grep 'ADD' | sort -gk12 > p.NY_gcase_COJO_$line.assoc.linear

echo $line > marker.txt
cp p.NY_gcase_COJO_$line.assoc.linear assoc
cp ~/runs/emsom/projects/GCase_GWAS/Columbia/allChrs_NY_addsam.info infos
R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
paste lambda.txt >> Columbia/results/plots/NY_lambdas.tab
mv QQ.tiff Columbia/results/plots/QQ_NY_gcase_COJO_$line.tiff    
mv Metal.tab Columbia/results/sumstats/Metal_NY_gcase_COJO_$line.tab
mv COJO.tab Columbia/results/sumstats/COJO_NY_gcase_COJO_$line.tab
mv fullSTATS.tab Columbia/results/sumstats/fullSTATS_NY_gcase_COJO_$line.tab
mv ManH.tiff Columbia/results/plots/ManH_NY_gcase_COJO_$line.tiff
