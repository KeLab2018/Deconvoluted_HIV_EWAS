
#test Gran CD4T CD8T NK Bcell Mono  ... methylome EWAS, 20 resi resi_850k_PC
#################
mkdir -p ~/temp/GSE217633/TCA/hiv_ewas_preART
cd ~/temp/GSE217633/TCA/hiv_ewas_preART
# 20 resiPC
#ewas Gran
split_file_with_header.pl ../Gran_tcaMethy.txt 208 Y
pbs_infile_glm.pl Gran_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p Gran_methylome_hiv_ewas -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd Gran_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas CD8T
split_file_with_header.pl ../CD8T_tcaMethy.txt 208 Y
pbs_infile_glm.pl CD8T_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p CD8T_methylome_hiv_ewas -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd CD8T_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas NK
split_file_with_header.pl ../NK_tcaMethy.txt 208 Y 
pbs_infile_glm.pl NK_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p NK_methylome_hiv_ewas --r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd NK_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas Mono
split_file_with_header.pl ../Mono_tcaMethy.txt 208 Y
pbs_infile_glm.pl Mono_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p Mono_methylome_hiv_ewas -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd Mono_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas Bcell
split_file_with_header.pl ../Bcell_tcaMethy.txt 208 Y
pbs_infile_glm.pl Bcell_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p Bcell_methylome_hiv_ewas -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd Bcell_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas CD4T
split_file_with_header.pl ../CD4T_tcaMethy.txt 208 Y
pbs_infile_glm.pl CD4T_tcaMethy_0*.txt -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -p CD4T_methylome_hiv_ewas -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd CD4T_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
#ewas PBMC
cd /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/split_files
pbs_infile_glm.pl betaVarLT00002_X_*.txt -p ~/temp/GSE217633/TCA/hiv_ewas_preART/PBMC_methylome_hiv_ewas -ifbeta T -e /gpfs/ycga/work/xu_ke/xz345/soft/git/McCleary/mybiotools/r/docker/hiv_tca_ewas_addt/hiv_ewas_addt_preART.R -pherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/PCA/residual_C30contPCA/resiPCA.txt -addpherowid id -addt /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/TCA/PBMC_beta_ct.txt -addpherowid id -t /vast/palmer/scratch/xu_ke/xz345/mytemp/GSE217633/geo_phenotype.txt -r exp_value -c HIV age gender CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30 resi_850k_PC1 resi_850k_PC2 resi_850k_PC3 resi_850k_PC4 resi_850k_PC5 resi_850k_PC6 resi_850k_PC7 resi_850k_PC8 resi_850k_PC9 resi_850k_PC10 resi_850k_PC11 resi_850k_PC12 resi_850k_PC13 resi_850k_PC14 resi_850k_PC15 resi_850k_PC16 resi_850k_PC17 resi_850k_PC18 resi_850k_PC19 resi_850k_PC20 resi_850k_PC21 resi_850k_PC22 resi_850k_PC23 resi_850k_PC24 resi_850k_PC25 resi_850k_PC26 resi_850k_PC27 resi_850k_PC28 resi_850k_PC29 resi_850k_PC30 -f HIV gender -v F -am binomial
cd ~/temp/GSE217633/TCA/hiv_ewas_preART/PBMC_methylome_hiv_ewas && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 6:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..
