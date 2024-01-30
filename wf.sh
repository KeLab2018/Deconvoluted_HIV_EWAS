#!/bin/bash
#for GSE217633 hiv ewas analysis
root=/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
GSE=GSE217633
#fro mcc
module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0"

mkdir -p /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633

download_geo.R GSE217633
#produce phe file
Rscript ~/mybiotools/r/produce_phe_from_geo.R geo_accession ch1$ geo_phenotype.txt ~/mybiotools/docker/hiv_tca_ewas_addt/get_batch.R

#download GSE217633_RAW.tar mannually from web page
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217633

#raw data processing
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
mkdir idatfiles qc rawdata TCA
cd idatfiles
Rscript ~/mybiotools/r/produce_idat_csv.R .
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0" && runR.sh  ~/mybiotools/r/minfi_QC_850k.R
#trim sample id
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/rawdata
Rscript ~/mybiotools/r/rem_batch_barcode.R ssNood_beta.txt 
Rscript ~/mybiotools/r/rem_batch_barcode.R ssNood_beta_nona.txt
Rscript ~/mybiotools/r/rem_batch_barcode.R ssNood_m.txt
Rscript ~/mybiotools/r/rem_batch_barcode.R ssNood_m_nona.txt

#cell type estimation & TCA
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA
glint.py --datafile ../rawdata/ssNood_beta_nona.txt --houseman --out PBMC_beta
Rscript ~/mybiotools/r/TCA_scripts/change_glint_celltype_res.R PBMC_beta.houseman_estimates.txt PBMC_beta_ct.txt
Rscript ~/mybiotools/r/MCseq_TCA/01HPC_tca_pre_selecting_probe_varLT0002.R ../rawdata/ssNood_beta_nona.txt > 01HPC_tca_pre_selecting_prbe_varLT0002.R.log
#split data
mkdir split_files
cd split_files
split_file_with_header.pl ../betaVarLT00002_X.txt 250 Y
cd ..
Rscript ~/mybiotools/r/TCA_scripts/01produce_slurm.R split_files 02HPC_tca_ctmethy_common.R MC tca_methyCT_slurm.sh
pbsv2.pl -q scavenge -i tca_methyCT_slurm.sh -ppn 8 -pmem 48G
pbsv2.pl -q scavenge -pmem 300G -ppn 10 -wt 3:00:00 "Rscript ~/mybiotools/r/TCA_scripts/03merge_ctmethy.R split_files MC"
#deleting useless figures
rm *_METHY_* -f

#bug fix, can be omit in future
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA
tf=`ls *y.txt`
for i in ${tf[@]};do
	Rscript ~/mybiotools/r/rem_batch_barcode.R ${i}
done
Rscript ~/mybiotools/r/rem_batch_barcode.R betaVarLT00002_X.txt

#PositiveControlPCA
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
mkdir -p PCA/PositiveControlPCA split_files
#cd PCA/PositiveControlPCA
cp ~/mybiotools/parameter_files/cont_pca_para_geo.R cont_pca_para.R
#! if need, vi cont_pca_para.R
runR.sh ~/mybiotools/r/PositiveControlPCA_850k.R

#making residual dataset
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/split_files
split_file_with_header.pl ../rawdata/ssNood_beta.txt 20 Y
residual_pbs_glm_PCA.pl ssNood_beta_*.txt  -t ../geo_phenotype.txt -pherowid id -addt ../PCA/PositiveControlPCA/cont_PC_850k.txt -addpherowid id  -addt ../TCA/PBMC_beta_ct.txt -addpherowid id -p ../PCA/residual_C30contPCA -r exp_value -v F -m gaussian -f gender -c gender age CD8T CD4T Gran NK Bcell Mono Cont_Pr_850k_PC1 Cont_Pr_850k_PC2 Cont_Pr_850k_PC3 Cont_Pr_850k_PC4 Cont_Pr_850k_PC5 Cont_Pr_850k_PC6 Cont_Pr_850k_PC7 Cont_Pr_850k_PC8 Cont_Pr_850k_PC9 Cont_Pr_850k_PC10 Cont_Pr_850k_PC11 Cont_Pr_850k_PC12 Cont_Pr_850k_PC13 Cont_Pr_850k_PC14 Cont_Pr_850k_PC15 Cont_Pr_850k_PC16 Cont_Pr_850k_PC17 Cont_Pr_850k_PC18 Cont_Pr_850k_PC19 Cont_Pr_850k_PC20 Cont_Pr_850k_PC21 Cont_Pr_850k_PC22 Cont_Pr_850k_PC23 Cont_Pr_850k_PC24 Cont_Pr_850k_PC25 Cont_Pr_850k_PC26 Cont_Pr_850k_PC27 Cont_Pr_850k_PC28 Cont_Pr_850k_PC29 Cont_Pr_850k_PC30
cd ../residual/residual_C30contPCA
check_glm_results.pl dm
pbs_result_collect.pl *_dm.txt -o residual_dm.txt
cp *000000_residual_dm.R backup_Rscript.R
rm -f *efile *ofile *\.sh\.?*
rm -f *_residual_dm.txt *_dm.R

#residual PCA
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
cp ~/mybiotools/parameter_files/resipca_para_geo.R resipca_para.R
runR.sh ~/mybiotools/r/resi_pca2phe2pheatmap_850k.R # for 850k

cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633
#post-ART vs. control
bash ~/mybiotools/docker/hiv_tca_ewas_addt/batch_geo_hiv_ewas.sh
#pre-ART vs. control
bash ~/mybiotools/docker/hiv_tca_ewas_addt/batch_geo_hiv_ewas_preART.sh
# for removing t.value=inf probes (pvalue = 0)
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/filter_inf_glm.R
# post ewas
module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0" && Rscript ~/mybiotools/r/TCA_scripts/04ctmethylomeEWAS_results_check_v2.R _methylome_hiv_ewas PBMC_methylome_hiv_ewas 0.05 topright 2.2 glm_pv_adj_anno.xls F
si *png *pdf */*png */*pdf */top5000.xls */*log */*R */inflation_lambda.txt

# mv all TCA results from /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633 to /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA
#meta analysis
mkdir -p /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/batch_meta_hiv_3datasets.R
#paper meta mht
doit_Specific_dir_pbs.pl * -p F -c "head -n 5001 EWAS_Meta_analysis_sampleSize_adj.xls > top5000.xls"
#original
#module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0" &&  Rscript ~/mybiotools/r/TCA_scripts/meta_results_mht_paper.R "" PBMC 0.05 top 2.2 EWAS_Meta_analysis_sampleSize_adj.xls
#adding heterog pv > 0.05 filter
module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0" &&  Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/meta_results_mht_heterogPvFilter_paper.R "" PBMC 0.05 top 2.2 EWAS_Meta_analysis_sampleSize_adj.xls
#paper circos
Rscript ~/mybiotools/r/TCA_scripts/meta_sig_Rcircos_paper.R 
#merging clusterProfile gmt enrich resutls for cell types
vi merge_gsea_enrich_hcluster_para.R
Rscript ~/mybiotools/r/merge_clusterPrile_enrich_hcluster.R
#new 3D effectsize plot for paper
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/effectsize_plot_paper_3cohort.R

#venn plot
module load R/4.2.0-foss-2020b
mkdir -p /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/venn_plot_vs_2cohorts_meta
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/venn_plot_vs_2cohorts_meta
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/gse_vs_2cohortMeta_venn_plot.R
mkdir -p /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/venn_plot_vs_vacs
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/venn_plot_vs_vacs
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/gse_vs_vacs_venn_plot.R

#cell type proportion plot
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/ct_jitter_plot.R
#bar plot for meta
module load R/4.2.0-foss-2020b
Rscript /gpfs/ycga/project/xu_ke/xz345/soft/git/mybiotools/r/TCA_scripts/TSS_barplot_meta.R "" /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets/barplot
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/promoter_test.R

# WIHS whole 850K HIV EWAS
#ewas PBMC
#cd /gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/splitfile
mkdir -p /gpfs/ycga/work/xu_ke/xz345/work/other/Brad_022020_Methy/ewas
cd /gpfs/ycga/work/xu_ke/xz345/work/other/Brad_022020_Methy/ewas
split_file_with_header.pl ../rawdata/ssNoob_Dsnp_beta.txt 100 Y
pbs_infile_glm.pl ssNoob_Dsnp_beta_*.txt -wt 05:00:00 \
	-p /gpfs/ycga/work/xu_ke/xz345/work/other/Brad_022020_Methy/ewas/hiv \
	-ifbeta T -pherowid microarrayID \
	-addt /gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/PCA/PositiveControlPCA/cont_PC_Brad.txt \
	-addpherowid microarrayID \
	-addt /gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/residual/residualPC_noVL.txt \
	-addpherowid microarrayID -e /gpfs/ycga/project/xu_ke/xz345/soft/git/mybiotools/r/TCA_scripts/brad_ext_v2_orig.R \
	-t /gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/phenotype/Brad_phe_v3.txt \
	-r exp_value -c HIV smoker drinker race AGEATVIS CD8T CD4T Gran NK Bcell Mono \
	Cont_Pr_Brad_PC1 Cont_Pr_Brad_PC2 Cont_Pr_Brad_PC3 Cont_Pr_Brad_PC4 Cont_Pr_Brad_PC5 Cont_Pr_Brad_PC6 Cont_Pr_Brad_PC7 Cont_Pr_Brad_PC8 Cont_Pr_Brad_PC9 Cont_Pr_Brad_PC10 Cont_Pr_Brad_PC11 Cont_Pr_Brad_PC12 Cont_Pr_Brad_PC13 Cont_Pr_Brad_PC14 Cont_Pr_Brad_PC15 Cont_Pr_Brad_PC16 Cont_Pr_Brad_PC17 Cont_Pr_Brad_PC18 Cont_Pr_Brad_PC19 Cont_Pr_Brad_PC20 Cont_Pr_Brad_PC21 Cont_Pr_Brad_PC22 Cont_Pr_Brad_PC23 Cont_Pr_Brad_PC24 Cont_Pr_Brad_PC25 Cont_Pr_Brad_PC26 Cont_Pr_Brad_PC27 Cont_Pr_Brad_PC28 Cont_Pr_Brad_PC29 Cont_Pr_Brad_PC30 PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20 PC21 PC22 PC23 PC24 PC25 PC26 PC27 PC28 PC29 PC30 \
	-f smoker drinker HIV race -v F -am binomial
cd /gpfs/ycga/work/xu_ke/xz345/work/other/Brad_022020_Methy/ewas/hiv && pbsv2.pl -q scavenge -pmem 12000 -ppn 2 -wt 24:00:00 "850k_methy_glm_annot_mht.sh 0.1" && cd ..

#PreART vs PostART
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA
#post-ART vs. control
bash ~/mybiotools/docker/hiv_tca_ewas_addt/batch_geo_hiv_ewas_preART_vs_postART.sh
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/preART_vs_postART
rm -f *_tcaMethy_*txt

#boxplot for top 100 CpGs in all cohors
cd /gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets
module load "R-bundle-Bioconductor/3.15-foss-2020b-R-4.2.0"
Rscript ~/mybiotools/docker/hiv_tca_ewas_addt/topCpG_boxplot.R

