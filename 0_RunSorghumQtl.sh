#! /bin/bash

TASSEL5="perl /home/jgwall/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"


# Not the best pipeline because of writing the first part on my travel Windows laptop, but it works

pedfile=0_sorghum_pedigree.fsfhap.txt
imputedir=1c_imputed

if [ ! -e $imputedir ]; then mkdir $imputedir; fi

# # Do quick checks on genotypes and phenotypes to know what sort of quality control I should do
# Rscript 1_CheckGenos.r
# Rscript 1_CheckPhenos.r

# # Filter the genotypes and transform the phenotypes so they actually look right
# Rscript 1a_FilterGenos.r
# Rscript 1a_TransformPhenos.r

# # Check LD among the markers to make sure chromosomes look okay
# Rscript 1b_ConvertGenosToHapmap.r
# $TASSEL5 -h 1b_genos_filtered.hmp.txt -ld -ldType All -ldHetTreatment Homozygous -ldd png -o 1b_genos_filtered.ld.jpg -ldplotsize 2000


# # Impute with FILLIN
# chroms=`tail -n +2 1b_genos_filtered.hmp.txt | cut -f 3 | sort | uniq | grep "^Chr" | sed s/Chr// `
# for chr in $chroms; do
#     $TASSEL5 -h 1b_genos_filtered.hmp.txt -separate $chr -export $imputedir/1c_genos.$chr.hmp.txt.gz
#     $TASSEL5 -h $imputedir/1c_genos.$chr.hmp.txt.gz -FSFHapImputationPlugin -pedigrees $pedfile -logfile $imputedir/1d_fsfhap.$chr.log -phet 0.5 -outParents false -endPlugin \
#         -export $imputedir/1d_imputed${chr}_.hmp.txt
#     mv $imputedir/1d_imputed${chr}_1.hmp.txt $imputedir/1d_imputed$chr.hmp.txt
#     rm $imputedir/1d_imputed${chr}_2.txt # don't care about keeping haplotype file
# #     break
# done

# TODO: Figure out why Chroms 5 & 8 failed and fix?
# # Build up final genotype file (have to test because some chromosomes failed for some reason 
# commands=""
# inputs=""
# runs=""
# for chr in $chroms; do
#     if [ -e $imputedir/1d_imputed$chr.hmp.txt ]; then
#         commands="$commands -fork$chr -h $imputedir/1d_imputed$chr.hmp.txt"
#     else
#         commands="$commands -fork$chr -h $imputedir/1c_genos.$chr.hmp.txt.gz"
#     fi
#     
#     inputs="$inputs -input$chr"
#     runs="$runs -runfork$chr"
# done
# $TASSEL5 $commands -combine999 $inputs -mergeGenotypeTables -export 1d_genos.imputed.hmp.txt $runs

# Remove redundant markers
# $TASSEL5 -h 1d_genos.imputed.hmp.txt -NumericalGenotypePlugin -endPlugin -export 1d_genos.imputed.numeric.txt -exportType ReferenceProbability
# Rscript 1e_FindRedundantMarkers.r -i 1d_genos.imputed.numeric.txt -o 1e_redundant_sites.txt 
# $TASSEL5 -h 1d_genos.imputed.hmp.txt -excludeSiteNamesInFile 1e_redundant_sites.txt  -export 1f_genos.imputed.unique.hmp.txt

# # Calculate heritability of each phenotype
# $TASSEL5 -h 1f_genos.imputed.unique.hmp.txt -ck -export 1f_genos.imputed.unique.kinship.txt
# head -n 2 1d_genos.imputed.hmp.txt > 1g_dummy_genos.txt
# scriptdir=0_HeritabilityScripts 
# blups="1a_phenos_transformed.tassel.txt"
# dummy_genos=1g_dummy_genos.txt
# kinship=1f_genos.imputed.unique.kinship.txt
# nperms=10000
# cutoff=0.01 # Considered a "good" phenotype if empirical p-value is below this
# targetdir=1g_heritability
# outprefix=1g_heritability
# if [ ! -e $targetdir ] ; then mkdir $targetdir; fi
# bash 1g_RunTasselHeritabilityPerms.sh $scriptdir $blups $dummy_genos $kinship $nperms $cutoff $targetdir $outprefix "$TASSEL5"

# # Filter phenotypes and run a quick GWAS in TASSEL.
# best_phenos="noeggs_categories noeggs_log eggperroot_categories eggperroot_log Rootwt BRIX"
# python3 1h_FilterGoodTraits.py -i 1a_phenos_transformed.tassel.txt -o 1h_good_phenos.tassel.txt --traitnames $best_phenos
# $TASSEL5 -fork1 -h 1d_genos.imputed.hmp.txt -fork2 -t 1h_good_phenos.tassel.txt -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -permute true -nperm 10000 -endPlugin -export 1i_glm   
# mv 1i_glm1.txt 1i_glm.sitefile.txt
# mv 1i_glm2.txt 1i_glm.allelefile.txt
# # Plot Manhatten plots - Modify existing script to do this
# python3 1j_PlotGwasData.py -i 1i_glm.sitefile.txt -o 1j_gwas_results --chromlengths 0_sorghum_chrom_lengths.txt #--debug

# Do the same pre-imputation to make sure imputation didn't mess up the data
# $TASSEL5 -fork1 -h 1b_genos_filtered.hmp.txt -fork2 -t 1h_good_phenos.tassel.txt -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -permute true -nperm 10000 -endPlugin -export 1i_glm_noimpute   
# mv 1i_glm_noimpute1.txt 1i_glm_noimpute.sitefile.txt
# mv 1i_glm_noimpute2.txt 1i_glm_noimpute.allelefile.txt
# grep -v "SUPER_" 1i_glm_noimpute.sitefile.txt > 1i_glm_noimpute.sitefile.core_chroms.txt
# python3 1j_PlotGwasData.py -i 1i_glm_noimpute.sitefile.core_chroms.txt -o 1j_gwas_results_noimpute --chromlengths 0_sorghum_chrom_lengths.txt # --zoom 05 #--debug

# # Now add the best marker as a covariate and test it again. Note that marker is included as both additive and dominant so as to capture the "genotype" model TASSEL uses
# Rscript 1k_MakeCovariateFromBestMarker.r -i 1i_glm.sitefile.txt --genos 1d_genos.imputed.numeric.txt -o 1k_best_marker 
# # All the eggs-related ones have the identical marker as best, so just need to use that
# $TASSEL5 -fork1 -h 1d_genos.imputed.hmp.txt -fork2 -t 1h_good_phenos.tassel.txt -fork3 -t 1k_best_marker.noeggs_log.txt -combine4 -input1 -input2 -input3 -intersect \
#     -FixedEffectLMPlugin -permute true -nperm 10000 -endPlugin -export 1l_glm_with_covariate
# mv 1l_glm_with_covariate1.txt 1l_glm_with_covariate.sitefile.txt
# mv 1l_glm_with_covariate2.txt 1l_glm_with_covariate.allelefile.txt
# python3 1j_PlotGwasData.py -i 1l_glm_with_covariate.sitefile.txt -o 1m_gwas_results_with_covariate --chromlengths 0_sorghum_chrom_lengths.txt #--debug

# No more signal, so there appears to be just a single QTL

# Finally, get LD of the target SNP. It's s1221, so site #1220 (java starts numbering at 0)
# mysite=1220
# $TASSEL5 -h 1b_genos_filtered.hmp.txt -ld -ldType SiteByAll -ldTestSite $mysite -ldHetTreatment Homozygous -export 1n_best_snp_ld.txt
# python3 1n_PlotLd.py -i 1n_best_snp_ld.txt -o 1n_ld_plot --chromlengths 0_sorghum_chrom_lengths.txt 

# # # Check genetic map in R/QTL just to be sure - TODO - didn't finish this one
# echo "ENTRY107" > 1o_parentA.txt
# echo "ENTRY57" > 1o_parentB.txt
# $TASSEL5 -h 1d_genos.imputed.hmp.txt -GenosToABHPlugin -o 1o_genos_imputed.abh.csv -parentA 1o_parentA.txt -parentB 1o_parentB.txt -outputFormat c
# markers_to_drop="s15 s85 s499 s562 s563 s2149 s2150 s2151 s2152 s2454 s2455"
# Rscript 1p_CheckMapInRqtl.r -i 1o_genos_imputed.abh.csv -o 1p_rqtl_map -p 1h_good_phenos.tassel.txt --badmarkers $markers_to_drop
# Rscript 1q_RunQtlMappingInRqtl.r -i 1p_rqtl_map.csv -o 1q_rqtl_analysis --perms 1000 --zoom 05 --num-cores 7

# Don't have this setup in the pipeline yet, but basically it gives a 3-lod interval from s1220 to s1222 (so, the target SNP and the two to either side)


# # # Publication graphics
# # # GWAS
# python3 1j_PlotGwasData.py -i 1i_glm.sitefile.txt -o 2_gwas --chromlengths 0_sorghum_chrom_lengths.txt 
# python3 1j_PlotGwasData.py -i 1i_glm.sitefile.txt -o 2_gwas_zoom --chromlengths 0_sorghum_chrom_lengths.txt --zoom 05
# # # LD
# python3 1n_PlotLd.py -i 1n_best_snp_ld.txt -o 2_ld --chromlengths 0_sorghum_chrom_lengths.txt
# python3 1n_PlotLd.py -i 1n_best_snp_ld.txt -o 2_ld_zoom --chromlengths 0_sorghum_chrom_lengths.txt --zoom 05

# Manual inspection of the LD output with a cutoff of 0.8 impluies the QTL is from 5:2868700 to about 5:10382873, but there's a lot of wiggle room on the right slope

# # Publication SNPs stats for text
# target_snp=s1221
# Rscript 2a_SnpStats.r -i 1d_genos.imputed.hmp.txt --phenofile 1h_good_phenos.tassel.txt --snp $target_snp -o 2a_snp_stats --reverse-log
# Rscript 2a_SnpStats.r -i 1d_genos.imputed.hmp.txt --phenofile 1a_phenos_transformed.tassel.txt --snp $target_snp -o 2a_snp_stats.orig --reverse-log # Original phenotypes; easier to work with for paper main text

