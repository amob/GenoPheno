# GenoPheno
Scripts with files used and made by scripts, packages, and notes
*indicates files must be provided by authors
**indicates files downloaded and from other published work

Script 000_gxeonly.R
	purpose: fit models for plastic respones, GxE of plant population or family to soil biota
	input:	phenomat.csv*
			pedmat.csv*
			covmat.csv*
			AlphabeticalPopEnvDat.csv*
	output:	popdamtraits_soils.pdf
			traitgxemodels.Rdata
			traitgxemodelsELEV.Rdata
			traitgxemodelsTAnn.Rdata
			traitgxemodelsPann.Rdata
			traitgxemodelsSWC.Rdata
			04_02_gxeresults.txt
	packages: MCMCglmm

Script: 001_formatgenos.sh (uses tassel version 5 to handle h5 file)
	purpose: formats genotype file, filters to sites with data in 95% of samples
	input:	AOMexGBS2_ZeaGBSv27impV5.h5*
	output: AOMexGBS2_ZeaGBSv27impV5_95filter.plk.ped 
			AOMexGBS2_ZeaGBSv27impV5_95filter.plk.map
	packages: gcc jdk/1.8 tassel/5
			 
Script: 002_RAFM_taxa95.R  
	purpose: fits coancestry matrix, makes coancestry figure
	of note: one pop is alternately labeled ML vs M. Factors in R sort by alphabet, but RAFM sorts by order it encounters a population in the input dataframes
	input:	AOMexGBS2_ZeaGBSv27impV5_95filter.plk.ped
			popmomlabel.csv*
			AlphabeticalPopEnvDat.csv*
	output: AOMexGBS2_ZeaGBSv27impV5_95filter_Rthin.ped
			allcoan_Rthin_f95_LG.txt
			ALLrafmcoan_f95_LG2MXANNA.pdf
	packages: RAFM 1.2 and dependencies
			
Script 003_drifsel_med95_scaledtraits_allsoilsandhigh.R 
	purpose: performs driftsel analysis, S and H tests
	of note: different amounts of thinning for some of the highland samples only runs
	input:	allcoan_Rthin_f95_LG.txt
			phenomat.csv*
			pedmat.csv*
			covmat.csv*
			AlphabeticalPopEnvDat.csv*
	output:	scaledtraits.samp.mt_med_LG95
			scaledtraits.samp.mc_med_LG95.Rdata
			scaledtraits.samp.tc_med_LG95.Rdata
			scaledtraits.samp.mt.h_med_LG95_440k.Rdata
			scaledtraits.samp.mc.h_med_LG95_440k.Rdata
			scaledtraits.samp.tc.h_med_LG95_440k.Rdata
			scaled.results_mLG95.txt
			scaled.h.results_HmLG95.txt
			resultH_med_LG95.txt
			resultH_med_LG95.h.txt
	packages: RAFM 1.2, driftsel 2.1.2, and dependencies

Script 004_Driftsel_results_figures.R
	purpose: figures from driftsel analysis
	input:	scaledtraits.samp.mt_med_LG95
			scaledtraits.samp.mc_med_LG95.Rdata
			scaledtraits.samp.tc_med_LG95.Rdata
			scaledtraits.samp.mt.h_med_LG95_440k.Rdata
			scaledtraits.samp.mc.h_med_LG95_440k.Rdata
			scaledtraits.samp.tc.h_med_LG95_440k.Rdata
			AlphabeticalPopEnvDat.csv*
	output:	dtf and stw 440k.pdf
			sigdiffs by pop and trait CI 95 names means torder 440.pdf
			bivariate sigdiffs by trait CI 95 names means torder 440.pdf
	packages: driftsel 2.1.2, SDMTools, and dependencies

Script 005_5traitAnimalModelComparison.R
	purpose: fit G matrices
	input:	phenomat.csv*
			pedmat.csv*
			covmat.csv*
	output:	popsoil.animal.mods.5trait.Rdata
			popsoil.animal.mods.nearlyimp.prior.5trait.Rdata
			popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata
			Popsoil_Peds.Rdata
	packages: MCMCglmm

Script 006_checkconv.R
	purpose: inspect traces of G matrix model fitting
	input:	popsoil.animal.mods.5trait.Rdata
			popsoil.animal.mods.nearlyimp.prior.5trait.Rdata
			popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata
	output:	trace plots of MCMC chains for each of 30 G matrices for each of the three priors
	packages: MCMCglmm

Script 007_heritability.R
	input:	popsoil.animal.mods.5trait.Rdata
	output:	heritabilities.pdf
	packages: MCMCglmm
	
Script 008_eigentensors.R
	purpose: eigentensor analysis of G matrices
	input:	popsoil.animal.mods.5trait.Rdata
			popsoil.animal.mods.nearlyimp.prior.5trait.Rdata
			popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata
			Popsoil_Peds.Rdata
			AlphabeticalPopEnvDat.csv*
	output:	EigenValS_AnPopSoilG5trait.pdf
			TraitsPopSoilEigenVofEigenT5traits.pdf
			PopsoilTensor1.5trait.pdf
			[variance explained by eigentensors and eigenvalues of eigentensors]
	packages: gdata, matrixcalc, MCMCglmm, dependencies
	
Script 009_selproj.R
	purpose: Project responses to selection, generate figures
	input:	popsoil.animal.mods.5trait.Rdata
			popsoil.animal.mods.nearlyimp.prior.5trait.Rdata
			popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata
			Popsoil_Peds.Rdata
			AlphabeticalPopEnvDat.csv*
	output:	cont2SelSkews3_5Trait.pdf
			soil_sel_all.pdf
	packages: ks, MCMCglmm, dependencies
			
Script 010_Gmatplots.R
	purpose: Plot of G matrices
	input:	popsoil.animal.mods.5trait.Rdata
			popsoil.animal.mods.nearlyimp.prior.5trait.Rdata
			popsoil.animal.mods.nearlynoninf.prior.5trait.Rdata
			Popsoil_Peds.Rdata
			AlphabeticalPopEnvDat.csv*
	output: RawGAnPopsoil5trait.pdf
	packages: MCMCglmm


Script 011_swartssnps1.sh
	purpose: generate missingness data, filter to taxa to use in PCA
	input:	AOMexGBS2_ZeaGBSv27impV5.h5*
			MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf**
			012_keeptaxa.R*
	output:	AOMexGBS2_ZeaGBSv27impV5.vcf
			missing_swarts2017.imiss
			missing_AOMexGBS2.imiss
	packages: gcc jdk/1.8 tassel/5 vcftools/0.1.13

Script 012_keeptaxa.R
	purpose: get a set of samples with enough SNP calls from both datasets to use in PCA
	input:	swarts_taxatokeep_firstpass.csv*
			missing_swarts2017.imiss
			missing_AOMexGBS2.imiss
	output:	SwartsKeepTaxa_indFilt70.txt
			AOMexGBS2_indFilt70.txt
	packages: base R only

Script 013_matchrefgen_SwartsAnnaSNPs1.R
	purpose: export as a bed file
	input: MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf**
	output: MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed
	packages: base R
	
Manual action on ENSEMBL required for coordinate conversion of swarts snp files, and identification of non-convertible SNPs
	url: http://plants.ensembl.org/Oryza_sativa/Tools/AssemblyConverter?db=core
	input: MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed
	action: convert coodinates of that file from AGPv3 to AGPv2, then, use the output as input and convert back to AGPv3
	output:	output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed
			output_output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed

Script 014_matchrefgen_SwartsAnnaSNPs2.R
	purpose: set up files to get coordinates for two GBS datasets in PCA analysis
	input: 	MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf**
			AOMexGBS2_ZeaGBSv27impV5.vcf
			output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed
			output_output_MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_forconversion.bed
	output: MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2.vcf
	packages: base R

Script 015_swartssnps2.sh
	purpose: merge two Zea GBS datasets, filter SNPs, run PCA analysis
	input:	MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes.vcf**
			MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2.vcf
	output:	header.txt
			MinimallyFilteredWithNAM_InbredLandraces_InbredTeosintes_agpv3_ConvertedTo_agpv2_newheader.vcf
			SwartsKeepTaxa_indFilt70.txt
			AOMexGBS2_indFilt70.txt
			MergeSwartsAnna_80.plk.ped
			MergeSwartsAnna_80.plk.map
			pca_MergeSwartsAnna_801.txt
			pca_MergeSwartsAnna_802.txt
			pca_MergeSwartsAnna_803.txt
	packages: gcc jdk/1.8 tassel/5.2.14 plink/1.90 vcftools/0.1.13

Script 016_plotPCA_swartsanna.R
	purpose: Plot of PCA of current and other Zea genotype data
	input:	pca_MergeSwartsAnna_801.txt
			pca_MergeSwartsAnna_802.txt
	output:	swartsannaPCA.pdf
	packages: base R only
	
	
AOMexGBS2_ZeaGBSv27impV5.h5 is available in an alternate format on figshare, available upon publication 10.6084/m9.figshare.4714030

Data from Swarts et al 2017 (Science, DOI: 10.1126/science.aam9425) obtained from: http://datacommons.cyverse.org/browse/iplant/home/shared/panzea/dataFromPubs/Swarts2017Science
