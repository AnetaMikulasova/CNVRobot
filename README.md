# CNVRobot

<p align="left">
  <img src="https://github.com/AnetaMikulasova/CNVRobot/blob/main/CNVRobot_logo.png" alt="CNVRobot logo" width="200" height="207"/>
</p>

Welcome to CNVRobot v4.1!

IMPORTANT UPDATES:

Please note the following UPDATES compared to the v3 version:
1) Analysis is now available in CHM13v2.0. 
    - The FASTA file can be downloaded at [T2T aws](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/analysis_set/) (`chm13v2.0.fa.gz`)
    - To view profiles of whole chromosomes, read counts must also be collected in regions of low mapping quality. To enable this (by disabling "MappingQualityReadFilter"), `READ_COUNT_MODE` should be set to `DF` in `master_project`. Please note, this is a new column in `master_project`.
2) A new R package is required: [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html).
3) Checkpoints have been introduced: Several control points are implemented that halt the process by giving an error message when something unexpected is detected (e.g., unknown or empty value in master table, BAM file not found, database file missing). These checkpoints are still in development and may be error-prone. Therefore, there's an option to disable them by changing `CHECKPOINTS="yes"` to `CHECKPOINTS="no"` at the beginning of the `run.sh` script. If an error is discovered in the checkpoints or if any other checkpoint can be suggested, I would appreciate it if you [contact me](mailto:Aneta.Mikulasova@newcastle.ac.uk).

Work in progress: I have developed a module that allows outputs from CNVRobot to be directly processed by [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software) for matching tumor-germline analyses. This provides information about ploidy and absolute copy-number analyses. If you are interested, please do not hesitate to [contact me](mailto:Aneta.Mikulasova@newcastle.ac.uk) for the in-progress code.

## 1. Description
CNVRobot is an integrated pipeline designed to detect rare germline and somatic copy-number variants (CNVs) as well as loss of heterozygosity (LOH) regions in the human genome using short-read DNA sequencing from any NGS platform. This includes targeted sequencing (TS), whole-exome sequencing (WES), and whole-genome sequencing (WGS). It integrates sequencing depth with SNP zygosity across the genome, incorporates advanced data denoising, segmentation and variant prioritization options, and provides detailed visualization for manual inspection of detected CNVs.

## 2. Download
The code, databases, and example data can be downloaded from [HERE](https://newcastle-my.sharepoint.com/:f:/g/personal/nam320_newcastle_ac_uk/Eu2m1nidAKVDoXxuDONbqToBJfbYPwpbjeRuU6ERhcqHrg?e=yshGbA). The latest version is CNVRobot v4.1.

## 3. Dependencies (License)
- [GATK v4.2+](https://github.com/broadinstitute/gatk/releases) (Apache License 2.0)
- [PICARD](https://broadinstitute.github.io/picard/) (MIT License)
- [R](https://www.r-project.org) (GPL-2 | GPL-3)
- [Samtools, HTSlib (bgzip and tabix), BCFtools](http://www.htslib.org) (MIT/Expat License)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (GNU General Public License v2.0)
- [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/) (UCSC Genome Browser License)
- [dos2unix](http://dos2unix.sourceforge.net) (BSD License)
- [java](https://www.java.com/) (Oracle Binary Code License)
- [gzip/gunzip](https://www.gnu.org/software/gzip/) (GNU General Public License)
- R packages: 
    - [tidyverse](https://www.tidyverse.org) (MIT License)
    - [expss](https://cran.r-project.org/web/packages/expss/index.html) (GPL-2 | GPL-3)
    - [splitstackshape](https://cran.r-project.org/web/packages/splitstackshape/index.html) (GPL-3)
    - [karyoploteR](http://bioconductor.org/packages/release/bioc/html/karyoploteR.html) (Artistic-2.0)
    - [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) (Artistic-2.0)

## 4. Inputs
- BAM files
- FASTA file - the reference file used for alignment
- capture BED file (for WES or TS)
  - first three columns required - contig, start, and end; no header
  - requires the same genome assembly as the reference file (Contigs not found in the reference will be automatically excluded.)
- master files - spreadsheets of the project, controls, samples, etc.; see below for details

## 5. Initiate run
1. Fill `./CNVRobot_vX.X/Masters/setting_base.sh` to specify paths for software.
2. Fill `./CNVRobot_vX.X/Masters/master_projects` to specify conditions for your project.
3. Fill `./CNVRobot_vX.X/Masters/master_samples.txt` and `./CNVRobot_vX.X/Masters/master_controls.txt` to specify samples and control sets for the project.
4. Run `./CNVRobot_vX.X/Scripts/run.sh`.

How CNVRobot regulates what project(s), sample(s) and control(s) are executed?
- The `INCLUDE` column in each master file serves as a switch to determine which projects to execute, which samples to analyze, and which controls to employ. This way, it's possible to maintain multiple projects within the project master file, designating `no` in the `INCLUDE` column for the projects you do not wish to activate.
- The same principle applies to the sample master file; by marking `INCLUDE` as `yes` for a particular sample, you're indicating your intention to analyze it. Conversely, if you want to exclude certain samples from the current analysis, simply assign `no` to the `INCLUDE` field.
- For the control master file, you might have multiple controls listed. If you need to exclude some controls, perhaps due to poor quality control, assign `no` to their respective `INCLUDE` fields. CNVRobot will register this modification and execute the analysis excluding the specified controls.
- CNVRobot will link samples and controls for each project by matching `PROJECT_ID`, `CAPTURE_ID`, and `GENOME_VERSION` columns.

- Details about each master file are provided below.
- Please avoid using space characters in any identifier.

## 6. Example data
To run this test data, please fill `./CNVRobot_vX.X/Masters/setting_base.sh` and run `./CNVRobot_vX.X/Scripts/run.sh`.  
Simulated example data contains one sample, two controls, and a reference file.  
Platform: targeted sequencing  
Chromosome locus: 21q21  
Genomic position: chr21:15.0-16.0Mb  
Genome assembly: hg38  
CNVs in the sample: loss of chr21:15.2-15.3Mb, gain of chr21:15.7-15.8Mb  

## 7. Control set selection
GATK recommends that the control set consists of at least 40 samples without pathogenic CNVs, which were prepared, sequenced, and analyzed using the same laboratory protocols and bioinformatic processing. However, any simultaneously processed control is better than no control, and the use of just 2–5 controls can already be highly effective in reducing noise from technical artifacts. The control set can overlap with the sample set if no pathological CNVs are expected. For example, the control set can consist of unaffected parents in a pedigree dataset or germline samples in a cancer dataset. Any related samples are automatically excluded from controls during analysis using the main identifier stated in the sample and control spreadsheet. 

## 8. Sex test
The sex test is a part of the preparation script. It compares local coverage of two genes, the tested gene and the control gene, to verify the sex of each control. This helps to recognize potential sample swaps that could disrupt normalization, denoising, and analysis of gonosomes. Genomic coordinates for collection are provided by the sex master file, and the test on/off is regulated in the project's master file. The default setting involves the *SRY* and *GAPDH* genes, but requires both genes to be covered. If the default genes are not covered by the sequencing, alternative loci can be specified in the sex master file. Pairing customized loci with the project is defined by `PROJECT_ID`. If a customized setting is used, the tested loci (on chromosome Y) and the control loci (on any autosome) must be of the same size (bp). The sex test may not be possible to perform in the case of some targeted panels where no target at chromosome Y is included.

## 9. Data segmentation
The segmentation of denoised coverage and SNP zygosity data is an important step to recognize abnormal segments (losses, gains, and LOH).
### 9.1. Segmentation conditions
The pipeline offers a default segmentation setting. However, results can be reviewed and segmentation conditions adapted by the user based on data type, quality, and expectation. There are two options for customizing segmentation conditions using the project's master file `SEGMENTATION_ID` column:  
a) **smart segmentation**: The identifier pattern `smart-XX` is provided within `SEGMENTATION_ID` (in the project's master file), where `XX` represents a segmentation coefficient, a number between 0.5 (for high sensitivity) and 1 (for high specificity). If smart segmentation is used, segmentation conditions are not loaded from the segmentation master file, but are instead calculated automatically by the given segmentation coefficient and detected sample quality. For subclonal analysis, the pattern is `smart-XX-sub-YY`, with `YY` representing a number between 0.1 and 1.0 (for example, 0.5 indicates a subclonal analysis to detect monoallelic CNV in 50% of cells).  
b) **custom segmentation**: A short identifier (for example, *my_segm*, *somatic_segmentation*, etc.) is provided within `SEGMENTATION_ID` (in the project's master file). This identifier matches a unique column name within the segmentation master file, where all segmentation conditions are manually defined by the user.
### 9.2. Sex-flip 
For calling gonosomal abnormalities, controls of the same sex are required for data denoising. The pipeline operates in a sex-flip recognizing mode. If the opposite sex is recognized for denoising (as defined in the project's master file), the pipeline will automatically edit segmentation conditions for gonosomes to prevent the calling of loss/gain due to sex-flip. It should be noted that if female controls are used for male samples, chromosome Y will be excluded from the analysis.


## 10. Master files

### 10.1. Base master file
`/Masters/setting_base.sh`  
Paths to software, main output folders, and basic parameters for the identification of noisy and rare SNPs in the gnomAD database

### 10.2. Projects master file
`/Masters/master_projects.txt`  
Main information for each project

**(!)** = requires explicit values  
**(!-)** = requires explicit values in specific situations  
**(D)** = accepting `default`

Columns:
- `INCLUDE` **(!)** - `yes`, `no` - helps to regulate what project(s) will be executed
  - `yes` - project will be executed
  - `no` - project will not be executed
  - It is recommended to run one project per time.
- `PROJECT_ID` - short and unique identifier for the project
- `PROJECT_TYPE` **(!)** - `germline`, `tumor`, `other` - important value for variant origin prediction in the final report; also used to define abnormal segments smoothing during segmentation
  - `germline` - final report is generated for a pedigree project as heredity prediction; sample 1 is a proband and samples 2 and 3 are parents (works also if only one parent is available)
  - `tumor` - final report is generated for a cancer project as germline/somatic prediction (paired germline has to be available)
  - `other` - no report is generated; anything else 
- `SEQ_TYPE` **(!)** - `WGS`, `WES` , `TS` - type of the sequencing capture
  - `WGS` - whole genome sequencing
  - `WES` - whole exome sequencing
  - `TS` - targeted sequencing
- `CAPTURE_ID` - short identifier for the capture (examples: nextera, twist, sureselect_V4, wgs, wgs_10X etc.)
- `CAPTURE_FILE` - path to the capture file
  - `path/to/file` - required for any sequencing with targets (TS or WES)
  - for WGS, this column is not used and can be left blank or contains anything like not_used, na etc.
- `GENOME_VERSION` **(!)** - `GRCh37-hg19`, `GRCh38-hg38`, `CHM13v2.0`
  - `GRCh37-hg19` - for any version of human genome assembly such as GRCh37, hg19 and b37
  - `GRCh38-hg38` - for any version of human genome assembly such as GRCh38 and hg38
  - `CHM13v2.0` - for T2T-CHM13v2.0 genome assembly
- `REF` - path to fasta file
- `READ_COUNT_MODE` **(D)** **(!)** = `F`, `DF` - introduced with T2T-CHM13v2.0 genome assembly to enable collect low mapping quality reads
  - `default` = `F`
  - `F` - GATK recommended setting when all read filters are enabled
  - `DF` - enables collect reads with "-DF "MappingQualityReadFilter" argument. Low quality reads are collected which enable us to see regions like centromeres.
- `GENE_ANNOTATION` **(D)** - path to gene database that will be used for annotation and plots
  - `default` - refseq database will be used, provided in `/Databases/GENE_ANNOTATION/`
  - `/path/to/file` - custom file can be prepared, see `Custom gene annotation` section below
- `BIN` **(D)** - numeric value - length (bp) of intervals for processing, used in WGS or when capture targets are too long; if no bins are required, value is `0`
  - `default` - automatic use of `1000`
- `PADDING` **(D)** - numeric value - extension (bp) of each interval for depth collection
  - `default` - automatic use of `0` for WGS and read length for WES and TS (read length is calculated from the first control sample)
- `WAY_TO_BAM` **(D)** **(!)** - `absolute`, `find_in_dir`
  - `default` = `absolute`
  - `absolute` - BAM file for each control/sample is provided as an absolute path. The final absolute path consists of `CTRL_BAM_DIR` (projects master file) + `CTRL_PATH_TO_BAM` (controls master file) for controls and `SMPL_BAM_DIR` (projects master file) + `SAMPLE_PATH_TO_BAM` (samples master file) for samples.
  - `find_in_dir` - BAM file for each control/sample is found recursively in folder `CTRL_BAM_DIR`/`SMPL_BAM_DIR` using control/sample identifier (controls/samples master file) and shared pattern `CTRL_BAM_PATTERN`/`SMPL_BAM_PATTERN` (projects master file). Importantly, this pattern has to fit only one BAM file in the provided folder. This setting is not primarily recommended, but it can be very useful when BAM files are in folders with many subfolders and absolute path would be very long and not very enjoyable to specify for each control/sample. The final file is found as `CTRL_ID*CTRL_BAM_PATTERN`/`SAMPLE1_ID*SMPL_BAM_PATTERN`. Patterns can be, for example, `*.bam`, `*_final.bam`, etc.
- `CTRL` **(!)** - `yes`, `no`
  - `yes` - controls will be included in the analysis
  - `no` - controls will not be included in the analysis; coverage profile will be denoised only by GC correction
- `CTRL_BAM_DIR` - path to the folder with BAM file of controls
- `CTRL_BAM_PATTERN` **(D)** - pattern (file suffix) for the BAM files of controls (examples: `*.bam`, `*_final.bam`, etc.); used only when `WAY_TO_BAM` is `find_in_dir`
  - `default` = `na` (which means it is not used as `WAY_TO_BAM` `default` is `absolute`)
- `CTRL_SEX_TEST` **(D)** **(!)** - `custom`, `no` - see Sex test for more information
  - `default` - sex test of controls will be performed using the default setting (*SRY* and *GAPDH* genes coverage)
  - `custom` - sex test of controls will be performed using custom genomic locations defined within the sex master file under `PROJECT_ID`
  - `no` - sex test of controls will not be performed (for example because chromosome Y has no targets, sex of controls is unknown or user just does not require this test to be done)
- `CTRL_PON_SEX_SELECT` **(D)** **(!)** - `M`, `F`, `matched`, `mixed` - During pre-processing, controls are denoised by each other and segmented. This allows us to look at controls qc and the number of segments with a given segmentation setting, and possibly exclude controls with extreme values if wanted. Denoised data is also further used to determine CNVs in control in annotation and plots. Controls itself and controls with the same `MAIN_ID` are excluded from denoising. Setting `matched` is recommended, others can be used in case of a limited number of controls and/or controls with each sex.
  - `default` = `matched`
  - `M` - only male controls are used for denoising of controls
  - `F` - only female controls are used for denoising of controls
  - `matched` - male controls are used for denoising of male controls and female controls are used for denoising of female controls
  - `mixed` - male and female controls are mixed together for denoising of samples; NOT RECOMMENDED setting for two reasons: 1) gonozomes cannot be analyzed, and 2) mixing sex was found to increase false positivity rate for autosomes during pipeline testing (one single-sex control provided better sensitivity/specificity than six mixed controls). This setting is recommended to be used only if the sex is unknown and cannot be determined from Y capture.
- `SMPL_BAM_DIR` - path to the folder with BAM file of samples
- `SMPL_BAM_PATTERN` **(D)** - pattern (file suffix) for the BAM files of samples (examples: `*.bam`, `*_final.bam`, etc.); used only when `WAY_TO_BAM` is `find_in_dir`
  - `default` = `na` (which means it is not used as `WAY_TO_BAM` `default` is `absolute`)
- `SMPL_PON_SEX_SELECT` **(D)** **(!)** - `M`, `F`, `both_separated`, `matched_main`, `matched_each`, `mixed`, `all_options`. This value allows specifying what controls are used for denoising of samples. Controls with the same `MAIN_ID` as the sample are automatically excluded from denoising. Selection depends on user's requirements and the number of controls with each sex available.
  - `default` = `matched_each`
  - `M` - only male controls are used for denoising of samples
  - `F` - only female controls are used for denoising of samples
  - `both_separated` - `M` and `F` will be performed separately
  - `matched_main` - male or female controls are selected based on the sex of the main sample (proband or tumor sample) and used for all related samples (family members or germline sample of the tumor)
  - `matched_each` - male controls are used for denoising of male samples and female controls are used for denoising of female samples
  - `mixed` - male and female controls are mixed together for denoising of samples; NOT RECOMMENDED setting for two reasons: 1) gonozomes cannot be analyzed, and 2) mixing sex was found to increase false positivity rate for autosomes during pipeline testing (one single-sex control provided better sensitivity/specificity than six mixed controls). This setting is recommended to be used only if the sex is unknown and cannot be determined from Y capture. 
  - `all_option` - `M`, `F`, `mixed`, `matched_main`, and `matched_each` will be performed separately
- `CN_FREQ_IN_CTRLS` **(!)** - `yes`, `no`
  - `yes` - CNV analysis in controls will be included in the annotation and plots
  - `no` - CNV analysis in controls will not be included in the annotation and plots, for example because of no or a very small number of controls.
- `GNOMAD_SELECTION` **(D)** **(!)** - `capture_filter`, `capture_filter_and_ctrl_denois`, `af_filter`, `af_filter_and_ctrl_denois` - the original gnomAD database of SNP that is provided has to be filtered for efficient performing; usually `capture_filter`/`capture_filter_and_ctrl_denois` for TS or WES and `af_filter`/`af_filter_and_ctrl_denois` for WGS. Option `_and_ctrl_denois` requires a reasonable number of controls.
  - `default` - automatically selected `af_filter_and_ctrl_denois` for WGS and `capture_filter_and_ctrl_denois` for WES and TS
  - `capture_filter` - SNPs are filtered for those being in the covered regions (for TS and WES)
  - `capture_filter_and_ctrl_denois` - `capture_filter` + noisy SNPs (recognized in controls) are excluded; noise in controls is found by default parameters in the base master file
  - `af_filter` - SNPs with very rare alternative alleles are excluded; allelic frequency is defined by default parameter in the base master file
  - `af_filter_and_ctrl_denois` - `af_filter` + noisy SNPs (recognized in controls) are excluded; noise in controls is found by default parameters in the base master file
- `SEGMENTATION_ID` **(D)** **(!-)** (see `Data segmentation` section above for more information)
  - `default` = `smart-0.65` for single-sex denoising and `smart-0.85` for mixed-sex denoising
  - smart segmentation - `smart-XX` (without sub-clonal analysis) or `smart-XX-sub-YY` (with subclonal analysis); XX = number between 0.5 (high sensitivity) and 1.0 (high specificity), YY = number between 0.1 and 1.0 (for example, 0.5 indicates that subclonal analysis can detect monoallelic CNV in 50% of cells)
  - custom segmentation - short identifier (for example, `my_segm`, `somatic_segmentation`, etc.) to select required segmentation conditions from the segmentation master file
- `SEGMENTATION_ID_USE` **(D)** **(!)** - `segm_id`, `full` - This controls what identifier of segmentation will be used in filenames.
  - `default` = `segm_id`
  - `segm_id` - RECOMMENDED; segmentation condition will be used in filenames by its short identifier (full version is still kept inside of the segmentation table)
  - `full` - all segmentation parameters in numeric version will be used in filenames; not primarily recommended due to creating long file names that can be difficult in some systems
  - If smart segmentation is used, `segm_id` is selected automatically.
- `NOTE` - any additional information that the user wants to keep with the project

### 10.3. Controls master file
`/Masters/master_controls.txt`  
Controls spreadsheet

**(!)** = requires explicit values

Columns:
- `INCLUDE` **(!)** - `yes`, `no` - helps to regulate which controls will be included and at the same time to keep excluded controls with excluding reason in a note column
  - `yes` - control will be included
  - `no` - control will not be included
- `PROJECT_ID` - short and unique identifier for the project, matching `PROJECT_ID` in the other master files
- `CAPTURE_ID` - short identifier for the capture, matching `CAPTURE_ID` in the other master files
- `GENOME_VERSION` **(!)** - `GRCh37-hg19`, `GRCh38-hg38`, matching `GENOME_VERSION` in the other master files
  - `GRCh37-hg19` - for any version of human genome assembly such as GRCh37, hg19, and b37
  - `GRCh38-hg38` - for any version of human genome assembly such as GRCh38 and hg38
- `MAIN_ID` - controls/samples group identifier; This identifier is the one that decides which controls to exclude for denoising because of relatedness. Therefore, controls/samples that are related have to share the same identifier in the controls and the samples master file.
- `CTRL_ID` - control identifier; If `WAY_TO_BAM` is `find_in_dir`, this identifier has to be part of the BAM file name as `CTRL_ID*CTRL_BAM_PATTERN` (pattern determined in the projects master file)
- `CTRL_SEX` **(!)** - `M`, `F`, `unk` - control sex
  - `M` - control is male
  - `F` - control is female
  - `unk` - sex of control is unknown
- `CTRL_PATH_TO_BAM` - control BAM file path when `WAY_TO_BAM` is `absolute`, the final path for each BAM file is a combination of `CTRL_BAM_DIR` (projects master file) and `CTRL_PATH_TO_BAM`. If `WAY_TO_BAM` is `find_in_dir`, this column is not used.
- `NOTE` - any additional information that the user wants to keep with the control

### 10.4. Samples master file
`/Masters/master_samples.txt`  
Samples spreadsheet

**(!)** = requires explicit values

Columns:
- `INCLUDE` **(!)** - `yes`, `no` - helps to regulate what samples are analyzed and skipped during run
- `PROJECT_ID` - short and unique identifier for the project, matching `PROJECT_ID` in the other master files
- `CAPTURE_ID` - short identifier for the capture, matching `CAPTURE_ID` in the other master files
- `GENOME_VERSION` **(!)** - `GRCh37-hg19`, `GRCh38-hg38`, matching `GENOME_VERSION` in the other master files
  - `GRCh37-hg19` - for any version of human genome assembly such as GRCh37, hg19, and b37
  - `GRCh38-hg38` - for any version of human genome assembly such as GRCh38 and hg38
- `MAIN_ID` - controls/samples group identifier; This identifier is the one that decides what controls to exclude from the denoising because of relatedness. Therefore, controls/samples that are related have to share the same identifier in the controls and the samples master files.
- `SAMPLE1_ID` - sample_1 (main sample) identifier; If `WAY_TO_BAM` is `find_in_dir`, this identifier has to be part of BAM file name as `SAMPLE1_ID*SMPL_BAM_PATTERN` (pattern determined in the projects master file)
- `SAMPLE1_TYPE` - sample_1 (main sample) type; for example *proband*, *child-affected*, *child-unaffected*, *tumor*, etc.
- `SAMPLE1_SEX` - **(!)** - `M`, `F`, `unk` - sample_1 (main sample) sex
  - `M` - sample is male
  - `F` - sample is female
  - `unk` - sex of the sample is unknown
- `SAMPLE1_PATH_TO_BAM` - sample_1 (main sample) BAM file path when `WAY_TO_BAM` is `absolute`, the final path for each BAM file is a combination of `SMPL_BAM_DIR` (projects master file) and `SAMPLE1_PATH_TO_BAM`. If `WAY_TO_BAM` is `find_in_dir`, this column is not used.
- `SAMPLE2_ID` - sample_2 identifier; If `WAY_TO_BAM` is `find_in_dir`, this identifier has to be part of BAM file name as `SAMPLE2_ID*SMPL_BAM_PATTERN` (pattern determined in the projects master file)
- `SAMPLE2_TYPE` - sample_2 (main sample) type; for example *father*, *mother*, *germline*, etc.
- `SAMPLE2_SEX` - **(!)** - `M`, `F`, `unk` - sample_2 sex
  - `M` - sample is male
  - `F` - sample is female
  - `unk` - sex of the sample is unknown
- `SAMPLE2_PATH_TO_BAM` - sample_2 BAM file path when `WAY_TO_BAM` is `absolute`, the final path for each BAM file is a combination of `SMPL_BAM_DIR` (projects master file) and `SAMPLE2_PATH_TO_BAM`. If `WAY_TO_BAM` is `find_in_dir`, this column is not used.
- `SAMPLE3_ID` - sample_3 identifier; If `WAY_TO_BAM` is `find_in_dir`, this identifier has to be part of BAM file name as `SAMPLE3_ID*SMPL_BAM_PATTERN` (pattern determined in the projects master file)
- `SAMPLE3_TYPE` - sample_3 (main sample) type; for example *father*, *mother*, etc.
- `SAMPLE3_SEX` - **(!)** - `M`, `F`, `unk` - sample_3 sex
  - `M` - sample is male
  - `F` - sample is female
  - `unk` - sex of the sample is unknown
- `SAMPLE3_PATH_TO_BAM` - sample_3 BAM file path when `WAY_TO_BAM` is `absolute`, the final path for each BAM file is a combination of `SMPL_BAM_DIR` (projects master file) and `SAMPLE3_PATH_TO_BAM`. If `WAY_TO_BAM` is `find_in_dir`, this column is not used.
- `NOTE` - any additional information that user wants to keep with the sample



### 10.5. Sample and regions of interest to plot master file
`/Masters/master_reg_of_interest_to_plot.txt`  
Spreadsheet that allows to plot required region for required sample; if `PROJECT_ID` is recognized to be running and at least one of the lines of this master file has `INCLUDE` = `yes`, a separate detailed plot for the required sample within `/detail_regions_of_interest/` folder will be produced. It can be useful when a user wants to print a suspicious region that was not recognized as abnormal and therefore not printed as a detailed plot.

**(!)** = requires explicit values

Columns:
- First 17 columns are the same as columns in the samples master file
- `CONTIG` - chromosome of the region to plot
- `START` - starting position of the region to plot
- `END` - ending position of the region to plot
- `ZOOM` - numeric value, zoom for the defined region to be printed to the plot, examples:
  - `0`: `CONTIG`: `START` - `END` will be printed to the plot
  - `1`: `CONTIG`: [`START`-1✗(`END`-`START`)] - [`END`+1✗(`END`-`START`)] will be printed to the plot
  - `2`: `CONTIG`: [`START`-2✗(`END`-`START`)] - [`END`+2✗(`END`-`START`)] will be printed to the plot

### 10.6. Segmentation master file
`/Masters/setting_segmentation.txt`  
Parameters for CNVs and LOH segmentation, used when custom segmentation is selected (see column `SEGMENTATION_ID` in the projects master) and paired by a short identifier.

Columns:
- `FEATURE` - Segmentation parameter
- `SEGM-EXAMPLE` - Example of segmentation setting
- `...` - Additional columns can be added by the user; a short and unique ID of segmentation is required

The segmentation identifier defined in the column header (e.g., *my_segm*, *somatic_segmentation*, etc.) is used to pair the segmentation setting and project as `SEGMENTATION_ID` within the projects master file. This identifier is also used in filenames if `SEGMENTATION_ID_USE` = `segm_id`.

List of parameters:
- `SEGM_DIFFERENCE`
- `SEGM_MINSIZE`
- `SEGM_MINSIZE_SUB`
- `SEGM_MINPROBE`
- `SEGM_MINPROBE_SUB`
- `SEGM_MINKEEP`
- `SEGM_MINLOSS`
- `SEGM_MINLOSS_SUB`
- `SEGM_MINLOSS_BIAL`
- `SEGM_MINGAIN`
- `SEGM_MINGAIN_SUB`
- `SEGM_GAP`
- `SEGM_SMOOTHPERC`
- `SEGM_AFDIF`
- `SEGM_AFSIZE`
- `SEGM_AFPROBE`

Subclonal analysis is defined by parameters ending with `_SUB`. If no subclonal analysis is required, these variables are set to `none`. If `SEGM_GAP` is set to `none`, segments continue independently regardless of their distance from each other in the capture.

### 10.7. Sex master file
`/Masters/sex_test_regions.txt`  
Genomic locations for sex test.

**(!)** = requires explicit values

Columns:
- `CONTIG`
- `START`
- `END`
- `GENE_TYPE` **(!)** - `TEST_GENE`, `CTRL_GENE`
  - `TEST_GENE` - Genomic position of the tested gene 
  - `CTRL_GENE` - Genomic position of the control gene
- `PROJECT_ID` - Short identifier for the project (must match all `PROJECT_ID` values in the other master files)
- `GENE` - Gene symbol

## 11. Outputs

### 11.1. Data outputs
- `...allelicCounts...` - RDS file with SNP zygosity collection
- `...counts...` - HDF5 and TSV files with coverage collection
- `...denoisedCR...` - TSV file with denoised coverage data
- `...standardizedCR...` - TSV file with standardized coverage data
- `...SEGMENTS...segmentation...` - TSV file with segmentation data before annotation, contains all normal and abnormal segments
- `...SEGMENTS...` - TSV file with segmentation and annotation, contains all normal and abnormal segments
- `...SEGMENTS...REPORT...` - TSV file with segmentation, annotation, and variant origin prediction, contains ONLY abnormal segments


### 11.2. Plot outputs
- **genome** - Figure showing all chromosomes (contigs with centromere) extracted from the reference (usually chr1-chr22, chrX, and chrY)
- **chromosome** - Figure is generated for each chromosome
- **detail** - Figure is generated for each abnormal segment (does not include sub-clonal findings)

### 11.3. IGV outputs
- `mainID_sampletype_ID_sex_abnormal_segments.bed` - BED file with abnormal segments

- `mainID_sampletype_ID_sex_segments.bw` - BigWig with segments, `mainID_sampletype_ID_sex_CN.bw` - BigWig with denoised coverage data
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min -2.5, mid 0.0, max 2.5
  - Note: log2 ratio values smaller than -2.25 are shifted to -2.25 and log2 ratio bigger than 2.25 are shifted to 2.25 in IGV output (as well as in the PNG plots)
  - Tip: when these two BigWig files in the same setting, mark them and by right click select "Overlay Tracks"

- `mainID_sampletype_ID_sex_SNP.bw` - BigWig with SNP zygosity
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min 0.0, mid 0.5, max 1.0

## 12. Custom gene annotation
How to prepare gene annotation other than default RefSeq provided:
- Download the required gene annotation from [UCSC website](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=968304787_9UHI0R5esZuD1z4azUHoV9dEdsoU).
- Columns have to be named as follows (order does not need to be the same):
  - `gene_id`: Identifier that will be displayed at plots and used for annotation (usually gene symbol)
  - `chrom`: Chromosome
  - `strand`: + or - for strand
  - `txStart`: Transcription start position
  - `txEnd`: Transcription end position
  - `cdsStart`: Coding region start
  - `cdsEnd`: Coding region end
  - `exonStarts`: Exon start positions
  - `exonEnds`: Exon end positions
 Most of the columns will match automatically, only `gene_id` needs to be picked and renamed by the user, for example, by using the command `sed -i 's/name2/gene_id/g' file.txt`. The table should never be edited in EXCEL due to renaming some gene symbols to dates. Contig IDs have to match IDs in the reference file. 

## 13. Project-specific plot
This is an advanced option that allows the production of personalized detail plots for additional regions with personalized annotation.
A customized R script is required, provided in the example. This script has to be placed in the `/CNVRobot_vX.X/Scripts/R/` folder and named as `plot_CN_PROJECT_ID.R`. It is then recognized as a source in `plot_CN.R` and produces plots within the `/project_specific_regions/` folder.

## 14. Databases source
- Processing of the databases from the original source is available upon request.
- Original sources of databases:
  - gnomAD SNP:
    - hg19/GRCh37 and hg38/GRCh38: [gnomAD SNP v2.1](https://gnomad.broadinstitute.org/downloads) (exome and genome data was downloaded, combined, and processed into a simple VCF file in hg19/GRCh37 version; hg38/GRCh38 version prepared by liftover)
    - CHM13v2.0: [T2T aws](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/liftover/) - `chm13v2.0_dbSNPv155.vcf.gz` (filtered for gnomAD SNP only and processed into a simple VCF file)
  - [Database of Genomic Variants (DGV, TCAG) v2020-02-25](http://dgv.tcag.ca/dgv/app/home) (hg19/GRCh37 and hg38/GRCh38 prepared using files from the website, CHM13v2.0 prepared by liftover from hg38/GRCh38 version, unlifted regions excluded)
  - Chromosome bands:
    - hg19/GRCh37 and hg38/GRCh38: [UCSC chromosome bands and centromeres](https://hgdownload.soe.ucsc.edu/downloads.html)
    - CHM13v2.0: [T2T aws](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/) - `chm13v2.0_cytobands_allchrs.bed`
  - Chromosome centromeres: Border position between p arm and q arm of each chromosome, generated from chromosome bands
  - Gene annotation:
    - hg19/GRCh37 and hg38/GRCh38: [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) - `RefSeq All (ncbiRefSeq)` (under Genes and Gene Predictions, NCBI RefSeq)
    - CHM13v2.0: [T2T UCSC](https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.catLiftOffGenesV1/) - `catLiftOffGenesV1.bb` (processed to UCSC-looking format)
  - Mappability:
    - GRCh37/hg19: [ENCODE from UCSC](https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) - `wgEncodeCrgMapabilityAlign100mer.bigWig`
    - GRCh38/hg38: [Hoffman Lab from UCSC](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=mappability) `k100.Umap.MultiTrackMappability.bw`
    - CHM13v2.0: [T2T aws](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/mappability/) - `chm13v2.0_100mer_mappability.bb` (note this data is binary and not quantitative as the others)
  - ENCODE black list:
    - GRCh37/hg19: [UCSC Table Browser](https://genome-euro.ucsc.edu/cgi-bin/hgTables) - `ENCODE Blacklist (encBlackList)` (under Mapping and Sequencing, Problematic Regions)
    - GRCh38/hg38 and CHM13v2.0: liftover from GRCh37/hg19, unlifted regions excluded

## 15. Limitations
The CNVRobot pipeline is not suitable for the detection of common population CNVs and CNVs in low mappability regions. It should be noted that it is a depth of coverage analysis within certain bins and thus starts and ends of abnormalities are not precisely mapped. CNVRobot can detect only unbalanced changes. Similarly to DNA microarrays, coverage-based analysis can be problematic if the ploidy of the sample is changed, as it is based on median centering normalization.

## 16. Contact
Aneta.Mikulasova@newcastle.ac.uk
