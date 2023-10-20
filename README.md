# CNVRobot

<p align="left">
  <img src="https://github.com/AnetaMikulasova/CNVRobot/blob/main/CNVRobot_logo.png" alt="CNVRobot logo" width="200" height="207"/>
</p>

Welcome to CNVRobot v4.2!

## 1. Description
CNVRobot is an integrated pipeline designed to detect rare germline and somatic copy-number variants (CNVs) as well as loss of heterozygosity (LOH) regions in the human genome using short-read DNA sequencing from any NGS platform. This includes targeted sequencing (TS), whole-exome sequencing (WES), and whole-genome sequencing (WGS). It integrates sequencing depth with SNP zygosity across the genome, incorporates advanced data denoising, segmentation and variant prioritization options, and provides detailed visualization for manual inspection of detected CNVs.

## 2. Download
The code, databases, and example data are available for download [HERE](https://newcastle-my.sharepoint.com/:f:/g/personal/nam320_newcastle_ac_uk/Eu2m1nidAKVDoXxuDONbqToBJfbYPwpbjeRuU6ERhcqHrg?e=yshGbA). The latest version is CNVRobot v4.2.
Please note that the databases are not available here on GitHub.

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
- only for [CNVkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html) analysis: see [CNVkit github](https://github.com/etal/cnvkit) for instalation via Conda and setting up a new Python environment 
- only for [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software) analysis: R packages [ASCAT](https://github.com/VanLoo-lab/ascat), [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) (Artistic-2.0) and [plyranges](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html) (Artistic-2.0)


## 4. Inputs
- BAM files
- FASTA file - the reference file used for alignment (possible assemblies: hg19, GRCh37, hg38, GRCh38, or CHM13v2.0)
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
You can browse `/Test_example/Test_example_output` to see the output from this example data, including figures and IGV files.


## 7. Control set selection
GATK recommends that the control set consists of at least 40 samples without pathogenic CNVs. These samples should be prepared, sequenced, and analyzed using the same laboratory protocols and bioinformatic processes. However, even if you only have 2–5 controls processed simultaneously, it is still better than having no control at all. Even a small number of controls can already be highly effective in reducing noise from technical artifacts. The control set can overlap with the sample set, as long as no pathological CNVs are expected. For instance, the control set can comprise unaffected parents in a pedigree dataset or germline samples in a cancer dataset. By default, any related samples are automatically excluded from the controls during analysis. This is done by excluding all controls that match the main identifier `MAIN_ID` provided in both the sample and control spreadsheet.

### 7.1 Denoising tumor sample by correspoinding germline only
CNVRobot offers an option to denoise a tumor sample using its corresponding germline sample. The resulting profile might be noisier (as only one control is technically used), but it can provide a better understanding of what is somatic and what is germline. This option can be configured under `SETTING_MODE` in the `/Masters/master_projects.txt`. See the `Specific project setting options` section for further details.

## 8. Sex and gonosomes
For the accurate analysis of CNV on chrX and chrY, denoising using controls and segmentation are carefully processed. To detect gonosomal CNVs correctly, controls of the same sex are required for GATK data denoising. This means that `SMPL_PON_SEX_SELECT` within `/Masters/master_projects.txt` should be `M`, `F`, `both_separated`, `matched_main`, or `matched_each`. The CNVRobot segmentation algorithm will always recognize which control sex was used against which sample sex and adjust log2 ratio cut-offs for abnormalities accordingly. For chrY analysis in males, the Panel of Normals should consist of male controls. If the Panel of Normals consists of female controls, chrY data will be excluded from the analysis by GATK algorithms. If `SMPL_PON_SEX_SELECT` is mixed, the analysis of chrX will still be correct (due to coverage doubling in males, see below), but chrY cannot be reported accurately as it will be influenced by the ratio of males and females in the GATK Panel of Normals.
Note that CNVRobot has a module for denoising data using CNVkit, which can be beneficial for the analysis of chrY if a small control set is available and/or includes female controls only. See the `Specific project setting options` section for further details.
Before including a control in the Panel of Normals, its sex is verified or predicted if unknown. Sex tests and predictions do not run for samples, as gonosomal aneuploidies can be present.

### 8.1 Sex test of controls
Each control is tested to ensure the sex is correctly determined by the `CTRL_SEX` value in `/Masters/master_controls.txt`. This test uses global coverage (normalized for genomic size) of gonosomes compared to autosomes. If there are no gonosomal sequencing targets, the test will be skipped. If there are no autosomal sequencing targets, but both gonosomes have sequencing targets, the test will use the ratio of chrY to chrX. If the predicted sex differs from what was provided by the user, the pipeline will halt and report this error.

### 8.2 Sex prediction for controls with unknown sex
If `CTRL_SEX` is defined as `unk` in `/Masters/master_controls.txt`, the same test will be run to predict the sex, and the control will be further processed using this predicted sex.

### 8.3. Gonosomal analysis
By default, the collected coverage of chrX and chrY is now doubled in male samples and controls. This step was introduced in response to observations of GATK "over-denoising" gonosomes relative to autosomes when the Panel of Normals is composed of males (haploid chrX and chrY). The normal (diploid) values for log2 ratios are as follows:
- F sample denoised by F controls: chrX = 0, chrY = no data
- M sample denoised by F controls: chrX = 0, chrY = no data
- F sample denoised by M controls: chrX = 0, chrY = deep minus
- M sample denoised by M controls: chrX = 0, chrY = 0

The segmentation algorithm automatically adjusts log2 ratio cut-offs for gonosomes based on the sample/control sex. Therefore, chrY will not be reported as deleted in a female sample if the Panel of Normals consists of male samples.

If, for any reason, a user prefers the original analysis without doubling the male gonosomal coverage, it is possible to specify this for your project run under `SETTING_MODE` in `/Masters/master_projects.txt`. Refer to the `Specific project setting options` section for details. Note that the segmentation algorithm will detect this and adjust log2 ratio cut-offs for gonosomes accordingly, avoiding reports of changes due to sex-flips. In this setting, the normal (diploid) values for log2 ratios are as follows:
- F sample denoised by F controls: chrX = 0, chrY = no data
- M sample denoised by F controls: chrX = -1, chrY = no data
- F sample denoised by M controls: chrX = 1, chrY = deep minus
- M sample denoised by M controls: chrX = 0, chrY = 0



## 9. Data segmentation
The segmentation of denoised coverage and SNP zygosity data is an important step to recognize abnormal segments (losses, gains, and LOH).
### 9.1. Segmentation conditions
The pipeline offers a default segmentation setting. However, results can be reviewed and segmentation conditions adapted by the user based on data type, quality, and expectation. There are two options for customizing segmentation conditions using the project's master file `SEGMENTATION_ID` column:  
a) **smart segmentation**: The identifier pattern `smart-XX` is provided within `SEGMENTATION_ID` (in the project's master file), where `XX` represents a segmentation coefficient, a number between 0.5 (for high sensitivity) and 1 (for high specificity). If smart segmentation is used, segmentation conditions are not loaded from the segmentation master file, but are instead calculated automatically by the given segmentation coefficient and detected sample quality. For subclonal analysis, the pattern is `smart-XX-sub-YY`, with `YY` representing a number between 0.1 and 1.0 (for example, 0.5 indicates a subclonal analysis to detect monoallelic CNV in 50% of cells).  
b) **custom segmentation**: A short identifier (for example, *my_segm*, *somatic_segmentation*, etc.) is provided within `SEGMENTATION_ID` (in the project's master file). This identifier matches a unique column name within the segmentation master file, where all segmentation conditions are manually defined by the user.



### 10. Specific project setting options
CNVRobot offers several run options beyond the default settings that can be useful for specific situations. These settings are found under `SETTING_MODE` in `/Masters/master_projects.txt`. Due to their complexity, they are explained in this specific section rather than in the `/Masters/master_projects.txt` description. The `SETTING_MODE` can simply be set to `default`, in which case CNVRobot will operate in its standard setting. 

However, `SETTING_MODE` consists of seven positions, each of which can switch the run to a specific project setting.

**Position 1 - Low mapping quality read count collection**: `F`/`-` (default) or `D`:
- Introduced with the T2T-CHM13v2.0 genome assembly to enable the collection of low mapping quality reads.
  - `F`/`-` (default) - This is the GATK-recommended setting when all read filters are enabled. While `F` and `-` are equivalent, if Microsoft Excel is used to fill the master files, `F` should be used since the `-` symbol can cause issues when placed at the beginning.
  - `D` - Enables the collection of reads with the "-DF "MappingQualityReadFilter" argument. This collects low-quality reads, allowing us to visualize regions like centromeres using GATK algorithms.

**Position 2 - Matching germline denoising of tumor sample**: `-` (default) or `G`
  - `-` (default) - The Panel of Normals is created based on the list of controls found within `/Masters/master_controls.txt` and by the defined sex (`SMPL_PON_SEX_SELECT` in `/Masters/master_projects.txt`). The matching germline is excluded from the denoising. If listed as `SAMPLE2` in `/Masters/master_samples.txt`, it will be processed alongside the tumor with the same Panel of Normals.
  - `G` - The tumor sample will be denoised using only the corresponding germline sample. This germline must be listed in `/Masters/master_controls.txt` with the same `MAIN_ID` identifier as the tumor sample in `/Masters/master_samples.txt`. If the germline is also listed as `SAMPLE2` in `/Masters/master_samples.txt`, it will be processed next to the tumor and denoised by the Panel of Normals as in the default setting.

**Position 3 - Not-doubling gonosomal coverage**: `-` (default) or `A` (refer to `Gonosomal analysis` section for context)
  - `-` (default) - The collected coverage of chrX and chrY is doubled in male samples and controls.
  - `A` - The collected coverage of chrX and chrY IS NOT doubled in male samples and controls (retains the original GATK setting).

**Position 4 - Speeding up CNVRobot by skipping certain steps**: `-` (default) or `Q`: 
  - `-` (default) - No steps are skipped.
  - `Q` - Skips the following: segmentation and QC metrics for controls, as well as annotation and reports for samples.

**Position 5 - Speeding up CNVRobot plotting by omitting detailed plots**: `-` (default) or `Q`
  - `-` (default) - A figure is generated for each abnormal segment (excluding sub-clonal findings).
  - `Q` - No figure is generated for abnormal segments. Genome and chromosome figures are still generated.

**Position 6 - Running [CNVkit](https://cnvkit.readthedocs.io/en/stable/quickstart.html)** (tested with v0.9.10): `-` (default) or `Y` 
  - `-` (default) - CNVkit is silent.
  - `Y` - CNVkit runs alongside GATK. Note that it needs to be correctly installed in a specific Python environment and requires paths to `CONDA` and `CNVKIT_ENV` (environment name) specified in `/Masters/setting_base.txt`. CNVkit is used for coverage collection and denoising and may offer advantages for gonosomal analysis. Data is segmented using the CNVRobot algorithm. Outputs are similar to those produced by GATK, including plots and IGV files.

  The normal (diploid) values for log2 ratios processed by CNVkit are as follows:
  - F sample denoised by F controls: chrX = 0, chrY = deep minus
  - M sample denoised by F controls: chrX = 0, chrY = 0
  - F sample denoised by M controls: chrX = 0, chrY = deep minus
  - M sample denoised by M controls: chrX = 0, chrY = 0

**Position 7 - Running [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software)** (tested with v2.5.2): `-` (default) or `Y`
  - `-` (default) - ASCAT is silent.
  - `Y` - ASCAT is running. This allows output from CNVRobot to be directly processed by ASCAT to determine ploidy and analyze absolute copy-number data in a tumor sample. Please note that this mode is still a work in progress and requires both a tumor and a matching germline sample. The germline sample must be specified as `SAMPLE2` in `/Masters/master_samples.txt` to be included in the analysis. Outputs are stored in a separate folder within each sample and include regular ASCAT outputs as well as IGV files with absolute CN for segments (`/IGV/*_absCN*`) from CNVRobot, and ASCAT CN segmentation for major and minor alleles (`/IGV/*_allele.bw`).


Examples of `SETTING_MODE`:
- `default`: All defaults apply - only reads with high mappability are collected by GATK, Panel of Normals is used for denoising, gonosomal coverage is doubled in males, no steps are skipped, and neither CNVkit nor ASCAT are run.
- `D` (equivalent to `D------`): Reads with low mapping quality will be collected by GATK. Other settings will remain as default.
- `-G` (equivalent to `-G-----`/`FG-----`): The tumor will be denoised by its matching germline. Other settings will remain as default.
- `---QQ` (equivalent to `---QQ--`/`F--QQ--`): CNVRobot will skip certain steps during processing and will not generate a plot for each abnormality. Other settings will remain as default.
- `D-A` (equivalent to `D-A----`): Reads with low mapping quality will be collected by GATK, and gonosomal coverage will not be doubled in males. Other settings will remain as default.
- `-----Y` (equivalent to `-----Y-`/`F----Y-`): CNVRobot runs CNVkit in addition to default settings. 
- `------Y` (equivalent to `F-----Y`): CNVRobot runs ASCAT in addition to default settings. 

Important: If Microsoft Excel is used to fill the master files, the pattern starting `F` (not `-`) should be used. Symbol `-` at the beggining can cause issues.
Please note that the `SETTING_MODE` component of CNVRobot is still under development. If any aspect is unclear, do not hesitate to [contact me](mailto:Aneta.Mikulasova@newcastle.ac.uk) for further details.



## 11. Master files

### 11.1. Base master file
`/Masters/setting_base.sh`  
Paths to software, main output folders, and basic parameters for the identification of noisy and rare SNPs in the gnomAD database

### 11.2. Projects master file
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
- `SETTING_MODE` **(D)** **(!)** = see `Specific project setting options section`
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

### 11.3. Controls master file
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

### 11.4. Samples master file
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

### 11.5. Sample and regions of interest to plot master file
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

### 11.6. Segmentation master file
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



## 12. Outputs

### 12.1. Data outputs
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/data/`  
- `...allelicCounts...rds` - RDS file with SNP zygosity collection
- `...counts...hdf` - HDF5 file with coverage collection
- `...denoisedCR...tsv` - TSV file with denoised coverage data
- `...SEGMENTS...segmentation...tsv` - TSV file with segmentation data before annotation, contains all normal and abnormal segments
- `...SEGMENTS...tsv` - TSV file with segmentation and annotation, contains all normal and abnormal segments
- `...SEGMENTS...REPORT...tsv` - TSV file with segmentation, annotation, and variant origin prediction, contains ONLY abnormal segments

### 12.2. Plot outputs
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/plot_SETTING/plot/`
- **genome** - Figure showing all chromosomes (contigs with centromere) extracted from the reference (usually chr1-chr22, chrX, and chrY)
- **chromosome** - Figure is generated for each chromosome
- **/detail/** - Folder with figures generated for each abnormal segment (does not include sub-clonal findings)

### 12.3. IGV outputs
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/plot_SETTING/IGV/`  

- `..._abnormal_CN_segments.bed` - BED file with segments for DNA losses and gains
- `..._abnormal_LOH_segments.bed` - BED file with segments for LOH (any regions with deviation from heterozygosity, which can be due to deletion, cnnLOH or gain); generated only if LOH region(s) were detected

- `..._segments.bw` - BigWig with all segments (normal and abnormal), `..._CN.bw` - BigWig with denoised coverage data
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min -2.5, mid 0.0, max 2.5
  - Note: log2 ratio values smaller than -2.25 are shifted to -2.25 and log2 ratio bigger than 2.25 are shifted to 2.25 in IGV output (as well as in the PNG plots)
  - Tip: when these two BigWig files in the same setting, mark them and by right click select "Overlay Tracks"

- `.._SNP.bw` - BigWig with SNP zygosity
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min 0.0, mid 0.5, max 1.0

### 12.4. CNVkit
#### 12.4.1. Data outputs
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/data/`  
- `...(anti)targetcoverate.cnn` - CNN file with coverage collection
- `...cnr` - denoised coverage data
- `...SEGMENTS...CNVkit.tsv` - TSV files with segmentation, same as above, just generated using CNVkit data
#### 12.4.2 Plot and IGV outputs
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/plot_SETTING/cnvkit_plot/`
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/plot_SETTING/cnvkit_IGV/`  
- same as above, just generated using CNVkit data

### 12.5. ASCAT
`/CNVRobot/Processing/Results/PROJECT_ID/PROJECT_ID_SETTING/ascat.../`
- contains data outputs from ASCAT (`/data/`) and IGV files (`/IGV/`)
- note, this is still a work in progress


## 13. Custom gene annotation
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

## 14. Project-specific plot
This is an advanced option that allows the production of personalized detail plots for additional regions with personalized annotation.
A customized R script is required, provided in the example. This script has to be placed in the `/CNVRobot_vX.X/Scripts/R/` folder and named as `plot_CN_PROJECT_ID.R`. It is then recognized as a source in `plot_CN.R` and produces plots within the `/project_specific_regions/` folder.

## 15. Databases source
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

## 16. Limitations
The CNVRobot pipeline is not designed for detecting common population CNVs and may produce errors in regions with low mappability. Users should be aware that CNVRobot performs a depth of coverage analysis within specified bins; therefore, the exact starting and ending points of abnormalities are not precisely determined. The tool is only capable of detecting unbalanced changes. Much like DNA microarrays, coverage-based analysis can pose challenges when there's a change in the sample's ploidy, given that it relies on median centering normalization.

## 17. Contact
Aneta.Mikulasova@newcastle.ac.uk
