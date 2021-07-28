# CNVRobot

<p align="left">
  <img src="https://github.com/AnetaMikulasova/CNVRobot/blob/main/CNVRobot_logo.png" alt="CNVRobot logo" width="200" height="207"/>
</p>

Welcome to CNVRobot v3.3!

## 1. Description
CNVRobot is an integrated pipeline designed to detect rare germline and somatic copy-number variants (CNVs) and loss of heterozygosity (LOH) regions in the human genome using short-read DNA sequencing from any NGS platform, including targeted sequencing (TS), whole-exome sequencing (WES) and whole-genome sequencing (WGS). It integrates sequencing depth with SNP zygosity across the genome, incorporates advanced data denoising, segmentation and variant prioritization options, as well as detailed visualization for manual inspection of detected CNVs.

## 2. Download
Code, databases and example data can be downloaded here: [CNVRobot_v3.3.zip](https://newcastle-my.sharepoint.com/:u:/g/personal/nam320_newcastle_ac_uk/EWjjB9vnEfJEvhKk1mwT6PYBNuB01tQ60-EJE-ML9H1RQw?e=SQaJzD)

## 3. Dependencies (License)
- [GATK v4.2+](https://github.com/broadinstitute/gatk/releases) (Apache License 2.0)
- [PICARD](https://broadinstitute.github.io/picard/) (MIT License)
- [R](https://www.r-project.org) (GPL-2 | GPL-3)
- [Samtools, HTSlib (bgzip and tabix), BCFtools](http://www.htslib.org) (MIT/Expat License)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (GNU General Public License v2.0)
- [bedGraphToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/) (UCSC Genome Browser Licence)
- [dos2unix](http://dos2unix.sourceforge.net) (BSD License)
- [java](https://www.java.com/) (Oracle Binary Code License)
- [gzip/gunzip](https://www.gnu.org/software/gzip/) (GNU General Public License)
- R packages: 
    - [tidyverse](https://www.tidyverse.org) (MIT License)
    - [expss](https://cran.r-project.org/web/packages/expss/index.html) (GPL-2 | GPL-3)
    - [splitstackshape](https://cran.r-project.org/web/packages/splitstackshape/index.html) (GPL-3)
    - [karyoploteR](http://bioconductor.org/packages/release/bioc/html/karyoploteR.html) (Artistic-2.0)

## 4. Inputs
- BAM files
- FASTA file - reference file used for alignment
- capture BED file (for WES or TS)
  - first three columns required - contig, start and end; no header
  - requires the same genome assembly as reference file
- master files - spreadsheets of project, control, samples etc.; see below for details

## 5. Initiate run
1. Fill ```./Masters/setting_base.sh``` to specify paths for software.
2. Fill ```./Masters/master_projects``` to specify condition for your project.
3. Fill ```./Masters/master_samples.txt``` and ```./Masters/master_controls.txt``` to specify samples and controls sets for the project.
4. Run ```./CNVRobot/Scripts/run.sh```.

- CNVRobot will link samples and controls for each project by matching ```PROJECT_ID```, ```CAPTURE_ID``` and ```GENOME_VERSION```.
- First column in each master file (```INCLUDE```)  helps to regulate what project to run, what samples to analyze and what control to use.
- Details about each master file are stated below.
- Please avoid space character in any identifier.

## 6. Example data
To run this testing data, please fill ```./Masters/setting_base.sh``` and run ```./CNVRobot/Scripts/run.sh```.  
Simulated example data contains one sample, two controls and reference file.  
Platform: targeted sequencing  
Chromosome locus: 21q21  
Genomic position: chr21:15.0-16.0Mb  
Genome assembly: hg38  
CNVs in the sample: loss of chr21:15.2-15.3Mb, gain of chr21:15.7-15.8Mb  

## 7. Control set selection
GATK recommends the control set consists of at least 40 samples without pathogenic CNVs that were prepared, sequenced and analyzed by the same laboratory protocols and bioinformatic processing. However, any simultaneously processed control is better than no control and the use of just 2–5 controls can already be highly effective in reducing the noise from technical artefacts. The control set can overlap with sample set if no pathological CNV are expected. For example, the control set can consist of unaffected parents in pedigree dataset or germline samples in a cancer dataset. Any related samples are automatically excluded from controls during analysis using the main identifier stated in the sample and control spreadsheet. 

## 8. Sex test
Sex test is a part of preparation script. It compares a local coverage of two genes, tested and control, to verify the sex of each control. This helps to recognize potential sample swop that could disturb normalization, denoising and analysis of gonozomes. Genomic coordinates for collection are provided by the sex master file and test on/off is regulated in the project master file. Default setting involves *SRY* and *GAPDH* genes, but requires both genes to be covered. If default genes are not covered by the sequencing, alternative loci can be specified in the sex master file. Pairing customized loci with the project is defined by ```PROJECT_ID```. If customized setting is used, tested and control loci have to be in the same size (bp). Sex test may not be possible to perform in case of some targeted panels where no target at chromosome Y is included.

## 9. Data segmentation
Segmentation of denoised coverage and SNP zygosity data is an important step to recognize abnormal segments (losses, gains and LOH).
### 9.1. Segmentation conditions
The pipeline offers default segmentation setting. However, results can be reviewed and segmentation conditions adapted by user based on data type, quality and expectation. There are two options how to customise segmentation conditions using the project master file ```SEGMENTATION_ID``` column:  
a) **smart segmentation**: ```smart-XX``` identifier pattern is provided within ```SEGMENTATION_ID``` (project master file); where ```XX``` stays for a segmentation coeficient, a number between 0.5 (high sensitivity) and 1 (high specificity). If smart segmentation is used, segmentation conditions are not loaded from the segmentation master file, but calculated automatically by given segmentation coeficient and detected sample quality. For subclonal analysis, pattern is ```smart-XX-sub-YY``` with ```YY``` for a number between 0.1 and 1.0 (for example, 0.5 indicates that subclonal analysis to detect monoallelic CNV in 50% of cells).  
b) **custom segmentation**: Short identifier (for example *my_segm*, *somatic_segmentation* etc.) is provided within ```SEGMENTATION_ID``` (project master file). This identifier is matching with an unique column name within the segmentation master file where all segmentation conditions are defined manually by user.
### 9.2. Sex-flip 
For calling of gonozomal abnormalities, controls with the same sex are required for data denoising. The pipeline runs in sex-flip recognizing mode. If opposite sex is used for denoising (defined in the project master file), pipeline will automatically edit segmentation conditions for gonozomes accordingly, not to call loss/gain due to sex-flip. It should be noted that if female controls are used for male sample, chromosome Y will be excluded from the analysis.

## 10. Master files

### 10.1. Base master file
```/Masters/setting_base.sh```  
paths to software, main output folders and basic parameters for identification of noisy and rare SNPs in gnomAD database

### 10.2. Projects master file
```/Masters/master_projects.txt```  
main information for each project

**(!)** = requires explicit values  
**(!-)** = requires explicit values in specific situation  
**(D)** = accepting ```default```

Columns:
- ```INCLUDE``` **(!)** - ```yes```, ```no``` - helps to regulate what project(s) will be executed
  - ```yes``` - project will be executed
  - ```no``` - project will not be executed
- ```PROJECT_ID``` - short and unique identifier for the project
- ```PROJECT_TYPE``` **(!)** - ```germline```, ```tumor```, ```compare_independent``` - important value for variants origin prediction in the final report
  - ```germline``` - final report is generated for pedigree project as heredity prediction; sample 1 is a proband and samples 2 and 3 are parents (works also if only one parent is available)
  - ```tumor``` - final report is generated for cancer project as germline/somatic prediction (paired germline has to available)
  - ```compare_independent``` - no report is generated; anything else 
- ```SEQ_TYPE``` **(!)** - ```WGS```, ```WES``` , ```TS``` - type of the sequencing capture
  - ```WGS``` - whole genome sequencing
  - ```WES``` - whole exome sequencing
  - ```TS``` - targeted sequencing
- ```CAPTURE_ID``` - short identifier for the capture (examples: nextera, twist, sureselect_V4, wgs, wgs_10X etc.)
- ```CAPTURE_FILE``` - path to the capture file
  - ```path/to/file``` - required for any sequencing with targets (TS or WES)
  - for WGS, this column is not used and can be left blank or containg anything like not_used, na etc.
- ```GENOME_VERSION``` **(!)** - ```GRCh37-hg19```, ```GRCh38-hg38```
  - ```GRCh37-hg19``` - for any version of human genome assembly such as GRCh37, hg19 and b37
  - ```GRCh38-hg38``` - for any version of human genome assembly such as GRCh38 and hg38
- ```REF``` - path to fasta file
- ```GENE_ANNOTATION``` **(D)** - path to gene database that will be used for annotation and plots
  - ```default``` - refseq database will be used, provided in ```/Databases/GENE_ANNOTATION/```
  - ```/path/to/file``` - custom file can be prepared, see ```Custom gene annotation``` section below
- ```BIN``` **(D)** - numeric value - length (bp) of intervals for processing, used in WGS or when capture targets are too long; if no bins are required, value is ```0```
  - ```default``` -  automatic use of ```1000```
- ```PADDING``` **(D)** - numeric value - extension (bp) of each interval for depth collection
  - ```default``` - automatic use of ```0``` for WGS and read length for WES and TS (read length is calculated from first control sample)
- ```WAY_TO_BAM``` **(D)** **(!)** - ```absolute```, ```find_in_dir```
  - ```default``` = ```absolute```
  - ```absolute``` - BAM file for each control/sample is provides as absolute path. Final absolute path consists of ```CTRL_BAM_DIR```(projects master file) + ```CTRL_PATH_TO_BAM``` (controls master file) for controls and ```SMPL_BAM_DIR``` (projects master file) + ```SAMPLE_PATH_TO_BAM``` (samples master file) for samples.
  - ```find_in_dir``` - BAM file for each control/sample is found recursively in folder ```CTRL_BAM_DIR```/```SMPL_BAM_DIR``` using control/sample identifier (controls/samples master file) and shared pattern ```CTRL_BAM_PATTERN```/```SMPL_BAM_PATTERN``` (projects master file). Importantly, this pattern has to fit to only one BAM file in the provided folder. This setting is not primarily recommended, but it can be very useful when BAM files are in folders with many subfolders and absolute path would be very long and not very enjoyable to specify for each control/sample. Final file is found as ```CTRL_ID*CTRL_BAM_PATTERN```/```SAMPLE1_ID*SMPL_BAM_PATTERN```. Pattern can be for example *.bam*, *_final.bam* etc.
- ```CTRL``` **(!)** - ```yes```, ```no```
  - ```yes``` - controls will be included in the analysis
  - ```no``` - controls will not be included in the analysis; coverage profile will be denoise only by GC correction
- ```CTRL_BAM_DIR``` - path to folder with BAM file of controls
- ```CTRL_BAM_PATTERN``` **(D)** - pattern (file suffix) for the BAM files of controls (examples: *.bam*, *_final.bam* etc.); used only when ```WAY_TO_BAM``` is ```find_in_dir```
  - ```default``` = ```na``` (which means it is not used as ```WAY_TO_BAM``` ```default``` is ```absolute```)
- ```CTRL_SEX_TEST``` **(D)** **(!)** - ```custom```, ```no``` - see Sex test for more information.
  - ```default``` - sex test of controls will be performed using default setting (*SRY* and *GAPDH* genes coverage)
  - ```custom``` - sex test of controls will be performed using custom genomic locations defined within the sex master file under ```PROJECT_ID```
  - ```no``` - sex test of controls will not be performed (for example because chromosome Y has no targets, sex of controls is unknown or user just does not require this test to be done)
- ```CTRL_PON_SEX_SELECT``` **(D)** **(!)** - ```M```, ```F```, ```matched```, ```mixed``` - During pre-processing, controls are denoised by each other and segmented. This allows to look at controls qc and number of segments with given segmentation setting, and possible exclude controls with extreme values if wanted. Denoised data are also further used to determine CNVs in control in annotation and plots. Controls itself and controls with the same ```MAIN_ID``` are excluded from denoising. Setting ```matched``` is recommended, others can be used in case of limited number of controls and/or controls with each sex.
  - ```default``` = ```matched```
  - ```M``` - only male controls are used for denoising of controls
  - ```F``` - only female controls are used for denoising of controls
  - ```matched``` - male controls are used for denoising of male controls and female controls are used for denoising of female controls
  - ```mixed``` - male and female controls are mixed together for denoising of samples; NOT RECOMMENDED setting for two reasons: 1) gonozomes cannot be analyzed, and 2) mixing sex was found to increase false positivity rate for autosomes during pipeline testing (one single sex control provided better sensitivity/specificity than six mixed controls). This setting is recommended to be used only if the sex is unknown and cannot be determined from Y capture.
- ```SMPL_BAM_DIR``` - path to folder with BAM file of samples
- ```SMPL_BAM_PATTERN``` **(D)** - pattern (file suffix) for the BAM files of samples (examples: *.bam*, *_final.bam* etc.); used only when ```WAY_TO_BAM``` is ```find_in_dir```
  - ```default``` = ```na``` (which means it is not used as ```WAY_TO_BAM``` ```default``` is ```absolute```)
- ```SMPL_PON_SEX_SELECT``` **(D)** **(!)** - ```matched```, ```M```, ```F```, ```both_separated```, ```matched_main```, ```matched_each```, ```mixed```, ```all_options```. This value allows to specify what control are used for denoising of samples. Controls with the same ```MAIN_ID``` as sample are automatically excluded from denoising. Selection depends on user's requirements and number of controls with each sex available.
  - ```default``` = ```matched_each```
  - ```M``` - only male controls are used for denoising of samples
  - ```F``` - only female controls are used for denoising of samples
  - ```both_separated``` - ```M``` and ```F``` will be performed separately
  - ```matched_main``` - male or female controls are selected based on sex of the main sample (proband or tumor sample) and used for all related samples (family members or germline sample of the tumor)
  - ```matched_each``` - male controls are used for denoising of male samples and female controls are used for denoising of female samples
  - ```mixed``` - male and female controls are mixed together for denoising of samples; NOT RECOMMENDED setting for two reasons: 1) gonozomes cannot be analyzed, and 2) mixing sex was found to increase false positivity rate for autosomes during pipeline testing (one single sex control provided better sensitivity/specificity than six mixed controls). This setting is recommended to be used only if the sex is unknown and cannot be determined from Y capture. 
  - ```all_option``` - ```M```, ```F```, ```mixed```, ```matched_main``` and ```matched_each``` will be performed separately
- ```CN_FREQ_IN_CTRLS``` **(!)** - ```yes```, ```no```
  - ```yes``` - CNV analysis in controls will be included in the annotation and plots. 
  - ```no``` - CNV analysis in controls will not be included in the annotation and plots, for example because of no or very small number of controls.
- ```GNOMAD_SELECTION``` **(D)** **(!)** - ```capture_filter```, ```capture_filter_and_ctrl_denois```, ```af_filter```, ```af_filter_and_ctrl_denois``` - original gnomAD database of SNP that is provided has to be filtered for efficient performing; usually ```capture_filter```/```capture_filter_and_ctrl_denois``` for TS or WES and ```af_filter```/```af_filter_and_ctrl_denois``` for WGS. Option ```_and_ctrl_denois``` requires reasonable number of controls.
  - ```default``` - automatically selected ```af_filter_and_ctrl_denois``` for WGS and ```capture_filter_and_ctrl_denois``` for WES and TS
  - ```capture_filter``` - SNPs are filtered for those being in the covered regions (for TS and WES)
  - ```capture_filter_and_ctrl_denois``` - ```capture_filter``` + noisy SNPs (recognized in controls) are excluded; noise in controls is found by default parameters in the base master file
  - ```af_filter``` - SNPs with very rare alternative alleles are excluded; allelic frequency is defined by default parameter in the base master file
  - ```af_filter_and_ctrl_denois``` - ```af_filter``` + noisy SNPs (recognized in controls) are excluded; noise in controls is found by default parameters in the base master file
- ```SEGMENTATION_ID``` **(D)** **(!-)** (see ```Data segmentation``` section above for more information)
  - ```default``` = ```smart-0.65``` for single sex denoising and ```smart-0.85``` for mixed sex denoising
  - smart segmentation - ```smart-XX``` (without sub-clonal analysis) or ```smart-XX-sub-YY``` (with subclonal analysis); XX = number between 0.5 (high sensitivity) and 1.0 (high specificity), YY = number between 0.1 and 1.0 (for example, 0.5 indicates that subclonal analysis to detect monoallelic CNV in 50% of cells)
  - custom segmentation - short identifier (for example *my_segm*, *somatic_segmentation* etc.) to select required segmentation conditions from the segmentation master file
- ```SEGMENTATION_ID_USE``` **(D)** **(!)** - ```segm_id```, ```full``` - This controls what identifier of segmentation will be used in files names.
  - ```default``` = ```segm_id```
  - ```segm_id``` - RECOMMENDED; segmentation condition will be used in files names by its short identifier (full version is still kept inside of the segmentation table)
  - ```full``` - all segmetation parameters in numeric version will be used in files names; not primarily recommended due to creating long file names that can be difficult in some systems
- ```NOTE``` - any additional information that user wants to keep with the project

### 10.3. Controls master file
```/Masters/master_controls.txt```  
controls spreadsheet

**(!)** = requires explicit values

Columns:
- ```INCLUDE``` **(!)** - ```yes```, ```no``` - helps to regulate what controls will be included and at the same time to keep excluded controls with excluding reason in a note column
  - ```yes``` - control will be included
  - ```no``` - control will not be included
- ```PROJECT_ID``` - short and unique identifier for the project, matching ```PROJECT_ID``` in the other master files
- ```CAPTURE_ID``` - short identifier for the capture, matching ```CAPTURE_ID``` in the other master files
- ```GENOME_VERSION``` **(!)** - ```GRCh37-hg19```, ```GRCh38-hg38```, matching ```GENOME_VERSION``` in the other master files
  - ```GRCh37-hg19``` - for any version of human genome assembly such as GRCh37, hg19 and b37
  - ```GRCh38-hg38``` - for any version of human genome assembly such as GRCh38 and hg38
- ```MAIN_ID``` - controls/samples group identifier; This identifier is the one that decides what controls to exclude from for the denoising because of relatedness. Therefore, controls/samples that are related have to share the same identifier in the controls and the samples master file.
- ```CTRL_ID``` - control identifier; If ```WAY_TO_BAM``` is ```find_in_dir```, this identifier has to be part of BAM file name as ```CTRL_ID*CTRL_BAM_PATTERN``` (pattern determined in the projects master file)
- ```CTRL_SEX``` **(!)** - ```M```, ```F```, ```unk``` - control sex
  - ```M``` - control is male
  - ```F``` - control is female
  - ```unk``` - sex of control is unknown
- ```CTRL_PATH_TO_BAM``` - control BAM file path when ```WAY_TO_BAM``` is ```absolute```, final path for each BAM file is combination of ```CTRL_BAM_DIR``` (projects master file) and ```CTRL_PATH_TO_BAM```. If ```WAY_TO_BAM``` is ```find_in_dir```, this column is not used.
- ```NOTE``` - any additional information that user wants to keep with the control


### 10.4. Samples master file
```/Masters/master_samples.txt```  
samples spreadsheet

**(!)** = requires explicit values

Columns:
- ```INCLUDE``` **(!)** - ```yes```, ```no``` - helps to regulate what samples are analysed and skipped during run
- ```PROJECT_ID``` - short and unique identifier for the project, matching ```PROJECT_ID``` in the other master files
- ```CAPTURE_ID``` - short identifier for the capture, matching ```CAPTURE_ID``` in the other master files
- ```GENOME_VERSION``` **(!)** - ```GRCh37-hg19```, ```GRCh38-hg38```, matching ```GENOME_VERSION``` in the other master files
  - ```GRCh37-hg19``` - for any version of human genome assembly such as GRCh37, hg19 and b37
  - ```GRCh38-hg38``` - for any version of human genome assembly such as GRCh38 and hg38
- ```MAIN_ID``` - controls/samples group identifier; This identifier is the one that decides what controls to exclude from the denoising because of relatedness. Therefore, controls/samples that are related have to share the same identifier in the controls and the samples master files.
- ```SAMPLE1_ID``` - sample_1 (main sample) identifier; If ```WAY_TO_BAM``` is ```find_in_dir```, this identifier has to be part of BAM file name as ```SAMPLE1_ID*SMPL_BAM_PATTERN``` (pattern determined in the projects master file)
- ```SAMPLE1_TYPE``` - sample_1 (main sample) type; for example *proband*, *child-affected*, *child-unaffected*, *tumor* etc.
- ```SAMPLE1_SEX``` - **(!)** - ```M```, ```F```, ```unk``` - sample_1 (main sample) sex
  - ```M``` - sample is male
  - ```F``` - sample is female
  - ```unk``` - sex of sample is unknown
- ```SAMPLE1_PATH_TO_BAM``` - sample_1 (main sample) BAM file path when ```WAY_TO_BAM``` is ```absolute```, final path for each BAM file is combination of ```SMPL_BAM_DIR``` (projects master file) and ```SAMPLE1_PATH_TO_BAM```. If ```WAY_TO_BAM``` is ```find_in_dir```, this column is not used.
- ```SAMPLE2_ID``` - sample_2 identifier; If ```WAY_TO_BAM``` is ```find_in_dir```, this identifier has to be part of BAM file name as ```SAMPLE2_ID*SMPL_BAM_PATTERN``` (pattern determined in the projects master file)
- ```SAMPLE2_TYPE``` - sample_2 (main sample) type; for example *father*, *mother*, *germline* etc.
- ```SAMPLE2_SEX``` - **(!)** - ```M```, ```F```, ```unk``` - sample_2 sex
  - ```M``` - sample is male
  - ```F``` - sample is female
  - ```unk``` - sex of sample is unknown
- ```SAMPLE2_PATH_TO_BAM``` - sample_2 BAM file path when ```WAY_TO_BAM``` is ```absolute```, final path for each BAM file is combination of ```SMPL_BAM_DIR``` (projects master file) and ```SAMPLE2_PATH_TO_BAM```. If ```WAY_TO_BAM``` is ```find_in_dir```, this column is not used.
- ```SAMPLE3_ID``` - sample_3 identifier; If ```WAY_TO_BAM``` is ```find_in_dir```, this identifier has to be part of BAM file name as ```SAMPLE3_ID*SMPL_BAM_PATTERN``` (pattern determined in the projects master file)
- ```SAMPLE3_TYPE``` - sample_3 (main sample) type; for example *father*, *mother* etc.
- ```SAMPLE3_SEX``` - **(!)** - ```M```, ```F```, ```unk``` - sample_3 sex
  - ```M``` - sample is male
  - ```F``` - sample is female
  - ```unk``` - sex of sample is unknown
- ```SAMPLE3_PATH_TO_BAM``` - sample_3 BAM file path when ```WAY_TO_BAM``` is ```absolute```, final path for each BAM file is combination of ```SMPL_BAM_DIR``` (projects master file) and ```SAMPLE3_PATH_TO_BAM```. If ```WAY_TO_BAM``` is ```find_in_dir```, this column is not used.
- ```NOTE``` - any additional information that user wants to keep with the sample


### 10.5. Sample and regions of interest to plot master file
```/Masters/master_reg_of_interest_to_plot.txt```  
spreadsheet that allows to plot required region for required sample; if ```PROJECT_ID``` is recognized to be running and at least one of the lines of this master file has ```INCLUDE``` = ```yes```, separated detail plot for the required sample within ```/detail_regions_of_interest/``` folder will be produced. It can be useful when user wants to print a suspicious region that was not recognized abnormal and so not printed as a detail plot.

**(!)** = requires explicit values

Columns:
- first 17 columns are same as columns in the samples master file
- ```CONTIG``` - chromosome of the region to plot
- ```START``` - starting position of the region to plot
- ```END``` - ending position of the region to plot
- ```ZOOM``` - numeric value, zoom for the defined region to be printed to the plot, examples:
  - ```0```: ```CONTIG```: ```START``` - ```END``` will be printed to plot
  - ```1```: ```CONTIG```: [```START```-1✗(```END```-```START```)] - [```END```+1✗(```END```-```START```)] will be printed to plot
  - ```2```: ```CONTIG```: [```START```-2✗(```END```-```START```)] - [```END```+2✗(```END```-```START```)] will be printed to plot

### 10.6. Segmentation master file
```/Masters/setting_segmentation.txt```  
parameters for CNVs and LOH segmentation, paired with project by ```PROJECT_ID```

Columns:
- ```FEATURE``` - identifier of segmentation parameter
- ```SEGM-BASE``` - default segmentation setting
- ```...``` - additional columns can be added by user; short and unique ID of segmentation is required

Segmentation identifier defined in the column header (for example *my_segm*, *somatic_segmentation* etc.) is used to pair segmentation setting and project as ```SEGMENTATION_ID``` within the projects master file. This identifier is also used in files names if ```SEGMENTATION_ID_USE``` = ```segm_id```)

List of parameters:  
```SEGM_DIFFERENCE```  
```SEGM_MINSIZE```  
```SEGM_MINSIZE_SUB```  
```SEGM_MINPROBE```  
```SEGM_MINPROBE_SUB```  
```SEGM_MINKEEP```  
```SEGM_MINLOSS```  
```SEGM_MINLOSS_SUB```  
```SEGM_MINLOSS_BIAL```  
```SEGM_MINGAIN```  
```SEGM_MINGAIN_SUB```  
```SEGM_GAP```  
```SEGM_SMOOTHPERC```  
```SEGM_AFDIF```  
```SEGM_AFSIZE```  
```SEGM_AFPROBE```

Subclonal analysis is defined by parameters ending by ```_SUB```. If no subclonal analysis is required, these variables are ```none```.  
If ```SEGM_GAP``` is ```none```, segments continue independently on how far they are from each other in the capture.


### 10.7. Sex master file
```/Masters/sex_test_regions.txt```  
genomic locations for sex text

**(!)** = requires explicit values

Columns:
- ```CONTIG```
- ```START```
- ```END```
- ```GENE_TYPE``` **(!)** - ```TEST_GENE```, ```CTRL_GENE```
  - ```TEST_GENE``` - genomic position of tested gene 
  - ```CTRL_GENE``` - genomic position of control gene
- ```PROJECT_ID``` - short identifier for the project (has to match with all ```PROJECT_ID``` in the other master files)
- ```GENE``` - gene symbol

## 11. Outputs

### 11.1. Data outputs
- ```...allelicCounts...``` - rds file with SNP zygosity collection
- ```...counts...``` - hdf5 and tsv files with coverage collection
- ```...denoisedCR...``` - tsv file with denoised coverage data
- ```...standardizedCR...``` - tsv file with standardized coverage data
- ```...SEGMENTS...segmentation...``` - tsv file with segmentation data before annotation, contains all normal and abnormal segments
- ```...SEGMENTS...``` - tsv file with segmentation and annotation, contains all normal and abnormal segments
- ```...SEGMENTS...REPORT...``` - tsv file segmentation, annotation and variant origin prediction, contains ONLY abnormal segments


### 11.2. Plot outputs
- **genome** - figure showing all chromosomes together (chr1-chr22, chrX and chrY)
- **chromosome** - figure is generated for each chromosome
- **detail** - figure is generated for each abnormal segment (does not include sub-clonal findings)

### 11.3. IGV outputs
- ```mainID_sampletype_ID_sex_abnormal_segments.bed``` - BED file with abnormal segments

- ```mainID_sampletype_ID_sex_CN.bw``` - BigWig with denoised coverage data
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min -2.5, mid 0.0, max 2.5
  - Note: log2 ratio values smaller than -2.25 are shifted to -2.25 and log2 ration bigger than 2.25 are shifted to 2.25 in IGV output (as well as in the PNG plots)

- ```mainID_sampletype_ID_sex_SNP.bw``` - BigWig with SNP zygosity
  - IGV setting:
    - Type of Graph: Points
    - Windowing Function: None
    - Set Data Range: min 0.0, mid 0.5, max 1.0

## 12. Custom gene annotation
How to prepare gene annotation other then default RefSeq provided:
- download required gene annotation from [UCSC website](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=968304787_9UHI0R5esZuD1z4azUHoV9dEdsoU)
- columns have to be named as follows (order does not need to be the same):
  - ```gene_id```: identifier that will be displayed at plots and used for annotation (usually gene symbol)
  - ```chrom```: chromosome
  - ```strand```: + or - for strand
  - ```txStart```: transcription start position
  - ```txEnd```: transcription end position
  - ```cdsStart```: coding region start
  - ```cdsEnd```: coding region end
  - ```exonStarts```: exon start positions
  - ```exonEnds```: exon end positions
 Most of the columns will match automatically, only ```gene_id``` needs to picked and rename by user, for example by command ```sed -i 's/name2/gene_id/g' file.txt```. Table should never be edited in EXCEL due to renaming some gene symbols to dates.

## 13. Project-specific plot
This is an advance option that allows to produce personalized detail plots for additional regions with personalized annotation.
Customized R scripts is required, as provided example. This script has to be places to ```/Scripts/R/``` folder and named as ```plot_CN_PROJECT_ID.R```. It is then recognized as source in ```plot_CN.R``` and produces plots within ```/project_specific_regions/``` folder.

## 14. Databases source
- processing of the databases from original source is available under request
- original source of databases:
  - [gnomAD SNP v2.1](https://gnomad.broadinstitute.org/downloads) (exome and genome data was downloaded, combined and processed into a simple VCF file in hg19/GRCh37 version; hg38/GRCh38 version prepared by liftover)
  - [Database of Genomic Variants (DGV, TCAG) v2020-02-25](http://dgv.tcag.ca/dgv/app/home)
  - [UCSC chromosome bands and centromeres](https://hgdownload.soe.ucsc.edu/downloads.html)
  - [RefSeq](https://genome.ucsc.edu/cgi-bin/hgTables)


## 15. Limitations
CNVRobot pipeline is not suitable for detection of common population CNVs and CNVs in low mappability regions. It should be noted that it is a depth of coverage analysis within certains bins and so starts and ends of abnormalities are not precisily mapped. CNVRobot can detect only unbalanced changes. Similarly to DNA microarrays, coverage-based analysis can be problematic if the ploidy of the sample is changed, as it is based on median centering normalization.

## 16. Contact
Aneta.Mikulasova@newcastle.ac.uk
