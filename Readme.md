#  SVINGS: Structural Variants Identification in Next Generation Sequencing Data

Version -1 for CNV detection Pipeline
# RD based CNV detection tool

## Requirements
Make sure that these tools are added to your PATH variable 
* samtools ([link](https://sourceforge.net/projects/samtools/))
* bedtools (>=2.26) ([link](http://bedtools.readthedocs.org/en/latest/content/installation.html))
* openMPI  ([Installation Guide](http://lsi.ugr.es/~jmantas/pdp/ayuda/datos/instalaciones/Install_OpenMPI_en.pdf))
* R-packages - DNACopy, ParDNACopy, quantsmooth, GenomicAlignments, rtracklayers

## Installation
Download the source code from https://github.com/vogetihrsh/SVINGS/, extract the zip file

```
unzip SVINGS.zip
cd SVINGS
make 
```
Please download annotations from http://bioinf.iiit.ac.in/svings/ and add the folder to working directory.

## Usage

```
# CALCULATE OPTIMAL BIN/WINDOW SIZE
./calOptBinSize -c <config file> -i <input BAM>

# DATA PREPARATION
./prepareData -m <mappability file> -g <gc content file> --win <desired window size> --genome_file <Genome file> -o <Output file name prefix>

# PRETREATMENT STEP
./pretreatment -i <input BAM> -z <bed file containing bins> -m < file contaning mappability values> -o <output prefix> 

# SEGMENTATION
./runSegmentation -o <output prefix> <-t or -d or both>

# POSTPROCESSING AND CNV CALLING
./callCNV -o <output prefix> -z <bed file contaning bins> <genome flag (--hg18 or --hg19)>

# ANNOTATION 
./annotate -i <input bed file> -o <output prefix> <genome flag (--hg18 or --hg19)>

# PLOT
./plot -i <input bed file> -o <output file name> <genome flag (--hg18 or --hg19)>
./plot -i <input bed file> -d --start <start position> --end <end position> -o <output file name> <genome flag (--hg18 or --hg19)>



```
#### Note:
Pretreatment, runSegmentation, callCNV should have same value for the parameter -o.

### calOptBinSize
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-c |Yes/STRING |Config file contaning parameters to be used for calculations. Parameters similarity to those used in ReadDepth config file. For additional info, please refer to "config.txt" in the source directory. |Required |
|-i  | Yes/STRING | Input BAM File. | Required |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### prepareData
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-m |Yes/STRING | Mappability score file available for 100 bp | Required|
|-g |Yes/STRING | GC content score file availabled for 100 bp | Required|
|--win |Yes/INTEGER | Desired window size for CNV predictions | Required|
|--genome_file|Yes/STRING | A tab-separated genome file containing chromosome number and size of the chromosome | Required|
|-o | Yes/STRING | Output file name that is used as prefix for coordinate, GC content and mappability files generated| Required|

### pretreatment
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-i  | Yes/STRING | Input BAM File. | Required |
|-o |Yes/STRING | Output Prefix. | Required |
|-z  | Yes/STRING | Bed file contanning bins. | Required |
|--mapfile | Yes/STRING | File contaning mappability values ([0,1]) of bins present in the file passed as value to -z. | Required |
|-p | Yes/INTEGER | Default Value: 32. Number of processes to be created for parallelization. | Optional |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional |
|-y| No/- | When used Yoon et al. method for GC correction and mappability correction is applied to the data| Optional|
|--yoonGC| No/-|Use Yoon et al method for GC correction | Optional|
|--yoonMap| No/-|Use Yoon et al method for Mappability correction | Optional|
|--loessGC| No/-| Use Loess method to correct GC-Bias | Optional |
|--gcfile | No/- | File containing GC content values ([0,1]) generated in prepareData, required when using LoessGC/medianGC options| Optional|

### runSegmentation
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-o |Yes/STRING | Output Prefix. | Required |
|-p | Yes/INTEGER | Default Value: 32. Number of processes to be created for parallelization. | Optional |
|-t |No/- | When this flag is specified Total Variation Minimization algorithm is used for segmentation.  | Atleast one of -t and -d should be used. |
|-d |No/- | When this flag is specified Circular Binary Segmentation algorithm is used for segmentation.  | Atleast one of -t and -d should be used. |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional |

### callCNV
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-o |Yes/STRING | Output Prefix. | Required |
|-z  | Yes/STRING | Bed file contanning bins. | Required |
|--hg18 |No/- | For hg18 genome. This flag will be used for getting the gap info (hg18_gap.bed) from the annotations directory. | Not required when --hg19 or --genome is used. |
|--hg19 |No/- | For hg19 genome. This flag will be used for getting the gap info (hg19_gap.bed) from the annotations directory. | Not required when --hg18 or --genome is used. |
|--genome |No/- | For genomes other than hg18 and hg19. This flag will be used for getting the gap info (<genome>_gap.bed) from the annotations directory. | Not required when --hg19 or --hg18 is used. |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### annotate 
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
| -i  | Yes/STRING  |Input bed file contaning regions that needed to be annotated. | Required  |
| -o  | Yes/STRING  |Output prefix for the generated bed file.  | Required |
| --hg18  | No/-  | hg18 genome flag. When used checks for the files contaning gene info (hg18.bed), DGV annotations (hg18DGV.bed), location annotation (hg18LocInfo.bed and gap annotations (hg18_gap.bed) in the annotation directory |Not required when --hg19 or --genome is used. |
| --hg19  | No/-  | hg19 genome flag. When used checks for the files contaning gene info (hg19.bed), DGV annotations (hg19DGV.bed), location annotation (hg19LocInfo.bed and gap annotations (hg19_gap.bed) in the annotation directory |Not required when --hg18 or --genome is used |
| --genome  | Yes/STRING  | Any genome other than hg18 and hg19 can be specified through this flag. It requires <genome>.bed, <genome>DGV.bed, <genome>_gap.bed and <genome>LocInfo.bed  files to be present in annotation directory. | Not required when --hg18 or --hg19 is used. |
| --annDir | Yes/STRING  | Default: "SOURCEDIR/annotations". Path to the annotation directory. | Optional  |
| --help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### plot
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
| -i | Yes/STRING | Input bed file containing CNV cooridnates for visualization | Required |
| -o | Yes/STRING | Output file name prefix for the image file generated | Required |
| --hg18 |No/- | For hg18 genome. | Not requires when --hg19 is used |
| --hg19 |No/- | For hg19 genome. | Not required when --hg18 is used |
| -d | Yes/- | User may use this parameter for visualiing CNVs within specific user-defined coordinates | Optional, required when start and end coordinates are used |
|--start | Yes/INTEGER | Start coodinate | Optional |
|--end | Yes/INTEGER | Stop coordinate | Optional |

## Contact
  For any queries drop a mail @ vogetisri.harsha@research.iiit.ac.in
