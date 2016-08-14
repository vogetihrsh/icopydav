# pCNVD-1.0-dev

Version -1 for CNV detection Pipeline
# RD based CNV detection tool

## Requirements
Make sure that these tools are added to your PATH variable 
* samtools ([link](https://sourceforge.net/projects/samtools/))
* bedtools (>=2.26) ([link](http://bedtools.readthedocs.org/en/latest/content/installation.html))
* openMPI  ([Installation Guide](http://lsi.ugr.es/~jmantas/pdp/ayuda/datos/instalaciones/Install_OpenMPI_en.pdf))
* R-packages - DNACopy

## Installation
Download the source code from https://github.com/vogetihrsh/pCNVD-1.0-dev, extract the zip file

```
unzip pCNVD.zip
cd pCNVD
make 
```


## Usage
```
# PREPROCESSEING STEP
./preprocess -i <input BAM> -z <bed file containing bins> -m < file contaning mappability values> -o <output prefix> 

# SEGMENTATION
./runSegmentation -o <output prefix> <-t or -d or both>

# POSTPROCESSING AND CNV CALLING
./callCNV -o <output prefix> -z <bed file contaning bins> <genome flag (--hg18 or --hg19)>

# ANNOTATION 
./annotate -i <input bed file> -o <output prefix> <genome flag (--hg19 or --h19)>

# CALCULATE OPTIMAL BIN/WINDOW SIZE
./calOptBinSize -c <config file> -i <input BAM>
```
#### Note:
Preprocess, runSegmentation, callCNV should have same value for the parameter -o.

### annotate 
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
| -i  | Yes/STRING  |Input bed file contaning regions that needed to be annotated. | Required  |
| -o  | Yes/STRING  |Output prefix for the generated bed file.  | Required |
| - -hg18  | No/-  | hg18 genome flag. When used checks for the files contaning gene info (hg18.bed), DGV annotations (hg18DGV.bed), location annotation (hg18LocInfo.bed and gap annotations (hg18_gap.bed) in the annotation directory |Not required when --hg19 or --genome is used. |
| - -hg19  | No/-  | hg19 genome flag. When used checks for the files contaning gene info (hg19.bed), DGV annotations (hg19DGV.bed), location annotation (hg19LocInfo.bed and gap annotations (hg19_gap.bed) in the annotation directory |Not required when --hg18 or --genome is used |
| - -genome  | Yes/STRING  | Any genome other than hg18 and hg19 can be specified through this flag. It requires <genome>.bed, <genome>DGV.bed, <genome>_gap.bed and <genome>LocInfo.bed  files to be present in annotation directory. | Not required when --hg18 or --hg19 is used. |
| --annDir | Yes/STRING  | Default: "SOURCEDIR/annotations". Path to the annotation directory. | Optional  |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### callCNV
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-o |Yes/STRING | Output Prefix. | Required |
|-z  | Yes/STRING | Bed file contanning bins. | Required |
|- -hg18 |No/- | For hg18 genome. This flag will be used for getting the gap info (hg18_gap.bed) from the annotations directory. | Not required when --hg19 or --genome is used. |
|- -hg19 |No/- | For hg19 genome. This flag will be used for getting the gap info (hg19_gap.bed) from the annotations directory. | Not required when --hg18 or --genome is used. |
|- -genome |No/- | For genomes other than hg18 and hg19. This flag will be used for getting the gap info (<genome>_gap.bed) from the annotations directory. | Not required when --hg19 or --hg18 is used. |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### calOptBinSize
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-c |Yes/STRING |Config file contaning parameters to be used for calculations. Parameters similarity to those used in ReadDepth config file. For additional info, please refer to "config.txt" in the source directory. |Required |
|-i  | Yes/STRING | Input BAM File. | Required |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional | 

### preprocess
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-i  | Yes/STRING | Input BAM File. | Required |
|-o |Yes/STRING | Output Prefix. | Required |
|-z  | Yes/STRING | Bed file contanning bins. | Required |
|-m | Yes/STRING | File contaning mappability values ([0,1]) of bins present in the file passed as value to -z. | Required |
|-p | Yes/INTEGER | Default Value: 32. Number of processes to be created for parallelization. | Optional |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional |

### runSegmentation
| FLAG  |REQUIRES VALUE/TYPE  |DESCRIPTION  | OPTIONAL/REQUIRED |
| ---  | ---  | ---  | ---  |
|-o |Yes/STRING | Output Prefix. | Required |
|-p | Yes/INTEGER | Default Value: 32. Number of processes to be created for parallelization. | Optional |
|-t |No/- | When this flag is specified Total Variation Minimization algorithm is used for segmentation.  | Atleast one of -t and -d should be used. |
|-d |No/- | When this flag is specified Circular Binary Segmentation algorithm is used for segmentation.  | Atleast one of -t and -d should be used. |
| - -help | No/-  | Prints usage statement. Should be exclusively used. | Optional |





## Contact
  For any queries drop a mail @ vogetisri.harsha@research.iiit.ac.indd
