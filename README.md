# DeepSig
Deep Learning-based Mutational Signature Caller for MSK-IMPACT Data

## Overview
**DeepSig** contains single-base substitution (SBS) mutational signature inference tools intended for 
data sets derived from targeted sequencing (MSK-IMPACT) platforms. 
Models for major cancer types, pre-trained using synthetic data derived from whole-genome seqeuencing (WGS)
mutational catalogs and matched to corresponding MSK-IMPACT data, are used to predict presence of all
reference signatures in a new clinical sample based on estimated precision.

For a single sample or a small cohort derived from MSK-IMPACT, the main workflow comprises extracting signature
scores using a cancer type-specific model, making discrete ternary calls, and estimating exposures 
(number of mutations attributed to each signature).

## Input data
Input data are of the form of catalog matrix:

Sample_ID              | A[C>A]A | A[C>A]C  | A[C>A]G | A[C>A]T | C[C>A]A  
---------------------- | ------- | -------- | ------- | ------- | --------
Tumor_Sample_Barcode_1 |     0   |     3    |    5    |    0    |    0  
Tumor_Sample_Barcode_2 |     2   |     1    |    1    |    2    |    0 
Tumor_Sample_Barcode_3 |     5   |     0    |    0    |    1    |    1 

Column names are the the set of categories to which mutation data from sequencing experiments have been classified 
(`Sample_ID` column name is absent in the actual file). These trinucleotide contexts are the pyrimidine bases 
before and after mutation flanked by upstream and downstream nucleotides (96 in total). Each row corresponding 
to `Tumor_Sample_Barcode` contains non-negative counts of single nucleotide variants (SNVs). 
This trinucleotide matrix can be generated using the utility function [maf2cat3](https://github.com/mskcc/DeepSig/blob/main/R/utilities3.R). 
MAF files can be generated from VCF files using [vcf2maf](https://github.com/mskcc/vcf2maf). 
Note: the output from [maf2cat3] needs to be transposed so that rows contain samples.

## Main caller
The function

    DL.call(catalog, cancer.type = 'pancancer', model.path = './.DeepSig', min.attr = 1, ...)

will find the pre-trained model corresponding to `cancer.type` and perform signature fitting and filtering.
The `cancer.type` argument is used to choose among pre-trained models:

Cancer type                    | cancer.type  | OncoTree              | ICGC/PCAWG | TCGA
-------------------------------| ------------ | --------------------- | ---------- | ------
Breat cancer                   | breast       | Breast                | Breast     | BRCA
Ovarian cancer                 | ovary        | Ovary/Fallopian Tube  | Ovary      | OV
Prostate cancer                | prostate     | Prostate              | Prostate   | PRAD
Pancreatic cancer              | pancreas     | Pancreas              | Pancreas  | PAAD
Bladder cancer                 | bladder      | Bladder/Urinary Tract | Bladder    | BLCA
Colorectal cancer              | colorect     | Bowel                 | Colorectal | COAD
Melanoma                       | skin         | Skin                  | Skin       | SKCM
Glioma                         | cns          | CNS/Brain             | CNS        | GBM
Non-small cell lung cancer     | nsclc        | Lung                  | Lung       | LUAD/LUSC
Small cell lung cancer         | sclc         | CSCLC                 |            |
Head and neck cancer           | head_neck    | Head and Neck         | Head       | HNSC
Renal cell carcinoma           | kidney       | Kidney                | Kidney     | KICH/KIRC/KIRP
Endometrial cancer             | uterus       | Uterus                | Uterus     | UCEC
Esophagogastric cancer         | esophagus    | Esophagus/Stomach     | Esophagus  | ESCA/STAD
Germ cell tumor                | gct          |                       |            | TGCT
Pan-cancer model               | pancancer    |                       |            |


Input argument `cancer.type` will be checked against the second column above. If it does not match one, it will be 
thought of as an [oncoTree](https://oncotree.mskcc.org/) code, and an attempt will be made to match it to `cancer.type`
using [oncoTree()](https://github.com/mskcc/DeepSig/blob/main/R/oncotree.R). This matching of oncotree code to `cancer.type`
is mostly based on the top-level `tissue` name of oncotree hierarchy, except for germ cell tumor and `sclc`.
If `cancer.type` does not match an oncotree code, it will be thought of as a custom model (see below).

In general, the choice of which model to apply for a cohort should be based on the knowledge of how similar the cohort
is to one of the tissue of origin-based cancer types above. For rare cancer types, cancer of unknown primary, and samples
for which signatures not in the reference might be expected (e.g., lung cancer samples with history of treatment with 
temozolomide or hypermutation with POLE oncogenic mutations), `pancancer` model can be used. The drawback of `pancacer` model is
that detection precision and recall are generally lower than cancer type-specific models.

To reduce the size of package, default models are not included in the package. The `model.path` argument (by default `./.DeepSig`)
is the path where trained models can be found: directory containing SBS* subdirectories, similar
to those in [inst/extdata/dlsig/v0.95/breast](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast),
and other associated files (`refsig.txt` and `threshold_cut.txt`) will be looked for in
`model.path/[cancer.type]` subdirectory. 
The table of reference signatures (`refsig.txt`) in a model can be looked up from pre-trained model directories (e.g., for breast cancer,
see [inst/extdata/dlsig/v0.95/breast](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast)).
If this main subdirectory does not exist, an attempt will be made to 
download these files to the directory via queries to [GitHub REST API](https://docs.github.com/en/rest?api=&apiVersion=2022-11-28) with 
[modelFetch()](https://github.com/mskcc/DeepSig/blob/main/R/modelFetch.R).
This scheme also ensures that a second call for the same `cancer.type` will use the previously downloaded data. Runs without
internet access will require pre-downloaded models whose path is provided as `model.path`.

As an example, the following script will generate outputs for the TCGA-OV samples using `ovary` pre-trained model:

    > library(DeepSig)
    > xcat <- read.table(system.file('extdata', 'tcga_ov_xcat.txt', package = 'DeepSig'), header = TRUE, sep = '\t')
    > z <- DL.call(catalog = t(xcat), cancer.type = 'ovary')
    > names(z)
      [1] "score"        "ternary.call" "exposure"     "ref.sig"      "threshold"


### Output
The return value is a list of 5 data frames as shown in the above example. The component **score** contains
continuous scores for the presence of each sample-signature combination: 

sid            | M    | SBS1    | SBS2.13  | SBS3     | SBS5    | SBS6  
-------------- | ----:| ------: | -------: | -------: | ------: | -----:
TCGA.61.2092   | 141  | 5.24    | 2.72     | 9.03     | 2.67    | -4195  
TCGA.36.2540   | 10   | 2.93    | -1.39    | 1.57     | 1.57    | -73      
TCGA.24.1467   | 9    | -0.09   | -1.44    | -0.50    | -0.50   | -849 

where the columns show the repertore of reference signatures present in the **ovary** model being used.
The columnn labeled **M** shows the total number of mutations for each sample.

The threshold values accessbile with **threshold** are input data containing two threshold scores for minimum
estimated precision and minimum negative predictive value (NPV), both of 0.9 for most cases, shown in column
**min_value**:

S        | measure    | min_value | threshold | precision | recall | npv    | engine
-------- | ---------- | --------- | --------: | --------: | -----: | -----: | -------
SBS1     | precision  | 0.9       | 1.90      | 0.90      | 0.09   | 0.73   | mlp
SBS1     | npv        | 0.9       | -2.55     | 0.29      | 0.99   | 0.90   | mlp
SBS2.13  | precision  | 0.9       | 1.22      | 0.90      | 0.19   | 0.77   | mlp
SBS2.13  | npv        | 0.9       | -1.43     | 0.38      | 0.85   | 0.90   | mlp
SBS3     | precision  | 0.9       | 0.92      | 0.90      | 0.85   | 0.62   | mlp
SBS3     | npv        | 0.9       | -1.24     | 0.77      | 1.00   | 0.90   | mlp

The columns **precision**, **recall**, and **npv** show the values attained at the thresholds.
The column **engine** is informational only: either `mlp` or `cnn`, the flavor of DL engines that was used
for optimal performance in training. 

These threshold values are used to clip signature scores into discrete calls as shown in `ternary.call`:

sid            | M    | SBS1    | SBS2.13  | SBS3
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | P       | P        | P
TCGA.36.2540   | 10   | P       | I        | P
TCGA.24.1467   | 9    | I       | N        | I

There are three types of calls:
- `P` (positive): signature is present with minimum precision of 0.9
- `N` (negative): signature is not present with minimum NPV of 0.9
- `I` (indeterminate): signature is neither present with minimum precision of 0.9 nor absent with mininum NPV of 0.9

The component `exposure` is a data frame of exposure proportions for all signatures with either `P` or `I` calls:

sid            | M    | SBS1    | SBS2.13  | SBS3
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | 0.18    | 0.05     | 0.18
TCGA.36.2540   | 10   | 0.82    | 0        | 0.17
TCGA.24.1467   | 9    | 0.52    | 0        | 0

These proportions are obtained by fitting the subset of signatures with `P` or `I` calls to each catalog using
[extractSig](https://github.com/mskcc/DeepSig/blob/main/R/extract.R), which uses
maximum likelihood estimation.
Additionally, any proportions that lead
to attribution (proportion times `M`) less than `min.attr` parameter to `DL.call` (default 1) are set to zero.

The `ref.sig` is informational only since it has been used in training the model and refitting:

Mutation types | SBS1  | SBS3  | SBS2.13 | SBS8
-------------- | ----- | ----- | ------- | ----
A[C>A]A        | 0.027 | 0.014 | 0.004   | 0.02
A[C>A]C        | 0.013 | 0.010 | 0.002   | 0.01
A[C>A]G        | 0.004 | 0.002 | 0.001   | 0.003

The reference signatures are primarily based on existing databases of WGS-derived signatures including
[signal](https://signal.mutationalsignatures.com/).


## Standalone script
If you are not interested in interactive usages with more flexibility and functionality, or want to use **DeepSig** 
as a part of a pipeline, use the command-line script [deepsig.R](https://github.com/mskcc/deepsig/blob/master/exec/deepsig.R). 
If you installed **DeepSig** as an R package using `install_github`, find the path via
```Rscript
system.file('exec', 'deepsig.R', package = 'deepsig')
```

If you cloned the repository, the file is located at the `./exec` subdirectory of the github main directory. We denote this package directory path as `PKG_PATH`. The command syntax is
```shell
$ $PKG_PATH/exec/DeepSig.R -h
usage: ./deepsig.R [-h] -i CATALOG -o OUTPUT -c CANCER_TYPE [-q]

Extract mutational signatures using DeepSig algorithm

options:
  -h, --help            show this help message and exit
  -i CATALOG, --input CATALOG
                        input catalog data file
  -o OUTPUT, --output OUTPUT
                        output directory
  -c CANCER_TYPE, --cancer-type CANCER_TYPE
                        use a model for the specified cancer type
  -q, --quiet           Run quietly
```
The arguments `CATALOG`, `CANCER_TYPER`, and `OUTPUT`specify the paths of input catalog data, cancer type, and output directory.


## Installation

Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
On a Mac,

    $ brew install gsl

Also set 
```
PKG_CPPFLAGS = -I/opt/homebrew/Cellar/gsl/2.7.1/include
PKG_LIBS = -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
```
in `~/.R/MAKEVARS` on a Mac (check the links to make sure the correct version number is used).

DL-based filtering requires [tensorflow](https://tensorflow.org/install) and [pandas](https://pandas.pydata.org/getting_started.html). In a Mac,
  
    $ python3 -m pip install pandas h5py tensorflow==2.15.0


See the tensorflow install page for further information.

We currently recommend downloading a source tar.gz file and
installing **DeepSig** via

    > install.packages('DeepSig_v0.9.8.tar.gz', repos = NULL, dependencies = TRUE)

which will also install dependencies.

### Catalog matrix generation

The catalog file as input data can be generated from MAF file using [maf2cat3()](https://github.com/mskcc/DeepSig/blob/master/man/maf2cat3.Rd). It requires a reference genome package, either [BSgenome.Hsapiens.UCSC.hg19](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19) 
or [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg38) installed:

    > library(BSgenome.Hsapiens.UCSC.hg19)
    > maf <- system.file('extdata', 'brca.maf', package = 'DeepSig')
    > x <- maf2cat3(maf = maf, ref.genome = BSgenome.Hsapiens.UCSC.hg19)
    > write.table(x, file = 'brca_catalog.txt', row.names = TRUE, col.names = TRUE, sep = '\t', quote = F)
    
If you do not want to use R-interface, a command-line script is available, assuming that Bsgenome.Hsapiens.UCSC.hg19 package has been installed:

    $ ./maf2cat3.R -h
    usage: ./maf2cat3.R [-h] MAF CATALOG

    Construct mutational catalog from MAF file with Ref_Tri column

    positional arguments:
      MAF         input MAF file
      CATALOG     output catalog file

    optional arguments:
      -h, --help  show this help message and exit

### De novo inference

See the [vignettes](https://html-preview.github.io/?url=https://github.com/mskcc/DeepSig/blob/main/old/denovo.html)
