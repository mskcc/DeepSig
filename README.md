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

    DL.call(catalog, cancer.type = 'pancancer', model.path = NA, ref.sig = NA, threshold = NA, min.attr = 1, ...)

will find the pre-trained model corresponding to `cancer.type` and perform signature fitting and filtering.
Except for `catalog` and `cancer.type`, other arguments will default to those for pre-trained models if 
`cancer.type` is among those built-in:

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
Small cell lung cancer sclc    | sclc         | CSCLC                 |            |
Head and neck cancer head_neck | head_neck    | Head and Neck         | Head       | HNSC
Renal cell carcinoma kidney    | kidney       | Kidney                | Kidney     | KICH/KIRC/KIRP
Endometrial cancer uterus      | uterus       | Uterus                | Uterus     | UCEC
Germ cell tumor                | gct          |                       | TGCT       |
Pan-cancer model               | pancancer    |                       |            |


Input argument `cancer.type` will be checked against the second column above. If it does not match one, it will be 
thought of as an [oncoTree](https://oncotree.mskcc.org/) code, and an attempt will be made to match it to `cancer.type`
using [oncoTree()](https://github.com/mskcc/DeepSig/blob/main/R/oncotree.R). This matching of oncotree code to `cancer.type`
is mostly based on the top-level `tissue` name of oncotree hierarchy, except for germ cell tumor and `sclc`.
The table of reference signatures in a model can be looked up from pre-trained model directories (e.g., for breast cancer,
see [inst/extdata/dlsig/v0.95/breast](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast)).

The `model.path` argument is the path where trained models can be found (directory containing SBS* subdirectories similar
to those in [inst/extdata/dlsig/v0.95/breast](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast)), 
`ref.sig` is the path of the reference signature file 
(e.g., [refsig.txt](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast/refsig.txt)),
and `threshold` is the threshold file 
(e.g., [threshold_cut.txt](https://github.com/mskcc/DeepSig/tree/main/inst/extdata/dlsig/v0.95/breast/threshold_cut.txt)).

In general, the choice of which model to apply for a cohort should be based on the knowledge of how similar the cohort
is to one of the tissue of origin-based cancer types above. For rare cancer types, cancer of unknown primary, and samples
for which signatures not in the reference might be expected (e.g., lung cancer samples with history of treatment with 
temozolomide or hypermutation with POLE oncogenic mutations), `pancancer` model can be used. The drawback of `pancacer` model is
that detection precision and recall are generally lower than cancer type-specific models.
To reduce the size of package, default models are not included and are downloaded on demand via queries to 
[GitHub REST API](https://docs.github.com/en/rest?api=&apiVersion=2022-11-28) with 
[modelFetch()](https://github.com/mskcc/DeepSig/blob/main/R/modelFetch.R).

As an example, the following script will generate outputs for the TCGA-OV samples using `ovary` pre-trained model:

    > library(DeepSig)
    > token <- Sys.getenv('token')
    > xcat <- read.table(system.file('extdata', 'tcga_ov_xcat.txt', package = 'DeepSig'), header = TRUE, sep = '\t')
    > z <- DL.call(catalog = t(xcat), cancer.type = 'ovarian', token = token)
    > names(z)
      [1] "score"        "ternary.call" "exposure"     "ref.sig"      "threshold"

The argument `token` is the 
[personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)
in cases where the version of models are restricted-access. The above-use of `Sys.getev` assumes that the line

    token = '[your_access_token]'

has been placed in the file **~/.Renviron**.

### Output
The return value is a list of 6 data frames as shown in the above example. The component **score** contains
continuous scores for the presence of each sample-signature combination: 

sid            | M    | SBS1    | SBS2.13  | SBS3     | SBS5    | SBS6  
-------------- | ----:| ------- | -------- | -------- | ------- | -----
TCGA.61.2092   | 141  | 5.24    | 2.72     | 9.03     | 2.67    | -4195  
TCGA.36.2540   | 10   | 2.93    | -1.39    | 1.57     | 1.57    | -73      
TCGA.24.1467   | 9    | -0.09   | -1.44    | -0.50    | -0.50   | -849 

where the columns show the repertore of reference signatures present in the **ovary** model being used.
The columnn labeled **M** shows the total number of mutations for each sample.

The threshold values accessbile with **threshold** are input data containing threshold scores with minimum
estimated precision of 0.5 and 0.9:

S        | engine  | Threshold | Recall | Precision | Precision.cutoff
-------- | ------- | --------- | ------ | --------- | -----------------
SBS1     | mlp     | -0.83     | 0.62   | 0.5       | 0.5
SBS1     | mlp     | 1.90      | 0.09   | 0.9       | 0.9
SBS2.13  | mlp     | -0.93     | 0.65   | 0.5       | 0.5
SBS2.13  | mlp     | 1.22      | 0.19   | 0.9       | 0.9
SBS3     | mlp     | -Inf      | 1      | 0.75      | 0.5
SBS3     | mlp     | 0.92      | 0.85   | 0.9       | 0.9

The column `engine` is informational only: either `mlp` or `cnn`, the flavor of DL engines that was used
for optimal performance in training. `Recall` and `Precision` list estimated recall and precision values.

These threshold values are used to clip signature scores into discrete calls as shown in `ternary.call`:

sid            | M    | SBS1    | SBS2.13  | SBS3
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | P       | P        | P
TCGA.36.2540   | 10   | P       | N        | P
TCGA.24.1467   | 9    | I       | N        | I

There are three types of calls:
- `P` (positive): signature is present with minimum precision of 0.9
- `N` (negative): signature is not present even with low minimum precision of 0.5
- `I` (indeterminate): signature is present with minimum precision of 0.5 but not with 0.9

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
The arguments `CATALOG`, `CANCER_TYPER`, and `OUTPUT`specify the paths of input catalog data, cancer type, and output directory.\


## Installation

Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
In a Mac,

    $ brew install gsl

Also set 
```
PKG_CPPFLAGS = -I/opt/homebrew/Cellar/gsl/2.7.1/include
PKG_LIBS = -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
```
in `~/.R/MAKEVARS` on a Mac (check the links to make sure the correct version number is used).

DL-based filtering requires [tensorflow](https://tensorflow.org/install) and [pandas](https://pandas.pydata.org/getting_started.html). In a Mac,
  
    $ python3 -m pip install tensorflow pandas

See the tensorflow install page for further information.

We currently recommend downloading a source tar.gz file and
installing **DeepSig** via

    > install.packages('DeepSig_v0.9.8.tar.gz', repos = NULL, dependencies = TRUE)

which will also install dependencies.

### Catalog matrix generation

The catalog file as input data can be generated from MAF file using [maf2cat3()](https://github.com/mskcc/tempoSig/blob/master/man/maf2cat3.Rd). It requires a reference genome package, either [BSgenome.Hsapiens.UCSC.hg19](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19) 
or [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg38) installed:

    > library(BSgenome.Hsapiens.UCSC.hg19)
    > maf <- system.file('extdata', 'brca.maf', package = 'tempoSig')
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

See the [vignettes](https://github.com/mskcc/DeepSig/blob/main/old/denovo.html)

