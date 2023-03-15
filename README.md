# tempoSig
Single-Base Substitution Mutational Signature Inference for WES and MSK-IMPACT

## Overview
**tempoSig** contains single-base substitution (SBS) mutational signature inference tools intended for 
data sets derived from non-WGS (e.g., WES or MSK-IMPACT) platforms.
It combines a maximum likelihood-based algorithm, fitting mutational catalogs to proportions of a pre-defined 
reference signature set, with deep learning (DL)-based models making binary calls of the presence/absence of
reference signatures. The latter enables a quality control filtering of signature proportions parametrized
by the total mutation count of each sample.

## Single-sample workflow
For single samples or a small cohort derived from WES or MSK-IMPACT, the main workflow comprises refitting
and DL-filtering. Input data are of the form of catalog matrix:

Mutation context | Tumor Sample Barcode 1 | Tumor Sample Barcode 2
---------------- | ---------------------- | ----------------------
A[C>A]A          |                      0 |                      1
A[C>A]C          |                      3 |                      0
A[C>A]G          |                      2 |                      5
A[C>A]T          |                      0 |                      2
C[C>A]A          |                      0 |                      0
C[C>A]C          |                      1 |                      1
C[C>A]G          |                      0 |                      0
C[C>A]T          |                      1 |                      1

Mutation context is the set of categories to which mutation data from sequencing experiments have been classified. Trinucleotide contexts shown above are the pyrimidine bases before and after mutation flanked by upstream and downstream nucleotides (96 in total). Each column corresponding to **Tumor_Sample_Barcode** contains non-negative counts of single nucleotide variants (SNVs). This trinucleotide matrix can be generated using the utility function [maf2cat2] or [maf2cat3]. MAF files can be generated from VCF files using [vcf2maf](https://github.com/mskcc/vcf2maf).

The function

    DL.call(catalog, cancer.type, model.path, ref.sig, threshold)

will find the pre-trained model corresponding to **cancer.type** and perform signature fitting and filtering.
Except for **catalog** and **cancer.type**, other arguments will default to those for pre-trained models if 
**cancer.type** is among those built-in:

cancer.type    | ICGC/PCAWG | TCGA
-------------- | ---------- | ------
breast         | Breast     | BRCA
ovarian        | Ovary      | OV
prostate       | Prostate    | PRAD
pancreatic     | Pancraeas  | PAAD
bladder        | Bladder    | BLCA
colorectal     | Colorectal | COAD
melanoma       | Skin       | SKCM
cns            | CNS        | GBM
lung           | Lung       | LUAD/LUSC
head_neck      | Head_neck  | HNSC
renal_cell     | Kidney     | KICH/KIRC/KIRP
endometrial    | Uterus     | UCEC
germ_cell      |            | TGCT

As an example, the following script will generate outputs for the TCGA-OV samples using **ovarian** pre-trained model:

    > library(tempoSig)
    > xcat <- read.table(system.file('extdata', 'tcga_ov_xcat.txt', package = 'tempoSig'), header = TRUE, sep = '\t')
    > z <- DL.call(catalog = t(xcat), cancer.type = 'ovarian')
    > names(z)
      [1] "exposure.raw"   "exposure.fltrd" "binary.call" "ref.sig"  
      [5] "score"

### Outputs
The return value is a list of 5 data frames, which include the raw proportions **exposure.raw**

sid            | M    | SBS1    | SBS3     | SBS2.13  | SBS5   | SBS8  | SBS17  | SBS18  | SBS31  | SBS35
-------------- | -----| ------- | -------- | -------- | ------ | ----- | ------ | ------ | ------ | -------
TCGA.61.2092   | 141  | 5.e-8   | 0.222    | 0.014    | 0.53   | 9e-8  | 1e-10  | 1e-8   | 0.02   | 0.21
TCGA.36.2540   | 10   | 0.322   | 0.154    | 4e-12    | 0.49   | 9e-11 | 8e-12  | 1e-10  | 0.04   | 4e-10 
TCGA.24.1467   | 9    | 0.341   | 1e-12    | 7e-9     | 0.015  | 2e-10 | 0.19   | 2e-9   | 2e-9   | 0.45

where the columns show the repertore of reference signatures present in the **ovarian** model being used,
the binary calls for each signature **binary.call**

sid            | M    | SBS1    | SBS3     | SBS2.13
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | FALSE   | TRUE     | FALSE
TCGA.36.2540   | 10   | FALSE   | FALSE    | FALSE
TCGA.24.1467   | 9    | FALSE   | FALSE    | FALSE

the filtered exposures **exposure.fltrd**, where the raw proprtions with negative calls have been zero'ed out,

sid            | M    | SBS1    | SBS3     | SBS2.13
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | 0       | 0.222    | 0
TCGA.36.2540   | 10   | 0       | 0        | 0
TCGA.24.1467   | 9    | 0       | 0        | 0

and the DL-score **score** used for making binary calls via thresholds targeted for false positive rate of 0.05:

sid            | M    | SBS1    | SBS3     | SBS2.13
-------------- | -----| ------- | -------- | -----------
TCGA.61.2092   | 141  | -2.44   | 7.19     | -2.02
TCGA.36.2540   | 10   | -2.20   | -2.48    | 0.27
TCGA.24.1467   | 9    | -0.27   | 2.05     | -0.64

The **ref.sig** is informational only since it has been used in training the model and refitting:

Mutation types | SBS1  | SBS3  | SBS2.13 | SBS8
-------------- | ----- | ----- | ------- | ----
A[C>A]A        | 0.027 | 0.014 | 0.004   | 0.02
A[C>A]C        | 0.013 | 0.010 | 0.002   | 0.01
A[C>A]G        | 0.004 | 0.002 | 0.001   | 0.003

### Optional arguments
The following are optional arguments mainly for custom filtering:
**model.path** is the path of the directory where the trained models to be used can be found.
The argument **ref.sig** is the reference signature (see above) file path. The **threshold** is
the set of cutoffs for DL scores chosen to achieve false positive rate of 0.05 or less during training:

model  | Mmin  | Mmax  | epochs  | S       | tpr  | fpr  | thr
-------| ----- | ----- | ------- | ------- | ---- | ---- | -----
1a     | 1     | 10    | 10      | SBS1    | 0.11 | 0.05 | 0.26
1a     | 1     | 10    | 10      | SBS2.13 | 0.26 | 0.05 | 0.57      
1a     | 1     | 10    | 10      | SBS3    | 0.15 | 0.05 | 2.05


The raw exposure proportions are inferred by maximum likelihood estimation using the reference signature set,
which have been extracted from de novo inference and/or selected from publicly available repositories, including
[COSMIC](https://cancer.sanger.ac.uk/cosmic/signatures) and [signal](https://signal.mutationalsignatures.com/).

## Installation

Compilation requires GNU Scientific Library [(GSL)](https://www.gnu.org/software/gsl/). In Ubuntu Linux,

    $ sudo apt-get install libgsl-dev
    
In a Mac,

    $ brew install gsl

DL-based filtering requires [tensorflow](https://tensorflow.org/install). In a Mac,
  
    $ pip install tensorflow

See the tensorflow install page for further information.

Dependencies include [Rcpp](https://cran.r-project.org/package=Rcpp), [gtools](https://cran.r-project.org/package=gtools), [argparse](https://cran.r-project.org/package=argparse), [coneproj](https://cran.r-project.org/package=coneproj),
and [reticulate](https://cran.r-project.org/package=reticulate).

Installing **tempoSig** via

    > devtools::install_github("mskcc/tempoSig_DL")

will also install dependencies.

### Catalog matrix generation

If a MAF file contains the column `Ref_Tri`, the catalog matrix can be generated using [maf2cat()](https://github.com/mskcc/tempoSig/blob/master/man/maf2cat.Rd) or its command-line wrapper:

    $ ./maf2cat2.R -h
    usage: ./maf2cat2.R [-h] MAF CATALOG

    Construct mutational catalog from MAF file with Ref_Tri column

    positional arguments:
      MAF         input MAF file
      CATALOG     output catalog file

    optional arguments:
      -h, --help  show this help message and exit

If the MAF file does not contain the column `Ref_Tri`, use [maf2cat3()](https://github.com/mskcc/tempoSig/blob/master/man/maf2cat3.Rd). It requires a reference genome package, either [BSgenome.Hsapiens.UCSC.hg19](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19) 
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

Under construction.

