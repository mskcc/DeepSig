#' Convert OncoTree Code into Cancer Type
#' @param oncotree Oncotree code
#' @examples 
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt',package='DeepSig'))
#' z <- DL.call(catalog = t(data), cancer.type = 'nsclc', alpha = 0.05)
#' head(z$exposure.fltrd)
#' 
#' @export

oncoTree <- function(onco){
  
  tissue.map <- c(
    'CNS/Brain' = 'cns',
    'Lymphoid' = 'pancancer',
    'Head and Neck' = 'head_neck',
    'Bone' = 'pancancer',
    'Myeloid' = 'pancancer',
    'Soft Tissue' = 'pancancer',
    'Ovary/Fallopian Tube' = 'ovary',
    'Breast' = 'breast',
    'Bowel' = 'colorectal',
    'Kidney' = 'kidney',
    'Lung' = 'nsclc',
    'Cervix' = 'pancancer',
    'Eye' = 'pancancer',
    'Skin' = 'skin',
    'Uterus' = 'uterus',
    'Peritoneum' = 'pancancer',
    'Penis' = 'pancancer',
    'Biliary Tract' = 'pancancer',
    'Thyroid' = 'pancancer',
    'Pancreas' = 'pancreas',
    'Bladder/Urinary Tract' = 'bladder',
    'Peripheral Nervous System' = 'pancancer',
    'Volva/vignia' = 'pancancer',
    'Testis' = 'pancancer',
    'Liver' = 'liver',
    'Other' = 'pancancer',
    'Esophagus/Stomach' = 'pancancer',
    'Peura' = 'pancancer',
    'Thymus' = 'pancancer',
    'Ampulla of Vater' = 'pancancer',
    'Prostate' = 'prostate',
    'Adrenal Gland' = 'pancancer')
  
  codes <- system.file('extdata', 'oncotree.json', package = 'DeepSig')
  x <- jsonlite::fromJSON(codes)
  tissue <- unique(x$tissue)
  tissue <- tissue[!is.na(tissue)]
  if(!onco %in% x$code) stop(paste0(onco, ' not in Oncotree codes'))
  
  if(onco %in% c('BGCT','BMGCT','EGCT','OGCT','OMGCT','NSGCT','GCTSTM','MGCT',
                 'VGCT','VMGCT'))
  flag <- x$code==onco
  model <- x$tissue[flag]
  
  if(x$parent[flag] %in% c('BGCT','OGCT','NSGCT','VGCT') | onco == 'EGCT'){
    model <- 'germ_cell'
  } else if(onco == 'CSCLC'){
    model <- 'sclc'
  }
  
  return(model)
  
}