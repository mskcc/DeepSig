ctype <- 'cns'
dir <- paste0('/Users/wooh/work/sig/tf/paper/survey/',ctype,'/e3/impact/edist/esample3/theta2/mc2/n1e6/')
refsig <- paste0('/Users/wooh/work/sig/tf/paper/survey/',ctype,'/',ctype,'_gih_refsig.txt')
system(paste0('cp ',refsig,' ./refsig.txt'))
thr.fl <- paste0(dir,'threshold/threshold_cut.txt')
thr <- read.table(thr.fl,header=TRUE,sep='\t')
system(paste0('cp ',thr.fl,' .'))
S <- unique(thr$S)
engine <- thr[match(S,thr$S),'engine']
names(engine) <- S

for(s in names(engine)){
  fl <- paste0(dir,'/',engine[s],'/',s,'/',s)
  print(fl)
  if(!dir.exists(fl)) stop(paste0(fl,' does not exist'))
  system(paste0('cp -R ', fl, ' .'))
}
