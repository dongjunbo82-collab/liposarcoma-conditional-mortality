required<-c('data.table','survival','splines','ggplot2','svglite')
missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop('Install missing R packages: ',paste(missing,collapse=', '))
for(p in required)cat(p,as.character(packageVersion(p)),'\n')
