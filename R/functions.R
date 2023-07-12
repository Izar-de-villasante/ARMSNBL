
source("~/Projects/meth_functions.R")

# ¡¡¡¡¡ ATENTION !!!!!
# Modified (and uncommented) default bpplot function to add all labels and increase alpha.
# !!!!!!!!!!!!!!

#' # global opts:
#' options(future.globals.maxSize= 18912896000)
#' source("~/Projects/dmpdmrs.R")
#' # amke results folders:
#' make_results_dirs <- function(results_folder = "./results/", analysis_folder = "./analysis/",subf){
#'   params<-list(
#'     results_folder = results_folder,
#'     qc_folder  = paste(results_folder,subf,"/QC/",sep= .Platform$file.sep ),
#'     ss_clean_path = paste(analysis_folder,subf,sep= .Platform$file.sep),
#'     bplots_folder = paste(results_folder,subf,"plots/pca/bplots/",sep= .Platform$file.sep),
#'     corrplot_folder = paste(results_folder,subf,"plots/pca/corrplot/",sep= .Platform$file.sep),
#'     dmp_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
#'     dmpplots_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
#'     dmrs_folder = paste(results_folder,subf,"dmrs/",sep= .Platform$file.sep),
#'     pathway_folder = paste(results_folder,subf,"gopath/",sep= .Platform$file.sep)
#'     
#'   )
#'   # sapply(params,function(x)  dir.create(x,recursive=T,showWarnings = F))
#'   
#'   return(params)
#' }
#' 
#' #' @title Load data target factory.
#' #' @description Define 4 targets:
#' #' 1. Track the user-supplied data file.
#' #' 2. Read the data using `read_data()` (defined elsewhere).
#' #' 3. Strip category away.
#' #' 4. Generate category data.table.
#' #' @return A list of target objects.
#' #' @export
#' #' @param file Character, Sample sheet file path.
#' input_sample_sheet <- function(file,name) {
#'   name_ch<-deparse(substitute(name))
#'   
#'   
#'   list(
#'     tar_target_raw("file", file, format = "file", deployment = "main"),
#'     tar_target_raw("samps", quote(readRDS(file)), deployment = "main"),
#'     tar_target_raw("ss", quote(samps[-1,]), deployment = "main"),                                               # Strips category tags
#'     tar_target_raw("category", deployment = "main",
#'                    quote(data.table::as.data.table(t(samps[1,]),keep.rownames = T)))   # Category tags dict
#'   )
#' }
#' 
#' 
#' 
#' 
#' save_plot <-function(object,filename,path){
#'   dir.create(path,recursive = T)
#'   grDevices::png(file = paste0(path,"/",filename,".png"),
#'                  width = 480, # The width of the plot in inches
#'                  height = 620) # The height of the plot in inches
#'   object
#'   dev.off()
#'   
#' }
#' 
#' get_cols<-function(object,pal=NULL){
#'   g <- factor(object)
#'   if(is.null(pal))pal =  c(
#'     "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#FFE100", "#FF0010",
#'              "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
#'              "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
#'              "#740AFF", "#990000", "#FFFF80", "#5EF1F2", "#FF5005"
#'   )
#'   cols<-  pal[g]
#'   cols<-setNames(cols,g)
#'   return(cols)
#' }
#' name_rgset<-function(res,targets){
#'   
#'   res@colData <- methods::as(targets, "DataFrame")
#'   res@colData$barcode <- sapply(res@colData$barcode, function(x) substr(x,nchar(x)-18,nchar(x)))
#'   rownames(res@colData)<-sampleNames(res)<-targets$barcode
#'   # remove bad samples
#'   res[,!colnames(res) %in%"206702460034_R05C01"]->clean
#'   return(clean)
#' }
#' name_rgset2<-function(res,targets){
#'   cn<-colnames(res)
#'   class(res)
#'   data.table::setkey(targets,"barcode")
#'   pD <- data.frame(targets[cn,])
#'   pD$filenames <- targets$Basename
#'   #rownames(pD) <- colnames(res)
#'   res@colData <- methods::as(pD, "DataFrame")
#'   rownames(res@colData)<-cn
#'   
#'   sampleNames(res)<-rownames(res@colData)
#'   sampleNames(res)<-targets$barcode
#'   res[,!colnames(res) %in%"206702460034_R05C01"]->clean
#'   return(clean)
#' }
#' #Normalization functions:
#' 
#' noob <- function(rgSet){
#'   minfi::preprocessNoob(rgSet)
#' }
#' 
#' funn <- function(rgSet){
#'   minfi::preprocessFunnorm(rgSet)
#' }
#' noob_pq <- function(rgSet){
#'   minfi::preprocessNoob(rgSet)%>% minfi::preprocessQuantile()
#' }
#' pq <- function(rgSet){
#'   minfi::preprocessQuantile(rgSet)
#' }
#' 
#' Em2 <- function(rgSet, arraytype = NULL){
#'   pd <- rgSet@colData
#'   qc <- ENmix::QCinfo(rgSet,distplot = F)
#'   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#'   mSetSqn <- ENmix::mpreprocess(rgSet = rgSet,impute = T)
#'   if(is.null(arraytype)){
#'     arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
#'   }   
#'   
#'   if(arraytype=="EPIC"){
#'     mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
#'       mat = mSetSqn,
#'       array = "IlluminaHumanMethylationEPIC", 
#'       annotation = "ilm10b4.hg19"
#'     )
#'   }else mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
#'   
#'   mSetSqn@colData <- pd
#'   return(mSetSqn)
#' }
#' 
#' Em <- function(rgSet, arraytype = NULL){
#'   pd <- rgSet@colData
#'   qc <- ENmix::QCinfo(rgSet,distplot = F)
#'   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#'   mSetSqn <- minfi::mapToGenome(mdat)
#'   mSetSqn@colData <- pd
#'   return(mSetSqn)
#' }
#' 
#' filter<-function(targets, rgSet, sampGroups = NULL, sampNames = "Sample_Name",
#'                  frac=0.1,pval=0.01,remove_sex=TRUE,arraytype=NULL,cols=NULL,
#'                  qc_folder= "analysis/intermediate/QC"){
#'   # requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#'   # requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#'   if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
#'   n <- ncol(rgSet)
#'   #o <- rev(order(sampNames))
#'   #rgSet <- rgSet[, o]
#'   #sampNames <- sampNames[o]
#'   dir.create(qc_folder)
#'   if (is.null(sampGroups)) g <- rep(1, n) else g <- targets[[sampGroups]]
#'   if (is.null(g)) g<-1:n
#'   g<-factor(g)
#'   if(is.null(cols))cols<-get_cols(g)
#'   cols<-cols[!duplicated(cols)]
#'   # Check quality of the combined probe signal (testing against negative controls)
#'   # 1. Remove low quality samples. Defaults more than 10% we throw samples away.
#'   detP <- minfi::detectionP(rgSet, type = "m+u")
#'   grDevices::jpeg(file = paste0(qc_folder,"mean_detection_pvalues.jpeg"),width = 960,height = 1240)#,width = 480,height = 480)
#'   ylabels<-colnames(detP)
#'   par(mar=c(max(4.1,max(nchar(ylabels))/2.2) ,4.1 , 4.1, 2.1))
#'   
#'   barplot(colMeans(detP), col=cols, las=2,
#'           cex.names=0.8, ylim=c(0,max(0.002,max(colMeans(detP))*2)), main ="Mean detection p-values")
#'   graphics::abline(h=0.05,col="red")
#'   graphics::legend("topleft", legend=levels(g), fill=cols[1:length(levels(g))],
#'                    bg="white")
#'   grDevices::dev.off()
#'   
#'   # 2. Removing low-quality samples (with fraction of probes not passing pval)
#'   bad_samples <- colnames(detP)[colSums(detP >=pval)/nrow(detP) > frac]
#'   if(length(bad_samples)>0){
#'     warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
#'             frac*100,"% of the probes with p-val >", pval, "): \n ", paste(bad_samples,collapse = ", " ))
#'     rgSet <- rgSet[,setdiff(colnames(detP),bad_samples)]
#'   }else{
#'     cat("All samples passed detection P-value filter")
#'   }
#'   # 3. Removing low-quality probes (with p-value below pval)
#'   bad_probes<-which(rowSums(detP < pval) < ncol(rgSet)*(1-frac))
#'   rgSet <- rgSet[-c(bad_probes),]
#'   if(length(bad_samples)>0){
#'     warning("The following probes will be discarded since more than", frac*100,
#'             "% of the samples have detection p-values > ", pval, "): \n ", paste(bad_samples,collapse = ", " ))
#'   }else{
#'     cat("All samples passed detection P-value filter")
#'   }
#'   return(rgSet)
#' }
#' 
#' # comb <- function(results, x) {
#' #   i <- x$i
#' #   result <- x$result
#' #   if (x$error) {
#' #     cat(sprintf('master computing failed task %d\n', i))
#' #     # Could call function repeatedly until it succeeds,
#' #     # but that could hang the master
#' #     result <- try(fails_randomly(i))
#' #   }
#' #   results[i] <- list(result)  # guard against a NULL result
#' #   results
#' # }
#' # 
#' # comb <- function(a,b) {
#' #   minfi::combineArrays(object1 = a,object2=b)
#' # }
#' # 
#' read.metharray <- function(targets,folder,files=NULL,copy=FALSE, verbose = TRUE,
#'                            arraytype = NULL, ncores=NULL, extended = FALSE, force = TRUE){
#'   require(foreach)
#'   library(BiocGenerics)
#'   library(Biostrings)
#'   #Make cluster:
#'   if(is.null(files)) files<- targets$Basename
#'   if (copy ==TRUE){
#'     files<-paste0(folder,basename(targets$Basename))
#'   }
#'   
#'   
#'   ncores<-min(RcppParallel::defaultNumThreads(),nrow(ss))
#'   
#'   cl<- parallel::makePSOCKcluster(ncores,outfile="")
#'   parallel::clusterEvalQ(cl,{
#'     requireNamespace(c("minfi","S4Vectors","BiocGenerics"))
#'   })
#'   doParallel::registerDoParallel(cl)
#'   
#'   message("Reading multiple idat-files in parallel. Using ",ncores," cores.")
#'   res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
#'                         .combine='cbind',
#'                         .multicombine = F,
#'                         # .export = c("guessArrayTypes",".default.epic.annotation"),
#'                         .inorder=F,
#'                         .errorhandling = "pass",
#'                         .packages = c("Biostrings","BiocGenerics")
#'   )%dopar%{
#'     
#'     subdf<-as.data.frame(targets[it,])
#'     # handle idats path
#'     if (copy ==TRUE){
#'       requireNamespace("fs")
#'       fs::file_copy(paste0(subdf$Basename,"_Grn.idat"),new_path=folder,overwrite = T)
#'       fs::file_copy(paste0(subdf$Basename,"_Red.idat"),new_path=folder,overwrite=T)
#'     }
#'     
#'     # read idats:
#'     rgSet<-minfi::read.metharray(basenames = files[it], extended = extended, verbose = verbose, force =force)
#'     
#'     # arraytype:
#'     if (is.null(arraytype)){
#'       rgSet@annotation<-cnv.methyl:::guessArrayTypes(nrow(rgSet))
#'     }else{
#'       if (arraytype=="EPIC") {rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b4.hg19")}
#'       else if (arraytype=="450K"){rgSet@annotation <- c(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19")
#'       }else{rgSet@annotation<-cnv.methyl:::guessArrayTypes(nrow(rgSet))}
#'     }
#'     return(rgSet)
#'   }
#'   
#'   parallel::stopCluster(cl)
#'   # cn<-colnames(res)
#'   # class(res)
#'   # data.table::setkey(targets,"barcode")
#'   # pD <- data.frame(targets[cn,])
#'   # pD$filenames <- files
#'   # #rownames(pD) <- colnames(res)
#'   # res@colData <- methods::as(pD, "DataFrame")
#'   # rownames(res@colData)<-cn
#'   # colnames(res)<-cn
#'   #
#'   return(res)
#' }
#' 
#' 
#' densPlot <- function(dat, sampGroups = NULL, main = "", xlab = "Beta",
#'                      pal = brewer.pal(8, "Dark2"),cols=NULL,
#'                      xlim, ylim, add = TRUE, legend = TRUE, ...) {
#'   
#'   if (is(dat, "RGChannelSet") || is(dat, "MethylSet")) {
#'     b <- getBeta(dat)
#'   } else if (is(dat, "matrix")) {
#'     b <- dat
#'   } else {
#'     stop("argument 'dat' must be an 'RGChannelSet', a 'MethylSet' or ",
#'          "matrix.")
#'   }
#'   # NOTE: Have to ensure data are in memory as an ordinary vector before
#'   #       computing density
#'   d <- apply(b, 2, function(x) density(as.vector(x), na.rm = TRUE))
#'   if (missing(ylim)) ylim <- range(sapply(d, function(i) range(i$y)))
#'   if (missing(xlim)) xlim <- range(sapply(d, function(i) range(i$x)))
#'   if (is.null(sampGroups)) {
#'     sampGroups <- rep(1, ncol(b))
#'   } else if (length(sampGroups) == 1) {
#'     sampGroups <- rep(sampGroups, ncol(b))
#'   }
#'   sampGroups <- as.factor(sampGroups)
#'   if(is.null(cols))cols<-get_cols(sampGroups)
#'   
#'   # Save current graphical parameters
#'   opar <- par(no.readonly = TRUE)
#'   
#'   # Make sure d and SampGroups are in the same order:
#'   
#'   # Change the margins of the plot (the first is the bottom margin)
#'   par(mar = c(8, 4.1, 4.1, 2.1))
#'   
#'   if (add) {
#'     plot(x = 0,
#'          type = "n",
#'          ylim = ylim,
#'          xlim = xlim,
#'          ylab = "Density",
#'          xlab = xlab,
#'          main = main, ...)
#'     abline(h = 0, col = "grey80")
#'   }
#'   for (i in seq_along(d)) {
#'     lines(d[[i]], col = cols[i])
#'   }
#'   
#'   if (legend & length(levels(sampGroups)) > 1) {
#'     legend(x = "bottom",
#'            inset = c(0, -0.25), # You will need to fine-tune the second
#'            # value depending on the windows size
#'            legend = levels(sampGroups), 
#'            text.col = cols[!duplicated(cols)],
#'            # text.width=0.001,
#'            cex=0.7, pch=1, pt.cex = 1,
#'            lwd = 2,
#'            xpd = TRUE, # You need to specify this graphical parameter to add
#'            # the legend outside the plot area
#'            horiz = TRUE) # Horizontal legend. You can also set the number
#'     # of columns with the argument ncol
#'     # if horiz = FALSE
#'     
#'     # Back to the default graphical parameters
#'     
#'   }
#'   on.exit(par(opar))
#' }
#' 
#' 
#' #' Generate qc plots
#' #' @title generate qc plots for signal distribution prior to filtering
#' #' @param rgSet rgSet object containing channel intenisty values
#' #' @param sampNames variable containing barcodes or ids matching colnames of the rgsetdata
#' #' @param sampGroups variables to use for coloring groups
#' #' @param qc_folder path to the folder where plots will be saved
#' #' @return plots 
#' #' @author izar de Villasante
#' #' @export
#' #'
#' 
#' 
#' pathways_dict <- function(KEGG=T,GO=T ) {
#'   sym_KEGG = as.symbol(KEGG)
#'   sym_GO = as.symbol(GO)
#'   command_GO <- substitute(getGO(DO=GO),env = list(GO = GO))
#'   command_KEGG <- substitute(getKEGG(DO=KEGG),env = list(KEGG = KEGG))
#'   list(
#'     tar_target_raw("GO_dict",command_GO , deployment = "main"),
#'     tar_target_raw("KEGG_dict", command_KEGG, deployment = "main"),
#'     tar_target_raw("p_dict", quote(rbind(GO_dict,KEGG_dict)),
#'                    deployment = "main") 
#'   )
#' }
#' 
#' 
#' qc2 <- function(rgSet,sampGroups=NULL, sampNames= "Sample_Name", cols=NULL,
#'                 qc_folder="analysis/intermediate/QC/"){
#'   path=paste0(qc_folder,"Report.pdf")
#'   colData = rgSet@colData
#'   Sample_Group = colData[[sampGroups]]
#'   
#'   symrgSet = as.symbol(rgSet)
#'   symsG = as.symbol(Sample_Group)
#'   symsN = as.symbol(sampNames)
#'   symqc = as.symbol(qc_folder)
#'   symcols = as.symbol(cols)
#'   command_qcReport <- quote(
#'     minfi::qcReport(rgSet = symrgSet,pdf = path,
#'                     sampGroups = symsG))
#'   # sampNames = rgSet@colData[[sampNames]]),
#'   list(
#'     tar_target_raw("qc_report",command_qcReport )#,
#'     # tar_target_raw("KEGG_dict", command_KEGG, deployment = "main"),
#'     # tar_target_raw("p_dict", quote(rbind(GO_dict,KEGG_dict)),
#'     #                deployment = "main")
#'   )
#' }
#' qc<-function(rgSet,sampGroups=NULL, sampNames= "Sample_Name", cols=NULL,
#'              qc_folder="analysis/intermediate/QC/"){
#'   require(S4Vectors)
#'   require(Biostrings)
#'   require(Biobase)
#'   require(minfi)
#'   path=paste0(qc_folder,"Report.pdf")
#'   colData = rgSet@colData
#'   Sample_Group = colData[[sampGroups]]
#'   if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
#'   
#'   data.table::as.data.table(rgSet@colData)->ss
#'   ss[order(Sample_Group),barcode]->idx
#'   rgSet<-rgSet[,idx]
#'   if(is.null(cols))cols<-get_cols(factor(ss$Sample_Group))
#'   dir.create(qc_folder,recursive=T,showWarnings = F)
#'   minfi::qcReport(rgSet = rgSet,
#'                   pdf = paste0(qc_folder,"Report.pdf"),
#'                   sampGroups = rgSet@colData[[sampGroups]],
#'                   sampNames = rgSet@colData[[sampNames]])
#'   if(length(unique(rgSet@colData[[sampGroups]]))>1){
#'     grDevices::png(file = paste0(qc_folder,"density_plot.png"),
#'                    width = 960, # The width of the plot in inches
#'                    height = 1240) # The height of the plot in inches
#'     
#'     densPlot(
#'       rgSet, sampGroups = rgSet@colData[[sampGroups]],main = "Beta",
#'       pal =cols#[!duplicated(cols)]
#'     )
#'     grDevices::dev.off()
#'     grDevices::png(file = paste0(qc_folder,"bean_plot.png"),
#'                    width = 960, # The width of the plot in inches
#'                    height = 1240) # The height of the plot in inches
#'     minfi::densityBeanPlot(
#'       rgSet, sampGroups = rgSet@colData[[sampGroups]],
#'       pal=cols[!duplicated(cols)]
#'     )
#'     
#'     grDevices::dev.off()
#'     
#'   }
#'   mSet <- minfi::preprocessRaw(rgSet)
#'   qc   <- minfi::getQC(mSet)
#'   grDevices::png(file = paste0(qc_folder,"mean_qc.png"),   # The directory you want to save the file in
#'                  width = 960, # The width of the plot in inches
#'                  height = 960) # The height of the plot in inches
#'   minfi::plotQC(qc)
#'   
#'   grDevices::dev.off()
#' }
#' 
#' prep<-function(mSetSqn,remove_sex=TRUE,pval=0.01,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
#'   
#'   
#'   # 4. Removing probes with known SNPs at CpG site
#'   mSetSqn <-  minfi::mapToGenome(mSetSqn)
#'   mSetSqn <- minfi::dropLociWithSnps(mSetSqn)
#'   
#'   # 5. Removing cross reactive probes
#'   mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
#'   
#'   # 6. Sex. prediction & removal
#'   mSetSqn$predictedSex <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
#'   if(remove_sex){
#'     if(!is.null(arraytype)){anno<-get_anno(arraytype)
#'     }else{
#'       anno<-minfi::getAnnotation(mSetSqn)
#'       anno<-anno[!(anno$chr %in% c("chrX","chrY")),]
#'     }
#'   }
#'   return(mSetSqn)
#' }
#' 
#' 
#' addcol<-function(rgSet,col){
#'   require(data.table)
#'   dt<-data.table::as.data.table(rgSet@colData)
#'   droplevels.data.frame(dt)
#'   col<-unlist(col)
#'   data.table::setkey(dt, "barcode")
#'   dt[names(col),purity:=col]
#'   
#'   return(dt)
#' }
#' # Add info:
#' 
#' # Purity:
#' #cnv.methyl::purify()
#' 
#' # Celltype:
#' #1. cellCounts <- FlowSorted.Blood.450k::estimateCellCounts(rgSet)
#' #2. FlowSorted.Blood
#' 
#' # Copy Number Variation:
#' #cnv.methyl::Kc_get(ss )
#' 
#' # Sex:
#' # minfi::get_sex()
#' 
#' # Age:
#' # 
#' Enmix<-function(rgSet2){
#'   qcE<-ENmix::QCinfo(rgSet2)
#'   mdat<-ENmix::preprocessENmix(rgSet2, bgParaEst="oob", dyeCorr="RELIC",
#'                                QCinfo=qc, nCores=6)
#' }
#' 
#' # # obtaining the beta values
#' # beta_values <- getBeta(gmSet)
#' # colnames(beta_values) <- metadata$sample
#' # 
#' # saveRDS(beta_values, file = "results/beta_values.rds")
#' # 
#' # # PRINCIPAL COMPONENT ANALYSIS
#' # 
#' # # selecting the top 100 most variable CpG sites
#' top_beta <- function(beta_values, n=1000){
#'   sdv <- apply(beta_values, 1, sd)
#'   top100 <- names(head(sort(sdv,decreasing=T), n))
#'   beta_top100 <- beta_values[top100,]
#'   return(beta_top100)
#' }
#' pca_res <- function(beta_top100,scale=T, center=T){
#'   prcomp(t(beta_top100), scale=scale, center=center)
#' }
#' 
#' 
#' 
#' corpca <- function(beta_top100,metadata,vars=NULL,
#'                    path="./",filename="",title='PC1-6 clinical correlations'){
#'   requireNamespace("PCAtools")
#'   metadata<-data.frame(metadata,stringsAsFactors = T)
#'   rownames(metadata) <- metadata$barcode
#'   p<-PCAtools::pca(beta_top100[,rownames(metadata)],metadata = metadata, removeVar = 0.1)
#'   if(is.null(vars)){
#'     vars<-names(p$metadata)[sapply(p$metadata,function(x){
#'       !any(is.na(x))& length(unique(x))>1
#'     })
#'     ]
#'     
#'   }
#'   vars<-vars[sapply(p$metadata[,vars],function(x)length(unique(x)))>1]
#'   pcaplt<-PCAtools::eigencorplot(p,
#'                                  components = PCAtools::getComponents(p, 1:6),
#'                                  metavars = vars,
#'                                  col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
#'                                  cexCorval = 0.7,
#'                                  colCorval = 'white',
#'                                  fontCorval = 2,
#'                                  posLab = 'bottomleft',
#'                                  rotLabX = 45,
#'                                  posColKey = 'top',
#'                                  cexLabColKey = 1.5,
#'                                  scale = TRUE,
#'                                  main = title,
#'                                  colFrame = 'white',
#'                                  plotRsquared = FALSE)
#'   save_plot(pcaplt,path = path,filename = filename)
#'   return(pcaplt)
#' }
#' 
#' 
# bplot<-function(pca,ss,colgroup,s,combs=NULL,cols=NULL ,tit= NULL,folder = "analysis/pca/bplots/"){
# 
#   deps<-c("ggfortify","ggrepel","gplots","ggplot2")
#   sapply(deps,function(x){if(!require(x))renv::install(x)})
#   require(ggplot2)
#   require(gplots)
#   require(ggrepel)
#   require(ggfortify)
#   ss<-droplevels.data.frame(ss)
#   ss<-data.table::as.data.table(ss)
#   n<-nrow(pca$rotation)
#   lapply(colgroup, function(f) {
#     cols<-get_cols(object = as.character(unlist(ss[,.SD,.SDcols=f])),pal = cols)
#     cols <- cols[!duplicated(cols)]
# 
#     if(is.null(combs))combs<-combn(4,2)
#     # ss<-as.data.frame(ss)
#     #rownames(ss)<-ss$Sample_Name
#     for(i in 1:dim(combs)[2]){
#       if(is.null(tit))tit <- paste0("colored.by.",f, "_shape.", s)
#       ap<-ggplot2::autoplot(pca, x=combs[1,i], y=combs[2,i], data = ss, colour=f,shape=s,alpha=0.8,size=1)+
#         geom_text_repel(aes(label = Sample_Name, color = with(ss,get(f))),
#                         show.legend = FALSE, size = 1.5,max.overlaps = Inf,segment.size=0.2,min.segment.length = 0.8,point.size = 0.5)+
#         #scale_color_brewer(palette = "Paired")+
#         scale_color_manual(values = cols)+
#         #geom_point(aes(size=0.2))+
#         labs(colour=f)+
# 
#         ggtitle(tit)+
#         theme_bw(base_size = 7)+
#         theme(legend.key=element_blank(), legend.key.size=unit(1,"point"))
# 
# 
#       dir.create(folder)
#       #plot(ap)
#       ggsave(paste0(folder,tit,i,"_",n,".png"),plot=ap,width = 5.56, height = 2.80,units="in" )
#     }
#   })
#   return(folder)
# }
#' 
#' # Surrogate analysis:
#' surrogate<-function(grset,pheno,condition){
#'   mval<- getM(grset)
#'   pheno<-pData(grset)
#'   mod <- model.matrix(~as.factor(condition), data=pheno)
#'   mod0 <- model.matrix(~1, data=pheno)
#'   sva.results <- sva::sva(mval, mod, mod0)
#' }
#' 
#' 
#' # # PCA on the top 100 sites
#' # pca_res <- prcomp(t(beta_top100), scale=T, center=T)
#' # pca_all <- prcomp(t(beta_values), scale=T, center=T)
#' # 
#' # 
#' # 
#' # 
#' # 
#' # ## plotting PC1 and PC2 by condition
#' # autoplot(pca_res, x=1, y=2, data=metadata, colour="vascular_type", shape="type")+
#' #   geom_text_repel(aes(label=sample, color=vascular_type),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
#' #   labs(colour="Tissue", shape="Type")+
#' #   xlim(c(-0.5,0.3))+
#' #   theme_bw()+
#' #   ggtitle("PCA by tissue")
#' # 
#' 
#' #' Generate models
#' #' @title construct models and contrasts with limma
#' #' @param object your object containing beta values
#' #' @param group_var the variable used as independent variable
#' #' @param covs the set of variables to use as confounders
#' #' @param metadata the metadata or sample sheet
#' #' @param set a boolean vector to subset the observations
#' #' @param gr the group
#' #' @return fit2 ebayes model 
#' #' @author izar de Villasante
#' #' @export
#' #'
mod <- function(object,group_var,contrasts=NULL, covs=NULL, metadata,set = TRUE,gr=NULL,pairwise = T,
                singular=F,rename=NULL){
  data.table::setDT(as.data.frame(metadata))
  cont_sing=cont_pair=gr_cont_sing=gr_cont_pair=NULL
  metadata<-subset(metadata,set)
  metadata<-droplevels(metadata)
  object <- object[,metadata$barcode]
  covs_formula<-NULL
  if (!is.null(covs)&length(covs)>0)covs_formula<-paste0("+",paste0(covs,collapse=" + ",sep=""))
  design <- model.matrix(
    formula(
      paste("~ 0 +" , paste0(group_var),covs_formula,sep= " " )
    ),
    data = metadata
  )
  fit <- limma::lmFit(object,design)
  cols <- with(metadata,paste0(group_var, unique(get(group_var))))
  if(pairwise == T){
    cont_pair <- apply(combn(cols,2),2,function(x) paste(x,collapse = "-"))

    if(!is.null(gr)) {
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_pair <- apply(combn(gr_cols,2),2,function(x) paste(x,collapse = "-"))
    }
  }


  if(singular == T){
    cont_sing<-apply(combn(cols,length(cols)-1),2,function(x){
      var <- setdiff(cols,x)
      group <- contgroup(group_var,levels=x)
      contrast <- paste0(var,"-", group)
      return(contrast)
    } )
    if(!is.null(gr)) {
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_sing <- apply(combn(gr_cols,length(gr_cols)-1),2,function(x){
        var <- setdiff(gr_cols,x)
        group <- contgroup(group_var,levels=x)
        contrast <- paste0(var,"-", group)
        return(contrast)
      } )
    }
  }

  cont <- c(cont_sing,cont_pair)
  gr_cont <- c(gr_cont_sing,gr_cont_pair)

  if(is.null(contrasts)){
    contrasts<-c(cont,gr_cont)
  }else{
    contrasts=unlist(strsplit(contrasts,";"))
  }
  contMatrix <- limma::makeContrasts(
    contrasts=contrasts,
    levels=colnames(design)
  )
  # rename contrasts:
  # GR:
  if(!is.null(gr)) colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), gr_cols, gr, vectorize_all = FALSE)

  # Singular 1vsMean.
  large <- colnames(contMatrix) %in% c(cont_sing,gr_cont_sing)
  colnames(contMatrix)[large] <- sapply(
    colnames(contMatrix)[large], function(x) paste0("sing_",strsplit(x,"-")[[1]][1]))

  # remove group_var prefix:
  out<-tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      if(!is.null(rename)){colnames(contMatrix)<-rename
        # return("custom renamed")
      }else{colnames(contMatrix)<-stringi::stri_replace_all_fixed(colnames(contMatrix), group_var, "", vectorize_all = FALSE)
        # return("NULL rename, automatic names")
      }
    },error = function(cond) {
      colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), group_var, "", vectorize_all = FALSE)
      message(cond)
      # Choose a return value in case of error
      # return("some error, automatic names")
    }
  )
  print(out)
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit2 <- limma::eBayes(fit2)
  return(fit2)
}
#' 
#' 
#' contgroup<-function(name,levels){
#'   cols<-levels[grepl(name, levels, fixed = TRUE)]
#'   l<-length(cols)
#'   paste0("(",paste(cols,collapse = "+"),")/",l)
#' }
#' 
#' 
#' venns<-function(dt,gene_col="overlapping.genes",groupvar="Contrast",res="results/VENN/"){
#'   library(ggvenn)
#'   library(data.table)
#'   require(data.table)
#'   # Alternative install.packages("venn") -> Up to 7
#'   
#'   dir.create(res,recursive = T)
#'   # Hyper:
#'   
#'   glist<-dt[meandiff>0,unique(unlist(sapply(strsplit(overlapping.genes,","),"["))),by=groupvar]
#'   geneslist<-data.table::dcast.data.table(glist,V1~Contrast,fill = F)
#'   pdata<-lapply(geneslist[complete.cases(geneslist),] ,function(x)ifelse(x==F|is.na(x),F,T))|>as.data.table()
#'   
#'   # levels<-unique(as.character(unlist(glist[,.SD,.SDcols=groupvar])))
#'   
#'   g_venn_genes<-ggvenn(
#'     pdata, 
#'     setdiff(names(pdata),"V1"),#levels,
#'     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#'     stroke_size = 0.5, set_name_size = 4
#'   )
#'   
#'   ggsave(paste0("venn_genes_by_",groupvar,"_Hyper.png"),plot=g_venn_genes,path=res)
#'   
#'   # Hypo:
#'   
#'   glist<-dt[meandiff<0,unique(unlist(sapply(strsplit(overlapping.genes,","),"["))),by=groupvar]
#'   geneslist<-data.table::dcast.data.table(glist,V1~Contrast,fill = F)
#'   pdata<-lapply(geneslist[complete.cases(geneslist),] ,function(x)ifelse(x==F|is.na(x),F,T))|>as.data.table()
#'   
#'   # levels<-unique(as.character(unlist(glist[,.SD,.SDcols=groupvar])))
#'   
#'   g_venn_genes<-ggvenn(
#'     pdata, 
#'     setdiff(names(pdata),"V1"),#levels,
#'     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#'     stroke_size = 0.5, set_name_size = 4
#'   )
#'   
#'   ggsave(paste0("venn_genes_by_",groupvar,"_Hypo.png"),plot=g_venn_genes,path=res)
#'   
#'   
#'   # All:
#'   
#'   glist<-dt[,unique(unlist(sapply(strsplit(overlapping.genes,","),"["))),by=groupvar]
#'   geneslist<-data.table::dcast.data.table(glist,V1~Contrast,fill = F)
#'   pdata<-lapply(geneslist[complete.cases(geneslist),] ,function(x)ifelse(x==F|is.na(x),F,T))|>as.data.table()
#'   
#'   # levels<-unique(as.character(unlist(glist[,.SD,.SDcols=groupvar])))
#'   
#'   g_venn_genes<-ggvenn(
#'     pdata, 
#'     setdiff(names(pdata),"V1"),#levels,
#'     fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#'     stroke_size = 0.5, set_name_size = 4
#'   )
#'   
#'   ggsave(paste0("venn_genes_by_",groupvar,"_All.png"),plot=g_venn_genes,path=res)
#'   
#'   
#'   
#'   # ####
#'   # dt_genes[!is.na(gene_name1),.(sites=unique(rn)),by=condition]->sites
#'   # mutsites<-list(cond1=sites[condition==cond1,unique(sites)],
#'   #                cond2=sites[condition==cond2,unique(sites)]
#'   #                
#'   # )
#'   # names(mutsites)<-c(cond1,cond2)
#'   # 
#'   # 
#'   # g_venn_sites <- ggvenn(
#'   #   mutsites, 
#'   #   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#'   #   stroke_size = 0.5, set_name_size = 4,
#'   # )
#'   # g_venn_sites
#'   # ggsave(paste0(contrast,"_venn_sites.png"),plot=g_venn_sites,path=res)
#'   # # data.table::fwrite(genesets,paste0(res,"/",contrast,"_genes.txt"))
#'   # 
#'   # 
#'   # 
#'   # 
#'   # 
#'   # genesets_sites<-list(cond1=dt_genes[AI_cond==cond1,unique(gene_name1)],
#'   #                      cond2=dt_genes[AI_cond==cond2,unique(gene_name1)],
#'   #                      both=dt_genes[AI_cond=="both",unique(gene_name1)]
#'   #                      
#'   # )
#'   # names(genesets_sites)<-c(cond1,cond2,"both")
#'   # g_venn_comb<-ggvenn(
#'   #   genesets_sites, 
#'   #   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#'   #   stroke_size = 0.5, set_name_size = 4
#'   # )
#'   # g_venn_comb
#'   # ggsave(paste0(contrast,"_venn_comb.png"),plot=g_venn_comb,path=res)
#'   # 
#'   
#'   
#'   return(geneslist)
#' }
#' 
#' 
#' ##### PATHWAY
#' #gopath(dmrs_ANA,all.cpg=rownames(betas[,!is.na(ss_clean$ANA_dom)]),n=10,ann=ann)->pat
#' gopath <- function(object,all.cpg=NULL,n=Inf,ann=NULL,array.type = "EPIC",savepath=NULL){
#'   require(GenomicRanges)
#'   # object[[1]]
#'   #require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#'   object<-GRanges(object)
#'   requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
#'   future::plan("multisession")
#'   cont<-unique(object$Contrast)
#'   pathways <- future.apply::future_lapply(cont, function(x){
#'     
#'     results.ranges<- object[object$Contrast == x,]
#'     gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                     all.cpg = all.cpg,
#'                                                     collection = "GO", 
#'                                                     array.type = array.type,
#'                                                     anno = ann,
#'                                                     sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO"
#'     g
#'     },
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
#'     }
#'     )
#'     gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "GO", 
#'                                                       array.type = "EPIC",
#'                                                       genomic.features=c("TSS200","TSS1500","1stExon"),
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO_prom"
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
#'     }
#'     )
#'     gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "KEGG", 
#'                                                       array.type = "EPIC",
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "KEGG"
#'     g$ONTOLOGY <- NA
#'     g$TERM<-g$Description
#'     g$Description<-NULL
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
#'     }
#'     )
#'     result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
#'     return(result)
#'   },future.packages = c("S4Vectors","GenomicRanges"))
#'   
#'   
#'   out<-data.table::rbindlist(pathways)
#'   # if(!is.null(savepath))data.table::fwrite(out,path)
#'   return(out)
#'   
#' }
#' 
#' # cpgs <- GenomicRanges::GRanges(seqnames = anno$chr, 
#' #                                ranges = IRanges::IRanges(start = anno$pos, 
#' #                                                          end = anno$pos),
#' #                                strand = anno$strand,
#' #                                name = anno$Name)
#' # 
#' # overlaps <- GenomicRanges::findOverlaps(cpgs,regions)
#' # sig.cpg <- cpgs$name[overlaps@from]
#' # .getGO <- function(){
#' #   if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
#' #     stop("org.Hs.eg.db package required but not installed.")
#' #   egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
#' #   GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
#' #   d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
#' #   GeneID.PathID <- GeneID.PathID[d, ]
#' #   GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
#' #                                                       keys=unique(GeneID.PathID$go_id), 
#' #                                                       columns=c("GOID","ONTOLOGY","TERM"), 
#' #                                                       keytype="GOID"))
#' #   go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
#' #   
#' #   list(idList=go, idTable=GOID.TERM)
#' # }
#' # go <- .getGO()
#' # result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
#' #                   array.type=array.type, plot.bias=plot.bias, 
#' #                   prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
#' #                   fract.counts=fract.counts, 
#' #                   genomic.features = genomic.features,
#' #                   sig.genes = sig.genes)
#' # result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
#' # rownames(result) <- result$GOID
#' 
#' # subset(ss,!is.na(ss$ANA_dom1))
#' # technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
#' # 
#' # 
#' #  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
#' #  pheno<-ss[!is.na(ANA_dom1)]
#' #  mod <- model.matrix(
#' #    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
#' #    ,data = pheno
#' #  )
#' #  mod0 <- model.matrix(
#' #    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
#' #    ,data = pheno
#' #  )
#' #  sva.results <- sva(mval, mod, mod0)
#' # # 
#' # design <- model.matrix(
#' #   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
#' #   ,data = ss
#' #   )
#' # n.sv = sva::num.sv(betas,design,method="leek")
#' # 
#' # svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
#' # 
#' 
#' #
#' 
#' # From https://github.com/wjschne/spiro/blob/87f73ec37ceb0a7a9d09856ada8ae28d587a2ebd/R/spirograph.R
#' # Adapted under the CC0 1.0 Universal license: https://github.com/wjschne/spiro/blob/87f73ec37ceb0a7a9d09856ada8ae28d587a2ebd/LICENSE.md
#' 
#' #' gopath <- function(object,all.cpg=NULL,n=Inf,ann=NULL){
#' #'   require(GenomicRanges)
#' #'   # object[[1]]
#' #'   #require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#' #'   object<-GRanges(object)
#' #'   requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
#' #'   future::plan("multisession")
#' #'   cont<-unique(object$Contrast)
#' #'   pathways <- future.apply::future_lapply(cont, function(x){
#' #'     
#' #'     results.ranges<- object[object$Contrast == x,]
#' #'     gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#' #'                                                     all.cpg = all.cpg,
#' #'                                                     collection = "GO", 
#' #'                                                     array.type = "EPIC",
#' #'                                                     anno = ann,
#' #'                                                     sig.genes = T)
#' #'     g<-missMethyl::topGSA(gst, n=n)
#' #'     g$Contrast <- x
#' #'     g$method <- "GO"
#' #'     g
#' #'     },
#' #'     error = function(e) {
#' #'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
#' #'     }
#' #'     )
#' #'     gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#' #'                                                       all.cpg = all.cpg,
#' #'                                                       collection = "GO", 
#' #'                                                       array.type = "EPIC",
#' #'                                                       genomic.features=c("TSS200","TSS1500","1stExon"),
#' #'                                                       anno = ann,
#' #'                                                       sig.genes = T)
#' #'     g<-missMethyl::topGSA(gst, n=n)
#' #'     g$Contrast <- x
#' #'     g$method <- "GO_prom"
#' #'     g
#' #'     },
#' #'     
#' #'     error = function(e) {
#' #'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
#' #'     }
#' #'     )
#' #'     gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#' #'                                                       all.cpg = all.cpg,
#' #'                                                       collection = "KEGG", 
#' #'                                                       array.type = "EPIC",
#' #'                                                       anno = ann,
#' #'                                                       sig.genes = T)
#' #'     g<-missMethyl::topGSA(gst, n=n)
#' #'     g$Contrast <- x
#' #'     g$method <- "KEGG"
#' #'     g$ONTOLOGY <- NA
#' #'     g$TERM<-g$Description
#' #'     g$Description<-NULL
#' #'     g
#' #'     },
#' #'     
#' #'     error = function(e) {
#' #'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
#' #'     }
#' #'     )
#' #'     result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
#' #'     return(result)
#' #'   },future.packages = c("S4Vectors","GenomicRanges"))
#' #'   
#' #'   
#' #'   return(data.table::rbindlist(pathways))
#' #'   
#' #' }
#' 
#' 
#' gopath2 <- function(object,all.cpg=NULL,n=Inf,ann=NULL,array.type = "EPIC"){
#'   require(GenomicRanges)
#'   object<-GRanges(object)
#'   requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
#'   future::plan("multisession")
#'   cont<-unique(object$Contrast)
#'   pathways <- future.apply::future_lapply(cont, function(x){
#'     
#'     results.ranges<- object[object$Contrast == x,]
#'     gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                     all.cpg = all.cpg,
#'                                                     collection = "GO", 
#'                                                     array.type = array.type,
#'                                                     anno = ann,
#'                                                     sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO"
#'     g
#'     },
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
#'     }
#'     )
#'     gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "GO", 
#'                                                       array.type = "EPIC",
#'                                                       genomic.features=c("TSS200","TSS1500","1stExon"),
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO_prom"
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
#'     }
#'     )
#'     gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "KEGG", 
#'                                                       array.type = "EPIC",
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "KEGG"
#'     g$ONTOLOGY <- NA
#'     g$TERM<-g$Description
#'     g$Description<-NULL
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
#'     }
#'     )
#'     result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
#'     return(result)
#'   },future.packages = c("S4Vectors","GenomicRanges"))
#'   return(data.table::rbindlist(pathways))
#' }
#' 
#' # 
#' # pathways <- function(object,ann,array.type = NULL,contrast) {
#' #   require(GenomicRanges)
#' #   ann$start<-ann$end<-ann$pos
#' #   gr <- GRanges(ann)
#' # 
#' #   results.ranges<- GRanges(object[object$Contrast == contrast,])
#' #   sig <- subsetByOverlaps(gr,results.ranges)
#' # 
#' # 
#' # 
#' # 
#' # #
#' # #
#' #   gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges,
#' #                                                   all.cpg = ann$Name,
#' #                                                   collection = "GO",
#' #                                                   array.type = array.type,
#' #                                                   anno = ann,
#' #                                                   sig.genes = T)
#' #     {
#' #       go<-missMethyl:::.getGO()
#' #       collection<-go$idList
#' #       if(any(grepl("ALL", genomic.features)))genomic.features<-c("ALL", "TSS200","TSS1500",
#' #                                                                  "Body","1stExon","3'UTR",
#' #                                                                  "5'UTR","ExonBnd")
#' # 
#' #       result <- missMethyl::gsameth(sig.cpg=sig$Name, all.cpg=ann$Name, collection=collection,
#' #                                     array.type=array.type, plot.bias=F,
#' #                                     anno=ann, equiv.cpg=T,
#' #                                     genomic.features = "ALL",
#' #                                     sig.genes = T)
#' # 
#' #   }
#' #       {
#' #         {
#' # 
#' #         require(GenomicRanges)
#' # 
#' #         # Annotation:
#' #           if(any(grepl("ALL", genomic.features)))genomic.features<-c("ALL", "TSS200","TSS1500",
#' #                                                                      "Body","1stExon","3'UTR",
#' #                                                                      "5'UTR","ExonBnd")
#' # 
#' #           #remove probes not starting by cg:
#' #           ann <- ann[substr(row.names(ann),1,2) == "cg",]
#' # 
#' #           ann$start<-ann$end<-ann$pos
#' #         gr <- GRanges(ann)
#' # 
#' #         # DMRCate granges object:
#' #         results.ranges<- GRanges(object[object$Contrast == contrast,])
#' # 
#' #         # Significant genes annotation:
#' #         sig <- subsetByOverlaps(gr,results.ranges)
#' # 
#' #         # Pathways:
#' #         go<-missMethyl:::.getGO()
#' #         collection<-go$idList
#' # 
#' #         sig.cpg=sig$Name
#' #         dt_ann<-data.table::as.data.table(ann[ann$UCSC_RefGene_Name>0,])
#' #         dt_ann[,sig:=Name %in% sig.cpg]
#' #         # dt_sig <- dt_ann[Name %in% sig.cpg & ,]
#' # 
#' # 
#' #         dt_ann[,.(
#' #           GROUP=unlist(tstrsplit(UCSC_RefGene_Group, ";")),
#' #           SYMBOL=unlist(tstrsplit(UCSC_RefGene_Name, ";"))
#' #           ),by=c("Name","sig")]->dt_eg
#' # 
#' #         dt_eg<-dt_eg[GROUP %in% genomic.features,]
#' # 
#' #         dt_eg<-unique(dt_eg)
#' #         dt_eg$alias <- suppressWarnings(limma::alias2SymbolTable(dt_eg$SYMBOL))
#' # 
#' #         eg <- data.table::as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
#' #                                                      keys=AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db),
#' #                                                      columns=c("ENTREZID","SYMBOL"),
#' #                                                      keytype="ENTREZID"))
#' # 
#' #         # Add ENTREZID to annotation:
#' #         data.table::setkey(eg, SYMBOL)
#' #         dt_eg[,id:=paste0(Name,".",SYMBOL)]
#' #         merge(dt_eg,eg,by.x="SYMBOL",by.y="SYMBOL",all.x=T)->dtmp
#' #         merge(dt_eg,eg,by.x="alias",by.y="SYMBOL",all.x=T)->dtmp2
#' #         merge(dtmp,dtmp2,by=c("id","sig","alias","SYMBOL","Name","GROUP"))->dtmp3
#' #         dtmp3[,ENTREZID:=ifelse(is.na(ENTREZID.x),ENTREZID.y,ENTREZID.x)]
#' # 
#' #         # Remove trash vars:
#' #         dtmp3$ENTREZID.x<-dtmp3$ENTREZID.y<-dtmp3$Name<-NULL
#' #         dt_eg<-dtmp3
#' # 
#' # 
#' #         #
#' #         freq_genes <- table(dt_eg$ENTREZID)# frequence of each ENTREZID in annotation
#' #         eg.universe <- names(freq_genes) # All the ENTREZID in annotation
#' # 
#' # 
#' #         eg.sig <- unique(dt_eg[sig==T&!is.na(ENTREZID),ENTREZID]) # ENTREZID IN significative probes.
#' # 
#' #         test.de <- as.integer(eg.universe %in% eg.sig)
#' # 
#' #         ################################################# start here ############
#' #         multimap <- data.frame(table(dt_eg$Name))
#' #         multimap$Var1 <- as.character(multimap$Var1)
#' #         m3 <- match(flat.u$cpg, multimap$Var1)
#' #         flat.u$multimap <- multimap$Freq[m3]
#' # 
#' #         flat.u$inv.multimap <- 1/flat.u$multimap
#' # 
#' #         equivN <- tapply(flat.u$inv.multimap,flat.u$entrezid,sum)
#' #         mm <- match(eg.universe,names(equivN))
#' #         equivN <- equivN[mm]
#' # 
#' #         #sig.flat <- flat.u[!is.na(m1),]
#' #         if(any(grepl("ALL", genomic.features))){
#' #           sig.flat <- flat.u[!is.na(m1),]
#' #         } else {
#' #           # select only CpGs that map to certain genomic features
#' #           sig.flat <- flat.u[!is.na(m1) & flat.u$group %in% genomic.features, ]
#' #         }
#' # 
#' #         fract <- data.frame(weight=pmin(tapply(1/sig.flat$multimap,
#' #                                                sig.flat$entrezid,sum),
#' #                                         1))
#' # 
#' #         m4 <- match(sorted.eg.sig,rownames(fract))
#' #         fract.counts <- fract$weight[m4]
#' # 
#' #         out <- list(sig.eg = sorted.eg.sig, universe = eg.universe,
#' #                     freq = freq_genes, equiv =  equivN, de = test.de,
#' #                     fract.counts = data.frame(sigid=sorted.eg.sig,frac=fract.counts))
#' #         out
#' # 
#' #         ################################### end here ########################
#' #
#' #         out <- getMappedEntrezIDs(
#' #           sig.cpg=sig$Name, all.cpg=ann$Name,
#' #           array.type=array.type, ann, 
#' #           genomic.features = "ALL")
#' #         sorted.eg.sig <- out$sig.eg
#' #         eg.universe <- out$universe
#' #         freq_genes <- out$freq
#' #         test.de <- out$de
#' #         frac <- out$fract.counts
#' #         equiv <- out$equiv
#' #         
#' #         AnnotationDbi::select(
#' #           org.Hs.eg.db::org.Hs.eg.db,
#' #           keys = sorted.eg.sig,
#' #           columns = "SYMBOL")->gene.dict
#' #         dt_gene <- data.table::as.data.table(gene.dict)
#' #         data.table::setkey(dt_gene,ENTREZID)
#' #         # use "equivalent" no. of cpgs in odds calculation
#' #         pwf <- missMethyl:::.estimatePWF(D=test.de,bias=as.vector(equiv))
#' #         results <- data.table::data.table()#c("N","DE","P.DE","FDR")
#' #         results$cn <- names(collection)
#' #         results[,N:=unlist(lapply(collection,length))]
#' #         results[,DE:=unlist(lapply(collection, function(x) 
#' #           sum((sorted.eg.sig %in% x) * frac$frac)))]
#' #         Nuniverse <- length(eg.universe)
#' #         m <- length(sorted.eg.sig)
#' #          
#' #                                               
#' #         results[,InSet:=sapply(collection,function(x)sum(eg.universe %in% x,na.rm=T))]
#' #         
#' #         results[, odds:= (sum(pwf[InSet])/N) / (sum(pwf[!InSet])/(Nuniverse-N)),by=cn]
#' #         results[,P.DE:={
#' #           BiasedUrn::pWNCHypergeo(DE,
#' #                                   N,
#' #                                   Nuniverse-N,
#' #                                   m,odds,lower.tail=FALSE) + 
#' #             BiasedUrn::dWNCHypergeo(DE,
#' #                                     N,
#' #                                     Nuniverse-N,
#' #                                     m,odds)
#' #         },by=cn]
#' #         
#' #         ## This step needs better efficiency:
#' #         results[,
#' #                 SigGenesInSet:= paste(
#' #                   suppressMessages(AnnotationDbi::select(
#' #                     org.Hs.eg.db::org.Hs.eg.db,
#' #                     keys = intersect (sorted.eg.sig,collection[[cn]]),
#' #                     columns = "SYMBOL")
#' #                     ),collapse=",")
#' #                 ,by=cn]
#' # 
#' #           # Get gene symbols of significant genes
#' #           SigGenesEntrezID <- sorted.eg.sig[sorted.eg.sig %in% collection[[i]]]
#' #           SigGenesSymbol <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
#' #                                                                    keys = SigGenesEntrezID,
#' #                                                                    columns = "SYMBOL"))
#' #           results[i,SigGenesInSet:= paste(SigGenesSymbol$SYMBOL,collapse=",")]
#' #       }
#' #   
#' #   # gst$Contrast <- contrast
#' #   gst$method <- "GO"
#' #   gst
#' #   },
#' #   error = function(e) {
#' #     data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=contrast,method="GO")
#' #   }
#' #   )
#' #   gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#' #                                                     all.cpg = all.cpg,
#' #                                                     collection = "GO", 
#' #                                                     array.type = "EPIC",
#' #                                                     genomic.features=c("TSS200","TSS1500","1stExon"),
#' #                                                     anno = ann,
#' #                                                     sig.genes = T)
#' #   g<-missMethyl::topGSA(gst, n=n)
#' #   g$Contrast <- x
#' #   g$method <- "GO_prom"
#' #   g
#' #   },
#' #   
#' #   error = function(e) {
#' #     data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
#' #   }
#' #   )
#' #   gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#' #                                                     all.cpg = all.cpg,
#' #                                                     collection = "KEGG", 
#' #                                                     array.type = "EPIC",
#' #                                                     anno = ann,
#' #                                                     sig.genes = T)
#' #   g<-missMethyl::topGSA(gst, n=n)
#' #   g$Contrast <- x
#' #   g$method <- "KEGG"
#' #   g$ONTOLOGY <- NA
#' #   g$TERM<-g$Description
#' #   g$Description<-NULL
#' #   g
#' #   },
#' #   
#' #   error = function(e) {
#' #     data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
#' #   }
#' #   )
#' # }
#' # 
#' # plot_spirographs <- function(points) {
#' #   label <- "fixed_radius = %s, cycling_radius = %s"
#' #   points$parameters <- sprintf(label, points$fixed_radius, points$cycling_radius)
#' #   ggplot(points) +
#' #     geom_point(aes(x = x, y = y, color = parameters), size = 0.1) +
#' #     facet_wrap(~parameters) +
#' #     theme_gray(16) +
#' #     guides(color = "none")
#' # }
#' # 
#' # list(
#' #   tar_target(fixed_radius, sample.int(n = 10, size = 2)),
#' #   tar_target(cycling_radius, sample.int(n = 10, size = 2)),
#' #   tar_target(
#' #     points,
#' #     spirograph_points(fixed_radius, cycling_radius),
#' #     pattern = map(fixed_radius, cycling_radius)
#' #   ),
#' #   tar_target(
#' #     single_plot,
#' #     plot_spirographs(points),
#' #     pattern = map(points),
#' #     iteration = "list"
#' #   ),
#' #   tar_target(combined_plot, plot_spirographs(points))
#' # )
#' # 
#' # 
#' # library(missMethyl)
#' # result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
#' #                   array.type=array.type, plot.bias=plot.bias, 
#' #                   prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
#' #                   fract.counts=fract.counts, 
#' #                   genomic.features = genomic.features,
#' #                   sig.genes = T)
#' # result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
#' # rownames(result) <- result$GOID
#' # 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 










#' ############################## older stuff:
#' global opts:
#' options(future.globals.maxSize= 1891289600)
#' 
#' # amke results folders:
#' make_results_dirs <- function(results_folder = "./results/", analysis_folder = "./analysis/",subf){
#'   params<-list(
#'     results_folder = results_folder,
#'     qc_folder  = paste(analysis_folder,subf,"intermediate/QC/",sep= .Platform$file.sep ),
#'     ss_clean_path = paste(analysis_folder,subf,sep= .Platform$file.sep),
#'     bplots_folder = paste(results_folder,subf,"plots/pca/bplots/",sep= .Platform$file.sep),
#'     corrplot_folder = paste(results_folder,subf,"plots/pca/corrplot/",sep= .Platform$file.sep),
#'     dmp_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
#'     dmpplots_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
#'     dmrs_folder = paste(results_folder,subf,"dmrs/",sep= .Platform$file.sep),
#'     pathway_folder = paste(results_folder,subf,"gopath/",sep= .Platform$file.sep)
#'     
#'   )
#'   sapply(params,function(x)  dir.create(x,recursive=T,showWarnings = F))
#'   
#'   return(params)
#' }
#' name_rgset<-function(res,targets){
#'   require(minfi)
#'   # cn<-colnames(res)
#'   # class(res)
#'   # data.table::setkey(targets,"barcode")
#'   # pD <- data.frame(targets[cn,])
#'   # pD$filenames <- targets$Basename
#'   # #rownames(pD) <- colnames(res)
#'   # res@colData <- methods::as(pD, "DataFrame")
#'   # rownames(res@colData)<-cn
#'   
#'     # sampleNames(res)<-rownames(res@colData)
#'     sampleNames(res)<-targets$barcode
#'     return(res)
#' }
#' 
#' #Normalization functions:
#' 
#' noob <- function(rgSet){
#'   minfi::preprocessNoob(rgSet)
#' }
#' 
#' funn <- function(rgSet){
#'   minfi::preprocessFunnorm(rgSet)
#' }
#' noob_pq <- function(rgSet){
#'   minfi::preprocessNoob(rgSet)%>% minfi::preprocessQuantile()
#' }
#' pq <- function(rgSet){
#'   minfi::preprocessQuantile(rgSet)
#' }
#' 
#' Em2 <- function(rgSet, arraytype = NULL){
#'   pd <- rgSet@colData
#'   qc <- ENmix::QCinfo(rgSet,distplot = F)
#'   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#'   mSetSqn <- ENmix::mpreprocess(rgSet = rgSet,impute = T)
#'   if(is.null(arraytype)){
#'     arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
#'   }   
#'   
#'   if(arraytype=="EPIC"){
#'     mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
#'       mat = mSetSqn,
#'       array = "IlluminaHumanMethylationEPIC", 
#'       annotation = "ilm10b4.hg19"
#'     )
#'   }else mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
#'   
#'   mSetSqn@colData <- pd
#'   return(mSetSqn)
#' }
#' 
#' Em <- function(rgSet, arraytype = NULL){
#'   pd <- rgSet@colData
#'   qc <- ENmix::QCinfo(rgSet,distplot = F)
#'   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#'   mSetSqn <- minfi::mapToGenome(mdat)
#'   mSetSqn@colData <- pd
#'   return(mSetSqn)
#' }
#' 
#' filter<-function(targets, rgSet,sampGroups=NULL,sampNames="Sample_Name",frac=0.1,pval=0.01,remove_sex=TRUE,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
#'   # requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#'   # requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#'   if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
#'   n <- ncol(rgSet)
#'   #o <- rev(order(sampNames))
#'   #rgSet <- rgSet[, o]
#'   #sampNames <- sampNames[o]
#'   if (is.null(sampGroups)) g <- rep(1, n) else g <- targets[[sampGroups]]
#'   if (is.null(g)) g<-1:n
#'   g <- factor(g)
#'   
#'   pal =  c(
#'     "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
#'     "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
#'     "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
#'     "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
#'   
#'   
#'   # Check quality of the combined probe signal (testing against negative controls)
#'   # 1. Remove low quality samples. Defaults more than 10% we throw samples away.
#'   detP <- minfi::detectionP(rgSet, type = "m+u")
#'   grDevices::jpeg(file = paste0(qc_folder,"mean_detection_pvalues.jpeg"))#,width = 480,height = 480)
#'   ylabels<-colnames(detP)
#'   par(mar=c(max(4.1,max(nchar(ylabels))/2.2) ,4.1 , 4.1, 2.1))
#'   
#'   barplot(colMeans(detP), col=pal[g], las=2,
#'           cex.names=0.8, ylim=c(0,max(0.002,max(colMeans(detP))*2)), main ="Mean detection p-values")
#'   graphics::abline(h=0.05,col="red")
#'   graphics::legend("topleft", legend=levels(g), fill=pal[1:length(levels(g))],
#'                    bg="white")
#'   grDevices::dev.off()
#'   
#'   # 2. Removing low-quality samples (with fraction of probes not passing pval)
#'   bad_samples <- colnames(detP)[colSums(detP >=pval)/nrow(detP) > frac]
#'   if(length(bad_samples)>0){
#'     warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
#'             frac*100,"% of the probes with p-val >", pval, "): \n ", paste(bad_samples,collapse = ", " ))
#'     rgSet <- rgSet[,setdiff(colnames(detP),bad_samples)]
#'   }else{
#'     cat("All samples passed detection P-value filter")
#'   }
#'   # 3. Removing low-quality probes (with p-value below pval)
#'   bad_probes<-which(rowSums(detP < pval) < ncol(rgSet)*(1-frac))
#'   rgSet <- rgSet[-c(bad_probes),]
#'   if(length(bad_samples)>0){
#'     warning("The following probes will be discarded since more than", frac*100,
#'             "% of the samples have detection p-values > ", pval, "): \n ", paste(bad_samples,collapse = ", " ))
#'   }else{
#'     cat("All samples passed detection P-value filter")
#'   }
#'   return(rgSet)
#' }
#' 
#' # comb <- function(results, x) {
#' #   i <- x$i
#' #   result <- x$result
#' #   if (x$error) {
#' #     cat(sprintf('master computing failed task %d\n', i))
#' #     # Could call function repeatedly until it succeeds,
#' #     # but that could hang the master
#' #     result <- try(fails_randomly(i))
#' #   }
#' #   results[i] <- list(result)  # guard against a NULL result
#' #   results
#' # }
#' # 
#' # comb <- function(a,b) {
#' #   minfi::combineArrays(object1 = a,object2=b)
#' # }
#' # 
#' # read.metharray <- function(targets,folder,files=NULL,copy=FALSE, verbose = TRUE,
#' #                            arraytype = NULL, ncores=NULL, extended = FALSE, force = TRUE){
#' #   require(foreach)
#' #   library(BiocGenerics)
#' #   #Make cluster:
#' #   if(is.null(files)) files<- targets$Basename
#' #   if (copy ==TRUE){
#' #     files<-paste0(folder,basename(targets$Basename))
#' #   }
#' #   
#' #   
#' #   ncores<-min(RcppParallel::defaultNumThreads(),nrow(ss))
#' #   
#' #   cl<- parallel::makePSOCKcluster(ncores,outfile="")
#' #   parallel::clusterEvalQ(cl,{
#' #     requireNamespace(c("minfi","S4Vectors","BiocGenerics"))
#' #   })
#' #   doParallel::registerDoParallel(cl)
#' #   
#' #   message("Reading multiple idat-files in parallel. Using ",ncores," cores.")
#' #   res<-foreach::foreach(it=itertools::isplitIndices(nrow(targets), chunks=ncores),#isplitIndices(1400,chunks=ncores),
#' #                         .combine='comb',
#' #                         .multicombine = F,
#' #                         # .export = c("guessArrayTypes",".default.epic.annotation"),
#' #                         .inorder=F,
#' #                         .errorhandling = "pass"
#' #   )%dopar%{
#' #     
#' #     subdf<-as.data.frame(targets[it,])
#' #     # handle idats path
#' #     if (copy ==TRUE){
#' #       requireNamespace("fs")
#' #       fs::file_copy(paste0(subdf$Basename,"_Grn.idat"),new_path=folder,overwrite = T)
#' #       fs::file_copy(paste0(subdf$Basename,"_Red.idat"),new_path=folder,overwrite=T)
#' #     }
#' #     
#' #     # read idats:
#' #     rgSet<-minfi::read.metharray(basenames = files[it], extended = extended, verbose = verbose, force =force)
#' #     
#' #     # arraytype:
#' #     if (is.null(arraytype)){
#' #       rgSet@annotation<-cnv.methyl:::guessArrayTypes(nrow(rgSet))
#' #     }else{
#' #       if (arraytype=="EPIC") {rgSet@annotation <- c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b4.hg19")}
#' #       else if (arraytype=="450K"){rgSet@annotation <- c(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19")
#' #       }else{rgSet@annotation<-cnv.methyl:::guessArrayTypes(nrow(rgSet))}
#' #     }
#' #     return(rgSet)
#' #   }
#' #   
#' #   parallel::stopCluster(cl)
#' #   # cn<-colnames(res)
#' #   # class(res)
#' #   # data.table::setkey(targets,"barcode")
#' #   # pD <- data.frame(targets[cn,])
#' #   # pD$filenames <- files
#' #   # #rownames(pD) <- colnames(res)
#' #   # res@colData <- methods::as(pD, "DataFrame")
#' #   # rownames(res@colData)<-cn
#' #   # colnames(res)<-cn
#' #   # 
#' #   return(res)
#' # }
#' 
#' 
#' #' Generate qc plots
#' #' @title generate qc plots for signal distribution prior to filtering
#' #' @param rgSet rgSet object containing channel intenisty values
#' #' @param sampNames variable containing barcodes or ids matching colnames of the rgsetdata
#' #' @param sampGroups variables to use for coloring groups
#' #' @param qc_folder path to the folder where plots will be saved
#' #' @return plots 
#' #' @author izar de Villasante
#' #' @export
#' #'
#' qc <- function(rgSet,sampGroups=NULL, sampNames= "Sample_Name",qc_folder="analysis/intermediate/QC/"){
#'   if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
#'   # n <- ncol(rgSet)
#'   # o <- rev(order(sampNames))
#'   # rgSet <- rgSet[, o]
#'   # sampNames <- sampNames[o]
#'   # if (is.null(sampGroups))
#'   #   sampGroups <- rep(1, n)
#'   # sampGroups <- sampGroups[o]
#'   dir.create(qc_folder,recursive=T,showWarnings = F)
#'   minfi::qcReport(rgSet = rgSet,
#'                   pdf = paste0(qc_folder,"Report.pdf"),
#'                   sampGroups = rgSet@colData[[sampGroups]],
#'                   sampNames = rgSet@colData[[sampNames]])
#'   if(length(unique(rgSet@colData[[sampGroups]]))>1){
#'     grDevices::png(file = paste0(qc_folder,"density_plot.png"),
#'                    width = 480, # The width of the plot in inches
#'                    height = 620) # The height of the plot in inches
#'     pal =  c(
#'       "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
#'       "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
#'       "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
#'       "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
#'     minfi::densityPlot(
#'       rgSet, sampGroups = rgSet@colData[[sampGroups]],main = "Beta",
#'       pal =pal
#'     )
#'     grDevices::dev.off()
#'     grDevices::png(file = paste0(qc_folder,"bean_plot.png"),
#'                    width = 480, # The width of the plot in inches
#'                    height = 620) # The height of the plot in inches
#'     minfi::densityBeanPlot(
#'       rgSet, sampGroups = rgSet@colData[[sampGroups]],
#'       pal=pal
#'     )
#'     
#'     grDevices::dev.off()
#'     
#'   }
#'   mSet <- minfi::preprocessRaw(rgSet)
#'   qc   <- minfi::getQC(mSet)
#'   grDevices::png(file = paste0(qc_folder,"mean_qc.png"),   # The directory you want to save the file in
#'                  width = 480, # The width of the plot in inches
#'                  height = 480) # The height of the plot in inches
#'   minfi::plotQC(qc)
#'   
#'   grDevices::dev.off()
#' }
#' 
#' prep<-function(mSetSqn,remove_sex=TRUE,pval=0.01,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
#'   
#'   
#'   # 4. Removing probes with known SNPs at CpG site
#'   mSetSqn <-  minfi::mapToGenome(mSetSqn)
#'   mSetSqn <- minfi::dropLociWithSnps(mSetSqn)
#'   
#'   # 5. Removing cross reactive probes
#'   mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
#'   
#'   # 6. Sex. prediction & removal
#'   mSetSqn$predictedSex <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
#'   if(remove_sex){
#'     if(!is.null(arraytype)){anno<-get_anno(arraytype)
#'     }else{
#'       anno<-minfi::getAnnotation(mSetSqn)
#'       anno<-anno[!(anno$chr %in% c("chrX","chrY")),]
#'     }
#'   }
#'   return(mSetSqn)
#' }
#' 
#' 
#' # Add info:
#' 
#' # Purity:
#' #cnv.methyl::purify()
#' 
#' # Celltype:
#' #1. cellCounts <- FlowSorted.Blood.450k::estimateCellCounts(rgSet)
#' #2. FlowSorted.Blood
#' 
#' # Copy Number Variation:
#' #cnv.methyl::Kc_get(ss )
#' 
#' # Sex:
#' # minfi::get_sex()
#' 
#' # Age:
#' # 
#' Enmix<-function(rgSet2){
#'   qcE<-ENmix::QCinfo(rgSet2)
#'   mdat<-ENmix::preprocessENmix(rgSet2, bgParaEst="oob", dyeCorr="RELIC",
#'                                QCinfo=qc, nCores=6)
#' }
#' 
#' # # obtaining the beta values
#' # beta_values <- getBeta(gmSet)
#' # colnames(beta_values) <- metadata$sample
#' # 
#' # saveRDS(beta_values, file = "results/beta_values.rds")
#' # 
#' # # PRINCIPAL COMPONENT ANALYSIS
#' # 
#' # # selecting the top 100 most variable CpG sites
#' top_beta <- function(beta_values, n=1000){
#'   sdv <- apply(beta_values, 1, sd)
#'   top100 <- names(head(sort(sdv,decreasing=T), n))
#'   beta_top100 <- beta_values[top100,]
#'   return(beta_top100)
#' }
#' pca_res <- function(beta_top100,scale=T, center=T){
#'   prcomp(t(beta_top100), scale=scale, center=center)
#' }
#' 
#' 
#' 
#' corpca <- function(beta_top100,metadata,vars=NULL,title='PC1-6 clinical correlations'){
#'   requireNamespace("PCAtools")
#'   p<-PCAtools::pca(beta_top100,metadata = metadata, removeVar = 0.1)
#'   if(is.null(vars)){
#'     vars<-names(p$metadata)[sapply(p$metadata,function(x){
#'       !any(is.na(x))& length(unique(x))>1
#'     })
#'     ]
#'     vars[sapply(p$metadata,function(x)length(unique(x)))>1]
#'   }
#'   PCAtools::eigencorplot(p,
#'                components = PCAtools::getComponents(p, 1:6),
#'                metavars = vars,
#'                col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
#'                cexCorval = 0.7,
#'                colCorval = 'white',
#'                fontCorval = 2,
#'                posLab = 'bottomleft',
#'                rotLabX = 45,
#'                posColKey = 'top',
#'                cexLabColKey = 1.5,
#'                scale = TRUE,
#'                main = title,
#'                colFrame = 'white',
#'                plotRsquared = FALSE)
#'   
#' }
#' 
#' 
#' bplot<-function(pca,ss,colgroup,s,combs=NULL, tit= NULL,folder = "analysis/pca/bplots/"){
#'   library(ggplot2)
#'   library(gplots)
#'   library(ggrepel)
#'   library(ggfortify)
#'   ss<-droplevels.data.frame(ss)
#'   pal=c("#1b9e77","#d95f02", "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
#'         "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
#'         "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
#'         "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
#'   # pal =  c(
#'   #   "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
#'   #   "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
#'   #   "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
#'   #   "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
#'   lapply(colgroup, function(f) {
#'     f_len <- length(unique(with(ss,get(f))))
#'     cols <- pal[1:f_len]
#'     names(cols)<- unique(with(ss,get(f)))
#'     if(is.null(combs))combs<-combn(4,2)
#'     ss<-as.data.frame(ss)
#'     #rownames(ss)<-ss$Sample_Name
#'     for(i in 1:dim(combs)[2]){
#'       if(is.null(tit))tit <- paste0("colored.by.",f, "_shape.", s)
#'       ap<-ggplot2::autoplot(pca, x=combs[1,i], y=combs[2,i], data = ss, colour=f,shape=s,alpha=0.7,size=1)+
#'         geom_text_repel(aes(label = Sample_Name, color = with(ss,get(f))),
#'                         show.legend = FALSE, size = 1.5,max.overlaps = Inf,segment.size=0.2,min.segment.length = 0.8,point.size = 0.5)+
#'         #scale_color_brewer(palette = "Paired")+
#'         scale_color_manual(values = cols)+
#'         #geom_point(aes(size=0.2))+
#'         labs(colour=f)+ 
#'         
#'         ggtitle(tit)+
#'         theme_bw(base_size = 7)+
#'         theme(legend.key=element_blank(), legend.key.size=unit(1,"point"))
#'       
#'       
#'       dir.create(folder)
#'       #plot(ap)
#'       ggsave(paste0(folder,tit,i,".png"),plot=ap,width = 5.56, height = 2.80,units="in" )
#'     }
#'   })
#'   return(folder)
#' }
#' 
#' # Surrogate analysis:
#' surrogate<-function(grset,pheno,condition){
#'   mval<- getM(grset)
#'   pheno<-pData(grset)
#'   mod <- model.matrix(~as.factor(condition), data=pheno)
#'   mod0 <- model.matrix(~1, data=pheno)
#'   sva.results <- sva::sva(mval, mod, mod0)
#' }
#' 
#' 
#' # # PCA on the top 100 sites
#' # pca_res <- prcomp(t(beta_top100), scale=T, center=T)
#' # pca_all <- prcomp(t(beta_values), scale=T, center=T)
#' # 
#' # 
#' # 
#' # 
#' # 
#' # ## plotting PC1 and PC2 by condition
#' # autoplot(pca_res, x=1, y=2, data=metadata, colour="vascular_type", shape="type")+
#' #   geom_text_repel(aes(label=sample, color=vascular_type),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
#' #   labs(colour="Tissue", shape="Type")+
#' #   xlim(c(-0.5,0.3))+
#' #   theme_bw()+
#' #   ggtitle("PCA by tissue")
#' # 
#' 
#' #' Generate models
#' #' @title construct models and contrasts with limma
#' #' @param object your object containing beta values
#' #' @param group_var the variable used as independent variable
#' #' @param covs the set of variables to use as confounders
#' #' @param metadata the metadata or sample sheet
#' #' @param set a boolean vector to subset the observations
#' #' @param gr the group
#' #' @return fit2 ebayes model 
#' #' @author izar de Villasante
#' #' @export
#' #'
#' mod <- function(object, group_var, covs=NULL, metadata,set = TRUE,gr=NULL,pairwise = T,
#'                 singular=F){
#'   data.table::setDT(as.data.frame(metadata))
#'   cont_sing=cont_pair=gr_cont_sing=gr_cont_pair=NULL
#'   metadata<-subset(metadata,set)
#'   metadata<-droplevels(metadata)
#'   object <- object[,metadata$barcode]
#'   covs_formula<-NULL
#'   if (!is.null(covs)&length(covs)>0)covs_formula<-paste0("+",paste0(covs,collapse=" + ",sep=""))
#'   design <- model.matrix( 
#'     formula(
#'       paste("~ 0 +" , paste0(group_var),covs_formula,sep= " " )
#'     ),
#'     data = metadata
#'   )
#'   fit <- limma::lmFit(object,design)
#'   cols <- with(metadata,paste0(group_var, unique(get(group_var))))
#'   if(pairwise == T){
#'     cont_pair <- apply(combn(cols,2),2,function(x) paste(x,collapse = "-"))
#'     
#'     if(!is.null(gr)) { 
#'       gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
#'       gr_cont_pair <- apply(combn(gr_cols,2),2,function(x) paste(x,collapse = "-"))
#'     }
#'   }
#'   
#'   
#'   if(singular == T){
#'     cont_sing<-apply(combn(cols,length(cols)-1),2,function(x){
#'       var <- setdiff(cols,x)
#'       group <- contgroup(group_var,levels=x)
#'       contrast <- paste0(var,"-", group)
#'       return(contrast)
#'     } )
#'     if(!is.null(gr)) { 
#'       gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
#'       gr_cont_sing <- apply(combn(gr_cols,length(gr_cols)-1),2,function(x){
#'         var <- setdiff(gr_cols,x)
#'         group <- contgroup(group_var,levels=x)
#'         contrast <- paste0(var,"-", group)
#'         return(contrast)
#'       } )
#'     }
#'   }
#'   
#'   cont <- c(cont_sing,cont_pair)
#'   gr_cont <- c(gr_cont_sing,gr_cont_pair)
#'   contMatrix <- limma::makeContrasts(
#'     contrasts=c(cont,gr_cont),
#'     levels=colnames(design)
#'   )
#'   # rename contrasts:
#'   # GR: 
#'   if(!is.null(gr)) colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), gr_cols, gr, vectorize_all = FALSE)
#'   
#'   # Singular 1vsMean.
#'   large <- colnames(contMatrix) %in% c(cont_sing,gr_cont_sing)
#'   colnames(contMatrix)[large] <- sapply(
#'     colnames(contMatrix)[large], function(x) paste0("sing_",strsplit(x,"-")[[1]][1]))
#'   
#'   # remove group_var prefix:
#'   colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), group_var, "", vectorize_all = FALSE)
#'   fit2 <- limma::contrasts.fit(fit, contMatrix)
#'   fit2 <- limma::eBayes(fit2)
#'   return(fit2)
#' }
#' 
#' 
#' contgroup<-function(name,levels){
#'   cols<-levels[grepl(name, levels, fixed = TRUE)]
#'   l<-length(cols)
#'   paste0("(",paste(cols,collapse = "+"),")/",l)
#' }
#' 
#' 
#' betasdmps <- function(betas,dmps,rgSet,anno_cols=c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")){
#'   require(data.table)
#'   dmps<-data.table::setDT(as.data.frame(dmps))
#'   ss<-rgSet@colData
#'   pvals <- minfi::detectionP(rgSet)
#'   p<-pvals[dmps$Name,]
#'   colnames(p)<-paste0(ss[colnames(pvals),"Sample_Name"],".Detection_Pval")
#'   b<-betas[dmps$Name,]
#'   colnames(b)<-paste0(ss[colnames(betas),"Sample_Name"],".AVG_betas")
#'   d<-dmps[,.SD,.SDcols=anno_cols]
#'   df<-cbind(d,b,p)
#'   df<-data.table::setDT(as.data.frame(df))
#'   data.table::setorder(df,chr,pos)
#'   # df[order(chr,pos)]
#'   return(df)
#' }
#' 
#' plotDMP <- function(DMPann,names,path=NULL){
#'   library(ggplot2)
#'   DMPresults <- data.frame(table(DMPann[ ,c("Contrast","Type")]))
#'   # plot DMPs (hypo/hyper)
#'   g1<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
#'     geom_bar(position="dodge", stat= "identity")+
#'     theme_bw()+
#'     scale_fill_manual(values=c("red", "skyblue")) +
#'     theme(axis.text.x = element_text(angle = 45, hjust=1))+
#'     labs(x = "", y = "count", fill='Methylation')+
#'     ggtitle('Differently methylated probes')
#'   
#'   # ggplot2::ggsave(g1,paste0("analysis/DMPplots/",names,".png"))
#'   
#'   # plot with facets
#'   g2<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
#'     geom_bar(position="dodge", stat= "identity")+
#'     facet_wrap(.~Type, scales = "free_x") +
#'     theme_bw()+
#'     scale_fill_manual(values=c("red", "skyblue")) +
#'     theme(axis.text.x = element_text(angle = 45, hjust=1))+
#'     labs(x = "", y = "count", fill='Methylation')+
#'     ggtitle('Differently methylated probes')
#'   # ggplot2::ggsave("analysis/DNA_methylation/DMP_p0.01_m0.3.facet.png")
#'   
#'   # plot proportion of DMPs in CGI
#'   DMP_annCGI <- data.frame(DMPann[ ,c("Contrast","Type", "Relation_to_Island")])
#'   g3<-ggplot2::ggplot(DMP_annCGI, aes(Contrast, fill = Relation_to_Island)) +
#'     facet_wrap(.~Type, scales = "free_x") +
#'     geom_bar(position ="fill", width = 0.8) +
#'     theme_bw() +
#'     scale_fill_brewer(palette = "Set1") +
#'     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
#'     ylab("DMPs") +
#'     xlab("")
#'   # ggplot2::ggsave("analysis/DNA_methylation/DMP_annCGI.png")
#'   
#'   #table(DMPann[ ,c("Relation_to_Island","Type", "Contrast")])
#'   
#'   # plot proportion of DMPs in genomic elements
#'   DMPann$UCSC_RefGene_Group[which(DMPann$UCSC_RefGene_Group == "")] <- "."
#'   DMPann$UCSC_RefGene_Group_short <- unlist(lapply(strsplit(DMPann$UCSC_RefGene_Group, ";"),'[[', 1))
#'   
#'   DMP_annGenomic <- data.frame(DMPann[ ,c("Contrast","Type", "UCSC_RefGene_Group_short")])
#'   g4<-ggplot2::ggplot(DMP_annGenomic, aes(Contrast, fill = UCSC_RefGene_Group_short)) +
#'     facet_wrap(.~Type, scales = "free_x") +
#'     geom_bar(position = "fill", width = 0.8) +
#'     theme_bw() +
#'     scale_fill_brewer(palette = "Set1") +
#'     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
#'     labs(fill = "RefGene") +
#'     ylab("DMPs") +
#'     xlab("")
#'   # ggplot2::ggsave("analysis/DNA_methylation/DMP_annGenomic.png")
#'   
#'   # table(DMPann[ ,c("UCSC_RefGene_Group_short","Type", "Contrast")])
#'   plt_list<-list(g1,g2,g3,g4)
#'   n <- c("DMP_annCGI.png", "DMP_annGenomic.png", "DMP_count_facet.png", "DMP_count.png")
#'   
#'   if(!is.null(path)){
#'     sapply(1:length(plt_list),function(x){
#'       ggplot2::ggsave(
#'         filename = n[x],
#'         plot = plt_list[[x]],
#'         device = NULL,
#'         path = path,
#'         scale = 1,
#'         width = NA,
#'         height = NA,
#'         units = c("in", "cm", "mm", "px"),
#'         dpi = 300,
#'         limitsize = TRUE,
#'         bg = NULL
#'       )
#'     })
#'   } 
#'   return(plt_list)
#' }
#' 
#' save_plot <-function(object,filename,path){
#'   grDevices::png(file = paste(path,filename,sep="/"),
#'                  width = 480, # The width of the plot in inches
#'                  height = 620) # The height of the plot in inches
#'   object  
#'   dev.off()
#'   
#' } 
#' 
#' find_dmrs<-function(object, model, pcutoff = 0.01, betacutoff = 0.3, min.cpg=5){
#'   contrasts <- colnames(model$contrasts)
#'   future::plan("multisession")
#'   results <- future.apply::future_lapply(contrasts, function(x){
#' 
#'     myAnnotation <- DMRcate::cpg.annotate(
#'       object = object,                    
#'       datatype = "array", 
#'       what = "Beta", 
#'       analysis.type = "differential", 
#'       design = model$design, 
#'       contrasts = TRUE,
#'       cont.matrix = model$contrasts, 
#'       coef = x, 
#'       arraytype = "EPIC"
#'     )
#'     
#'     out <- tryCatch(
#'       {
#'         # Just to highlight: if you want to use more than one 
#'         # R expression in the "try" part then you'll have to 
#'         # use curly brackets.
#'         # 'tryCatch()' will return the last evaluated expression 
#'         # in case the "try" part was completed successfully
#'         
#'         # message("This is the 'try' part")
#'         DMRs <- DMRcate::dmrcate(myAnnotation, 
#'                                  pcutoff = pcutoff,
#'                                  betacutoff= betacutoff,
#'                                  min.cpgs =min.cpg
#'         )
#'      
#'       },
#'       error=function(cond) {
#'         
#'         message(cond)
#'         # Choose a return value in case of error
#'         return(NA)
#'       },
#'       warning=function(cond) {
#'         
#'         message(cond)
#'         # Choose a return value in case of warning
#'         return(NULL)
#'       },
#'       finally={
#'         # NOTE:
#'         # Here goes everything that should be executed at the end,
#'         # regardless of success or error.
#'         # If you want more than one expression to be executed, then you 
#'         # need to wrap them in curly brackets ({...}); otherwise you could
#'         # just have written 'finally=<expression>' 
#'         message(paste("Processed contrast:", x))
#'         message("Next.")
#'       }
#'     )
#'     
#'     
#'     if(!identical(out@no.cpgs, integer(0))){
#'       results.ranges <- DMRcate::extractRanges(out)
#'       results.ranges$Contrast = x
#'       # data.frame(Contrast = x, results.ranges)
#'     }else{
#'       results.ranges<-NULL
#'     }
#'     
#'     return(results.ranges)
#'   }
#'   )
#'   
#'   # return(data.table::rbindlist(results) )
#'   res<-suppressWarnings(do.call("c",results))
#'   return(res)
#' }
#' 
#' #gopath(dmrs_ANA,all.cpg=rownames(betas[,!is.na(ss_clean$ANA_dom)]),n=10,ann=ann)->pat
#' gopath <- function(object,all.cpg=NULL,n=20,ann=NULL){
#'   # object[[1]]
#'   #require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#'   requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
#'   future::plan("multisession")
#'   cont<-unique(object$Contrast)
#'   pathways <- future.apply::future_lapply(cont, function(x){
#'     
#'     results.ranges<- object[object$Contrast == x,]
#'     gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                     all.cpg = all.cpg,
#'                                                     collection = "GO", 
#'                                                     array.type = "EPIC",
#'                                                     anno = ann,
#'                                                     sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO"
#'     g
#'     },
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
#'     }
#'     )
#'     gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "GO", 
#'                                                       array.type = "EPIC",
#'                                                       genomic.features=c("TSS200","TSS1500","1stExon"),
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "GO_prom"
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
#'     }
#'     )
#'     gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
#'                                                       all.cpg = all.cpg,
#'                                                       collection = "KEGG", 
#'                                                       array.type = "EPIC",
#'                                                       anno = ann,
#'                                                       sig.genes = T)
#'     g<-missMethyl::topGSA(gst, n=n)
#'     g$Contrast <- x
#'     g$method <- "KEGG"
#'     g$ONTOLOGY <- NA
#'     g$TERM<-g$Description
#'     g$Description<-NULL
#'     g
#'     },
#'     
#'     error = function(e) {
#'       data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
#'     }
#'     )
#'     result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
#'     return(result)
#'   },future.packages = c("S4Vectors","GenomicRanges"))
#'   
#'   
#'   return(data.table::rbindlist(pathways))
#'   
#' }
#' 
#' # cpgs <- GenomicRanges::GRanges(seqnames = anno$chr, 
#' #                                ranges = IRanges::IRanges(start = anno$pos, 
#' #                                                          end = anno$pos),
#' #                                strand = anno$strand,
#' #                                name = anno$Name)
#' # 
#' # overlaps <- GenomicRanges::findOverlaps(cpgs,regions)
#' # sig.cpg <- cpgs$name[overlaps@from]
#' # .getGO <- function(){
#' #   if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
#' #     stop("org.Hs.eg.db package required but not installed.")
#' #   egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
#' #   GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
#' #   d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
#' #   GeneID.PathID <- GeneID.PathID[d, ]
#' #   GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
#' #                                                       keys=unique(GeneID.PathID$go_id), 
#' #                                                       columns=c("GOID","ONTOLOGY","TERM"), 
#' #                                                       keytype="GOID"))
#' #   go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
#' #   
#' #   list(idList=go, idTable=GOID.TERM)
#' # }
#' # go <- .getGO()
#' # result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
#' #                   array.type=array.type, plot.bias=plot.bias, 
#' #                   prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
#' #                   fract.counts=fract.counts, 
#' #                   genomic.features = genomic.features,
#' #                   sig.genes = sig.genes)
#' # result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
#' # rownames(result) <- result$GOID
#' 
#' # subset(ss,!is.na(ss$ANA_dom1))
#' # technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
#' # 
#' # 
#' #  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
#' #  pheno<-ss[!is.na(ANA_dom1)]
#' #  mod <- model.matrix(
#' #    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
#' #    ,data = pheno
#' #  )
#' #  mod0 <- model.matrix(
#' #    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
#' #    ,data = pheno
#' #  )
#' #  sva.results <- sva(mval, mod, mod0)
#' # # 
#' # design <- model.matrix(
#' #   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
#' #   ,data = ss
#' #   )
#' # n.sv = sva::num.sv(betas,design,method="leek")
#' # 
#' # svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
#' # 
#' 
#' #