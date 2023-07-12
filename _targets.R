# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint


ncores=RcppParallel::defaultNumThreads()

deps<-c("SummarizedExperiment", "S4Vectors","minfiData", "minfi", "maxprobes", "limma",
        "IRanges", "impute", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "IlluminaHumanMethylation450kanno.ilmn12.hg19", "GenomicRanges", 
        "conumee", "BiocGenerics", "Biobase","RcppParallel","HDF5Array","DMRcatedata")

deps_gh<-c("markgene/maxprobes")  

pkgs<-c("renv","devtools","targets","tarchetypes","tibble","foreach","S4Vectors",
        "minfi","limma","ggplot2","grDevices","graphics","ggfortify","PCAtools","gplots",
        "data.table","future","future.apply","stringi","DMRcate","missMethyl",
        "GenomicRanges"
)
pkgs_gh<-c("ijcBIT/cnv.methyl")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


install_pkgs <- function(pkg) {
  out <- tryCatch(
    {
      renv::install(pkgs)
      
    },
    error=function(cond) {
      BiocManager::install(pkg,Ncpus=ncores)
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      message("Installed package:",pkg)
    }
  )    
  return(out)
}


i_pkgs<-installed.packages()
pkgs_to_install <-unique(c(deps,pkgs)[! (c(deps,pkgs) %in% i_pkgs)])
sapply(pkgs_to_install,function(x) install_pkgs(x))

i_pkgs<-installed.packages()
if (! ("devtools" %in% i_pkgs)) install.packages("devtools")
gh_pkgs_to_install <- unique(c(deps_gh,pkgs_gh)[! (c(deps_gh,pkgs_gh) %in% i_pkgs)])
sapply(gh_pkgs_to_install,function(x) devtools::install_github(x))



# Load packages required to define the pipeline:
library(targets)

library(tarchetypes) # Load other packages as needed. # nolint



.libPaths("~/Projects/marato/renv/library/R-4.2/x86_64-pc-linux-gnu")
# Set target options:
tar_option_set(
  packages = c("tibble","foreach","S4Vectors","minfi"), # packages that your targets need to run
  #imports = "cnv.methyl",
  format = "qs", # default storage format
  storage = "worker", retrieval = "worker"
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore"# multiprocess or multicore, LSF, SGE, Slurm etc.
)

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Load the R scripts with your custom functions:
#lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint
source("R/functions.R")

# produce_data <- function() {
#   expand.grid(samplesheet = c("a", "b"), model = c("c", "d"), normalization = c(1, 2, 3))
# }
# list(
#   tar_group_by(data, produce_data(), samplesheet, model),
#   tar_target(group, data, pattern = map(data))
# )


results_folder = "./results/"
analysis_folder = "./analysis/"

samplesheets <- c(
  ARMS="data/ss_ARMS.rds",
  NBL="data/ss_NBL.rds",
  FULL="data/ss_FULL.rds"
  
)

values <- tibble::tibble( # Use all possible combinations of input settings.
  method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
  data_paths = samplesheets,
  data_names = names(samplesheets),
  group_var = "Sample_Group",
  Contrasts = c(paste(c("Sample_Grouptumour_ARMS-Sample_Groupmetastasis_ARMS", #T1 vs T2 -> mice
                        "Sample_Groupcell_line_ARMS-Sample_Groupprimary_ARMS", #T0 vs T3 -> start end
                        "Sample_Groupcell_line_ARMS-Sample_Grouptumour_ARMS",  #T0 vs T1 -> tumor
                        "Sample_Groupmetastasis_ARMS-Sample_Groupprimary_ARMS",  #T2 vs T3 -> metastasi
                        "Sample_Groupcell_line_ARMS-(Sample_Groupprimary_ARMS+Sample_Grouptumour_ARMS+Sample_Groupmetastasis_ARMS)/3",  #T0 vs T1&T2&T3
                        "(Sample_Groupprimary_ARMS+Sample_Groupcell_line_ARMS)/2-(Sample_Grouptumour_ARMS+Sample_Groupmetastasis_ARMS)/2"  #T0&T1 vs T2&T3
                        ),collapse = ";"
                      ),
                paste(c("Sample_Groupprimary_NBL-Sample_Grouptumour_NBL", #T0 vs T1 -> tumor
                        "Sample_Grouptumour_NBL-Sample_Groupmetastasis_NBL", #T1 vs T2 -> mice
                        "(Sample_Groupprimary_NBL+Sample_Grouptumour_NBL)/2-Sample_Groupmetastasis_NBL",  #T0&T1 vs T2 -> tumor vs metastasi
                        "Sample_Groupprimary_NBL-(Sample_Grouptumour_NBL+Sample_Groupmetastasis_NBL)/2" #T0 vs T1&T2 -> culture vs Tissue
                        ),collapse = ";"
                      ),
                paste(c("Sample_Grouptumour_ARMS-Sample_Grouptumour_NBL", # T1 -> tumors
                       "Sample_Groupmetastasis_ARMS-Sample_Groupmetastasis_NBL", #T2 -> metastasi
                       "Sample_Groupcell_line_ARMS-Sample_Groupprimary_NBL",  # T0 CL
                       "(Sample_Groupcell_line_ARMS+Sample_Groupprimary_NBL)/2-(Sample_Grouptumour_ARMS+Sample_Groupmetastasis_ARMS+Sample_Grouptumour_NBL+Sample_Groupmetastasis_NBL)/4" #T0 vs T1&T2 -> culture vs Tissue

                       ),collapse = ";"
                )
  ),
                
  Contrast_names=list(c("At1-At2:tissue","At0-At3:start/end","At0-At1:tumors","At2-At3:metastasis","At0-All:Cell_line","At0At1 - At2At3:TuMet"),
                   c("Nt0-Nt1:tumor","N1-Nt2:tissue","Nt0+Nt1-Nt2:tumor/metastasi","Nt0-All:culture/tissue"),
                   c("At1-Nt1:tumors","At2-Nt2:metastasi","At0-Nt0:Cultures","culture vs tissue")
                   )
                
                
  
  # model=c("Sample_Group",)
)


top_betas_N_target <- tar_target(top_betas_N,c(100,1000,5000,10000))

targets <- tar_map(
  values = values[1,],
  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(fparams,make_results_dirs(subf=data_names, results_folder = results_folder, analysis_folder = analysis_folder),deployment = "main"),
  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path),deployment = "main"),
  tar_target(ss,samplesheet[-1,],deployment = "main"),
  tar_target(category,samplesheet[1,],deployment = "main"),
  tar_target(nrgSet, minfi::read.metharray.exp(base = NULL, targets = ss, extended = T,
                                               recursive = FALSE, verbose = FALSE, force = T)),
  tar_target(rgSet,name_rgset(nrgSet,ss),deployment = "main"),
  # Qc report:
  tar_target(QC_plots, qc(rgSet,sampGroups = "Sample_Group",sampNames="barcode", qc_folder = fparams[["qc_folder"]])),

  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
  tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="Sample_Group",qc_folder = fparams[["qc_folder"]])),

  tar_target(normalize, method_function(filtered)),
  tar_target(clean, prep(normalize)),
  tar_target(ss_clean, droplevels.data.frame( cbind(clean@colData,purity)),deployment = "main"),
  tar_target(save_ss_clean,write.table(ss_clean,paste0(fparams[["ss_clean_path"]],"/","ss_clean.csv"),quote = F,sep = ",") ),

  tar_target(plotvars,                                                           # Variables to plot on correlation matrix
             c("predictedSex",colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch","covs")])),


  tar_target(ann, minfi::getAnnotation(clean)),
  tar_target(betas, minfi::getBeta(clean)),
  tar_target(top,top_beta(betas,n=top_betas_N),pattern=map(top_betas_N)),                                        # subsample for PCA & ploting
  tar_target(pca, pca_res(top),pattern=map(top)),                                                 # PCA

  tar_target(pca_corrplot,corpca(beta_top100 = top,                              # Correlation plot
                                 metadata=ss_clean,
                                 path=paste0(fparams[["corrplot_folder"]],"/",NROW(top),"/"),
                                 filename=paste0(data_names,"_pca_corrplot",NROW(top),".png"),
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars(top ",NROW(top)," )")
  ),pattern=map(top)
  ),


  tar_target(bplots, bplot(pca,                                                  # Bi plots for PCA components
                           ss=ss_clean,
                           colgroup=plotvars,
                           s="Type",
                           cols=NULL,
                           alfa=0.7,
                           folder = paste0(fparams$bplots_folder,"/",NROW(pca$rotation),"/")),
             packages = c("ggfortify","ggrepel","gplots","ggplot2"),
             pattern=map(pca)
  ),


  tar_target(model, mod(object = betas, group_var = group_var, contrasts = Contrasts,
                        covs= colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch")] ,
                        metadata = ss_clean, rename = Contrast_names)
  ),

  tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,
                                            ContrastsDM = colnames(model$contrasts),
                                            beta_normalized = betas,
                                            p.value = 0.1,
                                            mDiff = 0.1,
                                            ann = ann,
                                            writeOut = F
  )),
  tar_target(dmps_f , filter_dmps(dmps_mod1, p.value = 0.01, mDiff = 0.5)),
  # tar_target(save_dmps, writexl::write_xlsx(dmps_mod1,paste0(fparams$dmp_folder,data_names,"_dmps.xlsx"))),
  tar_target(save_dmps, data.table::fwrite(dmps_f,paste0(fparams$dmp_folder,data_names,"_dmps.txt"))),
  tar_target(betas_DIF,betasdmps(betas,dmps_f,rgSet)),
  tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet)),
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(fparams$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(fparams$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),


  tar_target(dmps_summary,{
    dt<-data.table::as.data.table(dmps_f)
    dt[,list(Hyper=sum(Type=="Hyper"),Hypo=sum(Type=="Hypo")),by=c("Contrast")]
  }),
  # tar_quarto(report,"1.Preprocess.qmd",execute_params = list(output_file=paste0(data_names,"_report.html"),ss = as.character(quote(samplesheet)))),
  tar_target(save_dmps_summary, writexl::write_xlsx(dmps_summary,paste0(fparams$dmp_folder,data_names,"_summary.xlsx"))),
  tar_target(dmpplot_mod1, plotDMP(dmps_f,path=fparams[["dmpplots_folder"]])),

  tar_target(dmrs1,                                                              # Finds DMRs with dmrcate
             find_dmrs(betas,model,
                       fdr = 0.25, p.value = 0.1,betacutoff = 0.1, min.cpg=3)),
  tar_target(dmrs_battery,  priority = 1, error ="continue",                                       # DMRs distribution along params
             apply_filter_dmrs(
               dmrs = dmrs1,path=paste0(fparams$dmrs_folder,data_names))),
  tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.05, mDiff = 0.2, min.cpg=3)),   # Filters DMRs
  tar_target(save_dmrs,                                                          # Saves DMRs
             writexl::write_xlsx(
               dmrs, paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(summary,summarize(dmps=dmps_f,dmrs=dmrs,path=paste0("./results/",data_names,"/summary.csv"))),
  # tar_target(dmrs1, find_dmrs(betas,model,fdr = 0.01, betacutoff = 0.25, min.cpg=5),deployment = "worker"),
  # tar_target(dmrs, filter_dmrs(dmrs1,p.value = 0.01, mDiff = 0.25, min.cpg=5)),
  #
  # tar_target(save_dmrs,
  #            writexl::write_xlsx(as.data.frame(dmrs),
  #                                paste0(fparams$dmrs_folder,"_",data_names,".xlsx"))),
  #
  # tar_target(allpathways, gopath(dmrs,all.cpg=rownames(betas),n=Inf,ann=ann)),   # Pathways with all
  tar_target(hyperpathways,                                                      # Pathways hyper
             gopath(dmrs[meandiff>0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hyper.csv"))),
  tar_target(hypopathways,                                                       # Pathways with hypo
             gopath(dmrs[meandiff<0,],
                    all.cpg=rownames(betas),n=Inf,ann=ann,savepath=paste0(fparams$pathway_folder,"pathways_hypo.csv"))),
  tar_target(results_gopathhypo,
             path_results(pathway=hypopathways,topN=50,group="method",pval=0.05,path=paste0(fparams$pathway_folder,data_names,"_hypoPathways.csv"))),
  tar_target(results_gopathhyper,
             path_results(pathway=hyperpathways,topN=50,group="method",pval=0.05,path=paste0(fparams$pathway_folder,data_names,"_hyperPathways.csv"))),

  tar_target(save_gopathhypo,
             writexl::write_xlsx(hypopathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"hypoPathways.xlsx"))),
  tar_target(save_gopathhyper,
             data.table::fwrite(hyperpathways[FDR<0.05,],
                                 paste0(fparams$pathway_folder,"hyperPathways.csv"))),

  
  NULL
)
#List all targets to run
list(
  top_betas_N_target,
  targets
  )
# # list(
# #   tar_target(samplesheet_path, "data/ss_H358.rds", format = "file"),
# #   tar_target(samplesheet, readRDS(samplesheet_path)),
# #   # tar_target(ss,samplesheet),
# #   tar_target(ss,samplesheet[-1,]),
# #   tar_target(category,samplesheet[1,]),
# #   tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T)),
# #   # Qc report:
# #   tar_target(QC_plots, cnv.methyl::qc(rgSet,sampGroups = "condition")),
# #
# #   # Calculate purity:
# #   tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
# #   tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="condition")),
# #   targets,
# #   # tar_target( alldmps,
# #   #             list(apply(tidyr::expand_grid(c("dmps_ANA","dmps_SANDRA","dmps_WT"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
# #   # )
# #
# #   #tar_combine(combined_gopath_ANA,
# #   #list(apply(tidyr::expand_grid(c("gopath_ANA"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
# #   NULL
# #   )
# #
# # subset(ss,!is.na(ss$ANA_dom1))
# # technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
# #
# #
# #  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
# #  pheno<-ss[!is.na(ANA_dom1)]
# #  mod <- model.matrix(
# #    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
# #    ,data = pheno
# #  )
# #  mod0 <- model.matrix(
# #    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
# #    ,data = pheno
# #  )
# #  sva.results <- sva(mval, mod, mod0)
# # #
# # design <- model.matrix(
# #   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
# #   ,data = ss
# #   )
# # n.sv = sva::num.sv(betas,design,method="leek")
# #
# # svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
# #
# # 
# # # subset ss
# # BiocManager::install(c("Biobase", "conumee", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "IRanges", "limma", "maxprobes", "minfi", "SummarizedExperiment"))
# 

