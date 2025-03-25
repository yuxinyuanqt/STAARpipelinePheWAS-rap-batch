Gene_Centric_Noncoding_PheWAS_preload <- function(chr,gene_name,category=c("all_categories","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                                  genofile,obj_nullmodel_list,
                                                  dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                                  dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                                  dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                                  dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                  rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                                  SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  category <- match.arg(category)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  if(category=="all_categories")
  {
    results <- noncoding_PheWAS_preload(chr,gene_name,genofile,obj_nullmodel_list,
                                        dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                        dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                        dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                        dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                        rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                        rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                        SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="downstream")
  {
    results <- downstream_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,
                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                 rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="upstream")
  {
    results <- upstream_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,
                               rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                               rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                               SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="UTR")
  {
    results <- UTR_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,
                          rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                          rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                          SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="promoter_CAGE")
  {
    results <- promoter_CAGE_PheWAS_preload(chr,gene_name,genofile,obj_nullmodel_list,
                                            dfPromCAGEVarGene.SNV=dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE=variant.id.SNV.PromCAGE,
                                            rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                            rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                            SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="promoter_DHS")
  {
    results <- promoter_DHS_PheWAS_preload(chr,gene_name,genofile,obj_nullmodel_list,
                                           dfPromrOCRsVarGene.SNV=dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs=variant.id.SNV.PromrOCRs,
                                           rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                           rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                           SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="enhancer_CAGE")
  {
    results <- enhancer_CAGE_PheWAS_preload(chr,gene_name,genofile,obj_nullmodel_list,
                                            dfHancerCAGEVarGene.SNV=dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE=variant.id.SNV.HancerCAGE,
                                            rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                            rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                            SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  if(category=="enhancer_DHS")
  {
    results <- enhancer_DHS_PheWAS_preload(chr,gene_name,genofile,obj_nullmodel_list,
                                           dfHancerrOCRsVarGene.SNV=dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs=variant.id.SNV.HancerrOCRs,
                                           rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                           rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                           SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
  }

  return(results)
}


noncoding_PheWAS_preload <- function(chr,gene_name,genofile,obj_nullmodel_list,
                                     dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                     dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                     dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                     dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                     QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                     Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                     Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                     SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
    as.character(x$id_include)
  }))

  #####################################
  #   Gene Info
  ## get SNV id
  filter <- seqGetData(genofile, QC_label)
  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }
  
  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  
  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  
  variant.id <- seqGetData(genofile, "variant.id")
  
  rm(filter)
  gc()
  
  ########################################
  #   Downstream
  
  GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
  is.in <- (GENCODE.Category=="downstream")&(SNVlist)
  variant.id.downstream <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.downstream,sample.id=phenotype.id)
  
  rm(variant.id.downstream)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
  
  rm(GENCODE.Info)
  gc()
  
  rm(variant_gene_num)
  gc()
  
  Gene <- as.character(unlist(GENCODE.Info.split))
  
  rm(GENCODE.Info.split)
  gc()
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_downstream <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_downstream <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- rep(NA,4)
        results_temp[3] <- "downstream"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_downstream <- rbind(results_downstream,results_temp)
      }
      
      if(!is.null(results_downstream))
      {
        if(!use_SPA)
        {
          colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
          colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
          colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_downstream)[dim(results_downstream)[2]] <- c("STAAR-B")
        }
        
        results_list_downstream[[num]] <- results_downstream
      }
    }
  }
  seqResetFilter(genofile)
  
  ########################################
  #   Upstream
  
  is.in <- (GENCODE.Category=="upstream")&(SNVlist)
  variant.id.upstream <- variant.id[is.in]
  
  seqSetFilter(genofile,variant.id=variant.id.upstream,sample.id=phenotype.id)
  
  rm(variant.id.upstream)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[,]")
  variant_gene_num <- sapply(GENCODE.Info.split,function(z) length(z))
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  variant.id.SNV <- rep(variant.id.SNV,variant_gene_num)
  
  rm(GENCODE.Info)
  gc()
  
  rm(variant_gene_num)
  gc()
  
  Gene <- as.character(unlist(GENCODE.Info.split))
  
  rm(GENCODE.Info.split)
  gc()
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_upstream <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_upstream <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- rep(NA,4)
        results_temp[3] <- "upstream"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_upstream <- rbind(results_upstream,results_temp)
      }
      
      if(!is.null(results_upstream))
      {
        if(!use_SPA)
        {
          colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
          colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
          colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_upstream)[dim(results_upstream)[2]] <- c("STAAR-B")
        }
        
        results_list_upstream[[num]] <- results_upstream
      }
    }
  }
  
  seqResetFilter(genofile)
  
  ########################################################
  #                UTR
  
  is.in <- ((GENCODE.Category=="UTR3")|(GENCODE.Category=="UTR5")|(GENCODE.Category=="UTR5;UTR3"))&(SNVlist)
  variant.id.UTR <- variant.id[is.in]
  
  rm(GENCODE.Category)
  gc()
  
  seqSetFilter(genofile,variant.id=variant.id.UTR,sample.id=phenotype.id)
  
  rm(variant.id.UTR)
  gc()
  
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
  
  rm(GENCODE.Info)
  gc()
  
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) z[1]))
  
  rm(GENCODE.Info.split)
  gc()
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  
  seqResetFilter(genofile)
  
  ### Gene
  is.in <- which(Gene==gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_UTR <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_UTR <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- rep(NA,4)
        results_temp[3] <- "UTR"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_UTR <- rbind(results_UTR,results_temp)
      }
      
      if(!is.null(results_UTR))
      {
        if(!use_SPA)
        {
          colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
          colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
          colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_UTR)[dim(results_UTR)[2]] <- c("STAAR-B")
        }
        
        results_list_UTR[[num]] <- results_UTR
      }
    }
  }
  
  seqResetFilter(genofile)

  #############################################
  #   Promoter-CAGE

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_promoter_CAGE <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_promoter_CAGE <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfPromCAGEVarGene.SNV[1, 1:4]
        results_temp[3] <- "promoter_CAGE"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
      }
      
      if(!is.null(results_promoter_CAGE))
      {
        if(!use_SPA)
        {
          colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_CAGE)[dim(results_promoter_CAGE)[2]] <- c("STAAR-B")
        }
        
        results_list_promoter_CAGE[[num]] <- results_promoter_CAGE
      }
    }
  }
  
  seqResetFilter(genofile)
  

  ##################################################
  #       Promoter-DHS

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_promoter_DHS <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_promoter_DHS <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfPromrOCRsVarGene.SNV[1, 1:4]
        results_temp[3] <- "promoter_DHS"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_promoter_DHS <- rbind(results_promoter_DHS,results_temp)
      }
      
      if(!is.null(results_promoter_DHS))
      {
        if(!use_SPA)
        {
          colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_promoter_DHS)[dim(results_promoter_DHS)[2]] <- c("STAAR-B")
        }
        
        results_list_promoter_DHS[[num]] <- results_promoter_DHS
      }
    }
  }
  
  seqResetFilter(genofile)

  ###########################################
  #        Enhancer-CAGE

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_enhancer_CAGE <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_enhancer_CAGE <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfHancerCAGEVarGene.SNV[1, 1:4]
        results_temp[3] <- "enhancer_CAGE"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
      }
      
      if(!is.null(results_enhancer_CAGE))
      {
        if(!use_SPA)
        {
          colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
          colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_enhancer_CAGE)[dim(results_enhancer_CAGE)[2]] <- c("STAAR-B")
        }
        
        results_list_enhancer_CAGE[[num]] <- results_enhancer_CAGE
      }
    }
  }
  
  seqResetFilter(genofile)

  ##################################################
  #       Enhancer-DHS

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  if(length(variant.is.in)>=1)
  {
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    seqResetFilter(genofile)
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate,
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
  } else
  {
    Geno <- NULL
  }
  
  results_list_enhancer_DHS <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results_enhancer_DHS <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfHancerrOCRsVarGene.SNV[1, 1:4]
        results_temp[3] <- "enhancer_DHS"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
      }
      
      if(!is.null(results_enhancer_DHS))
      {
        if(!use_SPA)
        {
          colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
          colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results_enhancer_DHS)[dim(results_enhancer_DHS)[2]] <- c("STAAR-B")
        }
        
        results_list_enhancer_DHS[[num]] <- results_enhancer_DHS
      }
    }
  }
  
  seqResetFilter(genofile)

  ############################################
  #           results
  
  results_noncoding_list <- list(upstream = results_list_upstream,
                                 downstream = results_list_downstream,
                                 UTR = results_list_UTR,
                                 promoter_CAGE = results_list_promoter_CAGE,
                                 promoter_DHS = results_list_promoter_DHS,
                                 enhancer_CAGE = results_list_enhancer_CAGE,
                                 enhancer_DHS = results_list_enhancer_DHS)
  
  seqResetFilter(genofile)
  
  return(results_noncoding_list)

}


promoter_CAGE_PheWAS_preload <- function(chr,gene_name,genofile,obj_nullmodel_list,
                                         dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                         rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                         QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                         Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                         Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                         SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
    as.character(x$id_include)
  }))

  ### Gene
  is.in <- which(dfPromCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromCAGE[is.in]
  
  if(length(variant.is.in)==0)
  {
    results_list <- rep(list(c()), length(obj_nullmodel_list))
    seqResetFilter(genofile)
    return(results_list)
  }

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  ## get AF, Missing rate
  AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
  REF_AF <- AF_AC_Missing$af
  Missing_rate <- AF_AC_Missing$miss
  variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
  seqResetFilter(genofile)
  
  Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                        sample.id=phenotype.id,
                                        REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                        Missing_rate=Missing_rate,
                                        rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        annotation_phred=Anno.Int.PHRED.sub)
  Geno <- Genotype_sp$Geno
  Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
  rm(Genotype_sp)
  gc()
  
  results_list <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfPromCAGEVarGene.SNV[1,1:4]
        results_temp[3] <- "promoter_CAGE"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results <- rbind(results,results_temp)
      }
      
      if(!is.null(results))
      {
        if(!use_SPA)
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[dim(results)[2]] <- c("STAAR-B")
        }
        
        results_list[[num]] <- results
      }
    }
  }
  
  seqResetFilter(genofile)
  
  return(results_list)
}


promoter_DHS_PheWAS_preload <- function(chr,gene_name,genofile,obj_nullmodel_list,
                                        dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                        rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                        QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                        Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                        Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                        SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
    as.character(x$id_include)
  }))

  ### Gene
  is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.PromrOCRs[is.in]
  
  if(length(variant.is.in)==0)
  {
    results_list <- rep(list(c()), length(obj_nullmodel_list))
    seqResetFilter(genofile)
    return(results_list)
  }

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  ## get AF, Missing rate
  AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
  REF_AF <- AF_AC_Missing$af
  Missing_rate <- AF_AC_Missing$miss
  variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
  seqResetFilter(genofile)
  
  Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                        sample.id=phenotype.id,
                                        REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                        Missing_rate=Missing_rate,
                                        rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        annotation_phred=Anno.Int.PHRED.sub)
  Geno <- Genotype_sp$Geno
  Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
  rm(Genotype_sp)
  gc()
  
  results_list <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
        results_temp[3] <- "promoter_DHS"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results <- rbind(results,results_temp)
      }
      
      
      if(!is.null(results))
      {
        if(!use_SPA)
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[dim(results)[2]] <- c("STAAR-B")
        }
        
        results_list[[num]] <- results
      }
    }
  }
  
  seqResetFilter(genofile)
  
  return(results_list)
}


enhancer_CAGE_PheWAS_preload <- function(chr,gene_name,genofile,obj_nullmodel_list,
                                         dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                         rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                         QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                         Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                         Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                         SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
    as.character(x$id_include)
  }))

  ### Gene
  is.in <- which(dfHancerCAGEVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerCAGE[is.in]
  
  if(length(variant.is.in)==0)
  {
    results_list <- rep(list(c()), length(obj_nullmodel_list))
    seqResetFilter(genofile)
    return(results_list)
  }

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  ## get AF, Missing rate
  AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
  REF_AF <- AF_AC_Missing$af
  Missing_rate <- AF_AC_Missing$miss
  variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
  seqResetFilter(genofile)
  
  Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                        sample.id=phenotype.id,
                                        REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                        Missing_rate=Missing_rate,
                                        rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        annotation_phred=Anno.Int.PHRED.sub)
  Geno <- Genotype_sp$Geno
  Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
  rm(Genotype_sp)
  gc()
  
  results_list <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfHancerCAGEVarGene.SNV[1,1:4]
        results_temp[3] <- "enhancer_CAGE"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results <- rbind(results,results_temp)
      }
      
      if(!is.null(results))
      {
        if(!use_SPA)
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[dim(results)[2]] <- c("STAAR-B")
        }
        
        results_list[[num]] <- results
      }
    }
  }
  
  seqResetFilter(genofile)
  
  return(results_list)
}


enhancer_DHS_PheWAS_preload <- function(chr,gene_name,genofile,obj_nullmodel_list,
                                        dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                        rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                        QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                        Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                        Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                        SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

  ## evaluate choices
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
    as.character(x$id_include)
  }))

  ### Gene
  is.in <- which(dfHancerrOCRsVarGene.SNV[,5]==gene_name)
  variant.is.in <- variant.id.SNV.HancerrOCRs[is.in]
  
  if(length(variant.is.in)==0)
  {
    results_list <- rep(list(c()), length(obj_nullmodel_list))
    seqResetFilter(genofile)
    return(results_list)
  }

  seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if(variant_type=="SNV")
  {
    if(Use_annotation_weights)
    {
      for(k in 1:length(Annotation_name))
      {
        if(Annotation_name[k]%in%Annotation_name_catalog$name)
        {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
          
          if(Annotation_name[k]=="CADD")
          {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          
          if(Annotation_name[k]=="aPC.LocalDiversity")
          {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  ## get AF, Missing rate
  AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
  REF_AF <- AF_AC_Missing$af
  Missing_rate <- AF_AC_Missing$miss
  variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
  seqResetFilter(genofile)
  
  Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                        sample.id=phenotype.id,
                                        REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                        Missing_rate=Missing_rate,
                                        rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        annotation_phred=Anno.Int.PHRED.sub)
  Geno <- Genotype_sp$Geno
  Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
  rm(Genotype_sp)
  gc()
  
  results_list <- rep(list(c()), length(obj_nullmodel_list))
  phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
  
  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
  {
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      ## Extract the genotype matrix and variant information of the current trait
      Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
      MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
      Missing_num.in <- Missing_num.sp(Geno.num)
      MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
      
      if (geno_missing_imputation == "mean")
      {
        Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
      }
      if (geno_missing_imputation == "minor")
      {
        Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
        MAF.in <- MAC.in/(2*samplesize.num)
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }else
        {
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      results <- c()
      if(inherits(pvalues, "list"))
      {
        results_temp <- dfHancerrOCRsVarGene.SNV[1,1:4]
        results_temp[3] <- "enhancer_DHS"
        results_temp[2] <- chr
        results_temp[1] <- as.character(gene_name)
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results <- rbind(results,results_temp)
      }
      
      if(!is.null(results))
      {
        if(!use_SPA)
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[dim(results)[2]] <- c("STAAR-B")
        }
        
        results_list[[num]] <- results
      }
    }
  }
  
  seqResetFilter(genofile)
  
  return(results_list)
}

