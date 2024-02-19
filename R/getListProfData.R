#' get a list of Profile Data of every available dimensions. This function load matrices of every dimension (Exp, CNA, Met, RPPA,miRNA,Mut) and save them in a list for every disease.
#' @usage getListProfData()
#' @return a list of data frame with Profiles Data
#' @export
#' @examples
#' readRDS(paste(path.package("canceR"),"/extdata/rdata/brca_tcga73genes.rds", sep=""))
#' \dontrun{
#' getListProfData()
#' head(ENV$ProfData$Expression)
#' }
getListProfData <- function(){
    ## Testing The existing  Cases and GenProf cgdsr references and getProfileData
    grepRef<-function(regex1, listRef1,regex2, listRef2, GeneList,Mut){
        if(length(grep(regex1,listRef1)) != 0){
            if(length(grep(regex2,listRef2))!= 0){
                if(Mut== 0){
                    print(paste("Getting Profile Data of ",regex2,"...",sep=""))
                    ProfData_X <- getProfileData(ENV$cgds,GeneList, regex2,regex1)
                    ########################
                    #ProfData_X <- ProfData_X[,as.factor(ENV$GeneList)]
                    #ProfData_X <- ProfData_X[ENV$GeneList,,drop=FALSE]


#                     > ENV$GeneList[c(179,306,338,400)]
#                     [1] "HBXIP"     "C16orf5"   "SELS"      "MAGI2-IT1"
#                     > colnames(ProfData)[c(82,198,205,404)]
#                     [1] "CDIP1"     "LAMTOR5"   "MAGI2.IT1" "VIMP"
#                     >
                    ####################
                    if(length(ProfData_X)== 0){
                        if(length(ENV$GeneList) <= 500){
                            ## built empty data frame with gene Symbol in colnames
                            ProfData_X <- as.data.frame(setNames(replicate(length(ENV$GeneList),numeric(1), simplify = FALSE), ENV$GeneList[order(ENV$GeneList)]))
                            return(ProfData_X)

                        }else{
                            ProfData_X <- as.data.frame(setNames(replicate(length(SubMegaGeneList),numeric(1), simplify = FALSE), SubMegaGeneList[order(SubMegaGeneList)]))
                            return(ProfData_X)
                        }
                    }else{
                        return(ProfData_X)
                    }

                }else if(Mut==1){
                    print(paste("Getting Mutation Data of ",ENV$checked_Studies[i],"...",sep=""))
                    MutData <- getMutationData(ENV$cgds,regex1, regex2, GeneList)
                    print(paste("MutData: ",dim(MutData)))
                    if(length(MutData)==0){
                        ## built emty data.frame as the same form of MutData
                        MutData <- data.frame("gene_symbol"=character(1),"mutation_type"=character(1), "amino_acid_change"=character(1))
                        return(MutData)
                    }else{
                        ## From Mut Data frame select only Gene_symbol, Mutation_Type, AA-Changes
                        #ENV$ListMutData[[StudiesRef[Si]]] <- MutData[,c(2,6,8)]# Gene Symbol, Mut type, AA change
                        return(MutData)
                    }
                    #return(MutData)

                }
            }else{tkmessageBox(message = paste("There is no genetic Profiles: ", regex2," for Study:",ENV$checked_Studies[i] ), icon="warning" )
                  print(paste("There is no genetic Profile: ", regex2," for Study:",ENV$checked_Studies[i],"..." ))
                  ## built empty data frame with gene Symbol in colnames
                  if(length(ENV$GeneList) <500){
                       ProfData_X <- as.data.frame(setNames(replicate(length(ENV$GeneList),numeric(1), simplify = FALSE), ENV$GeneList[order(ENV$GeneList)]))
                      return(ProfData_X)
                  }else{
                      ProfData_X <- as.data.frame(setNames(replicate(length(SubMegaGeneList),numeric(1), simplify = FALSE), SubMegaGeneList[order(SubMegaGeneList)]))
                      return(ProfData_X)
                  }
                  return( ProfData_X)
            }
        }else{
            tkmessageBox(message = paste("There is no Cases: ",regex1, " for Study:",ENV$checked_Studies[i] ), icon="warning" )
            print(paste("There is no Cases: ", regex1," for Study:",ENV$checked_Studies[i],"..." ))
            ## built empty data frame with gene Symbol in colnames
            if(length(ENV$GeneList) <500){
                ProfData_X <-as.data.frame(setNames(replicate(length(ENV$GeneList),numeric(1), simplify = FALSE), ENV$GeneList[order(ENV$GeneList)]))
                return(ProfData_X)
            }else{
                ProfData_X <-as.data.frame(setNames(replicate(length(SubMegaGeneList),numeric(1), simplify = FALSE), SubMegaGeneList[order(SubMegaGeneList)]))
                return(ProfData_X)
            }
            return(ProfData_X)


        }
    }

if(exists("ListProfData", envir = ENV)){
    rm(ListProfData, envir=ENV)

}
if(exists("ListMetData", envir = ENV)){
    rm(ListMetData, envir=ENV)

}
if(exists("ListMutData", envir = ENV)){
    rm(ListMutData, envir=ENV)

}

    Lchecked_Studies <- ENV$lchecked_Studies_forCases

    #get Study references
    StudiesRef <- getCancerStudies(ENV$cgds)[,1]


    LengthGenProfs <- 0
    LengthCases <- 0

    for (i in 1: Lchecked_Studies){
        Si = ENV$checked_StudyIndex[i]
        progressBar_ProfilesData <- tkProgressBar(title = ENV$Studies[Si], min = 0,
                                                  max = Lchecked_Studies, width = 400)

        Sys.sleep(0.1)
        setTkProgressBar(progressBar_ProfilesData, i, label=paste( round(i/Lchecked_Studies*100, 0),
                                                                   "% of Profiles Data"))


        ### get Cases and Genetic Profiles  with cgdsr references
        GenProf_CNA<- paste(ENV$checked_Studies[i],"_gistic", sep="")
        Case_CNA   <- paste(ENV$checked_Studies[i],"_cna", sep="")

        GenProf_Exp<- paste(ENV$checked_Studies[i],"_rna_seq_v2_mrna", sep="")
        Case_Exp   <- paste(ENV$checked_Studies[i],"_rna_seq_v2_mrna", sep="")

        GenProf_Met_HM450<- paste(ENV$checked_Studies[i],"_methylation_hm450", sep="")
        Case_Met_HM450   <- paste(ENV$checked_Studies[i],"_methylation_hm450", sep="")

        GenProf_Met_HM27<- paste(ENV$checked_Studies[i],"_methylation_hm27", sep="")
        Case_Met_HM27   <- paste(ENV$checked_Studies[i],"_methylation_hm27", sep="")

        GenProf_RPPA<- paste(ENV$checked_Studies[i],"_RPPA_protein_level", sep="")
        Case_RPPA   <- paste(ENV$checked_Studies[i],"_rppa", sep="")

        GenProf_miRNA<- paste(ENV$checked_Studies[i],"_mirna", sep="")
        Case_miRNA   <- paste(ENV$checked_Studies[i],"_microrna", sep="")

        GenProf_Mut<- paste(ENV$checked_Studies[i],"_mutations", sep="")
        Case_Mut   <- paste(ENV$checked_Studies[i],"_sequenced", sep="")


        ## Subsettint of Gene List if bigger than 500
        if(length(ENV$GeneList)>500){
            MegaGeneList <- ENV$GeneList
            if(is.integer(length(MegaGeneList)/500)){
                G <- lenght(MegaGeneList)/500
            }else{
                G <- as.integer(length(MegaGeneList)/500) + 1
            }

            MegaProfData_CNA <- 0
            MegaProfData_Exp <- 0
            MegaProfData_Met_HM450 <- 0
            MegaProfData_Met_HM27 <- 0
            MegaProfData_RPPA <- 0
            MegaProfData_miRNA <- 0
            MegaMutData <- 0
            SubMegaGeneList <- 0
            LastSubMegaGeneList <- 0
            for(g in 1: G){

                if (length(MegaGeneList) - LastSubMegaGeneList > 500){
                    SubMegaGeneList <- MegaGeneList[((g-1)*(500)+1):((g)*500)]
                    LastSubMegaGeneList <- LastSubMegaGeneList + length(SubMegaGeneList)
                    print(paste("SubMega",g,length(SubMegaGeneList), sep=":"))
                } else{
                    print(paste("SubMegaGeneList <-MegaGeneList[",LastSubMegaGeneList, ":", length(MegaGeneList),"]"))
                    SubMegaGeneList <- MegaGeneList[LastSubMegaGeneList:length(MegaGeneList)]
                    print(paste("SubMega else",g,length(SubMegaGeneList),sep=":"))
                }

                print("*************************************************")
                print(paste("Getting Profile Data of Genes from: ", (((g-1)*500)+1), "to",((g-1)*500)+length(SubMegaGeneList),"of", ENV$checked_Studies[i],sep= " "))
                print("*************************************************")

                ProfData_CNA<-grepRef(Case_CNA, ENV$CasesRefStudies, GenProf_CNA, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                ProfData_Exp<-grepRef(Case_Exp, ENV$CasesRefStudies, GenProf_Exp, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                ProfData_Met_HM450<-grepRef(Case_Met_HM450, ENV$CasesRefStudies, GenProf_Met_HM450, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                print(paste("SubMega:",length(SubMegaGeneList)))
                ProfData_Met_HM27<-grepRef(Case_Met_HM27, ENV$CasesRefStudies, GenProf_Met_HM27, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                ProfData_RPPA<-grepRef(Case_RPPA, ENV$CasesRefStudies, GenProf_RPPA, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                ProfData_miRNA<-grepRef(Case_miRNA, ENV$CasesRefStudies, GenProf_miRNA, ENV$GenProfsRefStudies, SubMegaGeneList, Mut=0)
                MutData <- grepRef(Case_Mut,ENV$CasesRefStudies ,GenProf_Mut, ENV$GenProfsRefStudies,SubMegaGeneList, Mut=1)

                print("1")
                MegaProfData_CNA <- cbind(MegaProfData_CNA, ProfData_CNA)
                print("2")
                MegaProfData_Exp <- cbind(MegaProfData_Exp, ProfData_Exp)
                print("3")
                MegaProfData_Met_HM450 <- cbind(MegaProfData_Met_HM450, ProfData_Met_HM450)
                print("4")
                MegaProfData_Met_HM27 <- cbind(MegaProfData_Met_HM27, ProfData_Met_HM27)
                print("5")
                MegaProfData_RPPA <- cbind(MegaProfData_RPPA, ProfData_RPPA)
                print("6")
                MegaProfData_miRNA <- cbind(MegaProfData_miRNA, ProfData_miRNA)
                MegaMutData <- rbind(MegaMutData, MutData)

            }
            ProfData_CNA <- MegaProfData_CNA[,-1]
            ProfData_Exp <- MegaProfData_Exp[,-1]
            ProfData_Met_HM450 <- MegaProfData_Met_HM450[,-1]
            ProfData_Met_HM27 <- MegaProfData_Met_HM27[,-1]
            ProfData_RPPA <- MegaProfData_RPPA[,-1]
            ProfData_miRNA <- MegaProfData_miRNA[,-1]
            MutData <- MegaMutData[-1,]


        } else if (length(ENV$GeneList) > 0){


            ProfData_CNA<- grepRef(Case_CNA, ENV$CasesRefStudies, GenProf_CNA, ENV$GenProfsRefStudies, ENV$GeneList, Mut=0)
            ProfData_Exp<- grepRef(Case_Exp, ENV$CasesRefStudies, GenProf_Exp, ENV$GenProfsRefStudies, ENV$GeneList, Mut=0)
            ProfData_Met_HM450 <- grepRef(Case_Met_HM450, ENV$CasesRefStudies, GenProf_Met_HM450, ENV$GenProfsRefStudies, ENV$GeneList, Mut=0)
            ProfData_Met_HM27 <- grepRef(Case_Met_HM27, ENV$CasesRefStudies, GenProf_Met_HM27, ENV$GenProfsRefStudies, ENV$GeneList,Mut=0)
            ProfData_RPPA<- grepRef(Case_RPPA, ENV$CasesRefStudies, GenProf_RPPA, ENV$GenProfsRefStudies, ENV$GeneList,Mut=0)
            ProfData_miRNA<- grepRef(Case_miRNA, ENV$CasesRefStudies, GenProf_miRNA, ENV$GenProfsRefStudies, ENV$GeneList,Mut=0)
            MutData <- grepRef(Case_Mut,ENV$CasesRefStudies ,GenProf_Mut, ENV$GenProfsRefStudies,ENV$GeneList, Mut=1)

        } else {
            tkmessageBox(message= "Load gene List", icon="warning")
            close(progressBar_ProfilesData)
            stop("Load Gene List")
        }


        print("Saving... ")

        ENV$ListProfData$CNA[[StudiesRef[Si]]] <- ProfData_CNA
        ENV$ListProfData$Expression[[StudiesRef[Si]]] <- ProfData_Exp
        ENV$ListMetData$HM450[[StudiesRef[Si]]] <- ProfData_Met_HM450
        ENV$ListMetData$HM27[[StudiesRef[Si]]] <- ProfData_Met_HM27
        ENV$ListProfData$RPPA[[StudiesRef[Si]]] <- ProfData_RPPA
        ENV$ListProfData$miRNA[[StudiesRef[Si]]] <- ProfData_miRNA
        ENV$ListMutData[[StudiesRef[Si]]] <- MutData

        print(" End Getting Profiles Data... ")
        close(progressBar_ProfilesData)
    }
#     print("Start Ordering ...")
    ## range matrices by the same order
#     ENV$ListProfData$CNA <- ENV$ListProfData$CNA[order(names(ENV$ListProfData$CNA))]
#     ENV$ListProfData$Expression <- ENV$ListProfData$Expression[order(names(ENV$ListProfData$Expression))]
#     ENV$ListMetData$HM450 <- ENV$ListMetData$HM450[order(names(ENV$ListMetData$HM450))]
#     ENV$ListMetData$HM27 <- ENV$ListMetData$HM27[order(names(ENV$ListMetData$HM27))]
#     ENV$ListProfData$RPPA <- ENV$ListProfData$RPPA[order(names(ENV$ListProfData$RPPA))]
#     ENV$ListMutData <- ENV$ListMutData[order(names(ENV$ListMutData))]

#    print("End Ordering ...")

    ## get Gene Mutation Frequency
    UnifyRowNames <- function(x){
        df_MutData <-as.data.frame(table(x$gene_symbol))
        rownames(df_MutData) <- df_MutData$Var1
        ## ordering genes in MutData as in GeneList

        df_GeneList <- as.data.frame(t(ENV$GeneList))
        #df_GeneList <- as.data.frame(ENV$GeneList)
        rownames(df_GeneList) <- df_GeneList[,1]
        df_merge <- merge(df_GeneList, df_MutData, by="row.names",all.x=TRUE)
        Freq_Mut <- df_merge[,c(-2,-3)]
        return(Freq_Mut)
    }
######
    print("Start getting Frequency of Mutation ...")
    #Freq_ListMutData <- plyr::laply(ENV$ListMutData,function(x) UnifyRowNames(x))
    Freq_ListMutData <- lapply(ENV$ListMutData,
                               function(x) UnifyRowNames(x))
    Freq_ArrayMutData <- array(unlist( Freq_ListMutData),
                               dim = c(nrow(Freq_ListMutData[[1]]),
                                       ncol( Freq_ListMutData[[1]]),
                                       length(Freq_ListMutData)))
    
    if (inherits(try(dimnames(Freq_ArrayMutData) <-
                     list(Freq_ListMutData[[1]][,1],
                          colnames(Freq_ListMutData[[1]]),
                          names(Freq_ListMutData)),
                     silent=TRUE),"try-error")){
       # p("There is a Study without Mutation Data.
        #  Use Mutation Panel to verify mutations data for selected studies.",
         # align="center", style = "color:blue")
    }else{
        dimnames(Freq_ArrayMutData) <-
            list(Freq_ListMutData[[1]][,1],
                 colnames(Freq_ListMutData[[1]]),
                 names(Freq_ListMutData))
    }
    #   ?getListProfData(Genes= empty)
    if(dim(Freq_ArrayMutData)[3]==1){
        Freq_DfMutData <- as.numeric(Freq_ArrayMutData[,2,])
        names(Freq_DfMutData) <- names(Freq_ArrayMutData[,2,])
        ## ordering gene list as in GeneList from MSigDB:
        ## grouping genes with the same biological process or gene Sets
        Freq_DfMutData <- Freq_DfMutData[ENV$GeneList]
        Freq_DfMutData <- data.frame(round(Freq_DfMutData,digits=2))
        names(Freq_DfMutData) <- names(Freq_ListMutData)
    }else{
        Freq_DfMutData <- apply(Freq_ArrayMutData[,2,],2,as.numeric)
        rownames(Freq_DfMutData) <- rownames(Freq_ArrayMutData[,2,])
        ## ordering gene list as in GeneList from MSigDB:
        ## grouping genes with the same biological process or gene Sets
        Freq_DfMutData <- Freq_DfMutData[ENV$GeneList,,drop=FALSE]
        Freq_DfMutData <- data.frame(round(Freq_DfMutData,digits=2))
    }
    
    
    
    
    ENV$Freq_DfMutData <- Freq_DfMutData
    
    print("End getting Mutation Frequency...")
######
 
    # 
    # #output1 <- adply(Freq_ListMutData,1)
    # 
    # ## convert the list of correlation matrices to Array
    # Freq_ArrayMutData <- array(unlist( Freq_ListMutData), dim = c(nrow(Freq_ListMutData[[1]]), ncol( Freq_ListMutData[[1]]), length(Freq_ListMutData)))
    # dimnames(Freq_ArrayMutData) <- list(Freq_ListMutData[[1]][,1], colnames(Freq_ListMutData[[1]]), names(Freq_ListMutData))
    # 
    # 
    # Freq_DfMutData <- apply(Freq_ArrayMutData[,2,],2,as.numeric)
    # rownames(Freq_DfMutData) <- rownames(Freq_ArrayMutData[,2,])
    # ## ordering gene list as in GeneList from MSigDB: grouping genes with the same biological process or gene Sets
    # Freq_DfMutData <- Freq_DfMutData[ENV$GeneList,,drop=FALSE]
    # 
    # ENV$Freq_DfMutData <- Freq_DfMutData
    # 
    # print("End getting Mutation Frequency...")
    # 

}
