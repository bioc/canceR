#' Built Expression Set (eSet) from profile data.
#' @usage geteSet()
#' @export
#' @return ExpressionSet
#' @examples 
#'  f <- 9
#'  \dontrun{
#'  readRDS(paste(path.package("canceR"),"/extdata/rdata/prad_michPhenoTest1021.rds", sep=""))
#'  geteSet()
#'  }
#' 
geteSet <- function(){
    
    #function to replace blanks with missing
    
    blank2na <- function(x){ 
        z <- gsub("\\s+", "", x)  #make sure it's "" and not " " etc
        x[z==""] <- NA 
        return(x)
    } 
        ##Test checked Cases and Genetic Profiles
        testCheckedCaseGenProf()
       
        
        Lchecked_Studies <- ENV$lchecked_Studies_forCases
        Lchecked_Cases <- length(ENV$curselectCases)
        Lchecked_GenProf <- length(ENV$curselectGenProfs)
        
        ProfDataAll=0
        ProfData=0
        LengthGenProfs=0
        LengthCases=0
        for (i in 1:Lchecked_Studies){
            Si = ENV$checked_StudyIndex[i]
            progressBar_ProfilesData <- tkProgressBar(title = ENV$Studies[Si], min = 0,
                                                      max = Lchecked_GenProf, width = 400)
            
            #tkfocus(progressBar_ProfilesData)
            LastLengthGenProfs = LengthGenProfs
            LengthGenProfs = LengthGenProfs + ENV$LGenProfs[i]+1
            LastLengthCases = LengthCases
            LengthCases= LengthCases + ENV$LCases[i]+1
            
            for (k in 1:length(ENV$curselectCases)){
                
                Sys.sleep(0.1)
                setTkProgressBar(progressBar_ProfilesData, k, label=paste( round(k/Lchecked_GenProf*100, 0),
                                                                           "% of Expression Set"))
                
                if (ENV$curselectGenProfs[k] <= LengthGenProfs && ENV$curselectGenProfs[k]>LastLengthGenProfs){    
                    
                    GenProf<- ENV$GenProfsRefStudies[ENV$curselectGenProfs[k]]
                    
                    Case<- ENV$CasesRefStudies[ENV$curselectCases[k]]
                    
                    
                    
                    if(length(ENV$GeneList)>500){
                        ProfData <- getMegaProfData(ENV$GeneList,k )
                    } else{
                        ProfData<-getProfileData(ENV$cgds,ENV$GeneList, GenProf,Case)
                        print(ncol(ProfData))
                    }
                    
                    #ProfData<- getProfileData(ENV$cgds,ENV$GeneList, GenProf,Case)
                    #ProfData <-rbind.na(colnames(ProfData), ProfData)
                    
                    print("getting Profile Data and removing all NAs rows...")
                    ##remove all NAs rows
                    ProfData<- ProfData[which( apply( !( apply(ProfData,1,is.na) ),2,sum)!=0 ),]
                    
                    
                    ## Display AssyData with Tcl Table
                    title <- paste(ENV$StudyRefGenProf[k],":",ENV$CasesStudies[ENV$curselectCases[k]+1])
                    getInTable(ProfData, title)
                    
                    #####nicData_MultipleCases function
                    Case<- ENV$CasesRefStudies[ENV$curselectCases[k]]
                    
                    ClinicalData<-getClinicalData(ENV$cgds,Case)

                    matrix <-rbind.na(colnames(ClinicalData), ClinicalData)
                    rnames <- rownames(ClinicalData)
                    cnames <- colnames(ClinicalData)
                    
                   
                    #apply blank2na function
                    ClinicalData <- data.frame(lapply(ClinicalData,  blank2na))
                    rownames(ClinicalData) <- rnames
                    names(ClinicalData) <- cnames
                    
                     
                    
                   
                    
                    ## getClinicalData generate CHARACTER class if is there "NA" value in any column
                    ## Convert character value to numeric if grep [0-9] != 0
                     for(i in 1:ncol(ClinicalData)){
                     ## substitute "Not Available" by "NA"
                     ClinicalData[,i]<- gsub("\\[Not Available\\]",NA, ClinicalData[,i], ignore.case=TRUE)
                     ClinicalData[,i]<- gsub("NA",NA, ClinicalData[,i], ignore.case=TRUE)
                     if(length(grep("[0-9]*\\.[0-9]*",ClinicalData[,i]))!=0){
                         ClinicalData[,i] <- as.numeric(ClinicalData[,i])
                       }
                     }
                   
                    if(length(ClinicalData[1,])==0){
                        msgNoClinData=paste("No Clinical Data are Available for\n", CasesStudies[curselectCases[k]+1])
                        tkmessageBox(message=msgNoClinData, title= paste("Study: ",ENV$StudyRefCase[k]))
                        close(progressBar_ProfilesData)
                        break
                    } 
                    
                    title <- paste(ENV$StudyRefCase[k],ENV$GenProfChoice[k], sep=": ")
                    getInTable(matrix,title)
                    
                  
                    ## Select only Cases (rownames) that exist in ClinicalDataSub and ProfData
                    merge <- merge(ClinicalData, ProfData, by="row.names")
                    print("merge Clinical and Profile Data")
                    ClinicalData<- merge[,1:(length(ClinicalData)+1)]
                    
                    
                    rownames(ClinicalData)<- ClinicalData[,1]
                    ClinicalData <- ClinicalData[-1]
                    ProfData<-merge[,!(merge %in% ClinicalData)]
                    
                    #row.names(ProfData)<- ProfData[,1]
                    
                    #ProfData <- ProfData[-1]
                    AssayData<- t(ProfData)
                    colnames(AssayData) <- AssayData[1,]
                    AssayData <- AssayData[-1,]
                    rnames <- rownames(AssayData) 
                    AssayData <- as.matrix(apply(AssayData,2 ,function(x) as.numeric(x)))
                    rownames(AssayData) <- rnames
                    
                    
                    ##Convert column with digital values from factor to numeric
                    for(i in 1:ncol(ClinicalData)){
                        ClinicalData[,i] <- sapply(ClinicalData[,i], function(x) if(length(grep("[a-z'-'+A-Z'/'' ']", as.character(ClinicalData[,i])))==0) { as.numeric(as.character(x)) } else {x})
                    }
                    
                

                    ENV$ClinicalData <- ClinicalData
                    ENV$ProfData <- ProfData
                    ENV$AssayData <- AssayData
                    
                    
                    #Test if the same length cases for phenoData and AssayData
                    if (all(rownames(ClinicalData)==colnames(AssayData))){
                        
                        
                        
                        ## create labelDescription for columns of phenoData. 
                        ## labeldescription is used by Biobase packages
                        ## In our case labelDescription is Equal to column names
                        metaData <- data.frame(labelDescription= colnames(ClinicalData), row.names=colnames(ClinicalData))        ## Bioconductor’s Biobase package provides a class called AnnotatedDataFrame   
                        ##that conveniently stores and manipulates 
                        ##the phenotypic data and its metadata in a coordinated fashion. 
                        phenoData<-new("AnnotatedDataFrame", data=ClinicalData, varMetadata=metaData)    
                        
                        ##Assembling an ExpressionSet  
                        ENV$eSet<-Biobase::ExpressionSet(assayData=AssayData, phenoData=phenoData, annotation="GO") 
                        print(paste("End of building eSet..."))
                        
                        #             for (i in 1:length(names(pData(eSet)))){
                        #             pData(eSet)[i] <- as.matrix(na.omit(pData(eSet)[i]))
                        #             }
                    }else {tkmessageBox( message= "The expression Gene Set and the Clinical Data do not have the same samples", icon="warning")}
                } 
            } 
            close(progressBar_ProfilesData)
        } 
        
}