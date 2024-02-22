#' get Methylation data for multiple genes
#' @usage
#' getMetDataMultipleGenes()
#'
#' @return a a dataframe with mean and median of methylation rate (threshold of silencing gene)
#' @export
#'
#' @examples
#' readRDS(paste(path.package("canceR"),"/extdata/rdata/ucec_tcga_pubGSEA1021.rds", sep=""))
#' \dontrun{
#' getMetDataMultipleGenes()
#' }
getMetDataMultipleGenes <-function(){
    
    tclRequire("BWidget")
    tclRequire("Tktable")
    
    testCheckedCaseGenProf()
    
    
    Lchecked_Studies <- length(ENV$checked_Studies_id)
    Lchecked_Cases <- length(ENV$curselectCases)
    Lchecked_GenProf <- length(ENV$curselectGenProfs)
    ######### Test if all cases were corresponded to appropriate Gen profs (same Study)
    for (i in 1:Lchecked_Cases){
        if(ENV$StudyRefCase[i]!=ENV$StudyRefGenProf[i]){
            
            msgBadChoice="Correspond the Genetic Profile to the Case for the same Study"
            tkmessageBox(message=msgBadChoice, icon="warning")
            tkfocus(ENV$ttCasesGenProfs)
            stop("Correspond the Genetic Profile to the Case for the same Study")
            
        }
    }
    
    ProfDataAll=0
    ProfData=0
    LengthGenProfs=0
    LengthCases=0
    for (i in 1:Lchecked_Studies){
        
        Si =ENV$checked_StudyIndex[i]
        progressBar_ProfilesData <- tkProgressBar(title = ENV$Studies$name[Si], min = 0,
                                                  max = Lchecked_GenProf, width = 400)
        
        #tkfocus(progressBar_ProfilesData)
        LastLengthGenProfs = LengthGenProfs
        LengthGenProfs = LengthGenProfs + ENV$n_GenProfs[i]+1
        LastLengthCases = LengthCases
        LengthCases = LengthCases + ENV$n_Cases[i]+1
        
        for (k in 1:Lchecked_GenProf){
            Sys.sleep(0.1)
            setTkProgressBar(progressBar_ProfilesData, k, label=paste( round(k/Lchecked_GenProf*100, 0),
                                                                       "% of Methylation Data"))
            
            
            if(ENV$curselectGenProfs[k] <= LengthGenProfs && 
               ENV$curselectGenProfs[k] > LastLengthGenProfs){    
                
                GenProf <- ENV$GenProfsRefStudies[ENV$curselectGenProfs[k]]
                
                if (length(grep("methylation", GenProf))==0){
                    msgNoMeth <- "Select Methylation data from Genetics Profiles"
                    tkmessageBox(message = msgNoMeth, icon='info')
                    break
                }
                
                Study_id <- ENV$CasesRefStudies[ENV$curselectCases[k]]
                
                # if (length(grep("methylation", Study_id))==0){
                #     msgNoMeth <- "Select Methylation data from Cases"
                #     tkmessageBox(message = msgNoMeth, icon='info')
                #     break
                # }
                
                
                #ProfData <- getProfileData(ENV$cgds,ENV$GeneList, GenProf,Study_id)
                ProfData <-  cBioPortalData::getDataByGenes(api = ENV$cgds,
                                               studyId = Study_id,
                                               genes = ENV$GeneList,
                                               by = "hugoGeneSymbol",
                                               molecularProfileIds = GenProf) |>
                    unname() |>
                    as.data.frame()|>
                    select("sampleId", "hugoGeneSymbol", "value") |>
                    tidyr::spread("hugoGeneSymbol", "value") |>
                   # `rownames<-`("sampleId")
                   tibble::column_to_rownames("sampleId")
                
                ##convert data frame to numeric structure
                #if( !is.numeric(ProfData[1,1])){
                #    for(i in 1:ncol(ProfData)){
                    
                 #       ProfData[,i] <- as.numeric(ProfData[,i])
                  #  }
                    
                #}
                ##More facter
                cidx <- !(sapply(ProfData, is.numeric) )
                ProfData[cidx] <- lapply(ProfData[cidx], as.numeric)
                
                ProfData <- round(ProfData, digits=3)
                
                dialogMetOption(ProfData,k)
                
            }
        } 
        close(progressBar_ProfilesData)
    }
    
}