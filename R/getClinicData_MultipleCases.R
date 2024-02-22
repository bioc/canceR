#' get Clinical Data for Multiple Cases. User needs to select at least one case to run this function. Get clinical data for more one or multiple cases.
#' @usage 
#' getClinicData_MultipleCases(getSummaryGSEAExists)
#' 
#' @param getSummaryGSEAExists  if equal to 0, the clinical data is displayed in table. if the argument is equal to 1, the clinical data is used to summarise GSEA analysis results.
#' 
#' @export
#' @return dataframe with clinical data
#' 
#' @examples 
#' ##Load Session
#' readRDS(paste(path.package("canceR"),"/extdata/rdata/brca_tcga73genes.rds", sep=""))
#' ## Select Case
#' ENV <- new.env(parent = emptyenv())
#' ENV$curselectCases <- 2
#' ## get Clinical data
#' \dontrun{
#' getClinicData_MultipleCases(getSummaryGSEAExists = 0)
#' }
#' 
getClinicData_MultipleCases<- function(getSummaryGSEAExists){
    
    ##getSummaryGSEAExists is an argument for the function getSummaryGSEA(). this function accept only one clinical dats.
    if(getSummaryGSEAExists ==1 && length(ENV$curselectCases)>1){
        
        tkmessageBox(message=paste("Multiple Cases are loaded. To get Summary of GSEA Results, Select only the Clinical data of study on which are you working."), icon="warning")
        stop("Multiple Cases are loaded. To get Summary of GSEA Results, Select only the Clinical data of study on which are you working.")
    }
    
    
    ClinicalData<-0
    ClinicalData_All<-NULL
    ClinicalDataSub <-NULL
    ClinicalDataSub_All <- NULL
    
    if(length(grep("curselectCases",ls(ENV)))==0){
        msgNoCases <- "Select at less one Case"
        tkmessageBox(message=msgNoCases)
        stop(msgNoCases)
    }
    
    for(c in 1:length(ENV$curselectCases)){
        Study_id <- ENV$CasesRefStudies[ENV$curselectCases[c]]
        
        #study_id <- gsub('_[a-z]*[0-9]*$', '', Study_id)
        ## extract only abbreviation of studies
        # unique(
        #     grep(paste(
        #         unique(unlist(lapply(selected_by_curser, 
        #                              function(x) unlist(strsplit(x,split = "_"))[1]))),
        #         collapse = "|"), checked_Studies, value=TRUE)
        # )
        
        ClinicalData <- cBioPortalData::clinicalData(ENV$cgds, Study_id) 

        
        if(length(ClinicalData[1,])==0){
            msgNoClinData=paste("No Clinical Data are Available for\n", ENV$CasesStudies[ENV$curselectCases[c]+1])
            tkmessageBox(message=msgNoClinData, title= paste("Study: ",ENV$StudyRefCase[c]))
        } else{
            
            ttClin<-tktoplevel()
            tktitle(ttClin) <- paste("Clinical Data of", ENV$CaseChoice[c], sep=" ")
            tkwm.geometry(ttClin, "430x420")
            yscr1 <- tkscrollbar(ttClin, repeatinterval=2,
                                 command=function(...)tkyview(ttc,...))
            xscr1 <- tkscrollbar(ttClin, repeatinterval=2,orient="horizontal",
                                 command=function(...)tkxview(ttc,...))
            ttc<-tklistbox(ttClin,height=20, width= 50 ,selectmode="multiple",
                           xscrollcommand=function(...)tkset(xscr1,...),
                           yscrollcommand=function(...)tkset(yscr1,...),background="white")
            
            tkgrid(ttc,yscr1, columnspan=1)
            tkgrid.configure(yscr1,rowspan=20, columnspan=1,sticky="nsw")
            tkgrid(xscr1,columnspan=2)
            tkgrid.configure(xscr1,rowspan=2,columnspan=2,sticky="ew")
            
            tkinsert(ttc,"end","All")
            cbIValue=0
            for(j in 1: length(names(ClinicalData))){
                tkinsert(ttc,"end",names(ClinicalData)[j])
                
            }
            
            OnOK <- function(){
                ENV$curselect <- as.numeric(tkcurselection(ttc))+1
                
                if(ENV$curselect[1]=="1"){
                    ENV$ClinicalData <- ClinicalData
                    
                    title=paste(ENV$StudyRefCase[c],ENV$CasesStudies[ENV$curselectCases[c]+1], sep=": ")
                    getInTable(ENV$ClinicalData, title)
                    tkdestroy(ttClin)
                } else{
                    #ENV$ClinicalData <- ClinicalData[!(is.na(ClinicalData[,ENV$curselect]) | ClinicalData[,ENV$curselect]==""), ] 
                    ENV$ClinicalData <- ClinicalData[ENV$curselect-1]
                    
                    
                    title=paste(ENV$StudyRefCase[c],ENV$CasesStudies[ENV$curselectCases[c]+1], sep=": ")
                    getInTable(ENV$ClinicalData, title)
                    
                    tkdestroy(ttClin)
                    
                }
                
                
            }
            OK.but <- tkbutton(ttClin,text="OK",command=OnOK)
            tkgrid(OK.but,columnspan=2)
            
            ##Waiting to checkbox before to access to the next clinical data
            tkwait.window(ttClin)
        }
        
        
        
        ##getSummaryGSEA funtion needs return(ClinicalSub_All)
        if(getSummaryGSEAExists==1){
            return(ENV$ClinicalData)
        }
    }
    
    
}