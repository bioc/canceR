myGlobalEnv <- new.env(parent = emptyenv())
#' main function
#' @usage
#' canceR()
#'
#' @return open the starting windows with cancer studies
#' @export
#'
#' @examples
#' myGlobalEnv <- new.env(parent = emptyenv())
#' \dontrun{
#' canceR()
#'}
#' @import tkrplot
#' @import tcltk
#' @import cBioPortalData
#' @import R.oo
#' 
#'@importFrom graphics axis image layout legend lines par plot points text
#'@importFrom stats as.formula cor density dist hclust median na.exclude p.adjust pnorm sd setNames window
#'@importFrom utils browseURL capture.output memory.limit read.table write.table read.delim
#'@importFrom grDevices colors dev.cur dev.off graphics.off jpeg pdf png rainbow dev.new savePlot
#'@importFrom R.methodsS3 setMethodS3
#'

canceR <- function(){
    
    #if (!require("cgdsr")) devtools::install_version("cgdsr",version="1.3")
    
    ## connect to cBioPortal API
    myGlobalEnv$cgds <- cBioPortal(
        hostname = "www.cbioportal.org",
        protocol = "https",
        api = "/api/v2/api-docs"
    )

    ## Get all Cancer Studies using column 2 (description)
    myGlobalEnv$Studies <- cBioPortalData::getStudies(myGlobalEnv$cgds) 
    myGlobalEnv$Studies_name <- myGlobalEnv$Studies |> pull('name')
    
    ## first dialog START or CANCEL
    myGlobalEnv$ttMain <- tktoplevel()
    tktitle(myGlobalEnv$ttMain) <- "Search for Cancer Studies"
    #tkwm.resizable(myGlobalEnv$ttMain, FALSE, TRUE)
    #tkwm.geometry(myGlobalEnv$ttMain, "300x170")
    
    topMenu <- tkmenu(myGlobalEnv$ttMain)           # Create a menu
    tkconfigure(myGlobalEnv$ttMain, menu = topMenu) # Add it to the 'tt' window
    fileMenu <- tkmenu(topMenu, tearoff = FALSE)
    HelpMenu<- tkmenu(topMenu, tearoff = FALSE)
    
    ## TopMenus
    tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
    tkadd(topMenu, "cascade", label = "Help", menu = HelpMenu)
    
    ## Second Menus
    # GeneMenu<- tkmenu(AnalysisMenu, tearoff = FALSE)
    
    quit <- function(){
        QuitMsg <- tkmessageBox(message="Are you sure want to quit?",
                                icon="question",type="yesnocancel",default="yes")   
        if ( tclvalue(QuitMsg )== "yes") {
            tkdestroy(myGlobalEnv$ttMain)
            Sys.chmod(getwd(), mode = "0777", use_umask = TRUE)
            rm(list=ls(), envir=.GlobalEnv)
        }
        else {tkfocus(myGlobalEnv$ttMain)}
    }
    
    
    tkadd(fileMenu, "command", label = "Set Workspace", command = function() setWorkspace())
    tkadd(fileMenu, "command", label = "Get Gene Exp", command = function() getGeneExpMatrix())
    tkadd(fileMenu, "command", label = "Get Clinical Data", command = function() getClinicalDataMatrix())
    tkadd(fileMenu, "command", label = "Quit", command = function() quit())
    tkadd(HelpMenu, "command", label= "Vignette", command= function() canceR_Vignette())
    tkadd(HelpMenu, "command", label= "Issue", command = function() canceR_Issue())
    tkadd(HelpMenu, "command", label = "About", command= function() about())
    
    
    
    loadAllStudies <- function(){
        #delete last search by key words
        if(exists("n_studies", envir = myGlobalEnv)){
            rm(n_studies, envir=myGlobalEnv)
            rm(Studies_name, envir=myGlobalEnv)
        }else if(exists("Matched_index", envir = myGlobalEnv)){
            rm(Matched_index, envir=myGlobalEnv)
            rm(matched_studies, envir = myGlobalEnv)
        }
        
        myGlobalEnv$n_studies <- nrow(myGlobalEnv$Studies)
        myGlobalEnv$Studies_name <- myGlobalEnv$Studies$name
        
        nbrStudies <- paste("Query result: ", myGlobalEnv$n_studies,
                            " Studies were loaded. Select one or more Studies to get Data Profiles.",
                            sep="")
        #tkgrid(tklabel(myGlobalEnv$ttMain,text= nbrStudies ))
        ##Insert ListBox
        tkdelete(tlMain,0,600)
        tkdelete(tlInfo,0,1)
        for(s in myGlobalEnv$Studies_name)
        {
            tkinsert(tlMain,"end", s )
        }
        # Default selection.  Indexing starts at zero.
        tkselection.set(tlMain,2) 
        
        
        tkinsert(tlInfo,"end",text= nbrStudies )
        tkdelete(tlInfo,0)
    }
    
    loadAllStudies.button <- tkbutton(myGlobalEnv$ttMain,
                                      text = "Load all Studies", command = loadAllStudies)
    
    loadMatchStudies <- function(Word){
        #delete last load of all Studies
        if(exists("n_studies", envir = myGlobalEnv)){
            rm(n_studies, envir=myGlobalEnv)
            rm(Studies_name, envir=myGlobalEnv)
        }else if(exists("Matched_index", envir = myGlobalEnv)){
            rm(Matched_index, envir=myGlobalEnv)
            rm(matched_studies, envir = myGlobalEnv)
        }
        # select raw with matched "string"
        myGlobalEnv$Matched_index <- grep(Word, myGlobalEnv$Studies$name, ignore.case=TRUE) 
        
        #myGlobalEnv$matched_studies <- myGlobalEnv$Studies[myGlobalEnv$Matched_index, "name"] |> pull()
        myGlobalEnv$matched_studies <- myGlobalEnv$Studies$name[myGlobalEnv$Matched_index]
        
        ##Count the nomber of Matched Studies and return the number.
       nbrMatchedStudies_msg <- paste("Query result: ",
                                   length(myGlobalEnv$Matched_index), 
                                   " Studies were Matched. Select one or more to get its features.",
                                   sep="")
        #tkgrid(tklabel(myGlobalEnv$ttMain,text= nbrMatchedStudies_msg ))
        tkdelete(tlMain,0,600)
        tkdelete(tlInfo,0,1)
        for( m in myGlobalEnv$matched_studies){ 
            tkinsert(tlMain,"end", m)
        }
        tkselection.set(tlMain,2)  # Default selection.  Indexing starts at zero.
        
        tkinsert(tlInfo,"end",text= nbrMatchedStudies_msg)
        tkdelete(tlInfo,0)
    }
    
    launchDialog <- function() {
        Word <- modalDialog("Search Studies", "Search by Key Word", "")
        if (Word == "ID_CANCEL") return()
        loadMatchStudies(Word)
        
    }
    
    
    loadMatchStudies.button <- tkbutton(myGlobalEnv$ttMain,
                                        text = "Search by key words", command = launchDialog)
    
    tkgrid(loadAllStudies.button,loadMatchStudies.button, column=0)
    tkgrid.configure(loadAllStudies.button, sticky="w")
    tkgrid.configure(loadMatchStudies.button, sticky="e")
    
    
    
    tlInfo<-tklistbox(myGlobalEnv$ttMain,height=1, width= 65 ,selectmode="multiple",xscrollcommand=function(...)tkset(xscrInfo,...),background="white", foreground="blue")
    
    xscrInfo <- tkscrollbar(myGlobalEnv$ttMain, repeatinterval=2,orient="horizontal",
                            command=function(...)tkxview(tlInfo,...))
    
    tkgrid(tlInfo)
    
    ##Define ListBox of Studies
    tlMain<-tklistbox(myGlobalEnv$ttMain,height=18, width= 65 ,selectmode="multiple",
                      xscrollcommand=function(...)tkset(xscr,...),
                      yscrollcommand=function(...)tkset(yscr,...),background="white")
    yscr <- tkscrollbar(myGlobalEnv$ttMain, repeatinterval=2,
                        command=function(...)tkyview(tlMain,...))
    xscr <- tkscrollbar(myGlobalEnv$ttMain, repeatinterval=2,orient="horizontal",
                        command=function(...)tkxview(tlMain,...))
    
    tkgrid(tlMain,yscr)
    tkgrid.configure(yscr,rowspan=20,sticky="nsw")
    tkgrid(xscr)
    tkgrid.configure(xscr,rowspan=1,sticky="ew")
    
    police <- tkfont.create(family="arial", size=11)
    tkconfigure(tlMain, foreground="black", font=police)
    
    ##listbox Info
    
    loadCasesGenProfs <- function(){ 
        if(exists("n_studies", envir = myGlobalEnv)){
            myGlobalEnv$StudyChoice <- myGlobalEnv$Studies_name[as.numeric(tkcurselection(tlMain))+1]
            myGlobalEnv$checked_StudyIndex <- as.numeric(tkcurselection(tlMain))+1
            
        } else if (exists("matched_studies", envir=myGlobalEnv)){
            myGlobalEnv$StudyChoice <- myGlobalEnv$matched_studies[as.numeric(tkcurselection(tlMain))+1]
            myGlobalEnv$checked_StudyIndex <- myGlobalEnv$StudyIndex[as.numeric(tkcurselection(tlMain))+1]
        }
        
        #Identify checked Studies and its Index in cgds
        myGlobalEnv$checked_Studies_id <- myGlobalEnv$Studies[myGlobalEnv$checked_StudyIndex,"studyId"] 

        ###tkmessageBox if No Study was selected
        if(length(myGlobalEnv$checked_Studies_id) ==0){
            msgStudies <- paste("You need to select at less one Study")
            tkmessageBox(title = "Greetings from R TclTk", message=msgStudies, icon = "Info" )
        } else {
            
            getCasesGenProfs()
            
        }
    }
    
    getCasesGenProfs.but <-tkbutton(myGlobalEnv$ttMain,
                                    text="Get Cases and Genetic Profiles for selected Studies",command=loadCasesGenProfs)
    tkgrid(getCasesGenProfs.but)
    tkfocus(myGlobalEnv$ttMain)
    
}

