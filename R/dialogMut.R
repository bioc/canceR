#' Dialog bos to set returned Mutation information
#' @usage
#' dialogMut(title, question, entryInit, entryWidth = 40,returnValOnCancel = "ID_CANCEL")
#' @param title title of the table
#' @param question question
#' @param entryInit entryInit
#' @param entryWidth 40
#' @param returnValOnCancel "ID_CANCEL"
#'
#' @return a check box with mutations variables
#' @export
#'
#' @examples
#' readRDS(paste(path.package("canceR"),"/extdata/rdata/ucec_tcga_pubGSEA1021.rds", sep=""))
#' \dontrun{
#' dialogMut("title", "question", "entryInit", entryWidth = 40, returnValOnCancel = "ID_CANCEL")
#' }
dialogMut <- function(title, question, entryInit, entryWidth = 40,
                      returnValOnCancel = "ID_CANCEL") {
    dlg <- tktoplevel()
    tkwm.deiconify(dlg)
    tkgrab.set(dlg)
    tkfocus(dlg)
    tkwm.title(dlg, title)
    
    textEntryVarTcl <- tclVar(paste(entryInit))
    textEntryWidget <- tkentry(dlg, width = paste(entryWidth),
                               textvariable = textEntryVarTcl)
    tkgrid(tklabel(dlg, text = "       "))
    tkgrid(tklabel(dlg, text = question), textEntryWidget)
    tkgrid(tklabel(dlg, text = " *\n .\n |",justify = "right"), column=0, row=2, sticky="n")
    tkgrid(tklabel(dlg, text = " several time the last character\n any character or number\n or", justify = "left"), column=1,row=2, sticky="w")
    ReturnVal <- returnValOnCancel
    
    onOK <- function() {
        ENV$ReturnVal <- tclvalue(textEntryVarTcl)
        tkgrab.release(dlg)
        tkdestroy(dlg)
        
        #tkfocus(ttMain)
    }
    onCancel <- function() {
        ENV$ReturnVal <- returnValOnCancel
        tkgrab.release(dlg)
        tkdestroy(dlg)
        return(ReturnVal)
        #tkfocus(ttMain)
    }
    OK.but <- tkbutton(dlg, text = "   OK   ", command = onOK)
    Cancel.but <- tkbutton(dlg, text = " Cancel ", command = onCancel)
    tkgrid(Cancel.but,OK.but, column=1, sticky="w")
    tkgrid.configure(Cancel.but, column=1, sticky="n")
    #tkgrid(tklabel(dlg, text = "    "))
    
    tkfocus(dlg)
    tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg); tkfocus(ENV$ttCasesGenProfs)})
    tkbind(textEntryWidget, "<Return>", onOK)
    tkwait.window(dlg)
    
    return(ENV$ReturnVal)
}