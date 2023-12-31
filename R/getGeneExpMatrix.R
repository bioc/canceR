#' get matrix with gene expression from file
#' @usage getGeneExpMatrix()
#' @export
#' @return dataframe of gene expression
#' @examples 
#' readRDS(paste(path.package("canceR"),"/extdata/rdata/brca_tcga73genes.rds"", sep=""))
#' \dontrun{
#' getGeneExpMatrix()
#' }
getGeneExpMatrix <- function(){
    
    Sys.chmod(getwd(), mode = "0777", use_umask = TRUE)
    if(exists("GeneExpMatrixfile", envir = myGlobalEnv)){
        rm(myGlobalEnv$GeneExpMatrixfile)
    }
    myGlobalEnv$GeneExpMatrixfile <- tclvalue(tkgetOpenFile(filetypes = "{{txt Files} {.txt}} {{All files} *}", title="Choose Gene Expression Matrix from File")) 
    
    if (!nchar(myGlobalEnv$GeneExpMatrixfile)) {
        tkmessageBox(message = "No file was selected!")
    } else {
        Sys.chmod(getwd(), mode = "0777", use_umask = TRUE)
        myGlobalEnv$ProfData<-read.table(myGlobalEnv$GeneExpMatrixfile)
        
        tkmessageBox(message = paste("The file selected was", myGlobalEnv$GeneExpMatrixfile))
        tkfocus(ttMain)
    }
    
}