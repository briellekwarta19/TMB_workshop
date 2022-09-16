setwd("C:\\courses\\FISH 559_22\\TMB Workshop\\Lecture Examples\\")
library(TMB)

MakeADFun <-
  function(..., DLL = TMB:::getUserDLL()){
    ## Set local working directory
    orig_dir <- getwd()
    setwd( tempdir() )
    on.exit( setwd(orig_dir) )
    ## Save inputs
    All_inputs <- list(..., DLL=DLL)
    save( All_inputs, file="All_inputs.RData")
    ## Copy DLL to tempdir
    DLL <- All_inputs$DLL
    DLLfull <- paste0(orig_dir,"/",DLL)
    ## Write file to source
    txt <- c("library( TMB )",
             paste0("dyn.load(dynlib('",DLLfull,"'))"),
             "load( 'All_inputs.RData' )",
             "Obj <- do.call(TMB::MakeADFun, All_inputs)"
    )
    writeLines(txt, paste0(DLL,".R"))
    ## Try running
    Bdg_output <- gdbsource(paste0(DLL, ".R"))
    # Sort out outcomes
    if(any(grepl("#0", Bdg_output))){
      message("Model has errors")
      print(Bdg_output)
      stop()
    }
    ## OK ==> Safe to run in main session:
    TMB::MakeADFun(...,DLL=DLL)
  }

compile <-
  function(...){
    args <- list(...)
    if (.Platform$OS.type == "windows") {
      args$flags <- "-O1 -g"
      args$DLLFLAGS <- ""
    } else {
      args$flags <- "-O0 -g"
    }
    do.call(TMB::compile, args)
  }

gdbsource("LectB3.R",interactive=TRUE)
