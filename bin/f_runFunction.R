runFunction <-
  function(funcfile='',outputfile='',NamedListOfAddArgs=NULL,save=TRUE){
    if (!file.exists(outputfile)) {
      if (!exists(funcfile)) {
        source(funcfile)
      }
      func <- gsub("\\.R","",gsub(".*f_","",funcfile))
      output <- gsub("\\.Rdata","",gsub(".*\\/","",outputfile))
      command <- paste(func,"(",sep='')
      acount <- 1
      for (a in names(NamedListOfAddArgs)) {
        if (acount == length(NamedListOfAddArgs) ) {
          command <- paste(command,a,"=",NamedListOfAddArgs[[a]],sep='')
        }  else {
          command <- paste(command,a,"=",NamedListOfAddArgs[[a]],',',sep='')
        }
        acount <- acount +1
    #    print(a)
    #    print(NamedListOfAddArgs[[a]])
      }
      command <- paste(command,")",sep='')
      cat("Running:",command,"\n")
      output <- eval(parse(text=command))
      if (save) {
        save(output,file=outputfile)
        cat("Made ",outputfile,"\n")
      }
    } else {
      if (save) {
        load(outputfile)
        cat("Loaded ",outputfile,"\n")
      }
    }
    invisible(output)
  }
