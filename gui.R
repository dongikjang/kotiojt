#####################################################################################
eval_gui <- function(guitoolkit="RGtk2"){
  guitoolkit="RGtk2"
  if(guitoolkit == "RGtk2"){
    options("guiToolkit"="RGtk2")
    require(RGtk2)
    require(gWidgetsRGtk2)
  }else{
    if(guitoolkit == "tcltk"){
      options("guiToolkit"="tcltk")
      require(tcltk2)
      require(gWidgetstcltk)
    }
  }
  require(gWidgets)	
  #require(gdata)
  library(xlsx)
  #library(fields)
  #library(ggmap)
  #library(ellipse)
  library(gstat)
  library(sp)
  library(automap)

  cat("\n")
  cat("#################################################################\n")
  cat("GUI of Traffic Volume Estimation.\n")
  cat("Version: 2013.02.13\n")
  cat(paste("Graphical user interface toolkit:", guitoolkit, "\n"))
  cat("#################################################################\n\n")
  
  os <- .Platform$OS.type
  
  ##################################################################	
  ## Set up top level
  
  wins <- gwindow("Prediction of Traffic Volume", width=700,  height=500, visible = FALSE, expand=TRUE)
  
  g <- ggroup(horizontal=FALSE, cont=wins, expand=TRUE, height=500, width=700)  # main group
  
  ##################################################################	
  ## Open button with icon
  
  widgets <- NULL
  tb <- list()
  tb$Open$handler = function(h,...){
    filesep <- ifelse(.Platform$OS.type =="windows", "\\\\", "/")
    datafile <- file.choose()
    #print(datafile)
    datastring <- strsplit(datafile, filesep)
    datastring <- datastring[[1]][length(datastring[[1]])]
    datastring <- gsub("-", "_", datastring)
    datastring <- strsplit(datastring, " ")[[1]][1]
    
    ## Open xls file: Main widndows
    winsheet <- gwindow("Infile the xls file", visible=FALSE, width=40)
    
    ## Open xls file: main group
    groupsheet <- ggroup(horizontal = FALSE, container=winsheet)
    
    ## Open xls file: Frame
    fblsheet <- gframe(cont=groupsheet, text="Option for sheet and header", fill="x", expand=TRUE)
    tblsheet <- glayout(cont=fblsheet, spacing = 12)
    
    ## Open xls file: Widgets
    widgetsheet <- list()
    
    ## Open xls file: Argument
    tblsheet[1,1:4, anchor=c(-1,0)] <- "Type the sheet number "
    tblsheet[2,2:4, anchor=c(-1,0)] <- (widgetsheet[['sheet']] <- gedit("1", container=tblsheet))
    size(widgetsheet[['sheet']]) <- c(40, 30)
    
    tblsheet[3,1:4, anchor=c(-1,0)] <- "Does the file contain the column names of the variables as its first line?"
    tblsheet[4,2:4, anchor=c(-1,0)] <- (widgetsheet[['cheader']] <- gdroplist(c("Yes", "NO"), container=tblsheet, selected=1))
    size(widgetsheet[['cheader']]) <- c(40, 30)
    
    tblsheet[5,1:4, anchor=c(-1,0)] <- "Does the file contain the row names of the variables as its first column?"
    tblsheet[6,2:4, anchor=c(-1,0)] <- (widgetsheet[['rheader']] <- gdroplist(c("Yes", "No"), container=tblsheet, selected=1))
    size(widgetsheet[['rheader']]) <- c(40, 30)
    
    widgetsheet[['geom']] <- gbutton("ok", cont=groupsheet)
    font(widgetsheet[['geom']]) <- c(color="red", style="bold")
    size(widgetsheet[['geom']]) <- c(30,40)
    
    widgetsheet[['datafile']] <- datafile 
    widgetsheet[['datastring']] <- datastring
    visible(winsheet) <- TRUE
    
    ## Hander of clicked and changed
    addHandlerClicked(widgetsheet[['geom']], handler = function(h,...) readxlssheet() )
    addHandlerChanged(widgetsheet[['sheet']], handler = function(h,...)  readxlssheet() )
    
    ## read xls handler
    readxlssheet <- function(){
      
      datafile <- widgetsheet[['datafile']]
      datastring <- widgetsheet[['datastring']]
      tmpsheet <-  as.numeric(lapply(widgetsheet, svalue)$sheet)
      cheader <- lapply(widgetsheet, svalue)$cheader
      cheader <- ifelse(cheader=="Yes", TRUE, FALSE)
      if(lapply(widgetsheet, svalue)$rheader == "Yes"){
        xx <- read.xlsx(datafile, sheetIndex=tmpsheet, header=cheader, startRow=1, encoding="UTF-8") #row.names=1)
      }else{
        xx <- read.xlsx(datafile, sheetIndex=tmpsheet, header=cheader, encoding="UTF-8")
      }
      #traffic <- read.xlsx(paste(path, "Data/140121_KRIGING모형.xlsx", sep=""), 1)
      
      #colnames(xx) <- paste("V", 1:ncol(xx), sep="")			
      ctmpsheet <- ifelse(tmpsheet < 10, paste("0", tmpsheet, sep=""), as.character(tmpsheet))
      
      datastring <- paste(strsplit(datastring, "[.]")[[1]][1], "_s", ctmpsheet, sep="")
      
      assign(datastring, xx, envir=.GlobalEnv)
      invisible()
      update(vb)
      w <- gwindow(datastring, width=800)
      datagp <- ggroup(horizontal=FALSE, cont=w)
      #tmpration <- NULL
      #for(i in 1:ncol(xx)){
      #  tmpvalr <- as.character(fractions(xx[,i]))
      #  tmpvalr[is.na(tmpvalr)] <- ""
      #  tmpration <- cbind(tmpration, tmpvalr)
      #}
      #colnames(tmpration) <- colnames(xx)
      #tmpration <- as.data.frame(tmpration)
      gtable(xx, chosencol=1, cont=datagp, expand=TRUE)
      
      ## close the data table 
      dispose(winsheet)
      
      ## update the drop-down list of variables 
      if(!is.null(widgets[['data']])){
        is.variables <- function(i) is.matrix(i) || is.data.frame(i)  || is.vector(i) 
        tmp0 <- sapply(ls(envir=.GlobalEnv), function(i) is.variables(get(i)))
        tmp0 <- names(tmp0)[tmp0]
        widgets[['data']][] <- tmp0
        widgets[['pdata']][] <- tmp0
        #mwidgets[['data']][] <- tmp0
        
#        if(options("guiToolkit")$guiToolkit == "RGtk2"){
#          (widgets[['data']]@widget@widget)$SetActive(which(tmp0 == datastring)-1)
#          #(mwidgets[['data']]@widget@widget)$SetActive(which(tmp0 == datastring)-1)
#        }
#        if(options("guiToolkit")$guiToolkit == "tcltk"){
#          if(length(tmp0) == 1) {
#            widgets[['data']][] <- c(tmp0, "")
#            #mwidgets[['data']][] <- c(tmp0, "")
#          }
#          svalue(widgets[['data']]@widget, index=TRUE) <- which(tmp0 == datastring)
#          #svalue(mwidgets[['data']]@widget, index=TRUE) <- which(tmp0 == datastring)
#        }
        
      }
    }
  }
  
  ## List of  data set	
  is.variables <- function(i) is.matrix(i) || is.data.frame(i)  || is.vector(i) 
  tmps <- sapply(ls(envir=.GlobalEnv), function(i) is.variables(get(i)) )
  tmps <- names(tmps)[tmps]
  
  if(length(tmps) == 0) tmps <- " "	
  
  tb$Open$lavel <- "Open"
  tb$Open$icon  <- "open"
  tb <- gtoolbar(tb, container = wins, icon="open", style="both-horiz")
  
  ##################################################################	
  ## Main layout with variable browser	
  
  pg <- gpanedgroup(cont=g, exand=FALSE, height=500, width=900)
  vb <- gvarbrowser(cont=pg, action="show", silent=FALSE, height=450, exand=TRUE)         # left varbrowser
  nb <- gnotebook(cont=pg,, height=500, width=300)				# right is notebook
  
  ##################################################################	
  ## K-means evaluation: main widgets
  widgets <- list()
  
  ## K-means evaluation: main group
  qpg <- ggroup(horizontal=FALSE, cont=nb, label="Spatial Prediction of Traffic Volume")
  
  ## K-means evaluation: Data Set frame	
  fbl <- gframe(cont=qpg, text="Modeling Data Set ", fill='x', expand=TRUE)
  font(fbl) <- c(weight="bold", size="large", style="italic")
  tbl <- glayout(cont=fbl)
  tbl[1, 1, anchor=c(1,0)] <- "Data frame: "
  tbl[1, 2:10] <- (widgets[['data']] <- gdroplist(tmps, cont=tbl, editable=TRUE, 
                                                selected=ifelse(length(which(tmps=="modelset")) != 0, which(tmps=="modelset"), 1) )) 
  size(widgets[['data']]) <- c(300, 30)


  addHandlerChanged(widgets[['data']], handler = function(h,...)  updatemodelbox() )


  ## K-means evaluation: Argument frame
  fbl2 <- gframe(cont=qpg, text="Model ", expand=TRUE)
  font(fbl2) <- c(weight="bold", size="large", style="italic")
  
  variableNames = "" #colnames(traffic)
  tbl2 = glayout(cont=fbl2)
  
  tbl2[1,1] <- (widgets[['response']] <- gdroplist(variableNames, cont=tbl2, editable=TRUE))
  tbl2[2,1] <- "response"
  size(widgets[['response']]) <- c(150, 30)
  
  ## add ~
  tbl2[1,2] <- " ~ "
  
  tbl2[1,3:10] <- (widgets[['predictor']] <- gdroplist(variableNames, cont=tbl2, editable=TRUE))
  tbl2[2,3:10] <- "predictor(s)"
  size(widgets[['predictor']]) <- c(600, 30)
  
  tbl2[1,11] <- (predictorEdit <- gbutton("edit", container=tbl2))
  

  #svalue(widgets[['predictor']])
  #svalue(widgets[['response']])
  
  

  #tbl2[1,1:3, anchor=c(-1,0)] <- "Response variable: "
  #tbl2[1,4:5] <- (widgets[['type1']] <- gdroplist(colnames(traffic), cont=tbl2))
  #size(widgets[['type1']]) <- c(200, 30)
  
  
  #tbl2[2,2:3, anchor=c(-1,0)] <- "Transform: "
  #tbl2[2,4:5] <- (widgets[['trim1']] <- gdroplist(c("None","log", "sqrt"), editable=TRUE, cont=tbl2, selected=2) )
  #size(widgets[['trim1']]) <- c(100, 30)
  
  
  #response = gedit(svalue(widgets[['response']]), cont=tbl2)
  #predictor = gedit(svalue(widgets[['predictor']]), cont=tbl2)

  #responsewidget =widgets[['response']]
  #predictorwidget = widgets[['predictor']]
  
  
  editPredictorHandler = function(h,...) {
    editFormulaDialog(data=widgets[['data']],
                      responsewidget=widgets[['response']],
                      predictorwidget=widgets[['predictor']])
  }
  addhandlerclicked(predictorEdit, handler=editPredictorHandler)  
    
#  addHandlerChanged(predictorEdit,
#                    handler=editPredictorHandler,
#                    action=list(
#                                responseEntry=responseEntry,
#                                predictorEntry=predictorEntry))
#




  addSpace(qpg, 30, horizontal=FALSE)


  fbl3 <- gframe(cont=qpg, text="Assign model output to: ")
  tbl3 <- glayout(cont=fbl3)
  tbl3[1,1:2] <- widgets[['assignto']]  <- gedit("outs", cont=tbl3) 
  size(widgets[['assignto']]) <- c(100, 30)
  
  ## K-means evaluation: Ok button
  widgets[['geom']] <- gbutton("ok", cont=qpg)
  font(widgets[['geom']]) <- c(color="red", style="bold")
  size(widgets[['geom']]) <- c(165,40)

  addHandlerClicked(widgets[['geom']], handler = function(h,...) fitmodel() )
  
  fitmodel <- function(){
    
    #print(svalue(widgets[['predictor']]))
    #print(Encoding(svalue(widgets[['predictor']])))
    if(.Platform$OS.type=="windows"){
      pred <- iconv(svalue(widgets[['predictor']]), "UTF-8", localeToCharset())
      resp <- iconv(svalue(widgets[['response']]), "UTF-8", localeToCharset())
    } else {
      pred <- svalue(widgets[['predictor']])
      resp <- svalue(widgets[['response']])
    }
    kmodel <- as.formula(paste(resp, "~", pred))
    dataset <- eval(as.symbol(svalue(widgets[['data']])))
    coordinates(dataset) <-  ~ LINK_MID_X + LINK_MID_Y
    set.seed(3245)
    kfit <- autoKrige.cv(kmodel, dataset, model = c("Exp"), nfold = 10, verbose=c(FALSE, FALSE))
    
    vm.uk.opt<- variogram(kmodel, dataset)  
    mopt <- autoKrige(kmodel, dataset, dataset, model="Exp")$var_model
    
    assign(svalue(widgets[['assignto']]), mopt, envir=.GlobalEnv)
  }


  addSpace(qpg, 30, horizontal=FALSE)

  


  ## K-means evaluation: Data Set frame  
  fbl4 <- gframe(cont=qpg, text="Prediction Data Set ", fill='x', expand=TRUE)
  font(fbl4) <- c(weight="bold", size="large", style="italic")
  tbl4 <- glayout(cont=fbl4)
  tbl4[1, 1, anchor=c(1,0)] <- "Data frame: "
  tbl4[1, 2:10] <- (widgets[['pdata']] <- gdroplist(tmps, cont=tbl4, editable=TRUE, 
                                                  selected=ifelse(length(which(tmps=="pdata")) != 0, which(tmps=="pdata"), 1) )) 
  size(widgets[['pdata']]) <- c(300, 30)

  addSpace(qpg, 30, horizontal=FALSE)


  fbl5 <- gframe(cont=qpg, text="Assign predicion output to: ")
  tbl5 <- glayout(cont=fbl5)
  tbl5[1,1:2] <- widgets[['assignto2']]  <- gedit("preds", cont=tbl5) 
  size(widgets[['assignto2']]) <- c(100, 30)
  
  ## K-means evaluation: Ok button
  widgets[['geom2']] <- gbutton("ok", cont=qpg)
  font(widgets[['geom2']]) <- c(color="red", style="bold")
  size(widgets[['geom2']]) <- c(165, 40)

  addHandlerClicked(widgets[['geom2']], handler = function(h,...) predmodel() )
  predmodel <- function(){
    if(.Platform$OS.type=="windows"){
    	pred <- iconv(svalue(widgets[['predictor']]), "UTF-8", localeToCharset())
      resp <- iconv(svalue(widgets[['response']]), "UTF-8", localeToCharset())
    } else {
      pred <- svalue(widgets[['predictor']])
      resp <- svalue(widgets[['response']])
    }
    kmodel <- as.formula(paste(resp, "~", pred))
    dataset <- eval(as.symbol(svalue(widgets[['data']])))
    pdataset <- eval(as.symbol(svalue(widgets[['pdata']])))
    mopt <- eval(as.symbol(svalue(widgets[['assignto']])))


    modelval <- gsub(" ", "", strsplit(pred, "[+]")[[1]])
    if(!all(modelval %in% names(pdataset))){
      
      confirmDialog <- function(message, handler=NULL) {
        window <- gwindow("Error", width=300, height=150)
        group <- ggroup(container = window)
        gimage("errors", size="dialog", container=group)
        ## A group for the message and buttons
        inner.group <- ggroup(horizontal=FALSE, container = group)
        glabel(message, container=inner.group, expand=TRUE)
        ## A group to organize the buttons
        button.group <- ggroup(container = inner.group)
        ## Push buttons to right
        addSpring(button.group)
        gbutton("close", , handler = function(h,...) dispose(window), container=button.group)
        #gbutton("cancel", handler = function(h,...) dispose(window), container=button.group)
        return()
      }
      
      confirmDialog(paste("The variable(s)", modelval[!modelval %in% names(pdataset)], "\n is(or are) not in prediction data set."), handler = function(h,...) {
        print("Please edit the column names of prediction dataset.")
        ## In this instance dispose finds its parent window and closes it
        dispose(h$obj)
      })
      
    }
    
    coordinates(dataset) <-  ~ LINK_MID_X + LINK_MID_Y
    coordinates(pdataset) <-  ~ LINK_MID_X + LINK_MID_Y
    
    
    
    predopt <- krige(kmodel, dataset, pdataset, model = mopt)
    if(.Platform$OS.type=="windows"){
      pred <- iconv(svalue(widgets[['predictor']]), "UTF-8", localeToCharset())
      resp <- iconv(svalue(widgets[['response']]), "UTF-8", localeToCharset())
    } else {
      pred <- svalue(widgets[['predictor']])
      resp <- svalue(widgets[['response']])
    }

    if(strsplit(resp, "[(]")[[1]][1] =="log"){
      assign(svalue(widgets[['assignto2']]), exp(predopt$var1.pred), envir=.GlobalEnv)
      print(exp(predopt$var1.pred))
    } else{
      assign(svalue(widgets[['assignto2']]), predopt$var1.pred, envir=.GlobalEnv)
      print(predopt$var1.pred)
    }
    
  }


  svalue(nb) <- 1
  visible(wins) <- TRUE
  




  ###################################################################################################################################################
  ###################################################################################################################################################
  
  updatemodelbox <- function(){
    tmp0 <- names(eval(as.symbol(svalue(widgets[['data']]))))
    widgets[['response']][] <- tmp0
    widgets[['predictor']][] <- tmp0
  }
  
  getObjectFromString = function(STRING="", envir=.GlobalEnv) {
    tmp = try(get(STRING, envir), silent = TRUE)
    if(!inherits(tmp, "try-error")) return(tmp)
    
    tmp = try(rpel(STRING,envir), silent=TRUE)
    if(!inherits(tmp, "try-error"))  return(tmp)
    
    ## out of chances
    return(NULL)
  }
  
  Paste = function(..., sep="", collapse="") {
    x = unlist(list(...))
    x = x[!is.na(x)]
    x = x[x != "NA"]
    out <- paste(x, sep=sep, collapse=collapse)
    #out <- iconv(out, , "UTF-8") 
    return(out)
  }
  
  
  editFormulaDialog = function(data=NULL, responsewidget = NULL,  predictorwidget = NULL ) {
    ## actually get data, not just a string
    if(is(data,"guiWidget") || is(data,"gComponentANY")) {
      data = svalue(data)
      dataName = deparse(substitute(data))
    }
    if(is.character(data) && length(data) == 1) {
      dataName = data 
      data = getObjectFromString(data)
    }
    
    if(is.na(data) || is.null(data)) {
      warning("Can't find data set")
      return(NA)
    }
    
    ## coerce data if possible
    if(!is.data.frame(data)) {
      tmp = try(as.data.frame(data), silent=TRUE)
      if(inherits(tmp,"try-error")) {
        warning("gtable shows data frames of vectors")
        return(NA)
      }
      data = tmp
    }
    varNames = names(data)
    
    ## define key widgets
    
    
    ## Set up the window
    ## main window
    win = gwindow("Edit model formula values")
    
    group = ggroup(horizontal=FALSE, container=win)
    
    datagroup = ggroup(container=group)
    size(datagroup) <- c(300,500)
    glabel("Dataset: ", container=datagroup)
    tmp = glabel(dataName, container=datagroup);
    font(tmp) <- c(weight="bold")
    addSpace(datagroup, 20)
    
    variables = gtable(varNames,multiple=TRUE, cont=datagroup, expand=TRUE)
    
    gseparator(cont=group)
    ## buttons
    
    
    buttonGroup = ggroup(container=group)
    glabel("Actions:", container=buttonGroup)
    tbl = glayout(cont=buttonGroup, expand=TRUE)
    tbl[1,1,anchor=c(-1,0)] = (addresponse <- gbutton("Response",cont=tbl))
    tbl[1,2,anchor=c(-1,0)] = (addlogresponse <- gbutton("log(Response)",cont=tbl))
    tbl[1,3,anchor=c(-1,0)] = (addterm <-  gbutton("+ (main effect)",cont=tbl))
    #tbl[2,1,anchor=c(-1,0)] = (addwithin <- gbutton(": (interaction)",cont=tbl))
    #tbl[2,2,anchor=c(-1,0)] = (addinteraction <- gbutton("* (main + interaction)",cont=tbl))
    #tbl[3,1,anchor=c(-1,0)] = (addsecondpowers <- gbutton("^2 (second-order)",cont=tbl))
    #tbl[3,2,anchor=c(-1,0)] = (subtractintercept <-  gbutton("remove intercept",cont=tbl))
    visible(tbl) <- TRUE
    
    gseparator(cont=group)
    tbl = glayout(container=group)
    
    response = gedit(svalue(responsewidget),cont=tbl)
    predictor = gedit(svalue(predictorwidget), cont=tbl)
    widgets2 <- list()
    tbl[1,1] <- response
    tbl[1,2] <- glabel(" ~ ",cont=tbl)
    tbl[1,3:6] <- (widgets2[['predictor']] <- predictor)
    tbl[2,1] <- "response"
    tbl[2,3:6] <- "predictor formula"
    visible(tbl) <- TRUE
    size(widgets2[['predictor']]) <- c(600, 30)
    
    
    
    buttonbox = ggroup(container=group)
    addSpring(buttonbox)
    okbutton = gbutton("ok", cont=buttonbox)
    addSpace(buttonbox,15)  
    clearbutton = gbutton("clear", cont=buttonbox)
    cancelbutton = gbutton("cancel", cont=buttonbox)
    
    ## Now add handlers
    addhandlerclicked(addresponse, handler=function(h,...) {
      vals = svalue(variables)
      if(!is.null(vals)) {
        svalue(response) <- vals[1]
      }
    })
    
    addhandlerclicked(addlogresponse, handler=function(h,...) {
      vals = svalue(variables)
      if(!is.null(vals)) {
        svalue(response) <- paste("log(", vals[1], ")", sep="")
      }
    })
    
  #  addhandlerclicked(addinteraction,handler=function(h,...) {
  #    vars = svalue(variables)
  #    if(!is.null(vars)) {
  #      oldval = svalue(predictor)
  #      if(!is.null(oldval) && oldval !="")
  #        oldval = Paste(oldval, " + ")
  #      else
  #        oldval = ""
  #      svalue(predictor) <- Paste(oldval, paste(vars, sep="", collapse=" * "))
  #    }
  #  })
    addhandlerclicked(addterm, handler=function(h,...) {
      vars = svalue(variables)
      #print("get encoding")
      #print(Encoding(vars))
      #print(vars)
      if(!is.null(vars)) {
        oldval = svalue(predictor)
        if(!is.null(oldval) && oldval !=""){
        	#oldval <- iconv(oldval, , "UTF-8")
        	#print("old encoding")
        	#print(Encoding(oldval))
        	#print(oldval)
          if(.Platform$OS.type=="windows"){
        	  oldval <- iconv(oldval, "UTF-8", localeToCharset())
          }
          oldval = Paste(oldval, " + ")
          #oldval <- iconv(oldval, localeToCharset(), "UTF-8")
          #print("new encoding")
          #print(Encoding(oldval))
          #print(oldval)
          	
        }
        else{
          oldval = ""
        }
        newval <- Paste(oldval, paste(vars, sep="", collapse=" + "))
        #print(Encoding(newval))
        #print(newval)
        svalue(predictor) <- newval
      }
    })
  #  addhandlerclicked(addsecondpowers,handler=function(h,...) {
  #    vars = svalue(variables)
  #    if(!is.null(vars)) {
  #      oldval = svalue(predictor)
  #      if(!is.null(oldval) && oldval !="")
  #        oldval = Paste(oldval, " + ")
  #      else
  #        oldval = ""
  #      svalue(predictor) <- Paste(oldval,
  #                                 " (",
  #                                 paste(vars, sep="", collapse=" + "),
  #                                 ")^2")
  #    }
  #  })
  #  addhandlerclicked(addwithin,handler=function(h,...) {
  #    vars = svalue(variables)
  #    if(!is.null(vars)) {
  #      oldval = svalue(predictor)
  #      if(!is.null(oldval) && oldval !="")
  #        oldval = Paste(oldval, " + ")
  #      else
  #        oldval = ""
  #      svalue(predictor) <- Paste(oldval,
  #                                 paste(vars, sep="", collapse=":")
  #      )
  #    }
  #    
  #  })
  #  addhandlerclicked(subtractintercept,handler=function(h,...) {
  #    svalue(predictor) <- Paste(svalue(predictor), " -1")
  #  })
    addhandlerclicked(okbutton, handler = function(h,...) {
      if(.Platform$OS.type=="windows"){
        svalue(responsewidget) <- iconv(svalue(response), "UTF-8", localeToCharset())
        #print(svalue(response))
        #print(Encoding(svalue(response)))
        svalue(predictorwidget) <- iconv(svalue(predictor), "UTF-8", localeToCharset())
        #print(svalue(predictor))
        #print(Encoding(svalue(predictor)))
      } else {
        svalue(responsewidget) <- svalue(response)
        svalue(predictorwidget) <- svalue(predictor)
      }
      dispose(win)
    })
    addhandlerclicked(clearbutton, handler=function(h,...) {
      svalue(response) <- ""
      svalue(predictor) <- ""
    })
    addhandlerclicked(cancelbutton,handler=function(h,...) {
      dispose(win)
    })
  }
}
