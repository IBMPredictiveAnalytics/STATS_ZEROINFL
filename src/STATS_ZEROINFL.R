#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.0"

# History
# 25-Jul-2014 Original Version
namap = list(omit=na.omit, fail=na.fail)
methodmap = list(neldermead="Nelder-Mead", bfgs="BFGS", 
                cg="CG", lbfgsb="L-BFGS-B", sann="SANN",
                 brent="Brent")

### MAIN ROUTINE ###
dozeroinfl = function(modelsource="estimate", dependent=NULL, zeromodel=NULL, countmodel=NULL,
    zerooffset=NULL, countoffset=NULL, sameregressors=FALSE, id=NULL,
    countdist="poisson", zerolink="logit", dataset=NULL,
    naaction="omit", startvalues="genlin", maxiter=1000, 
    tol=.Machine$double.eps^(1/1.6),
    optmethod="bfgs", modelfile=NULL,
    workspaceaction="clear", workspaceoutfile=NULL) {
    # Estimate or predict from a zero-inflated count model

    setuplocalization("STATS_ZEROINFL")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Zero Inflation Counts")
    warningsprocname = gtxt("Zero Inflation Counts: Warnings")
    omsid="STATSZEROINFL"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(pscl), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "pscl"),dostop=TRUE)
        }
    )

    if (modelsource == "estimate") {
        if (is.null(dependent) || is.null(countmodel)) {
            warns$warn(gtxt("In estimation mode, the dependent and count model variables are required"),
                dostop=TRUE)
        }
        if (sameregressors) {
            if (!is.null(zeromodel)) {
                warns$warn(gtxt("Cannot specify the zero model regressor list if same regressor option chosen"),
                    dostop=TRUE)
            }
            zeromodel = countmodel
        }
        frml = buildfrml(dependent, zeromodel, countmodel, zerooffset, countoffset, sameregressors)
        alldata=union(dependent, c(countmodel, zeromodel, zerooffset, countoffset))
        EM = startvalues == "em"
        allargsest = as.list(environment())
    } else {
        if (is.null(dataset)) {
            warns$warn(gtxt("A dataset name must be specified for predict mode"),
                dostop=TRUE)
        }
        if (modelsource == "file" && is.null(modelfile)) {
            warns$warn(gtxt("A file was specified as the model source, but no name was given"), dostop=TRUE)
        }
        if (modelsource == "workspace" && !is.null(modelfile)) {
            warns$warn(gtxt("The model source was given as workspace, but an input model file was specified"),
                       dostop=TRUE)
        }
    }
    
    if (!is.null(dataset)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
                dostop=TRUE)
        }
        if (dataset %in% alldatasets) {
            warns$warn(gtxt("The output dataset name must not already be in use"),
                dostop=TRUE)
        }
    }

    # ---------------

    thismodel = modelsource
    thismodelfile = modelfile
    if (modelsource != "estimate") {
        if (thismodel == "file") {
            load(thismodelfile)
        }
        if (!exists("reszeroinfl") || !("zeroinfl" %in% class(reszeroinfl)) ) {
            warns$warn(gtxt("The workspace does not contain a ZEROINFL model"),
                dostop=TRUE)
        }
        alldata=allargsest$alldata
        alldata[1] = NULL
    }
    allargsnow = as.list(environment())
    dta = tryCatch(spssdata.GetDataFromSPSS(alldata, row.label=id, missingValueToNA=TRUE,
        factorMode="levels"),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
        )

    if (modelsource == "estimate") {
        if ((!is.null(countoffset) && is.factor(dta[[countoffset]])) ||
            (!is.null(zerooffset) && is.factor(dta[[zerooffset]]))) {
            warns$warn(gtxt("Offset variables cannot be categorical variables"),
                dostop=TRUE)
        }

        reszeroinfl = tryCatch(zeroinfl(frml, data=dta, na.action=namap[[naaction]], dist=countdist,
                link=zerolink, y=FALSE, x=FALSE,
                control = zeroinfl.control(EM=EM, maxit=maxiter,
                    reltol=tol, method=methodmap[[optmethod]])),
                error = function(e) {warns$warn(e$message, dostop=TRUE)},
                warning = function(e) {warns$warn(e$message, dostop=TRUE)}
            )

        ressuminfl = summary(reszeroinfl)

        allargsest$estimationdate = date()
    }

    displayresults(allargsest, allargsnow, reszeroinfl, ressuminfl, warns)

    if (!is.null(dataset)) {
        savepred(allargsest, allargsnow, reszeroinfl, dta, warns)
    }
    if (!is.null(workspaceoutfile)) {
        if (modelsource != "estimate") {
            warns$warn(gtxt("Ignoring workspace save request because command is not estimating the model"),
                       dostop=FALSE)
        } else {
            #rm(thismodel, thismodelfile)
            save(reszeroinfl, ressuminfl, allargsest, allargsnow, file=workspaceoutfile)
        }
    }
    if (workspaceaction == "clear") {
        rm(list=ls())
        rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv)
    } else {
        if (modelsource == "estimate") {
            assign("reszeroinfl", reszeroinfl, envir=.GlobalEnv)
            assign("ressuminfl", ressuminfl, envir=.GlobalEnv)
            assign("allargsest", allargsest, envir=.GlobalEnv)
            
        }
    }
}

buildfrml = function(dep, zeromodel, countmodel, zerooffset, countoffset, sameregressors) {
    # Return formula expression as formula object
    # format is dep ~ countmodel | zeromodel
    # zeromodel can be just 1
    # offset() can be added to either part
    
    # dep is the name of dependent variable
    # warns is the error message object
    
    if (!sameregressors && is.null(zeromodel)) {
        zeromodel = "1"
    }
    countmodel = paste(countmodel, collapse="+")
    zeromodel = paste(zeromodel, collapse="+")
    if (!is.null(countoffset)) {
        countmodel = paste(countmodel, paste("+ offset(", countoffset, ")", collapse=""))
    }
    if (!is.null(zerooffset)) {
        zeromodel = paste(zeromodel, paste("+ offset(", zerooffset, ")", collapse=""))
    }
    frml = paste(dep, "~", countmodel, "|", zeromodel, sep=" ")

    return(as.formula(frml))
}
    
displayresults = function(allargsres, allargsnow, res, ressum, warns) {
    # display results
    # allargs is the parameter set (estimation or prediction)
    
    StartProcedure(allargsres[["procname"]], allargsres[["omsid"]])
    
    # summary results
    # input specifications
    lbls = c(gtxt("Dependent Variable"),
             gtxt("Count Model Distribution"),
             gtxt("Zero-Inflation Link Model"),
             gtxt("Count Model Offset"),
             gtxt("Zero-Inflation Model Offst"),
             gtxt("Missing Value Treatment"),
             gtxt("Starting Value Method"),
             gtxt("Convergence"),
             gtxt("Number of Cases"),
             gtxt("Log Likelihood"),
             gtxt("Log Likelihood D. F."),
             gtxt("AIC"),
             gtxt("Theta"),
             gtxt("SE log(theta)"),
             gtxt("Output Dataset"),
             gtxt("Computational Algorithm"),
             gtxt("Number of Iterations"),
             gtxt("Maximum Number of Iterations"),
             gtxt("Convergence Tolerance"),
             gtxt("Model Estimation Date")
    )

    vals = c(
            allargsres$dependent,
            allargsres$countdist,
            allargsres$zerolink,
            ifelse(is.null(allargsres$countoffset), gtxt("--NA--"), allargsres$countoffset),
            ifelse(is.null(allargsres$zerooffset), gtxt("--NA--"), allargsres$zerooffset),
            allargsres$naaction,
            allargsres$startvalues,
            ifelse(ressum$converged, gtxt("Yes"), gtxt("No")),
            ressum$n,
            round(ressum$loglik, 5),
            ressum$n - ressum$df.residual,
            round(AIC(res), 5),
            ifelse(is.null(ressum$theta), gtxt("--NA--"), round(ressum$theta, 4)),
            ifelse(is.null(ressum$SE.logtheta), gtxt("--NA--"), ressum$SE.logtheta),
            ifelse(is.null(allargsnow$dataset), gtxt("--NA--"), allargsnow$dataset),
            allargsres$optmethod,
            tail(na.omit(ressum$optim$count),1),
            ressum$control$maxit,
            format(ressum$control$reltol, digits=4),
            allargsres$estimationdate
    )

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="ZEROINFLSUMMARY", outline=gtxt("Summary"),
        caption = gtxt("Computations done by R package pscl")
    )
    if (allargsnow$modelsource == "estimate") {        
        # coefficients
        columnlabels = c(gtxt("Estimate"),
             gtxt("Std. Error"),
             gtxt("z Value"),
             gtxt("Significance")
        )
        df = data.frame(ressum$coef$count)
        names(df) = columnlabels
        spsspivottable.Display(df, 
            title=gtxt("Count Model Coefficients"), 
            templateName="ZEROINFLCOEFCOUNT",
            outline=gtxt("Count Coefficients"),
            caption=gtxtf("Dependent Variable: %s", allargsres$dependent)
        )
        if (!is.null(ressum$coef$zero)) {
            df = data.frame(ressum$coef$zero)
            names(df) = columnlabels
            spsspivottable.Display(df,
                title=gtxt("Zero-Inflation Model Coefficients"),
                templateName="ZEROINFLCOEFZERO",
                outline=gtxt("Zero-Inflation Coefficients"),
                caption=gtxtf("Dependent Variable: %s", allargsres$dependent)
            )
        }
    }
    
    spsspkg.EndProcedure()
}

savepred = function(allargsest, allargsnow, res, dta, warns) {
    # save fitted values
    # in estimate mode, these come from the fitted values
    # in predict mode, they are constructed here
    
    dict = list()
    
    depname = allargsest[["dependent"]]
    dataset = allargsnow[["dataset"]]  # must have  been specified
    idname = ifelse(is.null(allargsnow[["id"]]), "", allargsnow[["id"]])
    if (allargsnow$modelsource == "estimate") {
        pred = data.frame(
            fitted(res), 
            residuals(res, type="pearson"),
            residuals(res, type="response")
        )
        rnlength = max(nchar(row.names(pred))) * 3 #unicode worst-case expansion
    } else { # predicting on new data
        # check that all predictors are supplied
        predictors = allargsnow$alldata
        ###predictors[allargsest[depname]] = NULL
        missingpred = setdiff(predictors, names(dta))
        if (length(missingpred) > 0) {
            warns$warn(gtxtf("The following predictors are missing from the input data: %s",
                paste(missingpred, collapse=", " )), dostop=TRUE)
        }
        pred = tryCatch(data.frame(
            predict(res, type="response", newdata=dta),
            predict(res, type="count", newdata=dta),
            predict(res, type="zero")),
            error=function(e) {warns$warn(e$message, dostop=TRUE)}
            )
        if (nrow(pred) == 0) {
            warns$warn(gtxt("All cases have missing data.  No predictions are produced"),
                       dostop=TRUE)
        }
        rnlength = max(nchar(row.names(pred))) * 3
    }

    # row names will always be nominal strings
    dict[[1]] = c("ID", idname, rnlength, paste("A", rnlength, sep=""), "nominal")
    dict[[2]] = c("PredictedResponse", gtxtf("Predicted Response for %s", depname), 0, "F8.2", "scale")
    if (allargsnow$modelsource == "estimate") {
        dict[[3]] = c("PearsonResiduals", gtxt("Pearson Residuals"), 0, "F8.2", "scale")
        dict[[4]] = c("RawResiduals", gtxt("Raw Residuals"), 0, "F8.2", "scale")
    } else {
        dict[[3]] = c("PredictedCount", gtxt("Count Mean without Inflation"), 0, "F8.2", "scale")
        dict[[4]] = c("PredictedProb", gtxt("Predicted Probability of Zero"), 0, "F8.3", "scale")
    }

    dict = spssdictionary.CreateSPSSDictionary(dict)
    spssdictionary.SetDictionaryToSPSS(dataset, dict)
    tryCatch(spssdata.SetDataToSPSS(dataset, data.frame(row.names(pred), pred)),
        error=function(e) {warns$warn(e$message, dostop=TRUE)}
    )
    spssdictionary.EndDataStep()
    
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_ZEROINFL"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_ZEROINFL"))
}


Run = function(args) {
    #Execute the STATS ZEROINFL command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("MODELSOURCE", subc="", ktype="str", var="modelsource",
            vallist=list("estimate", "file", "workspace")),
        spsspkg.Template("MODELFILE", subc="", ktype="literal", var="modelfile"),
        spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", 
            var="dependent"),
        spsspkg.Template("ZEROMODEL", subc="", ktype="existingvarlist", var="zeromodel", islist=TRUE),
        spsspkg.Template("COUNTMODEL", subc="", ktype="existingvarlist", var="countmodel", islist=TRUE),
        spsspkg.Template("ZEROOFFSET", subc="", ktype="existingvarlist", var="zerooffset"),
        spsspkg.Template("COUNTOFFSET", subc="", ktype="existingvarlist", var="countoffset"),
        spsspkg.Template("SAMEREGRESSORS", subc="", ktype="bool", var="sameregressors"),
        spsspkg.Template("COUNTDIST", subc="", ktype="str", var="countdist",
            vallist=list("poisson", "negbin", "geometric")),
        spsspkg.Template("ZEROLINK", subc="", ktype="str", var="zerolink",
            vallist=list("logit", "probit", "cloglog", "cauchit", "log")),
        
        spsspkg.Template("STARTVALUES", subc="OPTIONS", ktype="str", var="startvalues",
            vallist=list("genlin", "em")),
        spsspkg.Template("OPTMETHOD", subc="OPTIONS", ktype="str", var="optmethod",
            vallist=list("neldermead", "bfgs", "cg",  "sann")),
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter"),
        spsspkg.Template("TOL", subc="OPTIONS", ktype="float", var="tol",
            vallist=list(1e-10)),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missing",
            vallist=list("omit", "fail")),
        
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
        spsspkg.Template("ID", subc="SAVE", ktype="existingvarlist", var="id"),
        spsspkg.Template("WORKSPACEACTION", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("WORKSPACEOUTFILE", subc="SAVE", ktype="literal", 
            var="workspaceoutfile")        
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dozeroinfl")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}