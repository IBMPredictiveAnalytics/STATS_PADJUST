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
# version__ = "1.0.1"

# History
# 17-apr-2014 Original Version


helptext='STATS PADJUST 
    DSNAMES=variable name DSPVALUES=variable name
    or
    NAMELIST = list of names PLIST = list of p values
    METHOD=FDR HOLM HOCHBERG HOMMEL BY BONFERRONI
    NTESTS=number
/OUTPUT PLOT=YES or NO.

Example:
STATS PADJUST NAMELIST="x" "y" "z" PLIST=.001 .005 .05.

Compute adjusted p values using one or more adjustment types
for multiple tests.

The input for this command is either from the active dataset or
specified on the command line.  For dataset input, specify the
name of a variable holding a label for each p value as DSNAMES
and the name of the variable holding the p values as DSPVALUES.

For input on the command line, specify the pvalues to PLIST
and, optionally, specify a list of names to NAMELIST.  The names
should be quoted.

Specify the methods to be used as METHOD.
FDR is Benjamini-Hochberg FDR and is the default
HOLM is Holm
HOMMEL is Hommel
HOCHBERG is Hochberg
BY is Benjamini-Yekutieli
BONFERRONI is Bonferroni.

FDR and BY are used to control the false discovery rate, which is
the expected proportion of false positives among the hypotheses
that are rejected.  The other methods are used to control the
familywise error rate and are more conservative.  The distinction
becomes more important as the number of tests grows.
Hochberg and Hommel assume that the tests are independent or
are not negatively correlated.

By default, the number of tests is assumed to be the
same as the number of p values supplied.  However, you
can use NTESTS to set a larger number.  This assumes that
the omitted p values are either larger than any of the 
included ones (Bonferroni and Holm) or are equal to 1.

PLOT specifies whether or not to produce a plot of the 
adjusted p values against the unadjusted ones.  The default
value is YES.

STATS GET R /HELP.  prints this information and does nothing else.
'

gtxt <- function(...) {
    return(gettext(...,domain="STATS_PADJUST"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_PADJUST"))
}
# method name translation mapping
# names should probably be left as these words but might be
# transliterated
transnames = list(p=gtxt("P"), fdr=gtxt("Benjamini-Hochberg"), holm=gtxt("Holm"), hochberg=gtxt("Hochberg"),
                hommel=gtxt("Hommel"), bonferroni=gtxt("Bonferroni"),
                BY=gtxt("Benjamini-Yekutieli"))

# main routine
dopadjust = function(dsnames=NULL, dspvalues=NULL, namelist=NULL, plist=NULL,
    method="fdr", ntests=NULL, plotit=TRUE) {
    
    setuplocalization("STATS_PADJUST")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Adjust P Values")
    warningsprocname = gtxt("Adjust P Values: Warnings")
    omsid="STATSPADJUST"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    sources = list(dsnames, dspvalues, namelist, plist)
    nulls = lapply(sources, is.null)
    datasetsource = sum(nulls[[1]], nulls[[2]]) < 2
    commandsource = sum(nulls[[3]], nulls[[4]]) < 2
    if (all(unlist(nulls)) || (datasetsource && commandsource) ||
        (commandsource && is.null(plist)) ||
        (datasetsource && is.null(dspvalues))) {
            warns$warn(gtxt("Invalid combination of probability and name arguments"), dostop=TRUE)
    }

    if (datasetsource) {
        dta = spssdata.GetDataFromSPSS(c(dsnames, dspvalues), missingValueToNA=TRUE)
        if (is.null(dsnames)) {
            dsnames="dsnames"
            dta[dsnames] = 1:nrow(dta)
        }
    } else {
        dta = data.frame(plist)
        if (!is.null(namelist)) {
            tryCatch({dta[["dsnames"]] = namelist}, 
                error=function(e) {warns$warn(gtxt("The lengths of the names and p-values lists are different"),
                    dostop=TRUE)}
            )
        } else {
            dta["dsnames"] = 1:nrow(dta)
        }
        dsnames = "dsnames"
        dspvalues = "plist"
    }
    hasby = match('by', method)
    if (!is.na(hasby)) {
        method[hasby] = "BY"
    }
    if (is.null(ntests)) {
        ntests = nrow(dta)
    }
    res = data.frame(dta[[dspvalues]])

    for (m in 1:length(method)) {
        res[, m+1] = p.adjust(p=dta[[dspvalues]], method=method[[m]], n=ntests)
    }
    names(res) = transnames[c("p", method)]
    
    # display
    StartProcedure(procname, omsid)
    res2 = res
    res2[is.na(res)] = gtxt("NA")  # for translation
    spsspivottable.Display(data.frame(res2), title=gtxt("Adjusted P Values"),
        templateName="PADJUSTPVALUES",
        outline=gtxt("Adjusted P Values"),
        rowlabels=as.character(dta[[dsnames]]),
        format=formatSpec.Significance,
        caption=gtxtf("Number of tests = %s\nResults computed by R p.adjust function", ntests)
    )

    if (plotit) {
        ordered = res[order(res[1]),]
        xlim = c(0, min(max(res[1], na.rm=TRUE) + .2,1))
        matplot(ordered[[1]], ordered, type="b", asp=1, lwd=2, cex=1,
            col=1:6, pch=0:5, xlim=xlim, xaxs="i",
            ylab=gtxt("Adjusted P values"),
            xlab=gtxt("Unadjusted P Values"),
            main = gtxt("P Values")
        )
        legend("bottomright", legend=names(res), col=1:6, 
            lty=1:6, inset=.02, cex=.8)
    }
    spsspkg.EndProcedure()
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
    lcl = list2env(lcl) # makes this list into an environment

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
            } else {
                procok = TRUE
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



Run = function(args) {
    #Execute the STATS DISAGG command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DSNAMES", subc="",  ktype="existingvarlist", var="dsnames"),
        spsspkg.Template("DSPVALUES", subc="", ktype="existingvarlist", var="dspvalues"),
        spsspkg.Template("NAMELIST", subc="", ktype="literal", var="namelist", islist=TRUE),
        spsspkg.Template("PLIST", subc="", ktype="float", var="plist", islist=TRUE),
        spsspkg.Template("METHOD", subc="", ktype="str", var="method", islist=TRUE,
            vallist=list("fdr", "holm", "hochberg", "hommel", "by", "bonferroni")),
        spsspkg.Template("NTESTS", subc="", ktype="int", var="ntests"),

        spsspkg.Template("PLOT", subc="OUTPUT",  ktype="bool", var="plotit")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dopadjust")
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