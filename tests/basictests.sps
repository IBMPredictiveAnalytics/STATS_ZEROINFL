* Encoding: UTF-8.
file handle data /name="c:/cc/misc2/extensions/r/stats_zeroinfl/tests".

get file="data/bioChemists.sav".
dataset name bio.
compute countoffset=1.
compute zerooffset=.3.
variable level countoffset zerooffset(scale).
compute randomcountoffset = rv.bernoulli(.5).
variable level randomcountoffset(scale).
compute id = -$casenum.

extension /SPECIFICATION    command="c:/extcommon/stats_zeroinfl.xml".

stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
  zeromodel=fem  mar  kid5 phd  ment.

* estimation dataset.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
  zeromodel=fem  mar  kid5 phd  ment
/save dataset=results id=id workspaceoutfile="c:/temp/zero.Rdata".

* equal 0/non prob for all.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment.

* explicit.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
countoffset=countoffset
sameregressors=no.

* factor offset.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
countoffset=randomcountoffset
sameregressors=no.

* different families.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
countdist=negbin zerolink=probit.

STATS ZEROINFL MODELSOURCE=ESTIMATE DEPENDENT=art COUNTMODEL=fem kid5 mar ment randomcountoffset
SAMEREGRESSORS=YES COUNTOFFSET=countoffset ZEROOFFSET=zerooffset
COUNTDIST=POISSON ZEROLINK=LOGIT
/OPTIONS STARTVALUES=GENLIN OPTMETHOD=BFGS. 



* predict from workspace file.
stats zeroinfl modelsource=file modelfile="c:/temp/zero.Rdata"
/save dataset=pred id=id.

* predict from retained workspace.
stats zeroinfl dependent=art countmodel=fem mar kid5  phd  ment
  zeromodel=fem  mar  kid5 phd  ment
/save workspaceaction=retain dataset=est.

stats zeroinfl modelsource=workspace
/save dataset=predfromworkspace workspaceaction=retain.

DATASET ACTIVATE bio.
STATS ZEROINFL MODELSOURCE=ESTIMATE DEPENDENT=art COUNTMODEL=fem kid5 mar ment phd
SAMEREGRESSORS=NO ZEROMODEL=mar 
COUNTDIST=POISSON ZEROLINK=CAUCHIT
/OPTIONS STARTVALUES=EM OPTMETHOD=CG MAXITER=1000 TOL=0.0000000001
/SAVE WORKSPACEACTION=CLEAR.

stats zeroinfl /help.

* another dataset.
* results match the reference very closely (but constant term is a little different).
get file="data/DebTrividi.sav".
dataset name DebTrividi.

DATASET ACTIVATE DebTrivedi.
STATS ZEROINFL MODELSOURCE=ESTIMATE DEPENDENT=ofp COUNTMODEL=hosp health numchron gender school 
    privins
SAMEREGRESSORS=YES 
COUNTDIST=NEGBIN ZEROLINK=LOGIT
/OPTIONS STARTVALUES=GENLIN OPTMETHOD=BFGS MAXITER=1000 TOL=0.0000000001
/SAVE WORKSPACEACTION=CLEAR.

STATS ZEROINFL MODELSOURCE=ESTIMATE DEPENDENT=ofp COUNTMODEL=hosp health numchron gender school 
    privins
SAMEREGRESSORS=NO 
COUNTDIST=NEGBIN ZEROLINK=LOGIT
/OPTIONS STARTVALUES=GENLIN OPTMETHOD=BFGS MAXITER=1000 TOL=0.0000000001
/SAVE WORKSPACEACTION=CLEAR.


* Generalized Linear Models with Tweedie distribution for comparison.
GENLIN ofp BY health gender privins (ORDER=ASCENDING) WITH hosp numchron school
  /MODEL health gender privins hosp numchron school INTERCEPT=YES
 DISTRIBUTION=TWEEDIE(1.5) LINK=LOG
  /CRITERIA METHOD=FISHER(1) SCALE=MLE COVB=MODEL MAXITERATIONS=100 MAXSTEPHALVING=5 
    PCONVERGE=1E-006(ABSOLUTE) SINGULAR=1E-012 ANALYSISTYPE=3(WALD) CILEVEL=95 CITYPE=WALD 
    LIKELIHOOD=FULL
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION.
