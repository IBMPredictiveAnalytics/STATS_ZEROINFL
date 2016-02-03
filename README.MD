# STATS ZEROINFL
## Estimate and predict a zero-inflated count model.
Estimate a mixture model with a Poisson or negative binomial count model and a point mass at zero. The predictors can be different for the two models. The estimated model can be saved and used for predictions on new data.

---
Requirements
----
- IBM SPSS Statistics 18 or later and the corresponding IBM SPSS Statistics-Integration Plug-in for R.

Note: For users with IBM SPSS Statistics version 23 or higher, the STATS ZEROINFL extension is installed as part of IBM SPSS Statistics-Essentials for R.

---
Installation intructions
----
1. Open IBM SPSS Statistics
2. Navigate to Utilities -> Extension Bundles -> Download and Install Extension Bundles
3. Search for the name of the extension and click Ok. Your extension will be available.

---
Tutorial
----

STATS ZEROINFL  
MODELSOURCE = ESTIMATE^&#42;&#42; or FILE or WORKSPACE  
MODELFILE = "file"  
DEPENDENT = variable  
COUNTMODEL = variables  
ZEROMODEL = variables  
SAMEREGRESSORS = YES or NO^&#42;&#42;  
COUNTOFFSET = variable  
ZEROOFFSET = variable  
COUNTDIST = POISSON^&#42;&#42; or NEGBIN or GEOMETRIC  
ZEROLINK = LOGIT^&#42;&#42; or PROBIT or CLOGLOG or CAUCHIT or LOG

/OPTIONS  
STARTVALUES = GENLIN^&#42;&#42; or EM  
OPTMETHOD = BFGS^&#42;&#42; or NELDERMEAD or CG or LBFGSB or SANN or BRENT  
MAXITER = number  
TOL = number  
MISSING = EXCLUDE^&#42;&#42; or FAIL

/SAVE  
DATASET=dataset  
ID = variable  
WORKSPACEACTION = CLEAR^&#42;&#42; or RETAIN  
WORKSPACEOUTFILE = "file"

/HELP

^&#42; Required  
^&#42;&#42; Default

Example:
```

STATS ZEROINFL DEPENDENT=y COUNTMODEL = x1 x2 x3
ZEROMODEL= w1 COUNTDIST= NEGBIN
/SAVE DATASET=results ID=id WORKSPACEOUTFILE="C:/myproject/zero.Rdata".
```
STATS ZEROINFL /HELP displays this help and does nothing else.

This command operates in two modes.  It can be used to estimate
a model, or it can be used with a model previously estimated
and saved in memory or a file to make predictions for new data.

**MODELSOURCE** specifies the model source.  The choices are
* ESTIMATE: estimate a model using data in the current active dataset
* WORKSPACE: use the model in the currently loaded workspace.
The model must come from a previously run command in
the current session with
MODELSOURCE=ESTIMATE and retained in memory with
WORKSPACEACTION=RETAIN.
* FILE: use a previously created model that was saved
to a file with WORKSPACEOUTFILE.

**MODELFILE** specifies the file containing the model
and is required if MODELSOURCE=FILE is used.

**DEPENDENT** specifies the dependent variable, which
is expected to contain counts.  It is required in
estimation mode.  It must have a scale measurement level.

**COUNTMODEL** specifies one or more independent variables
for the count and is required in estimation mode.

**ZEROMODEL** can specify a different set of predictors
for the zero/nonzero portion of the model.  If omitted,
the model depends on the setting for SAMEREGRESSORS.

**SAMEREGRESSORS** specifies whether the predictors
for the zero model are the same as for the count model.
If it is NO, and no ZEROMODEL variables are supplied,
the zero model is the same for all cases, i.e., the
zero model predictors consist just of a constant term.

**COUNTOFFSET** can specify an offset for the count model.

**ZEROOFFSET** can specify an offset for the zero model.
The SAMEREGRESSORS setting does not apply to the zerooffset
specification.

**COUNTDIST** specifies the distribution for the count model.
The choices are
* POISSON: Poisson
* NEGBIN: Negative binomial
* GEOMETRIC: Negative binomial with a size (number of failures) parameter of 1.

**ZEROLINK** specifies the link function for the binary
zero-inflation model.  The choices are
* LOGIT: logistic
* PROBIT: probit
* CLOGLOG: Complementary log log
* CAUCHIT: Cauchy quantiles
* LOG: log

OPTIONS
-------
**STARTVALUES** specifies whether the intial values for the
iterative estimation algorithm come from a generalized linear
model or from an EM (estimation maximization) model.  The EM
algorithm uses a generalized linear model first step, so it
is possible to get a failure message referring to glm when
EM is the setting.

**OPTMETHOD** specifies the iterative solving algorithm.
The choices are
* BFGS: Quasi Newton
* Nelder-Mead: A robust but slow method that does not
use gradients
* CG: Conjugate gradient
* SANN: variant of simulated annealing

**MAXITER** specifies the maximum number of iterations
and defaults to 1000.

**TOL** specifies the convergence tolerance and typically defaults to
1.64 &#42; 10^-10.

**MISSING** specifies whether cases with missing values are
omitted or cause the procedure to stop.

SAVE
----
**DATASET** specifies the output dataset name.  The name
must not already be in use.  For estimation, the dataset
will contain the fitted values and residuals.  For prediction,
it will contain the predicted values and is required.

**ID** optionally specifies an id variable to be
included in the dataset.  If one is not specified, the
cases are sequentially numbered.

**WORKSPACEACTION** specifies whether to clear the
in-memory workspace after the procedure (CLEAR) or
to retain it for following prediction calculations.

**WORKSPACEOUTFILE** specifies whether to save the
estimated model to a file.  It only applies to
estimation model.

REFERENCE
---------

---
License
----

- Apache 2.0
                              
Contributors
----

  - JKP, IBM SPSS
