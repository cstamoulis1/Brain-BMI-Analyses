Description of codes associated with publication:
Skylar J. Brooks, Calli Smith, Catherine Stamoulis
"Excess BMI in Early Adolescence Adversely Impacts Maturating Functional 
Circuits Supporting High-Level Cognition and Their Structural Correlates"
International Journal of Obesity (in press March 2023)
Codes last modified: March 22, 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILES:

run_statmodels.m: 
Script to load variables and run appropriate models

regression_mdls_multi_response.m: 
Function used to run linear (or logistic) regression for a given set of 
controls/covariates/predictors (including a single predictor of interest) 
across multiple response variables. Produces a table of statistics 
(FDR correction is performed as requested) and a map of models.

get_mdl_stats.m:
Function used to create a table of common statistics for a single model object 
and predictor of interest (called by regression_mdls_multi_response.m).

mediation_sobel_test.m:
Uses the Sobel test to estimate the significance of mediation, 
and consequently whether the relationship of interest is partially 
or fully mediated by the variable being evaluated. Applied to obesity sMRI fMRI 
models obtained using run_statmodels.m.

Node_Groups.mat:
In functional network analyses at the node level, 
this file is needed for adjusting p-values for the False Discovery Rate (FDR).
In the case of node properties, FDR corrections are used across nodes belonging to a 
particular network. Nodes are grouped according to the networks identified 
in Yeo et al, 2011

Structure_Groups.mat:
In structural region analyses, this file is needed for adjusting p-values for FDR.
FDR is applied across regions within a particular network, based on networks 
identified in Yeo et al, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VARIABLES REQUIRED:

Independent variables: 
1) Control and covariate variables: propensity weights, sex, age, race, ethnicity, 
family income, sleep duration (SDS), physical activity (YRB), screen time (STQ), 
and (for models with fMRI data) percent of frames censored for motion

2) Predictor of Interest: bmi (raw or zscored), or weight contrast

Dependent Variables: 

In primary models: rs-fMRI topological network properties for whole brain, 
network, and node scales. 
In this study, fMRI analyses were conducted using participants' a) best and 
second best run from release R2.0.1, b) median values
across all available runs from R2.0.1, c) best run from R4.0.

Whole brain and network level properties: Efficiency, global clustering, 
median connectivity, modularity (Newman), small-worldness, robustness, 
stability. 
Node level: degree, centrality, local clustering.

In mediation models: rs-fMRI median of node-level properties within 
structural regions, and structural MRI values (based on provided Destrieux 
and Desikan-Killiany values of cortical thickness, cortical volume, and 
white matter intensity)

In cognitive task models: cognitive task responses (NIH toolbox battery)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RUNNING MODELS: run_statmodels.m

After loading needed variables, run section-by-section to produce 
appropriate models and statistics.

Uses regression_mdls_multi_response and get_mdl_stats to run regression models
and calculate appropriate statistics for each set of models. 

Uses Node_Groups and Structure_Groups for appropriate FDR correction, based on node
or structure grouping according to the networks identified in Yeo et al, 2011. 




