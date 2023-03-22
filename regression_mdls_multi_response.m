function [stat_tbl, mdl_map] = regression_mdls_multi_response(pred_tbl, resp_tbl, pred_interest, store_mdls, FDR_arr)
%   Code associated with publication: 
%   Brooks, SJ, Smith, C, Stamoulis, C, Excess BMI in Early Adolescence Adversely 
%   Impacts Maturating Functional Circuits Supporting High-Level Cognition 
%   and Their Structural Correlates, Int J Obesity, 2023 (in press, March 2023)

%   PURPOSE: Code for running statistical models (and producing model statistics)
%   for a set of models with the same control and covariate variables and 
%   and a single predictor of interest (BMI/BMI status in this study) 
%   Models are run across a set of dependent variables 
%   (topological node properties, network or connectome properties, or cognitive 
%   task outcomes)

%   INPUT:
%   pred_tbl: table of control and covariate variables and predictor of interest 
%   resp_tbl: table of dependent variables (nodal properties) 
%   note: must have same row names (participants)
%   pred_interest: string name of single predictor of interest
%   store_mdls: optional boolean: 1 = save models in map container, 0 = do
%   not save (default)
%   FDR_arr: optional cell array: 1 column, each row contains a vector of 
%   indices across which FDR correction is applied (each index should appear 
%   once). If not included, will automatically correct across all response 
%   variables
%
%   OUTPUT: 
%   stat_tbl: table of statistics for each model (typically, one row per
%   node)
%   mdl_map: optional map container of models
%
%   Last modified: March 22, 2023

%% setup variable and output structures

if nargin<4 %by default, do not save models
    store_mdls=0;
end

if store_mdls %setup storage to save mdls, if needed
    responsenames = resp_tbl.Properties.VariableNames;   
    mdl_arr = cell(1, width(resp_tbl)); %array of mdls
end

%setup table of statistics
stat_cols = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
    'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
    'SE', 'Wald', 'pvalue','adjusted_pvalue','Intercept_pvalue','model_pvalue'};
stat_lst = setdiff(stat_cols, 'adjusted_pvalue', 'stable'); 

stat_tbl= array2table(NaN(width(resp_tbl), length(stat_cols)),...
    'RowNames', resp_tbl.Properties.VariableNames, 'VariableNames', stat_cols);

%find index of predictor of interest in model, +1 for indexing
%into coefficient table (the first row corresponds to the coefficient for the intercept)
ind = find(strcmp(pred_interest, pred_tbl.Properties.VariableNames)) +1;

%get which responses are binary: (binary=1)
num_values = varfun(@unique, rmmissing(resp_tbl), 'OutputFormat', 'cell');
bin_resp = cellfun(@length, num_values)==2;

%% run models, insert stats:

for r=1:width(resp_tbl)
    
    %model for this response variable
    if bin_resp(r) %if binary, use logistic regression
        mdl = fitglm([pred_tbl resp_tbl(:,r)],'Distribution', 'binomial','link','logit');
    else %use linear regression
        mdl = fitlm([pred_tbl resp_tbl(:,r)], 'linear'); 
    end
    
    %get statistics and store in table (not adjusted pvalue - calculated
    %later)
    stat_tbl(r,stat_lst) = get_mdl_stats(mdl, ind, stat_lst); 
    
    %store model, if needed
    if store_mdls 
        mdl_arr{1,r} = mdl;
    end
end

%% perform FDR correction
if nargin <5 %if FDR array is not provided, correct across all responses
    FDR_arr = {1:width(resp_tbl)};
end

for f=1:length(FDR_arr) %for each group of indices, apply FDR correction
    inds = FDR_arr{f,1};

    stat_tbl.adjusted_pvalue(inds) = mafdr(stat_tbl.pvalue(inds), 'BHFDR', true);
end

%% convert mdl_arr to map container, if needed
if store_mdls
    mdl_map = containers.Map(responsenames, mdl_arr); 
else
    mdl_map = NaN;
end

end

