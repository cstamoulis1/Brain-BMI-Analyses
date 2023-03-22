function mdl_stats = get_mdl_stats(mdl, ind, stat_lst)
%   PURPOSE: get a set of common statistics for a model/predictor of interest
%
%   INPUTS:
%   mdl: linear or logistic model object
%   ind: index of predictor of interest in coefficient table (+1 to account
%   for intercept row)
%   stat_lst: optional list of statistics you would like returned
%   (otherwise, will default to return all statistics)
%
%   OUTPUTS:
%   mdl_stats: single table row of (requested) statistics for this model
%
%   Last modified: March 22, 2023

%all available statistics
stats_all = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
    'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
    'SE', 'Wald', 'pvalue','Intercept_pvalue','model_pvalue'};

%setup stat storage
mdl_stats= array2table(NaN(1, length(stats_all)), ...
    'RowNames', {mdl.ResponseName}, 'VariableNames', stats_all);

%get statistics from coefficient table
mdl_stats{1, 'regression_coef'} = mdl.Coefficients.Estimate(ind);
cis = mdl.coefCI;
mdl_stats{1, {'beta_CI_low', 'beta_CI_up'}} = cis(ind, :);
mdl_stats{1, 'SE'} = mdl.Coefficients.SE(ind);
tstat_val = mdl.Coefficients.tStat(ind);
mdl_stats{1, 'Wald'} = tstat_val^2;
mdl_stats{1, 'pvalue'} = mdl.Coefficients.pValue(ind);
mdl_stats{1, 'Intercept_pvalue'} = mdl.Coefficients.pValue(1);
mdl_stats{1,'model_pvalue'} = coefTest(mdl); 

%standardize beta value and confidence interval if both predictor of
%interest and response are not binary:
bin_resp = length(unique(rmmissing(mdl.Variables{:, end})))==2;
bin_pred = length(unique(rmmissing(mdl.Variables{:, ind-1})))==2;
if ~bin_resp && ~bin_pred
    st_dev_resp = std(mdl.Variables{:, end}, 'omitnan'); %standard deviation of response variable
    st_dev_pred = std(mdl.Variables{:, ind-1}, 'omitnan'); %standard deviation of predictor of interest

    st_beta = mdl.Coefficients.Estimate(ind)*(st_dev_pred/st_dev_resp);
    mdl_stats{1, 'std_beta'} = st_beta;

    std_cis= cis(ind,:) .* (st_dev_pred/st_dev_resp);
    mdl_stats{1,{'std_beta_CI_low', 'std_beta_CI_up'}} = std_cis; 
end

%keep only what is requested
if nargin==3
    mdl_stats = mdl_stats(:, stat_lst);
end

end
