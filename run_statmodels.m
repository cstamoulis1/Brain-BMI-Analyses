%Code associated with publication: 
%Brooks, SJ, Smith, C, Stamoulis, C, Excess BMI in Early Adolescence Adversely Impacts
% Maturating Functional Circuits Supporting High-Level Cognition and Their Structural 
% Correlates, Int J Obesity, 2023 (in press, March 2023)
%
%
%   PURPOSE: run (fMRI) whole brain, network, and node level models (brain 
%   as outcome), obesity-sMRI-fMRI mediation models, and cognitive task models
%
%   REQUIRES: 
%   regression_multi_response.m (run models)
%   get_mdl_stats.m (get statistics)
%   Node_Groups.mat (FDR correction for node models)
%   Structure_Groups.mat (FDR correction for mediation models)
%
%   Last modified: March 22, 2023
%% SETUP

%load table of controls/covariates/predictor of interest as pred_tbl

%pred_interest = string of weight variable name: 'bmi', etc

%load response variables:
%whole brain: propstable = table of whole-brain properties

%network: netprops: struct, one field per network, each containing a
%table of network properties

%node: nodeprops: struct, one field per property, each contains a table of
%node values

%mediation set: 
%   props_fMRI_med = struct, one field per node-level property, each with 
%       table of median of node level property values within structural regions. 
%   props_sMRI = struct, one field per structural measure, each with table
%       of values for each structure
%   note: structures should be in same order across both sets of props

%cognitive task: cog_tasks = table of cognitive task results 

%Note: remove underweight subjects from models and make sure rows are the same
%across all tables

%% WHOLE BRAIN

[wb_stats, wb_mdls] = regression_mdls_multi_response(pred_tbl, propstable, pred_interest, 1);

%% NETWORK LEVEL

net_stats = struct;
net_lst = fieldnames(netprops); %list of networks

for n=1:length(net_lst) 
    [net_stats.(net_lst{n}),~] = regression_mdls_multi_response(pred_tbl, netprops.(net_lst{n}), pred_interest);
end

%% NODE LEVEL

load('Node_Groups.mat', 'node_groups'); %for FDR correction

node_stats = struct;
prop_lst = fieldnames(nodeprops); %list of node properties

for p=1:length(prop_lst) 
    [node_stats.(prop_lst{p}),~] = regression_mdls_multi_response(pred_tbl, nodeprops.(prop_lst{p}), pred_interest, 0, node_groups.group_inds);
end

%% MEDIATION SET
%general setup across paths:

%FDR correction:
%for desikan-killiany based results, use 'des_kil_groups' 
%for destrieux based results, use 'destrieux_groups' (used for path B results
%comparison only)
FDR_tbl = load('Structure_Groups.mat', 'des_kil_groups', 'destrieux_groups').des_kil_groups;
FDR_arr = FDR_tbl.group_inds;

%fMRI (node) and sMRI properties lists:
sMRI_prop_lst = fieldnames(props_sMRI);
fMRI_prop_lst = fieldnames(props_fMRI_med);

%structures list:
struct_lst = props_sMRI.(sMRI_prop_lst{1}).Properties.VariableNames;

%setup table of statistics (needed for paths C and D)
stat_cols = {'regression_coef', 'beta_CI_low', 'beta_CI_up',...
    'std_beta', 'std_beta_CI_low', 'std_beta_CI_up',...
    'SE', 'Wald', 'pvalue','adjusted_pvalue','Intercept_pvalue','model_pvalue'};
stat_lst = setdiff(stat_cols, 'adjusted_pvalue', 'stable'); 

%% obesity --> fMRI median properties (median within structure) (path A)
fMRI_med_stats = struct;
for f=1:length(fMRI_prop_lst)
    [fMRI_med_stats.(fMRI_prop_lst{f}),~] = regression_mdls_multi_response(pred_tbl, ...
        props_fMRI_med.(fMRI_prop_lst{f}), pred_interest, 0, FDR_arr);
end

%% obesity --> sMRI properties (path B)
%these models do not need percent censored:
pred_tbl_no_pc = pred_tbl;
pred_tbl_no_pc(:, 'percent_censored') = [];

sMRI_stats = struct;
sMRI_mdls = struct;
for s=1:length(sMRI_prop_lst)
    [sMRI_stats.(sMRI_prop_lst{s}), sMRI_mdls.(sMRI_prop_lst{s})] = regression_mdls_multi_response(pred_tbl_no_pc, ...
        props_sMRI.(sMRI_prop_lst{s}), pred_interest, 1, FDR_arr);
end

%% sMRI property --> fMRI median (used to subset results of sobel test)
% (path C)
%these models do not include obesity:
pred_tbl_no_bmi = pred_tbl;
pred_tbl_no_bmi(:, pred_interest) = [];

sMRI_fMRI_stats = struct;
for s=1:length(sMRI_prop_lst)
    for f=1:length(fMRI_prop_lst)
        %setup table
        stat_tbl= array2table(NaN(length(struct_lst), length(stat_cols)),...
            'RowNames', struct_lst, 'VariableNames', stat_cols);
        %run models for each structure and fill table
        for st=1:length(struct_lst)
            mdl = fitlm([pred_tbl_no_bmi, props_sMRI.(sMRI_prop_lst{s})(:, st), props_fMRI_med.(fMRI_prop_lst{f})(:, st)], 'linear'); 
            stat_tbl(st, stat_lst) = get_mdl_stats(mdl, width(pred_tbl_no_bmi)+2, stat_lst);
        end
        %FDR correction
        for i=1:length(FDR_arr) %for each group of indices, apply FDR correction
            inds = FDR_arr{i,1};
            stat_tbl.adjusted_pvalue(inds) = mafdr(stat_tbl.pvalue(inds), 'BHFDR', true);
        end
        %store results
        sMRI_fMRI_stats.(sMRI_prop_lst{s}).(fMRI_prop_lst{f}) = stat_tbl;
    end
end

%% obesity + sMRI property --> fMRI median (path D)
%indices of obesity variable and structural variable in predictors: 
ind_obesity = find(strcmp(pred_interest, pred_tbl.Properties.VariableNames)) +1;
ind_sMRI = width(pred_tbl) +2;

obesity_sMRI_fMRI_stats = struct; %will contain results relating to both 
%obesity variable and structural property 
for s=1:length(sMRI_prop_lst)
    for f=1:length(fMRI_prop_lst)
        %setup tables
        stat_tbl_obesity= array2table(NaN(length(struct_lst), length(stat_cols)),...
            'RowNames', struct_lst, 'VariableNames', stat_cols);
        stat_tbl_sMRI= array2table(NaN(length(struct_lst), length(stat_cols)),...
            'RowNames', struct_lst, 'VariableNames', stat_cols);
        %run models for each structure and fill tables
        for st=1:length(struct_lst)
            mdl = fitlm([pred_tbl, props_sMRI.(sMRI_prop_lst{s})(:, st), props_fMRI_med.(fMRI_prop_lst{f})(:, st)], 'linear'); 
            stat_tbl_obesity(st, stat_lst) = get_mdl_stats(mdl, ind_obesity, stat_lst);
            stat_tbl_sMRI(st, stat_lst) = get_mdl_stats(mdl, ind_sMRI, stat_lst);
        end
        %FDR corrections
        for i=1:length(FDR_arr) %for each group of indices, apply FDR correction
            inds = FDR_arr{i,1};
            stat_tbl_obesity.adjusted_pvalue(inds) = mafdr(stat_tbl_obesity.pvalue(inds), 'BHFDR', true);
            stat_tbl_sMRI.adjusted_pvalue(inds) = mafdr(stat_tbl_sMRI.pvalue(inds), 'BHFDR', true);
        end
        %store results
        obesity_sMRI_fMRI_tats.(sMRI_prop_lst{s}).(fMRI_prop_lst{f}).obesity = stat_tbl_obesity;
        obesity_sMRI_fMRI_stats.(sMRI_prop_lst{s}).(fMRI_prop_lst{f}).sMRI = stat_tbl_sMRI;
    end
end

%% Sobel test:
% reformat above models according to mediation_sobel_test.m:
%   one row per structure (labeled by structure name)
%   (nonsignificant results should be removed)
%   columns for relevant statistics, including 'adjusted_pvalue' and 
%   'regression_coef'

sobel_results_t = struct;
for s=1:length(sMRI_prop_lst)
    for f=1:length(fMRI_prop_lst)
        %for this combination of properties: 
        
        %load path A: obesity --> fMRI_med
        %load path B: obesity --> sMRI
        %load path C: sMRI --> fMRI_med
        %load path D_pred: *obesity* + sMRI --> fMRI_med
        %load path D_med: obesity + *sMRI* --> fMRI_med
        
        sobel_results_t.(sMRI_prop_lst{s}).(fMRI_prop_lst{f}) = mediation_sobel_test(pathA, pathB, pathC, pathD_pred, pathD_med, 1);
    end
end
 
%% COGNITIVE TASK OUTCOMES
%these models do not need to include percent of frames censored for motion
%since this variable is a relevant adjustment for models including brain parameters
pred_tbl_no_pc = pred_tbl;
pred_tbl_no_pc(:, 'percent_censored') = [];

[cog_stats, cog_mdls] = regression_mdls_multi_response(pred_tbl_no_pc, cog_tasks, pred_interest, 1);
