function t_tbl = mediation_sobel_test(pathA, pathB, pathC, pathD_pred, pathD_med, num_method)
%   PURPOSE: assess partial or complete mediation using the Sobel test 
%
%   INPUT: 
%   pathA: results of X--> Y (Models1) 
%   pathB: results of X --> M (Models2) 
%   pathC: results of M --> Y 
%   pathD_pred: results relating to X in X + M --> Y (Models3)
%   pathD_med: results relating to M in X + M -- Y (Models3)
%   num_method: numerator method (1 is more common)
%
%   Note on results table formatting:
%   one row per mediator (nonsignificant results should be removed)
%   columns for relevant statistics, including 'adjusted_pvalue' and 
%   'regression_coef'
%
%   OUTPUT: 
%   t_tbl: t value and partial (0) or full (1) mediation for each structure
%
%   SOBEL TEST:
%   X = independent variables 
%   M = mediator variables 
%   Y = dependent variables 
%   Models1: X-->Y (path A)
%   Models2: X-->M (path B)
%   Models3: X+M-->Y (path D)
%
%   t = numerator / SE_pooled 
%
%   Numerator method options: 
%       1 = a*B (product of coefficients)
%       2 = (T-T')
%   a = unstandardized coefficient of X in model2
%   B = unstandardized coefficient of M in model3
%   T = unstandardized coefficient of X in model1
%   T' = unstandardized coefficient of X in model3
%
%   methods 1 and 2 are equivalent only with equivalent samples used in
%   both models (method 1 is more common)
%   
%   SE_pooled = sqrt(a^2*SE_B^2 + B^2*SE_a^2)
%   
%   t = zscore (on normal distribution), so if abs(t) > 1.96, then p<0.05,
%   so there is a significant mediation effect
%
%   results are subset by whether all involved coefficients are significant:
%   X in models1 and models2, M in models3
%   
%   complete mediation (1) is when coefficients of X in models3 are 
%   non-significant, and partial mediation (0) is when coefficients of X in 
%   models3 are significant and different from coefficients of X in models1 
%
%   Notes:
%   See davidakenny.net/cm/mediate.htm for more information.
%
%   Last modified: March 22, 2023
    

%% run sobel test:
if ~(isempty(pathA) || isempty(pathB) || isempty(pathD_med) || isempty(pathC))
    % first, get the structures that have significant associations for X in paths A and B
    sig_resps = intersect(pathA.Properties.RowNames, pathB.Properties.RowNames);

    %this is not required by the sobel test, but we are requiring 
    %M significant in path C:
    sig_resps = intersect(sig_resps, pathC.Properties.RowNames); 

    %M must be significant in path D:
    sig_resps = intersect(sig_resps, pathD_med.Properties.RowNames);

    if ~isempty(sig_resps) %if the conditions for significance are met, 
        %determine partial or full mediation:

        %setup storage:
        t_tbl = array2table(NaN(length(sig_resps),2), ...
            'VariableNames', {'t_value', 'partial_full'},...
            'RowNames', sig_resps);

        for s = 1:length(sig_resps) %for each significant structure
            %get a, B, T, Tprime
            a = pathB.regression_coef(sig_resps{s});
            a_SE = pathB.SE(sig_resps{s});
            B = pathD_med.regression_coef(sig_resps{s});
            B_SE = pathD_med.SE(sig_resps{s});

            T = pathA.regression_coef(sig_resps{s});
            if ~isempty(pathD_pred) && ismember(sig_resps{s}, pathD_pred.Properties.RowNames) %X may or may not be significant in path D
                Tprime = pathD_pred.regression_coef(sig_resps{s});
            else
                warning('Tprime not available: num_method 1 will be used');
                num_method=1;
            end

            %calculate selected statistic and insert into table
            denominator = sqrt((a^2)*(B_SE^2 )+ (B^2)*(a_SE^2));

            if num_method==1  
                numerator = a* B;

            elseif num_method==2 
                %can't be used if path D X is not significant bc has already 
                %been subset for significance: use num_method 1 instead, or adjust 
                %to use unsubset data, and check significance as you go
                numerator = T - Tprime;
            end

            %insert t statistic
            t_stat = numerator/denominator;
            t_tbl{s,1} = t_stat;

            % if t is significant: determine partial or full mediation:
            %full mediation: X is not significant in path D
            %partial mediation: X is significant but different coefficient in
            %path D from path A
            if abs(t_stat)>1.96
                if isempty(pathD_pred) || ~ismember(sig_resps{s}, pathD_pred.Properties.RowNames)
                    t_tbl{s,2} = 1; %full
                else
                    if T ~= Tprime
                        t_tbl{s,2} = 0; %partial
                    end
                end   
            end %else, leave as NaN
        end
    else
        t_tbl = NaN; %no mediation
    end 
else
    t_tbl = NaN; %no mediation
end

    
end
    