clear;

%% Select dataset (CLHLS or HRS) and model (static, trend, or frailty)
% The '_resid' suffix refers to the model that includes the residence
% covariate -- only used for the CLHLS dataset

% --------------
% Select dataset
% --------------

% 1. CLHLS
dataName = 'clhls';

% % 2. HRS
% dataName = 'rndhrs';


% --------------
% Select model
% --------------

% 1. static / static_resid, 2. trend / trend_resid, 3. frailty / frailty_resid
model = 'frailty';

% For the frailty model, we consider two (2) cases for the latent process
% (1) Random walk and (2) AR(1)
if any(strcmp(model, {'frailty', 'frailty_resid'}))
    % Model 3.1. random walk
    latent = 'rw';

    % Model 3.2. auto regressive
%     latent = 'ar';
else
    latent = '';
end


%% Specify a model object

datatable = readtable([dataName, '_', 'transit.csv']);

Mdl = MLFIMdl;

Mdl.dataName = dataName;
Mdl.latent = latent;
Mdl.model = model;
Mdl.transitTable = readtable('transition.csv');
Mdl.S = datatable;
Mdl.t = datatable;
Mdl.N_H_STATE = datatable;

% % ---- No. of frailty paths used in estimation ----
% % - default is 1000
% % - set to a small number when testing the code
% Mdl.N_SIM = 5;


%% Estimate theta

% % Uncomment to estimate the parameters
% Mdl = estTheta(Mdl, datatable);

% If the model has been estimated and the estimation results have been
% saved, skip this section


%% Import estimation results from .mat file

Mdl = importResult(Mdl);


%% Print the estimated parameters 
% % with standard errors and significance levels in latex table format

% % Uncomment to print the LaTeX table
% printTable(Mdl);


%% Calculate infomation criteria, AIC and BIC

% % Uncomment to calculate AIC and BIC
% n = size(datatable, 1);
% [aic, bic] = calAIC(Mdl, n);


%% Recover the frailty process

% % Uncomment to recover the frailty process
% % the mean and variance will be saved in a .mat file
% [psi_hat, kf_V] = recoverPsi(Mdl, datatable);


%% Calculate estimated transition rates
% transition rate of the **static** and **trend** model

if ~any(strcmp({'frailty', 'frailty_resid'}, Mdl.model))
    Input = struct;
    Input.ageBeg = 65; % 65 or 75
    Input.ageEnd = 110;
    Input.residence = 'u';
    if ~MLFIMdl.isResid(Mdl.model)
        Input.residence = '';
    end
    Input.type = 'cohort'; % 'cohort' or 'period'
    
    Input.time_ageStart = 1998; 
    Input.iPsiPath = [];
    
    Input.gender = 'f'; 
    trsRateF = calTrsRate(Mdl, Input);
    
    Input.gender = 'm'; 
    trsRateM = calTrsRate(Mdl, Input);
    
    % save as a mat file --> to be used in main_plot.m
    matDir = './plotFiles/';
    if strcmp(model, 'static')
        tmp = [matDir, Mdl.dataName, '_', 'rate', '_', Mdl.matName, '.mat'];
    elseif strcmp(model, 'static_resid')
        tmp = [matDir, Mdl.dataName, '_', 'rate', '_', Mdl.matName, '_', Input.residence, '.mat'];
    elseif MLFIMdl.isResid(Mdl.model)
        tmp = [matDir, Mdl.dataName, '_', 'rate', '_', Mdl.matName, '_', Input.residence, '_', num2str(Input.time_ageStart), '.mat'];
    else
        tmp = [matDir, Mdl.dataName, '_', 'rate', '_', Mdl.matName, '_', num2str(Input.time_ageStart), '.mat'];
    end
    save(tmp, 'trsRateF', 'trsRateM')
end


%% Simulate future lifetime - micro-simulation

% Specify the cohort information
Input = struct;

% Maximum attainable age is 110
Input.ageEnd = 110;

% Residence: None or urban or rural 
% '' or 'r' or 'u'
Input.residence = 'r';
if ~MLFIMdl.isResid(Mdl.model)
    Input.residence = '';
end

% Simulate cohort life expectancy or period life expectancy 
% 'cohort' or 'period'
Input.type = 'cohort'; 

% No. of simulation
nSimObs = 10000;

genderList = {'f', 'm'};

% -----------------------------------------
% micro-simulation: statis + trend models
%   **skip if using the frailty model**
% -----------------------------------------
Input.iPsiPath = [];

% Starting year of the cohort
% 1998 or 2014
Input.time_ageStart = 1998;

% save the life expectancy results
out = table; 

% save the health distribution results
allProp = struct; aliveProp = struct; survProp = struct;
for ageBeg = [65, 75]
    for g = 1:numel(genderList)
        
        Input.ageBeg = ageBeg; % 65 or 75
        Input.gender = genderList{g}; % 'f' / 'm'
        
        % Initial health state: healthy, i.e. 1
        initHVector = ones(1, nSimObs); 

        hStateSim = simHState(Mdl, Input, initHVector);
        
        hStateCount = MLFIMdl.countHState(hStateSim);
        hStateProp = hStateCount ./ nSimObs;
        hAliveProp = hStateCount(:, 1:end-1) ./ sum(hStateCount(:, 1:end-1), 2);
        allProp.([Input.gender, num2str(ageBeg)]) = hStateProp;
        aliveProp.([Input.gender, num2str(ageBeg)]) = hAliveProp;
        survProp.([Input.gender, num2str(ageBeg)]) = 1 - hStateProp(:, end);
        
        % future life time
        lifeTime = MLFIMdl.calT(hStateSim, initHVector);  

        % onset of disability
        firstX = MLFIMdl.calOnsetX(Input.ageBeg, hStateSim, 2);

        % HLE / TLE
        ratio = MLFIMdl.calRatio(hStateSim, initHVector);

        % total LE; healthy LE; disabled LE; HLE/TLE; entry age (x)
        % '0' is used as separation lines
        out.([Input.gender, num2str(ageBeg)]) = ...
            [mean(lifeTime{end}); std(lifeTime{end})/sqrt(nSimObs); std(lifeTime{end}); 0
            mean(lifeTime{1}); std(lifeTime{1})/sqrt(nSimObs); std(lifeTime{1}); 0
            mean(lifeTime{2}); std(lifeTime{2})/sqrt(nSimObs); std(lifeTime{2}); 0
            mean(ratio); std(ratio)/sqrt(nSimObs); std(ratio); 0
            mean(firstX); std(firstX)/sqrt(numel(firstX)); std(firstX)];
    end
end

% save the results --> to be used in main_plot.m
if MLFIMdl.isResid(Mdl.model)
    % _resid models
    save(['./plotFiles/', Mdl.dataName, '_', Mdl.model, '_', Input.residence, '_', num2str(Input.time_ageStart), '.mat'], ...
        'allProp', 'aliveProp', 'survProp')
else
    save(['./plotFiles/', Mdl.dataName, '_', Mdl.model, '_', num2str(Input.time_ageStart), '.mat'], ...
        'allProp', 'aliveProp', 'survProp')
end


% ----end of micro-simulation for statis + trend models----



% ---------------------------------------
% micro-simulation: frailty model
% ---------------------------------------
nSimPsi = 1000; % no. of frailty paths in each simulation

% save the simulation results to a folder
matDir = './simFiles/';

% Note
% To make use of parallel computing, 
% simulate each starting age and each gender separately, 
% i.e. avoid using loop for ageBeg and gender

% Starting age: 65 or 75
Input.ageBeg = 65; % 65 / 75

% Gender: female or male 
Input.gender = 'f'; % 'f' or 'm'

time_ageStart = [1998 2014];

% ---- run the simulation ----
for tIndex = 1:numel(time_ageStart)
    Input.time_ageStart = time_ageStart(tIndex); % 1998 or 2014
    
    path = simulatePath(Mdl, Input, 2014, nSimPsi);
    
    initHVector = ones(1, nSimObs); % 10,000 individuals
    hStateCount = cell(nSimPsi, 1);
    lifeTime = cell(nSimPsi, Mdl.N_H_STATE);
    firstX = cell(nSimPsi, 1);
    ratio = cell(nSimPsi, 1);
    
    for iSim = 1:nSimPsi
        disp(num2str(iSim))
        tic
        Input.iPsiPath = path(:, iSim);
        
        hStateSim = simHState(Mdl, Input, initHVector);
        
        % no. of individuals in each state
        hStateCount{iSim} = MLFIMdl.countHState(hStateSim);      
        
        % future life time
        lifeTime(iSim, :) = MLFIMdl.calT(hStateSim, initHVector);  

        % onset of disability
        firstX{iSim} = MLFIMdl.calOnsetX(Input.ageBeg, hStateSim, 2);

        % HLE/TLE
        ratio{iSim} = MLFIMdl.calRatio(hStateSim, initHVector);
        toc
    end

    % save lifeTime and firstX
    save([matDir, Mdl.dataName, '_', Input.residence, Input.gender, '_', ...
        num2str(Input.ageBeg), '_', num2str(Input.ageEnd), '_', ...
        num2str(Input.time_ageStart), '.mat'], ...
        'hStateCount', 'lifeTime', 'firstX', 'ratio')

end


% ---- Analyse the simulation results of the frailty model ----
Input.time_ageStart = 1998; % 1998 or 2014

% store the life expectancy results
out = table; % 
lwrupr = table; % record 95% CI
Surv = struct; H1Prop = struct; H2Prop = struct;


for ageBeg = [65, 75]
    for gIndex = 1:2
        gender_i = genderList{gIndex};
        
        load([matDir, dataName, '_', Input.residence, gender_i, ...
            '_', num2str(ageBeg), '_', num2str(Input.ageEnd), '_', num2str(Input.time_ageStart)]);
        
        lifeTimeH1 = cell2mat(lifeTime(:, 1));
        lifeTimeH2 = cell2mat(lifeTime(:, 2));
        lifeTimeHT = cell2mat(lifeTime(:, 3));
        firstXAll = cell2mat(cellfun(@(x) x(:), firstX, 'un', 0));
        ratioAll = cell2mat(ratio); 
        
        
        meanLifeTimeH1 = mean(lifeTimeH1, 2);
        meanLifeTimeH2 = mean(lifeTimeH2, 2);
        meanLifeTimeHT = mean(lifeTimeHT, 2);
        meanRatioAll = mean(ratioAll, 2);
        
        varLifeTimeH1 = var(lifeTimeH1, [], 2);
        varLifeTimeH2 = var(lifeTimeH2, [], 2);
        varLifeTimeHT = var(lifeTimeHT, [], 2);
        varRatioAll = var(ratioAll, [], 2);        
        
        meanFirstXAll = cellfun(@mean, firstX);
        varFirstXAll = cellfun(@var, firstX);
        nFirstXAll = cellfun(@numel, firstX);

        % mean; standard error; standard deviation
        out.([num2str(ageBeg), gender_i]) = ...
            [mean(meanLifeTimeHT); sqrt(sum(varLifeTimeHT/nSimObs)) / nSimPsi; mean(sqrt(varLifeTimeHT)); 0
            mean(meanLifeTimeH1); sqrt(sum(varLifeTimeH1/nSimObs)) / nSimPsi; mean(sqrt(varLifeTimeH1)); 0
            mean(meanLifeTimeH2); sqrt(sum(varLifeTimeH2/nSimObs)) / nSimPsi; mean(sqrt(varLifeTimeH2)); 0
            mean(meanRatioAll); sqrt(sum(varRatioAll/nSimObs)) / nSimPsi; mean(sqrt(varRatioAll)); 0
            mean(meanFirstXAll); sqrt(sum(varFirstXAll./nFirstXAll)) / nSimPsi; mean(sqrt(varFirstXAll))];    
        
        % lwr and upr
        lwr = 0.025; upr = 0.975;
        lwrupr.([num2str(ageBeg), gender_i]) = ...
            [quantile(meanLifeTimeHT, [lwr, upr]); zeros(3, 2)
            quantile(meanLifeTimeH1, [lwr, upr]); zeros(3, 2)
            quantile(meanLifeTimeH2, [lwr, upr]); zeros(3, 2)
            quantile(meanRatioAll, [lwr, upr]); zeros(3, 2)
            quantile(meanFirstXAll, [lwr, upr]); zeros(3, 2)];
    end
end

%%

time_ageStart = [1998 2014];
nTime = numel(time_ageStart);

nSimPsi = 1000; % no. of frailty paths in each simulation



out = cell(nTime, Mdl.N_H_STATE-1);

N_H_STATE = Mdl.N_H_STATE;
for tIndex = 1:nTime
    Input.time_ageStart = time_ageStart(tIndex);
    
    path = simulatePath(Mdl, Input, 2014, nSimPsi);
    for initHState = 1:N_H_STATE-1
        
        avgYear = zeros(nSimPsi, N_H_STATE-1);
%         lifeExp = zeros(nSimPsi, 1);
        entryAge = zeros(nSimPsi, 1);
        
       
        for iSim = 1:nSimPsi
            
            if any(strcmp({'frailty', 'frailty_resid'}, model))
                Input.iPsiPath = path(:, iSim);
            else
                Input.iPsiPath = [];
            end             
            
            avgYear(iSim, :) = calAvgYear(Mdl, Input, initHState);
%             lifeExp(iSim) = calLifeExp(Mdl, Input, initHState);
            entryAge(iSim) = calFirstEntryAge(Mdl, Input, 2, initHState);            
        end
        
        lifeExp = sum(avgYear, 2);
        out{tIndex, initHState} = [avgYear, lifeExp, entryAge];
    end
end

if ~any(strcmp({'frailty', 'frailty_resid'}, Mdl.model))
    out = cell2mat(out);
elseif strcmp(Input.type, 'period')
    out = cell2mat(out);
end


%% Likelihood ratio test

% Perform the likelihood ratio test once all the models are estimated

clear

% --------------
% Select dataset
% --------------

% 1. CLHLS
dataName = 'clhls';

% % 2. HRS
% dataName = 'rndhrs';


% ----------------
% Setup model list
% ----------------

modelList = {'static', 'trend', 'rw', 'ar'};

% Model w/ the residence covariate -- only used for the CLHLS dataset
% modelList = {'static_resid', 'trend_resid', 'rw_resid', 'ar_resid'};


% ----------------
% Hypothesis test
% ----------------

% select the null model
iModel = 1; 

rModel = modelList{iModel};
rResult = load([dataName, '_', rModel, '.mat']);
rLogL = -rResult.fval;
rNPara = numel(rResult.theta);

% the alternative model
uModel = modelList{iModel+1};
uResult = load([dataName, '_', uModel, '.mat']);
uLogL = -uResult.fval;
uNPara = numel(uResult.theta);

dof = uNPara - rNPara;

[h, pValue, stat] = lratiotest(uLogL, rLogL, dof);


% --------------------
% Display the results
% --------------------
disp([num2str(stat), ' & ', MLFIMdl.calSignt(pValue)])
