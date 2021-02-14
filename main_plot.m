% plots

clear;
addpath('./plotFiles')

%% Read table and set Par

clhls = readtable(['clhls', '_', 'transit.csv']);
rndhrs = readtable(['rndhrs', '_', 'transit.csv']);

clhls_f = clhls(clhls.RAFEMALE==1, :);
rndhrs_f = rndhrs(rndhrs.RAFEMALE==1, :);

clhls_u = clhls(clhls.JOINURBAN==1, :);
clhls_r = clhls(clhls.JOINURBAN==0, :);

clhls_m = clhls(clhls.RAFEMALE==0, :);
rndhrs_m = rndhrs(rndhrs.RAFEMALE==0, :);

clhls_f_u = clhls_f(clhls_f.JOINURBAN==1, :);
clhls_f_r = clhls_f(clhls_f.JOINURBAN==0, :);
clhls_m_u = clhls_m(clhls_m.JOINURBAN==1, :);
clhls_m_r = clhls_m(clhls_m.JOINURBAN==0, :);

ParClhls = setPar(clhls);
ParHrs = setPar(rndhrs);

N_H_STATE = ParClhls.N_H_STATE;
S = ParClhls.S;
year_t1 = ParClhls.year_t1;

transitPair = getTransitPair(); %[fromState, toState];

fontsize = 12;
hStateList = {'Healthy', 'Disabled', 'Dead'};


%% Calculate crude transition rates

% CLHLS / HRS
CrudeAllClhls = calTransit(clhls, N_H_STATE);
CrudeAllHrs = calTransit(rndhrs, N_H_STATE);

% female
CrudeFClhls = calTransit(clhls_f, N_H_STATE);
CrudeFHrs = calTransit(rndhrs_f, N_H_STATE);

% male
CrudeMClhls = calTransit(clhls_m, N_H_STATE);
CrudeMHrs = calTransit(rndhrs_m, N_H_STATE);

% urban / rural in CLHLS
CrudeUClhls = calTransit(clhls_u, N_H_STATE);
CrudeRClhls = calTransit(clhls_r, N_H_STATE);


%% Exploratory Data Analysis 

% ---- plot exposure to risk ----
hState = 1; % <-- change here: 1 = healthy; 2 = disabled
lcn = {'northeast', 'northwest'};
figure
plot(CrudeAllClhls.age, CrudeAllClhls.transitETR{hState}, 'k-', 'LineWidth', 2)
hold on
plot(CrudeAllHrs.age, CrudeAllHrs.transitETR{hState}, 'k--', 'LineWidth', 2)
hold off
xlim([65, 110])
xlabel('Age')
ylabel('No. of exposure years')
legend('CLHLS', 'HRS', 'Location', lcn{hState})
legend boxoff
set(gca, 'fontsize', fontsize)

% ---- plot transition counts ----
Trs = CrudeAllClhls; % <-- change here: TrsAllClhls or TrsAllHrs
transitCount = cell(1, S);
for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);
    transitCount{s} = Trs.transitCount{fromState, toState};
end
transitCount = cell2mat(transitCount);
transitCountProp = transitCount ./ sum(transitCount, 2);

% use transition counts in proportion
yTransit = transitCountProp; 
yLabel = 'Proportion of transitions'; 

figure
plot(Trs.age, yTransit(:, 1), 'k-.', 'LineWidth', 2) % healthy to disabled
hold on
plot(Trs.age, yTransit(:, 2), 'k-', 'LineWidth', 2)
plot(Trs.age, yTransit(:, 3), 'k:', 'LineWidth', 2)
plot(Trs.age, yTransit(:, 4), 'k--', 'LineWidth', 2)
hold off
xlim([65, 110])
ylim([0, 1])
xlabel('Age')
ylabel(yLabel) 
legend(...
    'Healthy to disabled', ...
    'Disabled to healthy', ...
    'Healthy to dead', ...
    'Disabled to dead', ...
    'Location', 'best')
legend boxoff
set(gca, 'fontsize', fontsize)
% ---- end of transition counts ----


% --- plot crude transition rate by gender ---
for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);
    
    figure
    semilogy(CrudeFClhls.age, CrudeFClhls.transitRate{fromState, toState}, 'ro');
    hold on
    semilogy(CrudeMClhls.age, CrudeMClhls.transitRate{fromState, toState}, 'b*');
    semilogy(CrudeFHrs.age, CrudeFHrs.transitRate{fromState, toState}, 'md');
    semilogy(CrudeMHrs.age, CrudeMHrs.transitRate{fromState, toState}, 'c+');
    hold off
    xlim([65, 110])
    ylim([1e-3, 1])
    xlabel('Age')
    ylabel('Transition rate')
    title([hStateList{fromState}, ' to ', hStateList{toState}])
    set(gca, 'fontsize', fontsize)
    
    legend('Female (CLHLS)', 'Male (CLHLS)', ...
        'Female (HRS)', 'Male (HRS)', 'Location', 'best');
    legend boxoff    
end


% --- plot crude transition rate by residence ---
for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);

    figure
    semilogy(CrudeUClhls.age, CrudeUClhls.transitRate{fromState, toState}, 'ro');
    hold on
    semilogy(CrudeRClhls.age, CrudeRClhls.transitRate{fromState, toState}, 'kx');
    hold off
    xlim([65, 110])
    ylim([1e-3, 1])
    xlabel('Age')
    ylabel('Transition rate')
    title([hStateList{fromState}, ' to ', hStateList{toState}])
    set(gca, 'fontsize', fontsize)
    
    legend('Urban', 'Rural', 'Location', 'best');
    legend boxoff    
end


% ---- crude transition rates by time ----

CrudeClhls_t = getTransit_t(clhls, ParClhls);
CrudeHrs_t = getTransit_t(rndhrs, ParHrs);

ageMin = 65; 
ageMax = 106;
nAge = ageMax - ageMin + 1;
yMin = 10^(-3);

hStateNameList = {'h', 'd', 'dead'};

% CLHLS
for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);
    
    rate = CrudeClhls_t.transitRate{fromState, toState};
    
    begInd = ageMin - CrudeClhls_t.age(1) + 1;
    endInd = begInd + nAge - 1;
    rate = rate(begInd:endInd, :);
    
    nCol = size(rate, 2);

    figure
    h = semilogy(ParClhls.t, rate);
    set(h, {'color'}, num2cell(jet(size(rate, 1)), 2));
    set(gca, 'fontsize', fontsize)
    xlim([ParClhls.t(1), ParClhls.t(end)])
    ylim([yMin, 1])
    xlabel('Year')
    ylabel('Transition rate')
    title([hStateList{fromState}, ' to ', hStateList{toState}])
    yrList = getYrList(2, ParClhls, 2014);
    xticks(ParClhls.t)
    xticklabels(yrList)
    
    if s == S % add legend
        hl = legend(arrayfun(@num2str, ageMin:ageMax, 'UniformOutput', 0), ...
            'Orientation', 'horizontal', 'Location', 'best');
        hl.NumColumns = 6;
        legend boxoff
    end
    
%     % uncomment to save
%     tmp = ['trsRate_byWv_', hStateNameList{fromState}, '2', ...
%         hStateNameList{toState}, '_clhls'];
%     saveas(gca, tmp, 'epsc');
end


% HRS
for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);
    rate = CrudeHrs_t.transitRate{fromState, toState};
    
    begInd = ageMin - CrudeHrs_t.age(1) + 1;
    endInd = begInd + nAge - 1;
    rate = rate(begInd:endInd, :);
    
    nCol = size(rate, 2);

    figure
    h = semilogy(ParHrs.t, rate);
    set(h, {'color'}, num2cell(jet(size(rate, 1)), 2));
    set(gca, 'fontsize', fontsize)
    xlim([ParHrs.t(1), ParHrs.t(end)])
    ylim([yMin, 1])
    xlabel('Year')
    ylabel('Transition rate')
    title([hStateList{fromState}, ' to ', hStateList{toState}])
    yrList = getYrList(2, ParHrs, 2014);
    xticks(ParHrs.t)
    xticklabels(yrList)
    
    if s == S % add legend
        hl = legend(arrayfun(@num2str, ageMin:ageMax, 'UniformOutput', 0), ...
            'Orientation', 'horizontal', 'Location', 'best');
        hl.NumColumns = 6;
        legend boxoff
    end
    
%     % uncomment to save
%     tmp = ['trsRate_byWv_', hStateNameList{fromState}, '2', ...
%         hStateNameList{toState}, '_hrs'];
%     saveas(gca, tmp, 'epsc');
    
end


%% Crude transition rates overlaid by static results

model = 'static';
clhls_static = load(['clhls', '_', 'rate', '_', model]);
rndhrs_trend = load(['rndhrs', '_', 'rate', '_', model]);

clhls_f_static = getTrsMatrix(clhls_static.trsRateF, transitPair);
clhls_m_static = getTrsMatrix(clhls_static.trsRateM, transitPair);
hrs_f_static = getTrsMatrix(rndhrs_trend.trsRateF, transitPair);
hrs_m_static = getTrsMatrix(rndhrs_trend.trsRateM, transitPair);

% comment/uncomment to switch between CLHLS and HRS results
% HRS
TrsFemale = CrudeFHrs;
TrsMale = CrudeMHrs;
staticFemale = hrs_f_static;
staticMale = hrs_m_static;

% % CLHLS
% TrsFemale = CrudeFClhls;
% TrsMale = CrudeMClhls;
% staticFemale = clhls_f_static;
% staticMale = clhls_m_static;

for s = 1:S
    fromState = transitPair(s, 1);
    toState = transitPair(s, 2);
    
    figure
    semilogy(TrsFemale.age, TrsFemale.transitRate{fromState, toState}, 'ro');
    hold on
    semilogy(staticFemale(:, end), staticFemale(:, s), 'r-', 'LineWidth', 2);
    semilogy(TrsMale.age, TrsMale.transitRate{fromState, toState}, 'b*');
    semilogy(staticMale(:, end), staticMale(:, s), 'b--', 'LineWidth', 2);
    hold off
    xlim([65, 110])
    ylim([1e-3, 1])
    xlabel('Age')
    ylabel('Transition rate')
    title([hStateList{fromState}, ' to ', hStateList{toState}])
    set(gca, 'fontsize', fontsize)    
    legend('Female (crude)', 'Female (static)', ...
        'Male (crude)', 'Male (static)', 'Location', 'best');
    legend boxoff    
end


%% Trend model vs. Hanewald et al. (2019)

clear;

hStateList = {'Healthy', 'Disabled', 'Dead'};
transitPair = getTransitPair(); %[fromState, toState];
time_ageStart = 1998; 

% load Hanewald et al. (2019) results
h2019Mat = load(['h2019', '_', 'rate', '_', num2str(time_ageStart), '.mat']);
h2019Mat.tableVarName;
h2019_rate = h2019Mat.trsRate;

% load trend model results
model = 'trend_resid';
clhls_u_trend = load(['clhls', '_', 'rate', '_', model, '_', 'u_', num2str(time_ageStart), '.mat']);
clhls_r_trend = load(['clhls', '_', 'rate', '_', model, '_', 'r_', num2str(time_ageStart), '.mat']);

clhls_uf_trend = getTrsMatrix(clhls_u_trend.trsRateF, transitPair);
clhls_rf_trend = getTrsMatrix(clhls_r_trend.trsRateF, transitPair);
clhls_um_trend = getTrsMatrix(clhls_u_trend.trsRateM, transitPair);
clhls_rm_trend = getTrsMatrix(clhls_r_trend.trsRateM, transitPair);

N = size(h2019_rate, 2);

trend_rate = cell(1, N);
titleList = cell(1, N);
clhlsTrendCell = {clhls_um_trend; clhls_rm_trend; clhls_uf_trend; clhls_rf_trend};
typeList = {'Urban Male', 'Rural Male', 'Urban Female', 'Rural Female'};
sCell = {1, 3, 4}; % Hanewald et al. (2019) do not consider recovery, so s=2 is excluded
for i = 1:numel(sCell)
    s = sCell{i};
    for j = 1:4
        index = j + (i-1)*4;
        trend_rate{index} = clhlsTrendCell{j}(:, s);
        
        fromState = transitPair(s, 1); toState = transitPair(s, 2);
        tmp = [typeList{j}, ': ', hStateList{fromState}, ' to ', hStateList{toState}];
        titleList{index} = tmp;
    end
end
trend_rate = cell2mat(trend_rate);

fileNameList = h2019Mat.tableVarName;
for i = 1:N
    figure
    semilogy(h2019Mat.x_age+65, h2019_rate(:, i), 'k--', 'LineWidth', 2)
    hold on
    semilogy(clhls_uf_trend(:, end), trend_rate(:, i), 'k-', 'LineWidth', 2)
    hold off
    xticks(65:10:105)
    xlim([65, 105])
    ylim([1e-3, 1])
    title(titleList{i})
    xlabel('Age')
    ylabel('Transition rate')
    legend('Hanewald et al. (2019)', 'Trend model', 'Location', 'best')
    legend boxoff
    
%     % uncomment to save
%     set(gca, 'fontsize', 18)
%     tmp = ['trsRate_cf_', fileNameList{i}];
%     saveas(gca, tmp, 'epsc');
end


%% Estimated rates vs. Li et al. (2017) and Sherris and Wei (2020)
% **static model** and **trend model**

clear

hStateList = {'Healthy', 'Disabled', 'Dead'};
fontsize = 18;

transitPair = getTransitPair(); %[fromState, toState];
S = size(transitPair, 1);

% load results of Li et al. (2017) and Sherris and Wei (2020)
li_sw_mat = load('li2017_sw2020_rate.mat');
li_sw_trsRate = li_sw_mat.trsRate;

% load results of the static model or trend model
model = 'static'; % static / trend
time_ageStart = 2010;
if strcmp(model, 'static')
    tmp = ['rndhrs', '_', 'rate', '_', model, '.mat'];
else
    tmp = ['rndhrs', '_', 'rate', '_', model, '_', num2str(time_ageStart), '.mat'];
end
rndhrs_trend = load(tmp);

hrs_f = getTrsMatrix(rndhrs_trend.trsRateF, transitPair);
hrs_m = getTrsMatrix(rndhrs_trend.trsRateM, transitPair);

% female
for s = 1:S
    y2 = li_sw_trsRate.(['sw', '_' model, '_', 'f']);
    y3 = li_sw_trsRate.(['li', '_' model, '_', 'f']);
    
    figure
    semilogy(hrs_f(:, end), hrs_f(:, s), 'k-', 'LineWidth', 2)
    hold on
    semilogy(li_sw_mat.x_age, y2(:, s), 'b--', 'LineWidth', 2)
    semilogy(li_sw_mat.x_age, y3(:, s), 'r:', 'LineWidth', 2)
    hold off
    
    xticks(65:10:105)
    xlim([65, 105])
    ylim([1e-3, 1])
    
    xlabel('Age')
    ylabel('Transition rate')
    
    fromState = transitPair(s, 1); toState = transitPair(s, 2);
    title(['Female: ', hStateList{fromState}, ' to ', hStateList{toState}])
    legend([upper(model(1)), model(2:end), ' model'], 'Sherris and Wei (2020)', 'Li et al. (2017)', ...
        'Location', 'best')
    legend boxoff

    set(gca, 'fontsize', fontsize)
    
%     % uncomment to save
%     fileName = ['hRateCf', '_', 's', num2str(s), '_', 'f', '_', model, '_', num2str(time_ageStart)];
%     saveas(gca, fileName, 'epsc')
end

% male
for s = 1:S
    y2 = li_sw_trsRate.(['sw', '_' model, '_', 'm']);
    y3 = li_sw_trsRate.(['li', '_' model, '_', 'm']);
    
    figure
    semilogy(hrs_m(:, end), hrs_m(:, s), 'k-', 'LineWidth', 2)
    hold on
    semilogy(li_sw_mat.x_age, y2(:, s), 'b--', 'LineWidth', 2)
    semilogy(li_sw_mat.x_age, y3(:, s), 'r:', 'LineWidth', 2)
    hold off
    
    xticks(65:10:105)
    xlim([65, 105])
    ylim([1e-3, 1])
    
    xlabel('Age')
    ylabel('Transition rate')
    
    fromState = transitPair(s, 1); toState = transitPair(s, 2);
    title(['Male: ', hStateList{fromState}, ' to ', hStateList{toState}])
    legend([upper(model(1)), model(2:end), ' model'], 'Sherris and Wei (2020)', 'Li et al. (2017)', ...
        'Location', 'best')
    legend boxoff
    
    set(gca, 'fontsize', fontsize)
    
%     % uncomment to save
%     fileName = ['hRateCf', '_', 's', num2str(s), '_', 'm', '_', model, '_', num2str(time_ageStart)];
%     saveas(gca, fileName, 'epsc')
end




%% plot rural-urban transition

clear

Y = [0	13079	3725	4203	13462; ...
    0	1497	495	626	2038; ...
    0	409	109	196	479; ...
    1796	6788	1568	2071	5487; ...
    931	2525	570	1024	3034];

% modify Y to get 100% stacked
Y = bsxfun(@rdivide, Y, sum(Y,2)); Y = Y*100;

fontsize = 12;

hb = bar(Y, 'stacked', 'FaceColor', 'flat');
legend(fliplr(hb), ...
    fliplr({'Missing', ...
    'No change: Rural', 'Rural to Urban', ...
    'Urban to Rural', 'No change: Urban'}), 'Location', 'eastoutside')
legend boxoff
% set(gca, 'fontsize', fontsize)
xticklabels({'0','1','2','3','4'})
xlabel('Health transition type')
ylabel('Proportion of urban-rural transitions')
ytickformat('percentage')
% change color
for k = 1:size(Y, 2)
    hb(k).CData = [0,0,0] + k/size(Y, 2); 
end


%% plot survival curves and health distribution

clear;

time_ageStart = 1998;
ageEnd = 110;
dataName = {'clhls', 'rndhrs'};
gender = {'f', 'm'};

% used in figure title
dName = {'CLHLS', 'HRS'}; gName = {'Female', 'Male'};

fileDir = './simFiles/';

% compute survival curve and health distribution
Surv = struct; H1Prop = struct; H2Prop = struct;
for iData = 1:numel(dName)
    dataName_i = dataName{iData};
    
    for ageBeg = [65, 75]
        for gIndex = 1:numel(gender)
            gender_i = gender{gIndex};
            
            load([fileDir, dataName_i, '_', ...
                gender_i, '_', num2str(ageBeg), '_110_', num2str(time_ageStart)]);                
            
            nSimPsi = size(hStateCount, 1);
            h1Prop = cell(1, nSimPsi); h2Prop = h1Prop; survival = h1Prop;
            for iPsi = 1:nSimPsi
                hStateCount_i = hStateCount{iPsi};
                h1Count = hStateCount_i(:, 1);
                h2Count = hStateCount_i(:, 2);
                hAlive = h1Count + h2Count;
                hTotal = sum(hStateCount_i, 2);
                
                survival{iPsi} = (h1Count+h2Count) ./ hTotal;
                h1Prop{iPsi} = h1Count ./ hTotal;
                h2Prop{iPsi} = h2Count./ hTotal;
            end
            
            survMat = cell2mat(survival);
            h1PropMat = cell2mat(h1Prop);
            h2PropMat = cell2mat(h2Prop);
            
            Surv.([dataName{iData}, num2str(ageBeg), gender{gIndex}]) = survMat;
            H1Prop.([dataName{iData}, num2str(ageBeg), gender{gIndex}]) = h1PropMat;
            H2Prop.([dataName{iData}, num2str(ageBeg), gender{gIndex}]) = h2PropMat;
             
        end
    end
end


% load results of static + trend models
clhls_static = load('clhls_static_1998');
clhls_trend = load(['clhls_trend_', num2str(time_ageStart)]);
hrs_static = load('rndhrs_static_1998');
hrs_trend = load(['rndhrs_trend_', num2str(time_ageStart)]);

% plot survival curve CLHLS and HRS in the same figure
% frailty + static
for ageBeg = [65, 75]
    xAge = ageBeg:110;
    for gIndex = 1:2
        gender_i = gender{gIndex};
        
        % frailty
        y1 = Surv.(['rndhrs', num2str(ageBeg), gender{gIndex}]);
        y2 = Surv.(['clhls', num2str(ageBeg), gender{gIndex}]);
        
        % static
        y1Static = hrs_static.survProp.([gender_i, num2str(ageBeg)]); 
        y2Static = clhls_static.survProp.([gender_i, num2str(ageBeg)]);
        
        % trend
        y1Trend = hrs_trend.survProp.([gender_i, num2str(ageBeg)]); 
        y2Trend = clhls_trend.survProp.([gender_i, num2str(ageBeg)]);
        
        plotMeanCI_clhls_hrs(xAge, y1, y2, y1Static, y2Static)        
        title(gName{gIndex})
        
        legend('Frailty 95% CI (HRS)', 'Frailty Mean (HRS)', 'Static (HRS)', ...
            'Frailty 95% CI (CLHLS)', 'Frailty Mean (CLHLS)', 'Static (CLHLS)', 'Location', 'best')
        legend boxoff
        set(gca, 'fontsize', 12)

%         % uncomment to save
%         fileName = ['surv', '_', gender_i, '_', num2str(ageBeg), '_', ...
%                 'clhls_hrs', '_', num2str(time_ageStart), '_allModel'];
%         saveas(gca, fileName, 'epsc')
    end
end


%  plot disability distribution: static + trend + frailty
for iData = 1:2
    dataName_i = dataName{iData};    
        
    for ageBeg = [65, 75]
        xAge = ageBeg:ageEnd;
        
        for gIndex = 1:2
            gender_i = gender{gIndex};

            if iData == 1
                hProp_static = clhls_static.allProp.([gender_i,num2str(ageBeg)]);
                hProp_trend = clhls_trend.allProp.([gender_i,num2str(ageBeg)]);
            elseif iData == 2
                hProp_static = hrs_static.allProp.([gender_i,num2str(ageBeg)]);
                hProp_trend = hrs_trend.allProp.([gender_i,num2str(ageBeg)]);
            end
            
            h1prop_static = hProp_static(:, 1);
            h2prop_static = hProp_static(:, 2);
            
            h1prop_trend = hProp_trend(:, 1);
            h2prop_trend = hProp_trend(:, 2);

            h1prop_frailty = H1Prop.([dataName_i, num2str(ageBeg),gender_i]);
            h2prop_frailty = H2Prop.([dataName_i, num2str(ageBeg),gender_i]);

            figure
            hold on
            lower = quantile(h2prop_frailty, 0.025, 2);
            upper = quantile(h2prop_frailty, 0.975, 2);
            ciplot(lower, upper, xAge, rgb('light pink'))
            plot(xAge, mean(h2prop_frailty, 2), 'r-', 'LineWidth', 2)
            plot(xAge, h2prop_trend, 'b-.', 'LineWidth', 2)
            plot(xAge, h2prop_static, 'k:', 'LineWidth', 2)
                    
            hold off
            legend('Frailty 95% CI', 'Frailty mean', 'Trend', 'Static',...
                'Location', 'best')
            legend boxoff 

            xlim([xAge(1), xAge(end)])
            ylim([0, 0.2])

            xlabel('Age')
            ylabel('Probability')        
            titleName = [dName{iData}, ' ', gName{gIndex}];
            title(titleName)        
            set(gca, 'fontsize', 12)
            
%             % uncomment to save
%             fileName = ['disDist', '_', gender_i, '_', num2str(ageBeg), '_', ...
%                 dataName_i, '_', num2str(time_ageStart), '_allModel'];
%             saveas(gca, fileName, 'epsc')
        end
    end
end



%% urban-rural results

clear;

time_ageStart = 1998; % 1998 or 2014
residNameList = {'Urban', 'Rural'};

% load results of static + trend models
clhls_static = struct;
clhls_static.u = load('clhls_static_resid_u_1998');
clhls_static.r = load('clhls_static_resid_r_1998');

clhls_trend = struct;
clhls_trend.u = load(['clhls_trend_resid_u_', num2str(time_ageStart)]);
clhls_trend.r = load(['clhls_trend_resid_r_', num2str(time_ageStart)]);

Surv = struct; H1Prop = struct; H2Prop = struct;
resid = {'u', 'r'};
gender = {'f', 'm'};
rName = {'Urban', 'Rural'};  % used in figure title, match resid
gName = {'Female', 'Male'}; % used in figure title, match gender
fileDir = './simFiles/';
for iResid = 1:2
    resid_i = resid{iResid};
    
    for ageBeg = [65, 75]
        for gIndex = 1:2
            gender_i = gender{gIndex};
            
            load([fileDir, 'clhls', '_', ...
                resid_i, gender_i, '_', num2str(ageBeg), '_110_', num2str(time_ageStart)]);
         
            % survival curve and health distribution
            nSimPsi = size(hStateCount, 1);
            h1Prop = cell(1, nSimPsi); h2Prop = h1Prop; survival = h1Prop;
            for iPsi = 1:nSimPsi
                hStateCount_i = hStateCount{iPsi};
                h1Count = hStateCount_i(:, 1);
                h2Count = hStateCount_i(:, 2);
                hAlive = h1Count + h2Count;
                hTotal = sum(hStateCount_i, 2);
                
                survival{iPsi} = (h1Count+h2Count) ./ hTotal;
                h1Prop{iPsi} = h1Count ./ hTotal;
                h2Prop{iPsi} = h2Count./ hTotal;
            end
            survMat = cell2mat(survival);
            h1PropMat = cell2mat(h1Prop);
            h2PropMat = cell2mat(h2Prop);
            
            Surv.(['clhls', num2str(ageBeg), resid_i, gender_i]) = survMat;
            H1Prop.(['clhls', num2str(ageBeg), resid_i, gender_i]) = h1PropMat;
            H2Prop.(['clhls', num2str(ageBeg), resid_i, gender_i]) = h2PropMat;

        end
    end
end



% --- plot survival curve ---
for ageBeg = [65, 75]
    xAge = ageBeg:110;
    for rIndex = 1:2
        resid_i = resid{rIndex};
        
        % frailty
        y1 = Surv.(['clhls', num2str(ageBeg), resid_i, 'f']);
        y2 = Surv.(['clhls', num2str(ageBeg), resid_i, 'm']);
        
        % static
        y1Static = clhls_static.(resid_i).survProp.(['f', num2str(ageBeg)]); 
        y2Static = clhls_static.(resid_i).survProp.(['m', num2str(ageBeg)]); 
        
        % trend
        y1Trend = clhls_trend.(resid_i).survProp.(['f', num2str(ageBeg)]);
        y2Trend = clhls_trend.(resid_i).survProp.(['m', num2str(ageBeg)]);
        
        plotMeanCI_clhls_resid(xAge, y1, y2, y1Static, y2Static, y1Trend, y2Trend)
        
        title(residNameList{rIndex})
        legend('Frailty 95% CI (Female)', 'Frailty Mean (Female)', 'Static (Female)', ...
            'Frailty 95% CI (Male)', 'Frailty Mean (Male)', 'Static (Male)', 'Location', 'best')
        legend boxoff
        set(gca, 'fontsize', 12)

%         % uncomment to save
%         fileName = ['surv', '_', num2str(ageBeg), '_', ...
%                 'clhls_', resid_i, '_', num2str(time_ageStart), '_allModel'];
%         saveas(gca, fileName, 'epsc')
    end
end


% ---- plot disability distribution: static + trend + frailty ----
for rIndex = 1:2
    resid_i = resid{rIndex};
        
    for ageBeg = [65, 75]
        xAge = ageBeg:110;
        
        for gIndex = 1:2
            gender_i = gender{gIndex};

            hProp_static = clhls_static.(resid_i).allProp.([gender_i,num2str(ageBeg)]);
            hProp_trend = clhls_trend.(resid_i).allProp.([gender_i,num2str(ageBeg)]);

            
            h1prop_static = hProp_static(:, 1);
            h2prop_static = hProp_static(:, 2);
            
            h1prop_trend = hProp_trend(:, 1);
            h2prop_trend = hProp_trend(:, 2);

            h1prop_frailty = H1Prop.(['clhls', num2str(ageBeg), resid_i, gender_i]);
            h2prop_frailty = H2Prop.(['clhls', num2str(ageBeg), resid_i, gender_i]);

            figure
            hold on
            lower = quantile(h2prop_frailty, 0.025, 2);
            upper = quantile(h2prop_frailty, 0.975, 2);
            ciplot(lower, upper, xAge, rgb('light pink'))
            plot(xAge, mean(h2prop_frailty, 2), 'r-', 'LineWidth', 2)
            plot(xAge, h2prop_trend, 'b-.', 'LineWidth', 2)
            plot(xAge, h2prop_static, 'k:', 'LineWidth', 2)
                    
            hold off
            legend('Frailty 95% CI', 'Frailty mean', 'Trend', 'Static',...
                'Location', 'best')
            legend boxoff 

            xlim([xAge(1), xAge(end)])
            ylim([0, 0.2])

            xlabel('Age')
            ylabel('Probability')        
            titleName = [residNameList{rIndex}, ' ', gName{gIndex}];
            title(titleName)        
            set(gca, 'fontsize', 12)
            
%             % uncomment to save
%             fileName = ['disDist', '_', resid_i, gender_i, '_', num2str(ageBeg), '_', ...
%                 'clhls', '_', num2str(time_ageStart), '_allModel'];
%             saveas(gca, fileName, 'epsc')
        end
    end
end

