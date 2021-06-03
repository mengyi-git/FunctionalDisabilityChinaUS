
clear

addpath('./cfFiles')
addpath('./plotFiles')

%% Calculate transition rates in Hanewald et al. (2019)

% Model used in Hanewald et al. (2019)
% beta_0 + beta_1 * x + beta_2 * x^2 + beta_3 * t + beta_4 * t * x + beta_5 * t * x^2

h2019Table = readtable('h2019Table3.csv');

tableVarName = h2019Table.Properties.VariableNames;

coef_prep = h2019Table{:,:};
coef = transformBeta(coef_prep);


% -----------------------
% compute transition rate
% -----------------------

ageInterp = 65; % x = 0 at age 65 (see p. 152)
year_t0 = 1998; % t = 0 in year 1998 (see p. 152)

ageBeg = 65 - ageInterp;
ageEnd = 105 - ageInterp;
n = ageEnd - ageBeg + 1;


time_ageStart = 2002; % 1998 

x_age = ageBeg:ageEnd;
x_time = time_ageStart - year_t0 + (0:n-1);

X = [ones(1, n); x_age; x_age.^2; x_time; x_age .* x_time; x_time .* x_age.^2];

trsRate = exp(coef * X);
trsRate = trsRate';

notes = {'t = 0 in year 1998'; 'age variable is defined as x = age - 65'};

saveDir = './plotFiles/';
save([saveDir, 'h2019_rate', '_', num2str(time_ageStart)], ...
    'trsRate', 'x_age', 'x_time', 'tableVarName', 'time_ageStart', 'notes')


%% Calculate transition rates in Li et al. (2017) and Sherris and Wei (2020)

clear

transitPair = getTransitPair();
S = size(transitPair, 1); % type of transitions

ageBeg = 65; 
ageEnd = 109; % survival probability at age 110 is 0
n = ageEnd - ageBeg + 1;

% t = 1 in 1998, t = 2 in 2000, t = 3 in 2002, t = 4 in 2004 etc
year_t1 = 1998;

time_ageStart = 2010; 

genderList = {'f', 'm'};

% -----------------------
% static and trend model
% -----------------------

hatTable = readtable('li2017_sw2020_thetaHat.csv');
hatName = hatTable.Properties.VariableNames;
nName = numel(hatName);

notes = {'age last birthday'; 't = 1 in 1998, t = 2 in 2000, t = 3 in 2002, t = 4 in 2004 etc'};
trsRate = struct;
for iGender = 1:2
    gender_i = genderList{iGender};
    if strcmp(gender_i, 'f')
        x_gender = ones(1, n); % female
    elseif strcmp(gender_i, 'm')
        x_gender = zeros(1, n); % male
    end
   
    
    for i = 1:nName
        hatName_i = hatName{i};
        theta = hatTable.(hatName_i);
        params = reshape(theta, [S, numel(theta)/S]);
        
        if strcmp(hatName_i(1:2), 'fu')
            x_age = ((ageBeg:ageEnd) - 65)/10;
            x_time = (time_ageStart - (year_t1-1) + (0:n-1))/10;
        else
            x_age = ageBeg:ageEnd;
            x_time = 1 + (time_ageStart - year_t1 + (0:n-1))/2; 
        end
        
        X = [ones(1, n); x_age; x_gender; x_time];
        
        trsRate_i = exp(params * X); 

        % biannual rate in Li et al. (2017)
        if strcmp(hatName_i(1:2), 'li')
            trsRate_i = trsRate_i/2;
        end

        trsRate.([hatName_i, '_', gender_i]) = trsRate_i';
    end

end

saveDir = './plotFiles/';
save([saveDir, 'li2017_sw2020_rate.mat'], 'trsRate', 'time_ageStart', 'x_age', 'x_time', 'notes')



% -----------------------
% Frailty model
% -----------------------

% load estimated parameters
hatFrailtyTable = readtable('li2017_sw2020_thetaHat_frailty.csv');
hatFrailtyName = hatFrailtyTable.Properties.VariableNames;
nName = numel(hatFrailtyName);

% simulate path
N_SIM = 1000;
path = zeros(n, N_SIM);
path(1, :) = 0;  % initialise psi path

rng(2*n) % set seed
for tIndex = 2:n
    path(tIndex, :) = path(tIndex-1, :) + randn(1, N_SIM);
end

% compute transition rate
trsRateFrailty = struct;
for iGender = 1:2
    gender_i = genderList{iGender};
    if strcmp(gender_i, 'f')
        x_gender = ones(1, n); % female
    elseif strcmp(gender_i, 'm')
        x_gender = zeros(1, n); % male
    end
    
    trsRateCell = cell(N_SIM, 1);

    for i = 1:nName
        hatName_i = hatFrailtyName{i};
        theta = hatFrailtyTable.(hatName_i);
        params = reshape(theta, [S, numel(theta)/S]);

        if strcmp(hatName_i(1:2), 'fu')
            x_age = ((ageBeg:ageEnd) - 65)/10;
            x_time = (time_ageStart - (year_t1-1) + (0:n-1))/10;
        else
            x_age = ageBeg:ageEnd;
            x_time = 1 + (time_ageStart - year_t1 + (0:n-1))/2; 
        end
        
        for iSim = 1:N_SIM
            x_psi = path(:, iSim); x_psi = x_psi';
            X = [ones(1, n); x_age; x_gender; x_time; x_psi];

            trsRate_i = exp(params * X); 

            % biannual rate in Li et al. (2017)
            if strcmp(hatName_i(1:2), 'li')
                trsRate_i = trsRate_i/2;
            end

            trsRate_i = trsRate_i';

            trsRateCell{iSim, 1} = trsRate_i;

        end
        
        trsRateFrailty.([hatName_i, '_', gender_i]) = trsRateCell;
    end

end

nameT = fieldnames(trsRateFrailty);

% Rearrange the transition rates such that
% each element in the struct is a 1-by-S cell
trsRateFrailtyCopy = struct;
for i = 1:numel(nameT)
    name_i = nameT{i};
    trsRate_name_i = trsRateFrailty.(name_i);
    
    trsRateCopy = cell(1, S);
    for s = 1:S
        trsRateCopy{s} = zeros(n, N_SIM);
    end
       
    for iSim = 1:N_SIM
        trsRate_i = trsRate_name_i{iSim};
        for s = 1:S
            trsRateCopy{s}(:, iSim) = trsRate_i(:, s);
        end
    end
    
    trsRateFrailtyCopy.(name_i) = trsRateCopy;
end


% --------------------------------------
% plot transition rate - frailty model
% --------------------------------------

% ---- load crude rate ----

% Read table and set Par
rndhrs = readtable(['rndhrs', '_', 'transit.csv']);
rndhrs_f = rndhrs(rndhrs.RAFEMALE==1, :);
rndhrs_m = rndhrs(rndhrs.RAFEMALE==0, :);
ParHrs = setPar(rndhrs);

tIndex = find(ParHrs.t == time_ageStart - ParHrs.year_t1 + 1); 
crudeRateName = ['Crude rate: ', num2str(time_ageStart)];

% Calculate crude transition rates in each wave
CrudeFHrs_t = getTransit_t(rndhrs_f, ParHrs);
CrudeMHrs_t = getTransit_t(rndhrs_m, ParHrs);
Crude_t = {CrudeFHrs_t; CrudeMHrs_t};

% ---- load fitted rate ----
x_age = ageBeg:ageEnd;
trsRateLi = trsRateFrailtyCopy.li_frailty_f;    
trsRateSW = trsRateFrailtyCopy.sw_frailty_f;    
trsRateFu = trsRateFrailtyCopy.fu_frailty_f;

gName = {'Female', 'Male'}; % used in figure title
hStateList = {'Healthy', 'Disabled', 'Dead'};

for gIndex = 1:2
    gender_i = lower(gName{gIndex}(1));
    crude_i = Crude_t{gIndex};
    
    for s = 1:S    
        fromState = transitPair(s, 1); toState = transitPair(s, 2);    

        trsRateLi_s = trsRateLi{s};
        trsRateSW_s = trsRateSW{s};
        trsRateFu_s = trsRateFu{s};

        yMeanLi = mean(trsRateLi_s, 2);
        yMeanSW = mean(trsRateSW_s, 2);
        yMeanFu = mean(trsRateFu_s, 2);

        figure
        plot(crude_i.age, log10(crude_i.transitRate{fromState, toState}(:, tIndex)), 'kx');
        hold on
        plot(x_age, log10(yMeanFu), 'k-', 'LineWidth', 2)
        plot(x_age, log10(yMeanSW), 'b--', 'LineWidth', 2)
        plot(x_age, log10(yMeanLi), 'r:', 'LineWidth', 2)
        hold off
        xlim([65, 105])
        ylim([-3, 0])

        xticks(65:10:105)

        xlabel('Age')
        ylabel('log_{10} (transition rate)')
        title([gName{gIndex}, ': ', hStateList{fromState}, ' to ', hStateList{toState}])
        
        legend(crudeRateName, 'Frailty model', ...
            'Sherris and Wei (2020)', 'Li et al. (2017)', ...
            'Location', 'best')
        legend boxoff
        
        set(gca, 'fontsize', 18)

        % uncomment to save
        fileName = ['log_hRateCf', '_', 's', num2str(s), '_', gender_i, '_frailty', '_', num2str(time_ageStart)];
        saveas(gca, fileName, 'epsc')
    end
end

