function result = calTransit(dataset, N_H_STATE)
% Calculate crude transition rates

result = struct;

stateArray = dataset.RxHSTATE;
stateNextArray = dataset.RxHSTATE2;

ageInput = dataset.RxAGE;
ageNextInput = dataset.RxAGE2;

tauInput = dataset.TAU;

N_OBS = size(dataset, 1);

ageMin = floor(min(ageInput));
ageMax = ceil(max(ageNextInput));
N_AGE = ageMax - ageMin + 1; % total number of integer ages

%% Defining a matrix to deposit the central exposed to risk %
rawCETR = cell(N_H_STATE-1, 1);

for iHState = 1:N_H_STATE-1
    rawCETR{iHState} = zeros(N_OBS, N_AGE);
end

% Defining a matrix to deposit the number of transitions %
transitRawCell = cell(N_H_STATE - 1, N_H_STATE);
for iHState = 1:N_H_STATE-1
    for jHState = 1:N_H_STATE
        if iHState ~= jHState
            transitRawCell{iHState, jHState} = zeros(N_OBS, N_AGE);
        end
    end
end


for iObs = 1:N_OBS
    
    state = stateArray(iObs);
    stateNext = stateNextArray(iObs);
    
    age = ageInput(iObs);
    ageNext = ageNextInput(iObs);

    rawCETR{state}(iObs, :) = 0;
    rawCETR{state}(iObs, age - (ageMin-1)) = tauInput(iObs);

        
    % indicator of transitions
    if state ~= stateNext
        if ageNext < ageMax
            transitRawCell{state, stateNext}(iObs, floor(ageNext) - (ageMin-1)) = 1;
        else
            transitRawCell{state, stateNext}(iObs, N_AGE) = 1;
        end 
    end
    
end

% summarise central exposed to risk
centralETR = cellfun(@(x) sum(x, 1)', rawCETR, 'UniformOutput', false);

% summarise number of transitions
transitCell = cellfun(@(x) sum(x, 1)', transitRawCell, 'UniformOutput', false);

% compute raw transition rates
transitRateRawCell = cell(N_H_STATE - 1, N_H_STATE);
for iHState = 1:N_H_STATE-1
    for jHState = 1:N_H_STATE
        
        if iHState ~= jHState
            transitRateRawCell{iHState, jHState} = transitCell{iHState, jHState} ./ centralETR{iHState};
        end
        
    end
end


%% save results

result.age = (ageMin:ageMax)';
result.transitRate = transitRateRawCell;
result.transitCount = transitCell;
result.transitETR = centralETR;



