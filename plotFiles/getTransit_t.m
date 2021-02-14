function Result = getTransit_t(datatable, Par)
% Calculate crude transition rates by wave


ageMin = Par.ageMin;
ageMax = Par.ageMax;

N_H_STATE = Par.N_H_STATE;

t = Par.t; % unique waves

transitRate = cell(N_H_STATE-1, N_H_STATE);
transitCount = cell(N_H_STATE-1, N_H_STATE);
transitETR = cell(N_H_STATE-1, 1);

% exposure
for iHState = 1:N_H_STATE-1
    trsETR = nan(ageMax - ageMin + 1, Par.N);
    
    for tIndex = 1:Par.N

        iTime = t(tIndex);
        datatable_t = datatable(datatable.TIME == iTime, :);                

        Trs_t = calTransit(datatable_t, N_H_STATE);

        begInd = Trs_t.age(1) - Par.ageMin + 1;
        endInd = begInd + numel(Trs_t.age) - 1;

        trsETR(begInd:endInd, tIndex) = Trs_t.transitETR{iHState};

    end % end tIndex  
    
    transitETR{iHState} = trsETR;
end

for iHState = 1:N_H_STATE-1
    
    for jHState = 1:N_H_STATE
        
        if iHState ~= jHState % transition occurs
            rate = nan(ageMax - ageMin + 1, Par.N);
            trsCount = rate;
            
            for tIndex = 1:Par.N

                iTime = t(tIndex);
                datatable_t = datatable(datatable.TIME == iTime, :);                

                Trs_t = calTransit(datatable_t, N_H_STATE); 
                
                begInd = Trs_t.age(1) - Par.ageMin + 1;
                endInd = begInd + numel(Trs_t.age) - 1;

                rate_t = Trs_t.transitRate{iHState, jHState};
                trsCount_t = Trs_t.transitCount{iHState, jHState};
             
                rate(begInd:endInd, tIndex) = rate_t;
                trsCount(begInd:endInd, tIndex) = trsCount_t;
                
            end % end tIndex            
            
            transitRate{iHState, jHState} = rate;
            transitCount{iHState, jHState} = trsCount;
        end % end if
        
    end
end  


%%

Result = struct;
Result.age = (ageMin:ageMax)';
Result.transitRate = transitRate;
Result.transitCount = transitCount;
Result.transitETR = transitETR;
