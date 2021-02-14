function trs_matrix = getTrsMatrix(trs_cell, transitPair)
% Rearrange transition matrix: 
% from 
% cell (N_AGE-by-2): first column is age, second column is transition matrix
% to 
% matrix (N_AGE-by-(S+1)): first S columns correspond to each transition type, last column is age

nAge = size(trs_cell, 1);

S = size(transitPair, 1);

trs_matrix = zeros(nAge, S);

for i = 1:nAge
    for s = 1:S
        fromState = transitPair(s, 1);
        toState = transitPair(s, 2);
        
        trs_matrix(i, s) = trs_cell{i, 2}(fromState, toState);
    end
end

trs_matrix(:, S+1) = cell2mat(trs_cell(:, 1));
