function Par = setPar(datatable)

Par = MLFIPar;

T = readtable('transition.csv');
transition = table2array(T);
S = max(transition(:));

% t: unique(datatable.WAVE);
% N: numel(t);
[N, t] = setN(datatable.TIME);

ageMin = floor(min(datatable.RxAGE));
ageMax = ceil(max(datatable.RxAGE2));

nHState = numel(unique(datatable.RxHSTATE2));

%% set the properties in Par

Par.N = N;
Par.N_H_STATE = nHState;
Par.S = S;
Par.t = t;

Par.ageMin = ageMin;
Par.ageMax = ageMax;
