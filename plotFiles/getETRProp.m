function output = getETRProp(N_H_STATE, Trs)
% Get expose to risk in each health state in terms of percentage

yETR = cell(1, N_H_STATE-1);
for hState = 1:N_H_STATE-1
    yETR{hState} = Trs.transitETR{hState};
end
yETR = cell2mat(yETR);

output = yETR ./ sum(yETR, 2);
output(:, N_H_STATE) = Trs.age;
