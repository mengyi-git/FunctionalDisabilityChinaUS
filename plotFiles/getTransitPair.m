function transitPair = getTransitPair()

T = readtable('transition.csv');
transition = table2array(T);
S = max(transition(:));

fromState = zeros(S, 1);
toState = zeros(S, 1);

for i = 1:S
   [row, col] = find(transition == i);
   fromState(i) = row;
   toState(i) = col;
end
transitPair = [fromState, toState];
