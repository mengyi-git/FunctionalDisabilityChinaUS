function yrList = getYrList(digit, Par, lastyr)
% Concat string 'Year1 - Year2' where 
% Year1 and Year2 are interview years of two consecutive waves

N = Par.N;
yrArray = [Par.t + Par.year_t1 - 1; lastyr];

yrList = cell(1, N);
for tIndex = 1:N
    str1 = num2str(yrArray(tIndex));
    str1 = str1(end-digit+1:end);
    
    str2 = num2str(yrArray(tIndex+1));
    str2 = str2(end-digit+1:end);
    
    yrList{tIndex} = [str1, '-', str2];
 
end