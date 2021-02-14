function output = negloglik_frailty(obj, datatable, theta, thetaNo)
%Calculate the negative log likelihood function for frailty models
% The function is taken out of the MLFIMdl.m for faster computational speed

TAU = datatable.TAU;
N_TRS = obj.S;

Y = MLFIMdl.importY(datatable);
R = MLFIMdl.importR(datatable);            

[params, rho] = getCoef(obj, theta, thetaNo);

% simulate frailty path            
path = simFrailtyPath(obj, rho, obj.etaSim);

M = size(path, 2); 

lnpl = zeros(M, 1);

XB_tmp = getX(obj, datatable, 'B');
XA_tmp = getX(obj, datatable, 'A');

time_table = datatable(:, 'TIME');
uniqueTime = obj.t;

for iPath = 1:M % turn-on parfor to utilise parallel computing

    % load frailty; merge to datatable
    T = MLFIMdl.mergeFrailtyPath(uniqueTime, path(:, iPath), time_table); 
    FRAILTY = T.FRAILTY;

    XB = [XB_tmp; FRAILTY'];
    loglambdaB = params * XB;

    term1 = Y .* loglambdaB;

    XA = [XA_tmp; FRAILTY'];
    loglambdaA = params * XA;

    term2 = R .* exp(loglambdaA) .* repmat(TAU', [N_TRS, 1]);

    temp = term1 - term2;

    lnpl(iPath) = sum(temp(:)); % log likelihood for one observation

end

output0 = logsumexp(lnpl) - log(M);
output = -output0;


end % end function
