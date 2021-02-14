function y = artransform(params)

% Mapping real valued parameters to stationary region
% Source: https://github.com/helske/KFAS/blob/master/src/artransform.f95

phi = tanh(params);
p = length(params);

u = zeros(p);
for i = 1:p
    u(i, i) = phi(i);
end
for i = 2:p
    for j = 1:i-1
        u(i, j) = u(i-1, j) - phi(i) * u(i-1, i-j);
    end
end
y = u(p,:);
