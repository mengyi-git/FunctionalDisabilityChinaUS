classdef MLFIPar
    properties (Constant)
        N_FRAILTY_PATH = 1000;
        year_t1 = 1998;
    end
    properties
        % -----------------------------------------------------------------
        % Derived from data
        % -----------------------------------------------------------------
        ageMin
        ageMax
        deltaT   % t_i - t_{i-1}
        N    % = numel(t)
        N_H_STATE   % no. of health states
        Q    % variance used in the latent process
        S    % no. of transition types
        t    % unique time points (in year) at which interviews are held
             % t = 1 in 1998
%         transitPair % S x 2: Col 1 = from state, Col 2 = to state
        m    % no. of covariates, excluding the intercept
        % -----------------------------------------------------------------
        % Random numbers
        % -----------------------------------------------------------------       
%         etaSim   % eta ~ NIID(0, Q)
    end
    

end