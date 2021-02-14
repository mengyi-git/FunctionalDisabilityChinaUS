classdef MLFIMdl
    
    properties (Constant)
        N_FRAILTY_PATH = 1000; 
		year_t1 = 1998;		% t = 1 in year 1998
        default_rho00 = 2;
    end
    
    properties
%         datatable % don't include <-- it's too large --> slow down the code
        % -----------------------------------------------------------------
        % Model input - set by the user
        % -----------------------------------------------------------------
        dataName
        N_SIM = MLFIMdl.N_FRAILTY_PATH % no. of simulation
        latent = '' % rw or ar
        model
        transitTable %
        % -----------------------------------------------------------------
        % Initial value of the parameter - set by the user
        % -----------------------------------------------------------------  
        rho00  % initial value of rho        
        % -----------------------------------------------------------------
        % Derived from 'datatable'
        % -----------------------------------------------------------------        
        N_H_STATE   % no. of health states
        S           % no. of transition types, depends on 'datatable'
        t    		% unique time points (in year) at which interviews are held
					% t = 1 in 1998   
        % -----------------------------------------------------------------
        % Estimation results
        % -----------------------------------------------------------------
        thetaHat	% estimated theta
        se			% standard error of theta
        loglik		% maximised log likelihood function
    end
    
    properties (Dependent)
        deltaT      % t_i - t_{i-1}
        N           % = numel(t)
        Q           % variance used in the latent process
        % -----------------------------------------------------------------
        % Derived from 'model'
        % -----------------------------------------------------------------        
        matName
        m           % no. of covariates, excluding the intercept
        model_last	% to find theta_input
        thetaName
        thetaNo		% 2 elements: no. of coef and no. of rho's
        % -----------------------------------------------------------------
        % Derived from 'transitTable'
        % -----------------------------------------------------------------        
        transitPair
        % -----------------------------------------------------------------
        % Random numbers - depends on Q and N_SIM
        % -----------------------------------------------------------------        
        rho0
        etaSim      % eta ~ NIID(0, Q)
    end
    
    methods
        
% =========================================================================
% -------- 'set' methods --------------------------------------------------       
% =========================================================================
        function obj = set.dataName(obj, newName)
            % Control allowable Category values 
            nameList = {'clhls', 'rndhrs'};            
            switch newName 
                case nameList 
                    obj.dataName = newName;
                otherwise 
                    error(' Invalid data name! '); 
            end
        end
        
% -------------------------------------------------------------------------

        function obj = set.N_H_STATE(obj, datatable)
            state = [datatable.RxHSTATE, datatable.RxHSTATE2];
            obj.N_H_STATE = max(state(:));
        end
        
% -------------------------------------------------------------------------
       
        function obj = set.N_SIM(obj, nSim)
            if mod(nSim,1) == 0 && nSim > 0
                obj.N_SIM = nSim;
            else
                error('No. of simulations needs to be a positive integer!')
            end
        end

% -------------------------------------------------------------------------
        
        function obj = set.model(obj, newModel) 
            % Control allowable Category values 
            mList = {'static', 'trend', 'frailty', ...
                'static_resid', 'trend_resid', 'frailty_resid'};
            switch newModel 
                case mList 
                    obj.model = newModel;
                otherwise 
                    error('Invalid model name!'); 
            end 
        end
        
% -------------------------------------------------------------------------
       
        function obj = set.latent(obj, l)
            % Control allowable Category values 
            list = {'', 'rw', 'ar'};
            if any(strcmp(list, l))
                obj.latent = l;
            else
                error('Invalid latent name!'); 
            end
            
        end
        
% -------------------------------------------------------------------------
		
	function obj = set.rho00(obj, newRho)
		if isnumeric(newRho) && numel(newRho) == 1
			obj.rho00 = newRho;
		else
			error('Initial value of rho should be a number');
		end
	end
        
% -------------------------------------------------------------------------
        
        function obj = set.S(obj, datatable)
            state = [datatable.RxHSTATE, datatable.RxHSTATE2];
            transit = unique(state, 'rows');
            
            obj.S = sum(transit(:, 1) ~= transit(:, 2));
        end		
        
% -------------------------------------------------------------------------
        
        function obj = set.t(obj, datatable)
            obj.t = unique(datatable.TIME);
        end
        
% -------------------------------------------------------------------------
		
	function obj = set.thetaHat(obj, newTheta)
            obj.thetaHat = newTheta;
        end
        
% -------------------------------------------------------------------------
        
        function obj = set.transitTable(obj, newTable)
            if ~istable(newTable)
                error('transitTable needs to be a table!')
            else
                obj.transitTable = newTable;
            end         
        end
        
% -------------------------------------------------------------------------
        
% =========================================================================
% -------- 'get' methods --------------------------------------------------       
% =========================================================================
        
        function transitPair = get.transitPair(obj)
            transition = table2array(obj.transitTable);
            
            if max(transition(:)) ~= obj.S
                error('No. of transition types does not match the data.')
            elseif max(size(transition)) ~= obj.N_H_STATE
                error('No. of health states does not match the data.')
            end
            
            fromState = zeros(obj.S, 1);
            toState = zeros(obj.S, 1);

            for i = 1:obj.S
               [row, col] = find(transition == i);
               fromState(i) = row;
               toState(i) = col;
            end
            transitPair = [fromState, toState];
        end
        
% -------------------------------------------------------------------------
        
        function deltaT = get.deltaT(obj)
            deltaT = diff(obj.t);
            deltaT = [deltaT(1); deltaT];
        end
        
% -------------------------------------------------------------------------

        function etaSim = get.etaSim(obj)
	%Simulate the white noise in the latent process
			% eta ~ NIID(0, Q)
			% eta: N x N_SIM
			% Q: N x N
            switch obj.model
                
                case {'frailty', 'frailty_resid'}
                
                    if isempty(obj.N_SIM)
                        disp(['No. of simulations set to ', ...
                            num2str(obj.N_FRAILTY_PATH)])
                        obj = set.N_SIM(obj, obj.N_FRAILTY_PATH);
                    end
                    
                    rng(8) % set seed

                    eta = nan(obj.N, obj.N_SIM);
                    for iSim = 1:obj.N_SIM
                        eta(:, iSim) = diag(normrnd(0, obj.Q));
                    end

                    etaSim = eta;                    
                    
                otherwise                    
                    disp([obj.model, ' model. No simulation required.'])
            end            
        end
        
% -------------------------------------------------------------------------
        
        function m = get.m(obj)
            switch obj.model
                case 'static'
                    m = 2;
                case {'static_resid', 'trend'}
                    m = 3;
                case {'trend_resid', 'frailty'}
                    m = 4;
                case 'frailty_resid'
                    m = 5;
            end
        end

% -------------------------------------------------------------------------
        
        function matName = get.matName(obj)
            if ~any(strcmp({'frailty', 'frailty_resid'}, obj.model))
                matName = obj.model;
            else
                matName = obj.latent;
		if MLFIMdl.isResid(obj.model)
			matName = [matName, '_resid'];
		end

            end
        end
		
% -------------------------------------------------------------------------		
        
        function model_last = get.model_last(obj)
            switch obj.model
                case {'static', 'static_resid'}
                    model_last = [];
                case 'trend'
                    model_last = 'static';
                case 'trend_resid'
                    model_last = 'static_resid';
                case 'frailty'
                    model_last = 'trend';
                case 'frailty_resid'
                    model_last = 'trend_resid';
            end
        end
        
% -------------------------------------------------------------------------
        
        function N = get.N(obj)
            N = numel(obj.t);
        end
        
% -------------------------------------------------------------------------

        function rho0 = get.rho0(obj)
            rhoNo = obj.thetaNo(2);
            
            if any(strcmp(obj.model, {'frailty', 'frailty_resid'}))
                if rhoNo == 0
                    rho0 = ones(obj.N, 1);

                elseif rhoNo == 1
                    rho_tmp = obj.rho00;
					if isempty(rho_tmp)
						rho_tmp = obj.default_rho00;
					end
                    rho_tmp = artransform(rho_tmp); % <- ensure stationarity
                    rho0 = rho_tmp * ones(obj.N, 1);

                else
                    fprintf('No. of rho''s = %d\n', rhoNo)
                end
            else
                rho0 = [];
            end
        end
        
% -------------------------------------------------------------------------
 
        function Q = get.Q(obj)
            if obj.rho0(1) == 1
                Q = obj.deltaT;
            else
                Q = (1 - obj.rho0.^(2*obj.deltaT)) ./ (1 - obj.rho0.^2);
            end
            Q = diag(Q);
        end
        
% -------------------------------------------------------------------------
 
        function thetaName = get.thetaName(obj)
            if MLFIMdl.isResid(obj.model)
                thetaName = {'$\hat{\beta}_s$', '$\hat{\gamma}^\text{age}_s$', '$\hat{\gamma}^\text{female}_s$',...
                    '$\hat{\gamma}^\text{urban}_s$', '$\hat{\gamma}^\text{time}_s$', '$\hat{\alpha}_s$'};
            else
                thetaName = {'$\hat{\beta}_s$', '$\hat{\gamma}^\text{age}_s$', '$\hat{\gamma}^\text{female}_s$',...
                    '$\hat{\gamma}^\text{time}_s$', '$\hat{\alpha}_s$'};
            end
        end
        
% -------------------------------------------------------------------------
 
        function thetaNo = get.thetaNo(obj)
		if strcmp(obj.latent, 'ar')
			thetaNo = [(obj.m + 1)*obj.S, 1];
		else
			thetaNo = [(obj.m + 1)*obj.S, 0];
		end
        end

% -------------------------------------------------------------------------

% =========================================================================
% -------- other methods --------------------------------------------------       
% =========================================================================
        
        function [x, rho] = getCoef(obj, theta, thetaNo)
	%Get coef and rho: used in negloglik_frailty()
            if numel(theta) ~= sum(thetaNo)
                error('No. of parameters in theta does''t match!')
            end

            coefNo = obj.thetaNo(1);
            rhoNo = obj.thetaNo(2);


            coef = theta(1:coefNo);
            x = reshape(coef, [obj.S, numel(coef)/obj.S]);

            if any(strcmp(obj.model, {'frailty', 'frailty_resid'}))
                if rhoNo == 0
                    rho = ones(obj.N, 1);

                elseif rhoNo == 1
                    rho_tmp = theta(end);
                    rho_tmp = artransform(rho_tmp); % <- ensure stationarity
                    rho = rho_tmp * ones(obj.N, 1);

                else
                    fprintf('No. of rho''s = %d\n', rhoNo)
                end
            else
                rho = [];
            end
        end
        
% -------------------------------------------------------------------------        
        
        function X = getX(obj, datatable, type)
	%Get covariates, used in estimating parameters and recovering psi

            [AGE, AGETRS, FEMALE, ~, TIME, TIMETRS] = MLFIMdl.importAGT(datatable);
            
            if MLFIMdl.isResid(obj.model)
                URBAN = datatable.JOINURBAN;
            end

            N_OBS = size(datatable, 1);

            switch type
                case 'A' % A denotes the starting time
                    switch obj.model
                        case 'static'
                            X = [ones(1, N_OBS); AGE'; FEMALE'];

                        case {'trend', 'frailty'}
                            X = [ones(1, N_OBS); AGE'; FEMALE'; TIME'];

                        case 'static_resid'
                            X = [ones(1, N_OBS); AGE'; FEMALE'; URBAN'];

                        case {'trend_resid', 'frailty_resid'}
                            X = [ones(1, N_OBS); AGE'; FEMALE'; URBAN'; TIME'];

                        otherwise
                            error('Wrong model name!')
                    end

                case 'B' % B denotes time of transition
                    switch obj.model
                        case 'static'
                            X = [ones(1, N_OBS); AGETRS'; FEMALE'];

                        case {'trend', 'frailty'}
                            X = [ones(1, N_OBS); AGETRS'; FEMALE'; TIMETRS'];

                        case 'static_resid'
                            X = [ones(1, N_OBS); AGETRS'; FEMALE'; URBAN'];

                        case {'trend_resid', 'frailty_resid'}
                            X = [ones(1, N_OBS); AGETRS'; FEMALE'; URBAN'; TIMETRS'];

                        otherwise
                            error('Wrong model name!')
                    end

                otherwise
                    error('Type needs to be either A or B.')
            end
        end
        
% -------------------------------------------------------------------------        

        function loglambda = cal_loglambda(obj, theta, X)
	%Calculate log transition rate
            params = reshape(theta, [obj.S, numel(theta)/obj.S]);
            loglambda = params * X;
        end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------        
% -------- estimate parameters --------------------------------------------        
% -------------------------------------------------------------------------        

        function theta_input = findInitTheta_input(obj)
	%Obtain the input to determine theta0
			% empty for the static model
			% trend model <-- static model
			% frailty model <-- trend model
            if any(strcmp({'static', 'static_resid'}, obj.model))
                theta_input = [];
            else
                Result = load([obj.dataName, '_', obj.model_last, '.mat']);
                theta_input = Result.theta;    
            end
        end
        
% -------------------------------------------------------------------------        
       
        function thetaInit = findInitTheta_uniVar(obj, datatable)
        %Find theta0 to be used in fminunc() function
            
            %choose initial value
			theta_input = findInitTheta_input(obj);
 
            nPara = obj.thetaNo(1); % exclude rho

            start_index = numel(theta_input) + 1;

            thetaInit = zeros(nPara, 1);
            thetaInit(1:numel(theta_input)) = theta_input;

            % initialise next point
            next_pts = 0;

            next_pts_step_array = [4 2 1 0.5 0.1];
            for iParam = start_index:nPara
                tic
                fprintf('Parameter %d \n', iParam)

                param_index = iParam; % which parameter to check?    

                for iNextStep = 1:length(next_pts_step_array)

                    next_pts_step = next_pts_step_array(iNextStep);
                    grid_step = next_pts_step/10;

                    grid_low = next_pts - next_pts_step;
                    grid_high = next_pts + next_pts_step;

                    grid_pts = (grid_low:grid_step:grid_high);
                    nGrid = length(grid_pts);

                    % set the mid point = next_pts
                    % round to 2 decimal points to avoid floating error
                    grid_pts = round(grid_pts, 2);

                    tmp_theta = zeros(nPara, nGrid);
                    tmp_theta(1:numel(theta_input), :) = repmat(theta_input, [1, nGrid]);
                    for iGrid = 1:nGrid
                        tmp_theta(param_index, iGrid) = grid_pts(iGrid);
                    end

                    % negative log likelihood value
                    nll = zeros(1, nGrid);
                    for iGrid = 1:nGrid
                        nll(1, iGrid) = negloglik(obj, datatable, tmp_theta(:,iGrid), [nPara, 0]);
                    end

                    next_pts = grid_pts(nll==min(nll));

                    % multiple next_pts
                    if numel(next_pts) > 1
                        disp('Multiple minimum points found.')
                        fprintf('Minimum nll is %f \n', min(nll))
                        disp('Next points are:')
                        disp(next_pts)

                        next_pts = next_pts(1);
                    end

                end

                thetaInit(param_index) = next_pts;

                toc
            end
			
			if obj.thetaNo(2) ~= 0
				if isempty(obj.rho00)
					thetaInit = [thetaInit; obj.default_rho00];
				else
					thetaInit = [thetaInit; obj.rho00];
				end
			end
        end  
        
% -------------------------------------------------------------------------        

        function path = simFrailtyPath(obj, rho, wn)
	%Simulate latent process paths used for estimation
		
            phi = rho .^ obj.deltaT;
            
			% initialised at zero
            path_tmp = zeros(obj.N+1, obj.N_SIM);
            
            % forward iteration
            for iSim = 1:obj.N_SIM
                for tIndex = 2:obj.N+1
                    path_tmp(tIndex, iSim) = ...
                        phi(tIndex-1) * path_tmp(tIndex-1, iSim) + ...
                        wn(tIndex-1, iSim);
                end
            end

            path = path_tmp(2:end, :);
        end

% -------------------------------------------------------------------------        
        
        function output = negloglik(obj, datatable, theta, thetaNo)
	%Wrap negloglik_noFrailty and negloglik_frailty
        
            if any(strcmp(obj.model, {'frailty', 'frailty_resid'}))
                % frailty model
                output = negloglik_frailty(obj, datatable, theta, thetaNo);

            else    
                % non-frailty model
                output = negloglik_noFrailty(obj, datatable, theta);
            end
        end
        
% -------------------------------------------------------------------------        
      
        function output = negloglik_noFrailty(obj, datatable, theta)
	%Calculate the negative log likelihood function for non-frailty models

            TAU = datatable.TAU;

            XB = getX(obj, datatable, 'B');
            loglambdaB = cal_loglambda(obj, theta, XB);

            Y = MLFIMdl.importY(datatable);
            term1 = Y .* loglambdaB;

            XA = getX(obj, datatable, 'A');
            loglambdaA = cal_loglambda(obj, theta, XA);

            R = MLFIMdl.importR(datatable);
            term2 = R .* exp(loglambdaA) .* repmat(TAU', [obj.S, 1]);

            temp = term1 - term2;

            output = -sum(temp(:));
        
        end

% -------------------------------------------------------------------------

        function obj = estTheta(obj, datatable)
	%Estimate the parameters
		
			% 1. Obtain the initial values
            theta0 = findInitTheta_uniVar(obj, datatable);
            save([obj.dataName, '_', obj.matName, '_theta0', '.mat'], 'theta0')
                        
			% 2. Set the objective function
            nll = @(theta)negloglik(obj, datatable, theta, obj.thetaNo);
            
            % 3. Find the maximised log-likelihood function - use two steps
            % to improve the estimation accuracy 
            
            % 3.1. fminunc()
            options_unc = ...
                optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
                'MaxFunEvals', 20000, ...
                'Display', 'iter', ...
                'StepTolerance', 1e-10);

            [theta, fval, exitflag, output, grad, hessian] = fminunc(nll, theta0, options_unc);
            save([obj.dataName, '_', obj.matName, '_tmp', '.mat'], ...
                'theta', 'fval', 'exitflag', 'output', 'grad', 'hessian')

            % 3.2. fminunc(): 'central' to improve accuracy
            
            % reset theta0 to a random number drawn from a normal distribution
            theta0 = normrnd(theta, sqrt(diag(inv(hessian))));

            options_unc.FiniteDifferenceType = 'central';

            [theta, fval, exitflag, output, grad, hessian] = fminunc(nll, theta0, options_unc);
            save([obj.dataName, '_', obj.matName, '.mat'], ...
                'theta', 'fval', 'exitflag', 'output', 'grad', 'hessian')
        end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------- examine estimated parameters -----------------------------------
% -------------------------------------------------------------------------
        
        function obj = importResult(obj)
        %Import estimation result from .mat file
            Est = load([obj.dataName, '_', obj.matName, '.mat']);
			obj.thetaHat = Est.theta;
			obj.se = MLFIMdl.calStdErr(Est.hessian);
			obj.loglik = -Est.fval;            
        end
        
% -------------------------------------------------------------------------
		       
        function printTable(obj)
        %Print estimated parameters with se and significance levels

            theta = obj.thetaHat;

            % hypothesis test to get the p values
            pval = MLFIMdl.calPVal(theta, obj.se); % p-value

            % significance level
            star = arrayfun(@MLFIMdl.calSignt, pval, 'UniformOutput', false);
            
            theta_trs = reshape(theta, [obj.S, numel(theta)/obj.S])';
            se_trs = reshape(obj.se, [obj.S, numel(obj.se)/obj.S])';
            star_trs = reshape(star, [obj.S, numel(star)/obj.S])';

            nParamType = numel(theta)/obj.S;
            nRow = nParamType * 2;
            nCol = obj.S * 2;
            
            % print to LaTeX table
            for iRow = 1:nRow % nrow
                if mod(iRow, 2) == 1
                    fprintf('%30s', obj.thetaName{(iRow+1)/2})        
                end
                fprintf (' & ')

                for jCol = 1:nCol % ncol
                    if mod(iRow, 2) == 1 % coeff
                        if mod(jCol, 2) == 1 % coeff
                            fprintf('%8.4f', theta_trs((iRow+1)/2, (jCol+1)/2))
                        else
                            fprintf('%3s', star_trs{(iRow+1)/2, jCol/2})
                        end
                    else
                        if mod(jCol, 2) == 1 % std err
                            fprintf('(%6.4f)', se_trs(iRow/2, (jCol+1)/2))
                        else
                            fprintf('    ')  
                        end
                    end

                    if jCol ~= obj.S * 2
                        fprintf(' & ')
                    end
                end
                fprintf('\\\\ \n')
            end
            fprintf('Log likelihood & %8.0f ', obj.loglik)
            for jCol = 2:nCol % ncol
                fprintf(' & ')
            end
            fprintf('\\\\ \n')            
        end
        
% -------------------------------------------------------------------------
    
        function [a, b] = calAIC(obj, n)
        %Calculate AIC BIC
        
            obj = importResult(obj);
            k = sum(obj.thetaNo);

            a = 2*k - 2*obj.loglik; % aic
            b = log(n)*k - 2*obj.loglik; % bic
            
        end
        
% -------------------------------------------------------------------------

        function [psi_hat, kf_V] = recoverPsi(obj, datatable)
	%Recover the frailty process
	
            if strcmp(obj.latent, 'ar')
                error('The code is written for the random walk process')
            end
            
            timeArray = obj.t; %unique(datatable.TIME);
            N_TIME = obj.N; 
            
            % load estimation results: theta
            obj = importResult(obj);
            theta = obj.thetaHat;            
            
            % rearrange theta and extract the parameters of psi
            par_est = reshape(theta, [obj.S, numel(theta)/obj.S]);
            
            par_est_noPsi = par_est(:, 1:end-1); % exclude the parameters of psi            
            par_est_psi = par_est(:, end);
            
            % ---- load X, Y, R from 'datatable' ----
            % load covariates
            X_tmp = getX(obj, datatable, 'A');   
            X_tmp = X_tmp';
            
            % R_S1, ..., R_S9
            % datatable.RA_H RA_MD RA_SD
            R = MLFIMdl.importR(datatable);
            R = R';
            
            % y_S1, ..., y_S9 (9 types of transitions)
            Y_S = MLFIMdl.importY(datatable); 
            Y_S = Y_S';            
            % ---- end of load ----
            
            y_S = Y_S - X_tmp * par_est_noPsi';
            
            % initialise psi
            psi_frailty = zeros(N_TIME, 1);
            
            time_table = datatable(:, 'TIME');
            for i = 1:20 % iteration count: set to a larger number if covergence is slow
                
            	T = MLFIMdl.mergeFrailtyPath(obj.t, psi_frailty, time_table);                
            
            	% calculate lambda(t)
                X = [X_tmp, T.FRAILTY];
                ln_lambda_t = X * par_est';
            	lambda_t = exp(ln_lambda_t);
            
                
                % -------------------------------------------
                % calculate zeta and kappa2 given psi_frailty
                % -------------------------------------------            
          
                % allocate space
            	kappa2 = zeros(N_TIME, 1);
            	zeta = zeros(N_TIME, 1);
            
            	y_S_mean0 = zeros(size(y_S));
                
            	% Prepare: calculate kappa^2
            	kappa2_n_allWave = R * par_est_psi.^2; % kappa2 nominator all waves            
            	kappa2_d_allWave = R .* repmat(datatable.TAU, [1, obj.S]) .* lambda_t * par_est_psi.^2; % kappa2 denominator all waves
            
            	% Prepare: calculate zeta
            	zeta_n_term1 = Y_S .* R * par_est_psi;            
            	zeta_n_term2 = ln_lambda_t .* R * par_est_psi;            
                
                for iTime = 1:N_TIME                
                	time = timeArray(iTime);
                	wvIndex = datatable.TIME == time;
                
                	kappa2(iTime) = sum(kappa2_n_allWave(wvIndex)) / sum(kappa2_d_allWave(wvIndex));
                	
                	zeta_n_term3 = kappa2(iTime) * (Y_S * par_est_psi - ...
                		R .* repmat(datatable.TAU, [1, obj.S]) .* lambda_t * par_est_psi);
                
                	zeta_n_allWave = zeta_n_term1 - zeta_n_term2 - zeta_n_term3; % zeta nominator all waves
                	zeta_d_allWave = R * par_est_psi; % zeta denominator all waves
                
                	zeta(iTime) = sum(zeta_n_allWave(wvIndex)) / sum(zeta_d_allWave(wvIndex));
                	
                	% subtract mean 
                	y_S_mean0(wvIndex, :) = y_S(wvIndex, :) - zeta(iTime);
                end           
            
            
                % ----------------
            	% Kalman filtering
                % ----------------
            
            	% ---- initialise ----
            	kf_aCell = cell(N_TIME, 1);
            	kf_PCell = cell(N_TIME, 1);
            	kf_KCell = cell(N_TIME, 1);
            	kf_FCell = cell(N_TIME, 1);
            	kf_vCell = cell(N_TIME, 1);
                
                for iTime = 1:N_TIME            		
                	time = timeArray(iTime);
                	wvIndex = datatable.TIME == time;
                	
                	wvObs = sum(wvIndex);
                	
                	kf_aCell{iTime} = nan(wvObs*obj.S+1, 1);
                	kf_PCell{iTime} = nan(wvObs*obj.S+1, 1);
                	kf_KCell{iTime} = nan(wvObs*obj.S, 1);
                	kf_FCell{iTime} = nan(wvObs*obj.S, 1);
                	kf_vCell{iTime} = nan(wvObs*obj.S, 1);            		
                end
                
            	kf_aCell{1}(1) = 0;
            	kf_PCell{1}(1) = 1;
                
                % ---- end of initialise ----
            
                % ---- forward iteration ----
                diagQ = diag(obj.Q);
                for iTime = 1:N_TIME
                	
                	time = timeArray(iTime);
                	wvIndex = datatable.TIME == time;
                	
                	wvObs = sum(wvIndex);
                
                	kf_a = kf_aCell{iTime};
                	kf_P = kf_PCell{iTime};
                	kf_K = kf_KCell{iTime};
                	kf_F = kf_FCell{iTime};
                	kf_v = kf_vCell{iTime};
                	
                	y_S_mean0_wv = y_S_mean0(wvIndex, :);
                	R_wv = R(wvIndex, :);
                	
                	for iObs = 1:wvObs
                		
                		for iTransit = 1:obj.S
                			index = iTransit + (iObs - 1) * obj.S;
                
                			if R_wv(iObs, iTransit) == 1
                				kf_F(index) = kf_P(index) * par_est_psi(iTransit)^2 + kappa2(iTime);
                				kf_K(index) = kf_P(index) * par_est_psi(iTransit) / kf_F(index);
                				kf_v(index) = y_S_mean0_wv(iObs, iTransit) - par_est_psi(iTransit) * kf_a(index);
                				kf_a(index+1) = kf_a(index) + kf_K(index) * kf_v(index);
                				kf_P(index+1) = kf_P(index) - kf_K(index)^2 * kf_F(index);
                			else
                				kf_a(index+1) = kf_a(index);
                				kf_P(index+1) = kf_P(index);
                			end
                		end
                		
                	end
                	
                	kf_aCell{iTime} = kf_a;
                	kf_PCell{iTime} = kf_P;
                	kf_KCell{iTime} = kf_K;
                	kf_FCell{iTime} = kf_F;
                	kf_vCell{iTime} = kf_v;
                	
                	if iTime < N_TIME
                		kf_aCell{iTime+1}(1) = kf_aCell{iTime}(end);
                		kf_PCell{iTime+1}(1) = kf_PCell{iTime}(end) + diagQ(iTime);                
                	end
                end            
            
            	% cellfun(@(v)v(1), kf_aCell)
            	% cellfun(@(v)v(1), kf_PCell)
            
                % ----------------
            	% Kalman smoothing
                % ----------------
                
                % ---- initialise ----
            	kf_rCell = cell(N_TIME, 1);
            	kf_NCell = cell(N_TIME, 1);
            	kf_LCell = cell(N_TIME, 1);
            	for iTime = 1:N_TIME
            		
            		time = timeArray(iTime);
            		wvIndex = datatable.TIME == time;
            		
            		wvObs = sum(wvIndex);
            		
            		kf_rCell{iTime} = nan(wvObs*obj.S+1, 1);
            		kf_NCell{iTime} = nan(wvObs*obj.S+1, 1);
            		kf_LCell{iTime} = nan(wvObs*obj.S, 1);    
            	end
            
            	kf_rCell{end}(end) = 0;
            	kf_NCell{end}(end) = 0;
                % ---- end of initialise ----
                
            
                % ---- backward iteration ----
            	for iTime = N_TIME:-1:1
                    
            		time = timeArray(iTime);
            		wvIndex = datatable.TIME == time;
            		
            		wvObs = sum(wvIndex);
            		
            		R_wv = R(wvIndex, :);
            
            		kf_r = kf_rCell{iTime};
            		kf_N = kf_NCell{iTime};
            		kf_L = kf_LCell{iTime};
            
            		% from the filtering step
            		kf_K = kf_KCell{iTime};
            		kf_F = kf_FCell{iTime};
            		kf_v = kf_vCell{iTime};
            		
            		for iObs = wvObs:-1:1
            			
            			for iTransit = obj.S:-1:1
            				index = iTransit + (iObs - 1) * obj.S;
            
            				if R_wv(iObs, iTransit) == 1
            				   kf_L(index) = 1 - kf_K(index) * par_est_psi(iTransit);
            				   kf_r(index) = kf_v(index) * par_est_psi(iTransit) / kf_F(index) + ...
            					   kf_L(index) * kf_r(index+1);
            				   kf_N(index) = par_est_psi(iTransit)^2 / kf_F(index) + ...
            					   kf_L(index)^2 * kf_N(index+1);
            				else
            					kf_r(index) = kf_r(index+1);
            					kf_N(index) = kf_N(index+1);
            				end
            			end
            			
            		end
            		
            		kf_rCell{iTime} = kf_r;
            		kf_NCell{iTime} = kf_N;
            		kf_LCell{iTime} = kf_L;    
            		
            		if iTime > 1
            			kf_rCell{iTime-1}(end) = kf_rCell{iTime}(1);
            			kf_NCell{iTime-1}(end) = kf_NCell{iTime}(1);
            		end
            	end
            
            	kf_a1 = cellfun(@(v)v(1), kf_aCell);
            	kf_P1 = cellfun(@(v)v(1), kf_PCell);
            	kf_r1 = cellfun(@(v)v(1), kf_rCell);
            	kf_N1 = cellfun(@(v)v(1), kf_NCell);
            
            	psi_hat = kf_a1 + kf_P1 .* kf_r1;
            	kf_V = kf_P1 - kf_P1.^2 .* kf_N1;
            
            	fprintf('Iteration %d, ', i)
            	fprintf('Error = %f \n', norm(psi_frailty - psi_hat))
            	psi_frailty = psi_hat;
            
            	fprintf('%f\n', psi_hat')            
            end 
            
            save(['psi', '_', obj.dataName, '_', obj.matName, '.mat'], 'psi_hat', 'kf_V')
            
        end
                
% -------------------------------------------------------------------------
% ---- use estimated parameters to calculate transition rate/probability --
% -------------------------------------------------------------------------

        function X = createX(obj, Input)
        %Create covariate to calculate transition rate/probability
        
            ageBeg = Input.ageBeg;
            ageEnd = Input.ageEnd;
            gender = Input.gender;
            residence = Input.residence;
            time_ageStart = Input.time_ageStart;
            iPsiPath = Input.iPsiPath;
            type = Input.type;
            
            n = ageEnd - ageBeg + 1;
            
            % age
            x_age = MLFIMdl.transformAge(ageBeg:ageEnd);
            
            % gender
            switch gender
                case {'female', 'f'}
                    x_gender = ones(1, n);
                case {'male', 'm'}
                    x_gender = zeros(1, n);
                otherwise
                    error('Wrong gender input')
            end
            
            %
            if MLFIMdl.isResid(obj.model)
                switch residence
                    case {'urban', 'u'}
                        x_resid = ones(1, n);
                    case {'rural', 'r'}
                        x_resid = zeros(1, n);
                    otherwise
                        error('Wrong residence input')
                end
            else
                x_resid = [];
            end
            
            % time if not static model
            if ~any(strcmp({'static', 'static_resid'}, obj.model))
                if strcmp(type, 'period')
                    x_time = MLFIMdl.transformTime(time_ageStart*ones(1, n));
                elseif strcmp(type, 'cohort')
                    x_time = MLFIMdl.transformTime(time_ageStart + (0:n-1));
                else
                    error('Wrong type!')
                end
            else
                x_time = [];
            end
            
            % psi path if frailty model
            if any(strcmp({'frailty', 'frailty_resid'}, obj.model))
                % ensure iPsiPath is a row vector
                if ~isvector(iPsiPath)
                    error('Psi path should be a vector')
                elseif ~isrow(iPsiPath)
                    iPsiPath = iPsiPath';
                end
                
                % psi path
                if strcmp(type, 'period')
                    x_psi = iPsiPath(1)*ones(1, n);
                elseif strcmp(type, 'cohort')
                    x_psi = iPsiPath;
                else
                    error('Wrong type!')
                end
            else
                x_psi = [];
            end
            
            X = [ones(1, n); x_age; x_gender; x_resid; x_time; x_psi];            
        end

% -------------------------------------------------------------------------   
        
	function trsRateCell = calTrsRate(obj, Input)
        %Calculate transition rate matrix using estimated parameters
            % save in a cell
        
            ageBeg = Input.ageBeg;
            ageEnd = Input.ageEnd;
            
            n = ageEnd - ageBeg + 1;
            trsRateCell = cell(n, 2); % Col 1: age; Col 2: trs rate matrix
            
            X = createX(obj, Input);                               
            loglambda = cal_loglambda(obj, obj.thetaHat, X);            
            lambda = exp(loglambda);
            
            for i = 1:n
                trsRateCell{i, 1} = ageBeg + i - 1;
                
                % arrange it into matrix
                mat = zeros(obj.N_H_STATE);
                for s = 1:obj.S
                    row = obj.transitPair(s, 1); col = obj.transitPair(s, 2);
                    mat(row, col) = lambda(s, i);
                end
                mat = mat - diag(sum(mat, 2));
                
                trsRateCell{i, end} = mat;
            end

        end
        
% -------------------------------------------------------------------------

        function trsProbCell = calTrsProb(obj, Input)
        %Calculate transition probability using estimated parameters

	    ageBeg = Input.ageBeg;
	    ageEnd = Input.ageEnd;
            n = ageEnd - ageBeg + 1;
            trsProbCell = cell(n, 2); % Col 1: age; Col 2: trs rate matrix
            
            trsRateCell = calTrsRate(obj, Input);
            for i = 1:n
                trsProbCell{i, 1} = ageBeg + i - 1;
                trsProbCell{i, end} = expm(trsRateCell{i, end});
            end
        end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% --------  micro-simulation ----------------------------------------------
% -------------------------------------------------------------------------
        
        function path = simulatePath(obj, Input, lastyr, nSimPsi)
        %Simulate frailty paths in micro-simulation
        
            % load psi_hat
            PsiHat = load(['psi', '_', obj.dataName, '_', obj.matName, '.mat']);
            psiHat = PsiHat.psi_hat;
            
            psiHat = [0; psiHat];
            psiBeg = psiHat(Input.time_ageStart == [obj.t + obj.year_t1 - 1; lastyr]);
                        
            n = Input.ageEnd - Input.ageBeg + 1;            
            
            if strcmp(obj.latent, 'ar')
                rho_hat = obj.thetaHat(end);
            elseif strcmp(obj.latent, 'rw')
                rho_hat = 1;
            end        
            
            path = zeros(n, nSimPsi);
            path(1, :) = psiBeg;  % initialised at posterior mean
            
            rng(2*obj.N) % set seed
            for tIndex = 2:n
                path(tIndex, :) = rho_hat * path(tIndex-1, :) + randn(1, nSimPsi);
            end        
        
        end
        
% -------------------------------------------------------------------------    

        function hStateSim = simHState(obj, Input, initHVector)
        % Simulate health state path
            trsProbCell = calTrsProb(obj, Input);

	    ageBeg = Input.ageBeg;
	    ageEnd = Input.ageEnd;
            n = ageEnd - ageBeg + 1;
            
            if ~isvector(initHVector)
                error('Initial health state should be a vector!')            
            elseif iscolumn(initHVector)
                initHVector = initHVector';
            end
            
            nSim = numel(initHVector);            
            
            hStateSim = zeros(n, nSim);
            hStateSim(1, :) = initHVector;

            % set seed
            rng(2*obj.N-1) 

            for tIndex = 2:n
                transitProb_t = trsProbCell{tIndex-1, end};

                state = hStateSim(tIndex-1, :)';
                new = mnrnd(ones(nSim, 1), transitProb_t(state, :));

                % find first 1 in each row
                [~, newState] = max(new, [], 2);

                hStateSim(tIndex, :) = newState';
            end            
        end
        
    end % end of methods
    
    
    methods (Static)
        
	function output = isResid(model)
	% The model with residence always end with '_resid'

		newStr = split(model, '_');
		if strcmp(newStr{end}, 'resid')
			output = 1;
		else
			output = 0;
		end
        end		  
        
% -------------------------------------------------------------------------
        
% -------------------------------------------------------------------------
% --------  analyse micro-simulation results ------------------------------
% -------------------------------------------------------------------------

        function y = calT(hStateSim, initHVector)
        %Calculate future lifetime (random variable) spent in each state given 
        %  simulated health paths and initial health state
            if size(hStateSim, 2) ~= size(initHVector, 2)
                error('No. of rows do not match.')
            end
            
            hStateVec = unique(hStateSim);
            N_H_STATE = numel(hStateVec);
            
            y = cell(1, N_H_STATE);
            
            % future lifetime in each state
            for hState = 1:N_H_STATE-1
                % 1/2 means we assume the transition occurs at the middle of the year
                y{hState} = 1/2*(initHVector == hState) + ...
                    sum(hStateSim(2:end, :) == hState, 1);
            end

            tmp = y(1:end-1);
            tmp_array = cat(3, tmp{:});
            
            y{end} = sum(tmp_array, 3); % total future lifetime
        end
        
% -------------------------------------------------------------------------    
        
        function y = calOnsetX(ageBeg, hStateSim, enterState)
        %Calculate onset of disability conditional becoming disabled
            [~, B] = max(hStateSim==enterState, [], 1);
            y = B(B>1) + ageBeg - 1/2;
        end
        
% -------------------------------------------------------------------------    
       
        function y = calRatio(hStateSim, initHVector)
        %Calculate HLE/TLE               
            lifeTime = MLFIMdl.calT(hStateSim, initHVector);
            y = lifeTime{1} ./ lifeTime{end};
        end
        
% -------------------------------------------------------------------------    
        
        function y = countHState(hStateSim)
        %Count no. of people in each state at each age
        
            % no. of health states
            N_H_STATE = numel(unique(hStateSim));
            n = size(hStateSim, 1);     % age dimension
            y = zeros(n, N_H_STATE);    % sim dimension
            
            for i = 1:n
                for iHState = 1:N_H_STATE
                    y(i, iHState) = sum(hStateSim(i, :)==iHState);
                end
            end
        end       
        
% -------------------------------------------------------------------------    
        
% -------------------------------------------------------------------------
% --------  extract variables from datatable ------------------------------
% -------------------------------------------------------------------------
	
        function [AGE, AGETRS, FEMALE, TAU, TIME, TIMETRS] = importAGT(datatable)		
			
	    AGE = MLFIMdl.transformAge(datatable.RxAGE);
	    AGETRS = MLFIMdl.transformAge(datatable.RxAGETRS);

            FEMALE = datatable.RAFEMALE;

            TAU = datatable.TAU;

	    TIME = MLFIMdl.transformTime(datatable.TIME);
	    TIMETRS = MLFIMdl.transformTime(datatable.TIMETRS);
        end
        
% -------------------------------------------------------------------------    
        
        function R = importR(datatable)
            RA_S1 = datatable.R_H;
            RA_S2 = datatable.R_D;
            RA_S3 = RA_S1;
            RA_S4 = RA_S2;

            R = [RA_S1'; RA_S2'; RA_S3'; RA_S4'];            
        end
        
% -------------------------------------------------------------------------    
        
        function Y = importY(datatable)
            Y_S1 = datatable.Y_S1;
            Y_S2 = datatable.Y_S2;
            Y_S3 = datatable.Y_S3;
            Y_S4 = datatable.Y_S4;

            Y = [Y_S1'; Y_S2'; Y_S3'; Y_S4'];
        end
		
% -------------------------------------------------------------------------    
		
% -------------------------------------------------------------------------    
% -------- Analyse estimated parameters	-----------------------------------
% -------------------------------------------------------------------------    

	function se = calStdErr(hessian)
	%Calculate standard errors given the hessian matrix

		% check positive definite
		try 
			chol(hessian);
		%	disp('Calculate the standard error')
    		catch
			disp('Matrix is not symmetric positive definite')
		end

		fisher = inv(hessian);
		se = sqrt(diag(fisher));		
        end
        
% -------------------------------------------------------------------------    
		
        function pval = calPVal(theta, se)
	%Calculate p-value
            % hypothesis test
            wald = theta ./ se;
            pval = 1 - normcdf(abs(wald)); % p-value
        end
        
% -------------------------------------------------------------------------    
        
	function star = calSignt(pval)
	%Match p-value to stars

		if pval < 0.01 
			star = '$^{***}$';
		elseif pval < 0.05
			star = '$^{**}$';
		elseif pval < 0.1
			star = '$^{*}$';
		else
			star = '';
		end

        end
        
% -------------------------------------------------------------------------    

% -------------------------------------------------------------------------    
% -------- methods used in parameter estimation ---------------------------
% -------------------------------------------------------------------------    

        function newT = mergeFrailtyPath(t, iFrailty, oldT)
	%Merge simulated frailty path with oldT based on TIME

            frailtyTable = table;
            frailtyTable.TIME = t;
            frailtyTable.FRAILTY = iFrailty;

            newT = join(oldT, frailtyTable);
        end
		
% -------------------------------------------------------------------------    
		
	function newAge = transformAge(oldAge)
		newAge = (oldAge - 65)/10;
	end
		
% -------------------------------------------------------------------------

	function newTime = transformTime(oldTime)
		denom = 10;
		if numel(num2str(oldTime(1))) == 4 % 4-digit year format
			newTime = (oldTime - MLFIMdl.year_t1 + 1)/denom;
		else
			newTime = oldTime/denom;
		end
        end
		
% -------------------------------------------------------------------------
		
    end % end methods(Static)
    
end % end classdef
