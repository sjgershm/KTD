function param = KTD_defparam
    
    param.c = 1;                % prior variance
    param.s = 1;            % noise variance
    param.q = 0.01;              % diffusion variance
    param.g = 0.98;             % discount factor
%     param.K = 6;                % number of microstimuli
%     param.sigma = 0.08;         % width of microstimuli
%     param.decay = 0.985;        % trace decay
    param.lr = 0.3;             % learning rate for standard TD
    param.TD = 0;               % use standard TD?