function model = kalmanTD(X,r,param)
    
    if nargin < 3; param = KTD_defparam; end
    
    % initialization
    [N,D] = size(X);
    w = zeros(D,1);
    X = [X; zeros(1,D)];    % add buffer at end
    
    % parameters
    if nargin < 3 || isempty(param); param = KTD_defparam; end
    C = param.c*eye(D); % prior covariance
    s = param.s;        % noise variance
    g = param.g;        % discount factor
    
    if length(s)==1; s = zeros(N,1)+s; end
    if length(param.q) ==1; param.q = zeros(N,1)+param.q; end
    if param.TD && length(param.lr)==1; param.lr = zeros(N,1)+param.lr; end
    
    % run Kalman filter
    for n = 1:N
        
        Q = param.q(n)*eye(D); % transition covariance
        h = X(n,:) - g*X(n+1,:);    % temporal difference features
        rhat = h*w;
        dt = r(n) - rhat;            % prediction error
        C = C + Q;                  % a priori covariance
        P = h*C*h'+s(n);            % residual covariance
        K = C*h'/P;                 % Kalman gain
        w0 = w;
        if param.TD
            w = w + param.lr(n)*h'*dt;
        else
            w = w + K*dt;               % weight update
        end
        C = C - K*h*C;              % posterior covariance update
        
        % store results
        model(n).w = w;
        model(n).w0 = w0;
        model(n).C = C;
        model(n).K = K;
        model(n).dt = dt;
        model(n).rhat = rhat;
        
    end