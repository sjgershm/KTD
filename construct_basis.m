function F = construct_basis(basis,X,param)
    
    if nargin < 3; param = KTD_defparam; end
    
    [N, D] = size(X);
    
    switch basis
        
        case 'microstim'
            I = repmat((1:param.K)',1,D)/param.K;
            y = zeros(param.K,D);
            F = zeros(N,param.K*D);
            f = @(y) (1/sqrt(2*pi))*exp(-((y(:)-I(:)).^2)/(2*param.sigma^2));
            for t = 1:N
                y = bsxfun(@plus,y,X(t,:));
                F(t,:) = f(y).*y(:);
                y = y*param.decay;
            end
            
        case 'CSC'
            F = [];
            for d = 1:D
                F = [F diag(X(:,d))];
            end
            
    end