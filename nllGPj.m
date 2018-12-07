function [nll,dnlldhyp] = nllGPj(hyp,cov,X,y,p)
%NLLGPJ Computes the negative log likelihood of data set for GP
%
% In:
%    hyp           R*(E+1)+1 x 1    hyperparameter
%    cov           @fun             covariance function handle
%    X             E x N            Data points
%    y             1 x N (N x 1)    Data points
%    p             @{R x 1}         Cell array of function handles
% Out:
%    nll           1 x 1            Negative Log likelihood
%    dnlldhyp      1  x E+1         Derivative of nll wrt hyperparameter
% E: input dimensions
% R: Number of summands in cov function
% N: Number of training points

% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018


if nargin<5, p = []; end

[~,N] = size(X);
y = y(:); hyp = hyp(:);
% Verify Sizes
if size(y,1)~=N
    error('size mismatch');
end

if any(isinf(y)) || any(any(isinf(X))) ||...
    any(isnan(y)) || any(any(isnan(X)))
     error('Inf or NaN in training data');
end

if nargout >1
    [K,dKdhyp] = feval(cov,hyp(1:end-1),X,X,p);
else
    K = feval(cov,hyp(1:end-1),X,X,p);
end

sn2 = exp(2*hyp(end));
if isinf(sn2), error('sn2 is inf'); end

if sn2<1e-6
    Kn = K+1e-6*eye(N);
else
    Kn = K+sn2*eye(N);
end

% Work with Cholesky for numerical stability
L = chol(Kn);
alpha = L\(L'\y);

% From the GMPL book
% nll =  0.5*log(det(Kn)) +0.5*y'*(Kn\y) + (N/2)*log(2*pi);
% For comutational efficiency
nll =  sum(log(diag(L))) +0.5*y'*alpha + (N/2)*log(2*pi);


% Compute derivative if required
if nargout >1
    % Derivative wrt sn2 
    dKdhyp(:,:,end+1)=diag(ones(N,1)*sn2*2);
    
    % Preallocate 
    dnlldhyp = zeros(1,size(hyp,1));
    
    % Precompute for Speed
    Q = L\(L'\eye(N)) - alpha*alpha';     
    
    % Derivative wrt sf2 and ell
    for i = 1:length(dnlldhyp)
        % For comutational efficiency
        dnlldhyp(i) = 0.5*sum(sum(Q.*dKdhyp(:,:,i)));
        % From the GMPL book
        %dnlldhyp(i) = -0.5*trace(alpha*alpha'*dKdhyp(:,:,i) -  L\(L'\dKdhyp(:,:,i)));
    end
end

