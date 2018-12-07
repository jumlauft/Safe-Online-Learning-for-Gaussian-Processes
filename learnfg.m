function [muf,mug,hyps] = learnfg(xtr,ytr,mf,mg,U,hyps,uhyp)
%LEARNFG Estimates f and g from a cl measurement of a control affine system
% In:
%   xtr    E x Ntr          starting point of reference trajectory
%   ytr    E x Ntr          final point of reference trajectory
%   mf     @fun             prior mean function for f
%   mg     @fun             prior mean function for g
%   U      2 x Ntr          multiplier for for f and g
%   hyps   2*(E+1)+1 x 1    hyperparameters
%   uhyp   boolean          indicates if hyperparameters should be updated
% Out:
%   muf     @fun            GP  posterior mean estimate for f
%   mug     @fun            GP  posterior mean estimate for g
%   hyps    E+1 x N         hyperparameters
%
% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018

[E,Ntr] = size(xtr);

if Ntr ==0, muf = mf; mug = mg; return; end
sn = exp(hyps(end));
if uhyp
    hyps = minimize(hyps, @nllGPj, -30,'covSEaddmultj',xtr,ytr - mf(xtr).*U(1,:) - mg(xtr).*U(2,:) ,U);
end
Ktrtr = covSEaddmultj(hyps(1:end-1), xtr,xtr,U,U);
iKn = eye(Ntr)/(Ktrtr+eye(Ntr)*sn);
beta = iKn*(ytr - mf(xtr).*U(1,:) - mg(xtr).*U(2,:))';

muf = @(x) mf(x) + ((U(1,:).*covSEaddmultj(hyps(1:E+1), x, xtr))*beta)';
mug = @(x) mg(x) + ((U(2,:).*covSEaddmultj(hyps(E+2:end-1), x, xtr))*beta)';

end
