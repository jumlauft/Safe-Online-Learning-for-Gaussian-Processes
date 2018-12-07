% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 12/2018

clear; close all; clc;rng default; addpath(genpath('./gpml'));
%% Set Parameters
disp('Setting parameters...')

% Basic Parameters
Nstep = 40;         % Number of model updates to simulate
dt = 0.5;           % Time triggered update interval
sn = 1e-2;          % Observation noise (std deviation)
E = 2;              % State space dimension
x0 = [3; 2];        % initial State
Nuhyp = 1;          % Updates of hyperparameters ever Nuhyp-th update

% Reference Trajectory
reffun = @(t) refGeneral(t,E+1,@(tau) 1-sigmf(tau,[20 10]));  % smooth step

 % Controller gains
pFeLi.lam = ones(E-1,1); 
pFeLi.kc = 1;

% Define Systemdynamics
pdyn.f = @(x) 1-sin(x(1,:)) + 0.5*sigmf(x(2,:),[0.1 0]);
pdyn.g = @(x) 1+0.5*sin(0.5*x(2,:));

% GP prior mean functions and initial hyperparameters
mg = @(x) 2*ones(1,size(x,2));
mf = @(x) zeros(1,size(x,2));
hyps0 = [log(3); log(3); log(5); log(3); log(3); log(5); log(sn)];


% Visualization
setname = 'fgonline';
Nte = 1e4; XteMin = -4; XteMax = 4;


%% Learn fGP and gGP
disp('Learn fGP and gGP...')
Ntrf = 1e3;
Xtef = ndgridj(XteMin, XteMax,floor(nthroot(Ntrf,E))*ones(E,1)) ;
gprMdl = fitrgp(Xtef',pdyn.f(Xtef)); pdyn.f = @(x) predict(gprMdl,x')';
gprMdl = fitrgp(Xtef',pdyn.g(Xtef)); pdyn.g = @(x) predict(gprMdl,x')';


%% Simulating Adaptive Controller
disp('Simulating Controller...')
% Initialization
Xtr = []; Ytr = []; utr = [];pFeLi.f = mf; pFeLi.g = mg;hyps = hyps0;
X = x0; T = 0;

for nstep = 1:Nstep

    % Model update including hyperparameters if needed
    updatehyp = mod(nstep,Nuhyp)==0; if updatehyp, hyps=hyps0;
    disp(['Updated model with ' num2str(nstep) ' training points']);end
    [pFeLi.f,pFeLi.g,hyps] = learnfg(Xtr,Ytr,mf,mg,utr,hyps,updatehyp);
    
    % Set controller and dynamics
    ctrl = @(t,x) ctrlFeLi(t,x,pFeLi,reffun);
    dyn = @(t,x) dynAffine(t,x,ctrl,pdyn);
   
    % Add current point to training set
    Xtr = [Xtr X(:,end)];
    y  = dyn(T(end),X(:,end));
    Ytr = [Ytr y(2)+ sn*randn(1)];
    utr = [utr [1; ctrl(T(end),X(:,end))]];
    
     % Simulate for dt
    [t,x] = ode45(dyn,T(end)+[0 dt],X(:,end));
     T = [T t(2:end)'];  X = [X x(2:end,:)'];
end

%% Viualization
disp('Plotting results...'); 

Ndte = floor(nthroot(Nte,E)); % Nte = Ndte^E;
Xte = ndgridj(XteMin, XteMax,Ndte*ones(E,1)) ;
Xte1 = reshape(Xte(1,:),Ndte,Ndte); Xte2 = reshape(Xte(2,:),Ndte,Ndte);

figure; hold on,title('f vs fh'); view([10 10]); xlabel('x1');ylabel('x2');
fte = pdyn.f(Xte); mufte = pFeLi.f(Xte); ftr = pdyn.f(Xtr);
surf(Xte1,Xte2,reshape(mufte,Ndte,Ndte),'edgecolor','none','FaceColor','interp');
surf(Xte1,Xte2,reshape(fte,Ndte,Ndte),'edgecolor','none','FaceColor','interp');
plot3(Xtr(1,:),Xtr(2,:),ftr,'ro');

figure; hold on,title('g vs gh');  view([-80 10]);xlabel('x1');ylabel('x2');
gte = pdyn.g(Xte); mugte = pFeLi.g(Xte); gtr = pdyn.g(Xtr); 
surf(Xte1,Xte2,reshape(mugte,Ndte,Ndte),'edgecolor','none','FaceColor','interp');
surf(Xte1,Xte2,reshape(gte,Ndte,Ndte),'edgecolor','none','FaceColor','interp');
plot3(Xtr(1,:),Xtr(2,:),gtr,'ro');

figure; hold on, xlabel('x1');ylabel('x2'); axis equal; 
xd = reffun(T);
plot(X(1,:),X(2,:),'r');
plot(xd(1,:),xd(2,:),'k--');
plot(Xtr(1,:),Xtr(2,:),'ro');

figure; hold on; xlabel('t');ylabel('x1');
plot(T,xd(1,:),'--');
plot(T,X(1,:));

disp('Pau');