% Copyright (c) by Jonas Umlauft under BSD License
% Last modified: Jonas Umlauft 10/2018

clear;clc;  rng default; close all; 
%% Set Parameters
disp('Setting parameters...')
% Basic Parameters
Tsim = 30;         % Simulation time
sn = 1e-8;         % Observation noise (std deviation)
E = 2;             % State space dimension
x0 = [3; 2];       % initial State
r_min = 1e-5;      % precision bound

% Reference trajectory
reffun = @(t) refGeneral(t,E+1,@(x) sin(x));    % Cricle

% Controller gains and beta value
pFeLi.lam = ones(E-1,1);
pFeLi.kc = 1;
beta = 7;

% Define Systemdynamics
pdyn.f = @(x) 1-sin(x(1,:)) + 0.5*sigmf(x(2,:),[0.1 0]);
pdyn.g = @(x) 1+0.5*sin(0.5*x(2,:));
pFeLi.g = pdyn.g;

% GP learning parameters
optGPR = {'KernelFunction','ardsquaredexponential','fitMethod','none',...
    'OptimizeHyperparameters','none','KernelParameters',[5; 5; 5],...
    'ConstantSigma',true,'Sigma',sn};
odeopt = odeset('RelTol',1e-6,'AbsTol',1e-9);


% Visualization
Nte = 4e4; XteMin = -4; XteMax = 4;


%% Learn fGP
disp('Learn fGP...')
Ntrf = 1e3;
Xtef = ndgridj(XteMin, XteMax,floor(nthroot(Ntrf,E))*ones(E,1)) ;
gprMdl = fitrgp(Xtef',pdyn.f(Xtef)); pdyn.f = @(x) predict(gprMdl,x')';

%% Simulating Adaptive Controller
disp('Simulating Controller...')
% Initialization
X = x0; T = 0; Xtr = []; Ytr = [];  Tevent = [];Phi = 0;

while T(end) < Tsim
    % Add current point to training set
    Xtr = [Xtr X(:,end)];
    Ytr = [Ytr pdyn.f(X(:,end)) + sn*randn(1)];
    
    % Learn GP Model and update control law
    try
        gprModel = fitrgp(Xtr',Ytr',optGPR{:});
        pFeLi.f = @(x) predict(gprModel,x');
        pFeLi.varfun = @(x) (nth_output(2, @predict, gprModel,x')).^2;
    catch
        warning('Problem with fitrgp');
    end
    
    
    % Define control law and event trigger
    ctrl = @(t,x) ctrlFeLi(t,x,pFeLi,reffun);
    dyn = @(t,x) dynAffine(t,x,ctrl,pdyn);
    eventfun = @(t,x) eventPhi(t,x,reffun,beta,pFeLi,r_min);
    odeopt = odeset(odeopt,'Events',eventfun);
    
    % Check suspicious behavior
    if  eventfun(T(end),X(:,end))>0
        warning('Event was triggered at training point!');
    end
    
    
    % Simulate until event triggered or time is up
    [t,x,te,~,ie] = ode45(dyn,[T(end) Tsim],X(:,end),odeopt);
    
    if ~isempty(ie)
                disp(['Event triggered, added at t = ',num2str(t(end))]);
                Tevent = [Tevent t(end)];
    end    
    phi = getPhi(t,x',pFeLi,reffun,r_min);
    T = [T t(2:end)'];  X = [X x(2:end,:)'];  Phi = [Phi phi(2:end)'];
end

%% Viualization
disp('Plotting results...'); 

Ndte = floor(nthroot(Nte,E)); % Nte = Ndte^E;
Xte = ndgridj(XteMin, XteMax,Ndte*ones(E,1)) ;
Xte1 = reshape(Xte(1,:),Ndte,Ndte); Xte2 = reshape(Xte(2,:),Ndte,Ndte);

% Plot GP variance and trajectory
figure; hold on; axis tight; xlabel('x1'); ylabel('x2');axis equal;
title(' GP variance + trajectory'); xlim([-1 4]); ylim([-2 3]);
xd = reffun(T);
varte = pFeLi.varfun(Xte);sortvarte = sort(varte);
cutoff = sortvarte(length(varte) - floor(length(varte)*0.7));
varte(varte>cutoff) = cutoff;
surf(Xte1,Xte2,reshape(varte,Ndte,Ndte),'EdgeColor','none','FaceColor','interp'); colormap(flipud(parula));
plot(X(1,:),X(2,:),'r');
plot(xd(1,:),xd(2,:),'k--');
plot(Xtr(1,:),Xtr(2,:),'ro');

% Plot tracking error
figure;  
norme = sqrt(sum((X-xd(1:E,:)).^2,1));
semilogy(T,norme);
xlabel('t'); ylabel('|e|'); title('tracking error');
hold on; 
semilogy(Tevent,norme(any(Tevent'==T)),'ro');

% Plot Phi
figure; hold on,title('Events');xlabel('t'); ylabel('phi');
Philim = ones(size(T))*pFeLi.kc/beta;
Phievent = Phi(any(Tevent'==T));
plot(T,Phi,'k--');
plot(T,Philim,'b');
plot(Tevent,Phi(any(Tevent'==T)),'ro');

disp('Pau');