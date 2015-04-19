clear
close all
% 
% covvy=struct('covfn',@(hp,varargin)ndimsqdexp_cov_fn_withderivs(hp,1,varargin{:}));
%         
% covvy.hyperparams(1)=struct('name','logInputScale','priorMean',log(0.25),'priorSD',0.2,'NSamples',1,'type','real');
% covvy.hyperparams(2)=struct('name','logOutputScale','priorMean',log(16),'priorSD',0.6,'NSamples',1,'type','real');%30
% covvy.hyperparams(3)=struct('name','mean','priorMean',-6,'priorSD',1,'NSamples',1,'type','real');
% covvy.hyperparams(4)=struct('name','logNoiseSD','priorMean',-1,'priorSD',1,'NSamples',1,'type','real');
% 
% x0=0;
% fn=@(x) -humps(x);
% lowrbnd=-5;
% upprbnd=5;
% options.MaxFunEvals=5;
% options.Plot=true;

% covvy=struct('covfn',@(hp,varargin)ndimsqdexp_cov_fn_withderivs(hp,1,varargin{:}));
%         
% covvy.hyperparams(1)=struct('name','logInputScale1','priorMean',log(1),'priorSD',0.3,'NSamples',5,'type','real');
% covvy.hyperparams(2)=struct('name','logOutputScale','priorMean',log(0.10),'priorSD',0.2,'NSamples',5,'type','real');%30
% covvy.hyperparams(3)=struct('name','mean','priorMean',0,'priorSD',1,'NSamples',1,'type','real');
% covvy.hyperparams(4)=struct('name','logNoiseSD','priorMean',log(0.00001),'priorSD',0.6,'NSamples',1,'type','real');
% 
% x0=-7;
% fn=@(x) -normpdf(x,-2,2)-normpdf(x,5,1.5);
% lowrbnd=-15;
% upprbnd=15;
% options.MaxFunEvals=10;
% options.Plot=true;

fn=@(x) -normpdf(x,-2,2)-normpdf(x,5,1.5); %objective fn

ndims=1; % Number of dimensions

% Covariance functions
covvy=struct('covfn',@(hp,varargin)ndimsqdexp_cov_fn_withderivs(hp,ndims,varargin{:}));
        
covvy.hyperparams(1)=struct('name','logInputScale1','priorMean',log(1),'priorSD',0.3,'NSamples',5,'type','real');
covvy.hyperparams(2)=struct('name','logOutputScale','priorMean',log(0.10),'priorSD',0.2,'NSamples',5,'type','real');%30
covvy.hyperparams(3)=struct('name','mean','priorMean',0,'priorSD',1,'NSamples',1,'type','real');
covvy.hyperparams(4)=struct('name','logNoiseSD','priorMean',log(0.00001),'priorSD',0.6,'NSamples',1,'type','real');

x0=-7*ones(1,ndims); %starting point
lowrbnd=-15*ones(1,ndims);
upprbnd=15*ones(1,ndims);
options.MaxFunEvals=7;
options.Plot=true;



% fn=@Tom;
% covvy=struct('covfn',@(hp,varargin)ndimsqdexp_cov_fn_withderivs(hp,3,varargin{:}));
% covvy.meanfn=@(hp) ndimssqdexp_mean_fn(hp);
% covvy.hyperparams(1)=struct('name','logInputScale1','priorMean',log(0.4),'priorSD',0.2,'NSamples',3,'type','real');%0.54
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logInputScale2','priorMean',log(0.4),'priorSD',0.2,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logInputScale3','priorMean',log(0.4),'priorSD',0.2,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logOutputScale','priorMean',log(1),'priorSD',0.1,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','mean','priorMean',1,'priorSD',1,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logNoiseSD','priorMean',log(0.0001),'priorSD',1,'NSamples',1,'type','real');
% 
% x0=[0.5 0.5 0.5];
% lowrbnd=[0 0 0];
% upprbnd=[1 1 1];
% options.MaxFunEvals=20;

% fn=@(x) -Tom([x 0 1]);
% covvy=struct('covfn',@(hp,varargin)ndimsqdexp_cov_fn_withderivs(hp,1,varargin{:}));
% covvy.meanfn=@(hp) ndimssqdexp_mean_fn(hp);
% covvy.hyperparams(1)=struct('name','logInputScale1','priorMean',log(0.4),'priorSD',0.2,'NSamples',1,'type','real');%0.54
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logOutputScale','priorMean',log(1),'priorSD',0.1,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','mean','priorMean',1,'priorSD',1,'NSamples',1,'type','real');
% covvy.hyperparams(length(covvy.hyperparams)+1)=struct('name','logNoiseSD','priorMean',log(0.0001),'priorSD',1,'NSamples',1,'type','real');
% 
% x0=[0.5];
% lowrbnd=[0];
% upprbnd=[1];
% options.MaxFunEvals=5;
% options.Plot=true;

minimise


%[x,f,outputs] = maximise(fn,x0,lowrbnd,upprbnd,covvy)

% NX=20;
% NY=20;
% 
% X1s=linspace(lowrbnd(1),upprbnd(2),NX)';
% X2s=linspace(lowrbnd(1),upprbnd(2),NY)';
% 
% [X,Y]=meshgrid(X1s,X2s);
% surf(X,Y,-(1-X).^2-100*(Y-X.^2).^2)
    
%maximise(@(x) banana(x(1),x(2)),[0,0],covvy)