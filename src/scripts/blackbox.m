function output = blackbox(logTS,logLS,logNoise)
% output is the log-likelihood of a matern GP with timescale exp(logTS),
% noise SD exp(logNoise) and lengthscale exp(logLS) fitting a stored dataset.

% the dataset was generated according to the following:

ActualNoise=5;
ActualLS=10;
ActualTS=1;
% 
% NData=7;
% XData=(1:NData)'+2*(rand(NData,1)-0.5);
% K=@(a,b) matrify(@(s,t) covfn('matern',s,t,ActualTS,ActualLS),a,b);
% ZData=mvnrnd(0*XData,K(XData,XData)+ActualNoise^2*eye(NData))';
% plot(XData,ZData,'+')

XData = ...
[    0.1943
    2.6469
    3.3897
    3.6342
    5.9004
    5.0689
    6.8775];

ZData = ...
  [   14.4254
    3.3273
   10.5610
    4.6281
   -0.1083
   -5.1088
  -16.4388];

NData=length(XData);


Kfn= @(at,bt) covfn('matern',at,bt,exp(logTS),exp(logLS));
K=@(as,bs) matrify(Kfn,as,bs);

Mu=0;

Kchol=chol(K(XData,XData)+exp(logNoise)^2*eye(NData));
output=logmvnpdf(ZData,Mu*ones(NData,1),Kchol,'cholesky');       