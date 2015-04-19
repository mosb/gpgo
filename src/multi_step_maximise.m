%function [x,func,outputs] = maximise(fn,x0,lowerbnd,upperbnd,covvy,options)
%minimise function fn by sequentially greedily selecting the most valuable
%observation according to a GP fit to fn. 

NDims=length(x0);
defaultopt = struct(...
    'MaxFunEvals',100*NDims,...
    'ValueTol',1e-3...
    ); 

names=fieldnames(defaultopt);
for i=1:length(names);
    name=names{i};

    try
        optvalue = options.(name);
    catch
        optvalue = [];
        lasterr('');  % clean up last error
    end

    if isempty(optvalue)
        optvalue = defaultopt.(name);
    end
    
    options.(name)=optvalue;
end


optim_mean_options = optimset('GradObj','on','LargeScale','off','Hessian','off','Display','off');
optim_negval_options = optimset('Display','on');
%'GradConstr','on','MaxFunEvals',200,'TolFun',options.ValueTol,'TolCon',0.01^NDims);

negvalold=0;
x=x0;
XData=zeros(0,NDims);
YData=zeros(0,1);

NData=size(XData,1);


covvy=hyperparams(covvy);

personal_space=zeros(1,NDims);
for i=1:NDims
    personal_space(i)=max(covvy.hyperparams(cellfun(@(x) strcmp(x,['logInputScale',num2str(i)]),covvy.names)).samples);
end
personal_space=exp(personal_space);
length_scales=personal_space;
distances=[];

% We must know our final output to this level of accuracy (in the
% covariance)
thresh=max(covvy.hyperparams(cellfun(@(x) strcmp(x,'mean'),covvy.names)).samples)^2;

covvy=bmcparams(covvy);
% 
% figure(1);
% figure(2);
% figure(3);

for fn_eval=1:options.MaxFunEvals
    fn_eval
    x
    
    fx=fn(x);
    
    % If we do have conditioning probs, just get rid of one of the
    % redundant pts.
%     [mn,indredundant]=min(abs(x-XData));
%     % really this measure of distance should be informed by our TS - do
%     % this later
% %         [mn,indredundant]=min([fx,YData(indxclose)]);
% %         if indredundant==1
% 
% 
%     if mn<options.Tol
%         covvy=gpparams(XData,YData,covvy,'downdate',indredundant);
%         XData(indredundant,:)=[];
%         YData(indredundant,:)=[];
%     end


    % find the matrix of distances of every observation from every other
    % observation
    newdists=[((repmat(x,NData,1)-XData).^2*length_scales'.^2);nan];
    distances=[[distances;newdists(1:end-1,:)'],newdists];
    
    XData=[XData;x];
    YData=[YData;fx];
    NData=size(XData,1);
    
    covvy=gpparams(XData,YData,covvy,'update',length(YData));
   
    % move samples around?
    % covvy=hyperparams(covvy,'hyperactive');
    
    % efficiently update all our covariance terms and parameters given this
    %new observation

    rho=weights(covvy);
    
    [mxObserved,mxInd]=max(YData);
    mxLoc=XData(mxInd,:);
    negmean=@(XStar) weighted_neg_gpmean(rho,XStar,XData,covvy);
    nonlcon= @(x) certainabout(rho,x,XData,covvy,thresh);
    [mxLoc,mx]=fmincon(negmean,mxLoc,[],[],...
                        [],[],lowrbnd,upprbnd,nonlcon,optim_mean_options);
                    
%                     A*x <= b
%                     fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
     mx
     
    penalty=0.1*(mx-median(YData)); 

    negval = @(XStar) multi_step_negvalue(2,rho,XStar,XData,YData,covvy,mx,lowrbnd,upprbnd,7);
    %weightednegval(rho,XStar,XData,YData,covvy,sample,mx,personal_space,penalty);
    
    nonlcon= @(x) notnear(x,XData,personal_space);
     
    
    % because I do not trusy matlab's fminunc, I start it off at a number
    % of trial positions. 
    
    % find the loneliest observations - our variance is likely to be
    % largest in their vicinity
    mins=min(distances)';
    maxMins=max(mins);
    minMins=min(mins);
    XGuessInds=mins>=0.8*(maxMins-minMins)+minMins;
    
    % find the largest observations - our mean is likely to be largest in
    % their vicinity
    maxY=max(YData);
    minY=min(YData);
    XGuessInds2=YData>=0.8*(maxY-minY)+minY;
    
    % We guess that our utility will be largest where either the mean or
    % variance is largest
    XGuess=XData(or(XGuessInds,XGuessInds2),:);
    NGuess=size(XGuess,1);

    xtrial=nan(NGuess,NDims);
    negvalguess=nan(NGuess,1);
    sanitycheck=nan(NGuess,1);
    for guess=1:NGuess
        [xtrial(guess,:),negvalguess(guess),sanitycheck(guess)]=...
                fmincon(negval,XGuess(guess,:),[],[],[],[],...
                lowrbnd,upprbnd,[],optim_negval_options);
    end
    [negvaly,I]=min(negvalguess(sanitycheck>0));
    xtrial=xtrial(sanitycheck>0,:);
    x=xtrial(I,:);
        
        
%     % maybe put lengthscales in eyemat
%     eyemat=eye(NDims);
%     xtrial=nan(2*NDims*size(XGuess,1),NDims);
%     negvaltrial=nan(2*NDims*size(XGuess,1),1);
%     sanitycheck=nan(2*NDims*size(XGuess,1),1);
%     I=[];
%     while isempty(I)
%         for iter=1:2*NDims*size(XGuess,1)
%             guess=ceil(iter/(2*NDims));
%             trial=iter-(guess-1)*(2*NDims);
%             [xtrial(iter,:),negvaltrial(iter),sanitycheck(iter)]=...
%                 fmincon(negval,XGuess(guess,:)+(-1)^trial*(2*rand)*personal_space.*eyemat(ceil(trial/2),:),...
%                 [eye(NDims);-eye(NDims)],[upprbnd;-lowrbnd],...
%                 [],[],[],[],nonlcon,optim_negval_options);
%              %minimize(XGuess+(-1)^trial*options.Tol*eyemat(ceil(trial/2),:),negval,-10);
%         end
%         [negvaly,I]=min(negvaltrial(sanitycheck>0));
%         xtrial=xtrial(sanitycheck>0,:);
%         XGuess=XData(ceil((XData,1)*rand),:);
%     end
    
%     check=nonlcon(x)>0;
%     sanity=0;
%     while any(check)
%         sanity=sanity+1;
%         vec=x-XData(find(check,1));
%         mag=norm(vec);
%         if mag<eps^NDims || sanity>10
%             vec=rand(1,NDims)-0.5;
%             mag=norm(vec);
%         end
%         x=x+options.Tol*vec/mag;
%         check=nonlcon(x)>0;
%     end


%     if abs((negval-negvalold)/negvalold)<options.Tol
%         break
%     end

   
%     mess=sortrows([wf,XPlot],1);
%     maybeX=mess(:,2:end);
%     i=1;
%     x=maybeX(i,:);
%     while any(nonlcon(x)>0);
%         i=i+1;
%         x=maybeX(i,:);
%     end
% 
%     [c,I]=min(f);
%     x=XPlot(I,:);
%     

    %length_scales=exp(covvy.hypersamples(closestInd).hyperparameters(1:NDims));
    %convention is that input dimensions are always first
    
end

[rhomax,closestInd]=max(rho);
covvy.hypersamples(closestInd).hyperparameters

[maxY,maxInd]=max(YData);
x=XData(maxInd,:);
func=YData(maxInd,:);

outputs.xs=XData;
outputs.ys=YData;
outputs.covvy=covvy;