function negval=noisy_negvalue(thresh,rho,x1,XData,ZData, ...
																		covvy,mx,lowrbnd,upprbnd, ...
																		num_samples,guesses,varargin)
% Someday this should be upgraded to multi-step capability


if (~isfield(covvy, 'logNoiseSDPos'))
	logNoiseSDPos = cellfun(@(x) strcmp(x, 'logNoiseSD'), covvy.names);
	covvyout.logNoiseSDPos = logNoiseSDPos;
else
	logNoiseSDPos = covvy.logNoiseSDPos;
end

if nargin<9 || isempty(num_samples)
    num_samples=7;
end

optim_options = optimset('GradObj','on','Hessian','on','Display', ...
												 'off','MaxFunEvals',10);
N=@(xs,ys,covariance)matrify(@(x,y)normpdf(x,y,(covariance)^0.5),xs,ys);

%two step
negval=0;
for hs_ind=1:length(rho) %hs for 'hypersample'
    hyperparams=covvy.hypersample(hs_ind).hyperparameters;
    noisevar=exp(2*hyperparams(logNoiseSDPos));
    [m,C]=gpmeancov(x1,XData,covvy,hs_ind);
    sigma=sqrt(diag(C));
    % Add in noise
    sigma=sigma+noisevar;
    zs=linspacey(m-2*sigma,m+2*sigma,num_samples)';
    %width=sigma; % Any better ideas?
    width=0.5*sepn(zs);
    %w=N(m,zs,width^2+sigma^2)*N(zs,zs,width^2);
    w=N(m,zs,width^2+sigma^2)/N(zs,zs,width^2);
    XData1=[XData;x1];
    for z1_ind=1:num_samples
        
        
        negmean=@(XStar) weighted_neg_gpmean(rho,Augment(XStar),XData,covvy);
        nonlcon= @(XStar) certainabout(rho,Augment(XStar),XData,covvy,thresh);
      [mxLoc,mx] = fmincon(negmean,mxLocObserved,[],[],...
        [],[],lowrbnd,upprbnd,nonlcon,optim_mean_options);
      mx = -mx;
        
        
        z=zs(z1_ind);
        ZData1=[ZData;z];
        covvy1=gpparams(XData1,ZData1,covvy,'update');
        rhonew=weights(covvy1); % used purely for determining x
        mx=max(mx,y);
        NGuess = length(guesses);
        xtrial = nan(NGuess,size(x1,2));
        negvalguess=nan(NGuess,1);
        sanitycheck=nan(NGuess,1);
        for guess=1:NGuess
            [xtrial(guess,:),negvalguess(guess),sanitycheck(guess)]= ...
                            fmincon(@(x) weightednegval(rhonew,x,XData1,ZData1, ...
                                                                                    covvy1,mx,varargin{:}), ...
                                            guesses(guess,:),[],[],[],[],lowrbnd, ...
                                            upprbnd,[],optim_options);
        end
        [negval1,ind1]=min(negvalguess);
        x=xtrial(ind1,:);
        negval=negval+rho(hs_ind)*w(z1_ind)*negvalue(x,XData1,ZData1,covvy1,hs_ind,mx,varargin{:});
    end
end

