
num_cands=6; % the more the better, but increasing this above 6 or so is going to be pretty slow.

try
    candidate_combs_template=load(candidate_combs_cell);
catch
    candidate_combs_template{num_cands}=find_likelihood_samples(zeros(num_cands,1),ones(num_cands,1),100,300,false);
end

covvy.plots=false;

covvy.candidate_combs_template=candidate_combs_template;
% Initialises hypersamples
covvy=hyperparams(covvy);
% All hypersamples need to have their gpparams overwritten
covvy.lastHyperSampleMoved=1:numel(covvy.hypersamples);
covvy.explored=zeros(0,numel(covvy.hyperparams));%cat(1,covvy.hypersamples.hyperparameters);

% Initialises hyper2samples
covvy=hyper2params(covvy);

NData=size(XsFull,1);
NSamples=numel(covvy.hypersamples);

XData=[];
YData=[];

covvy=gpparams(XData, YData, covvy, 'overwrite', []);

count=0;

for ind=training_set_size+1:NData
    

    
    if params.print==1 && ind>0.01*(1+count)*NData
        count=count+1;
        ind
    elseif params.print==0 && ( rem(ind,100) == 0)
        fprintf('.');
    end
        
	X = % point added at last iteration;
	Y = % function value at X;
    
    mx = max(YData);
    % qs is a (NSamples x 1) vector of strictly positive numbers. 
    qs = nan(NSamples,1);
    % We choose to minimise the uncertainty of our integral around the
    % negvalue of the point added at the last iteration.
    for sample=1:NSamples   
        qs(sample) = -negvalue(X,XData,YData,covvy,sample,mx) - mx;
        % I want qs to be expected improvement so that it is strictly
        % positive, ROMAN can you check that I have done it right.
    end
    
    XData=[XData;X];
    YData=[YData;Y]; 

    covvy=gpparams(XData, YData, covvy, 'update', size(XData,1), ...
        setdiff(1:NSamples,covvy.lastHyperSampleMoved));
    covvy=gpparams(XData, YData, covvy, 'overwrite', [], ...
        covvy.lastHyperSampleMoved);
      

    % Use ML_ind from previous time-step -improve_bmc_conditioning could
    % potentially moved to after manage_hyper_samples

    [managed,hypermanaged,covvy] = improve_bmc_conditioning(XData,YData,covvy);   

    covvy = calculate_hyper2sample_likelihoods(covvy,qs);
    
    covvy = determine_candidates(covvy, num_cands);
    covvy = bmcparams_ahs(covvy);   
    covvy = sample_candidate_likelihoods(covvy);
    [rho,invRL_wderivs_wcandidates]=weights_ahs(covvy); 
    [hyperrho,covvy]=hyperweights(rho,qs,covvy);
    
    if ~hypermanaged 
        % if we manage hyper2samples before hypersamples, it ensures that
        % added hypersamples are not going to be poorly conditioned
        covvy2 = manage_hyper2_samples(covvy, 'all','move');
    	%hyperman1=[cat(1,covvy.hyper2samples.logL),cat(1,covvy.hyper2samples.hyper2parameters)]
    end
    if ~managed % don't do this if improve_bmc_conditioning has already moved samples

        % If we have been asked to make predictions about more than one point, we
        % arbitrarily choose to minimise our uncertainty around the first.
%         pt=wm(1);
%         %pt=Lvec'*q_mean_star(:,1)/sum(Lvec);
%         qs=normpdf(repmat(pt,NSamples,1),q_mean_star(:,1),q_SD_star(:,1));
        
        [no_change,no_change_uncertainty,covvy] = manage_hyper_samples(covvy,covvy2,qs,invRL_wderivs_wcandidates);
        monitor.t(ind).no_change_uncertainty=no_change_uncertainty;
        if no_change
            covvy = zoom_hyper_samples(covvy);
    		%zoom1=[cat(1,covvy.hyper2samples.logL),cat(1,covvy.hyper2samples.hyper2parameters)]

% If we allow the below, samples can be zoomed into the region around
% another sample that will be forbidden by the new hyper2samples
%             if ~hypermanaged
%                     covvy = manage_hyper2_samples(covvy, 'all','move');
%                 %hyperman2=[cat(1,covvy.hyper2samples.logL),cat(1,covvy.hyper2samples.hyper2parameters)]
% 
%             end
        end
    elseif ~hypermanaged % if hypermanaged, covvy2 dne
        covvy=covvy2;
    end
    if ~hypermanaged 
        covvy = manage_hyper2_samples(covvy, 'all','clear');
    end

    

end
