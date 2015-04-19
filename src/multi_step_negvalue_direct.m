function negval=multi_step_negvalue_direct(x1,XData,YData, ...
																		covvy,mx,lowrbnd,upprbnd, ...
																		num_samples,guesses,varargin)
% n = # of steps lookahead
XStar = XStar';
rho = Prob.user.rho;
XData = Prob.user.XData;
YData = Prob.user.YData;
covvy = Prob.user.covvy;
mx = Prob.user.mx;
num_samples=7;

N=@(xs,ys,covariance)matrify(@(x,y)normpdf(x,y,(covariance)^0.5),xs,ys);
switch lookahead
  case 2
    %two step
    negval=0;
    for hs_ind=1:length(rho) %hs for 'hypersample'
      try
        [m,C]=gpmeancov(x1,XData,covvy,hs_ind);
        sigma=sqrt(diag(C));
        ys=linspacey(m-2*sigma,m+2*sigma,num_samples)';
        %width=sigma; % Any better ideas?
        width=0.5*sepn(ys);
        %w=N(m,ys,width^2+sigma^2)*N(ys,ys,width^2);
        w=N(m,ys,width^2+sigma^2)/N(ys,ys,width^2);
        XData1=[XData;x1];
        for y1_ind=1:num_samples
          y=ys(y1_ind);
          YData1=[YData;y];
          covvy1=gpparams(XData1,YData1,covvy,'update');
          rhonew=weights(covvy1); % used purely for determining x
          mx=max(mx,y);
          NGuess = length(guesses);
          xtrial = nan(NGuess,size(x1,2));
          negvalguess=nan(NGuess,1);
          sanitycheck=nan(NGuess,1);
					for guess=1:NGuess
            [xtrial(guess,:),negvalguess(guess),sanitycheck(guess)]= ...
								fmincon(@(x) weightednegval(rhonew,x,XData1,YData1, ...
																						covvy1,mx,varargin{:}), ...
												guesses(guess,:),[],[],[],[],lowrbnd, ...
												upprbnd,[],optim_options);
          end
          [negval1,ind1]=min(negvalguess);
          x=xtrial(ind1,:);
          negval=negval+rho(hs_ind)*w(y1_ind)*negvalue(x,XData1,YData1,covvy1,hs_ind,mx,varargin{:});
        end
			catch
			end
    end
  case 3
    %three step 
    negval=0;

    for hs_ind=1:length(rho) %hs for 'hypersample'
      [m,C]=gpmeancov(x1,XData,covvy,hs_ind);
      sigma=sqrt(C);
      ys=linspacey(m-2*sigma,m+2*sigma,num_samples)';
      width=sigma; % Any better ideas?
      w1=N(m,ys,width^2+sigma^2)*N(ys,ys,width^2);
      XData1=[XData;x1];
      for y1_ind=1:num_samples
        y=ys(y1_ind);

        YData1=[YData;y];
        covvy1=gpparams(XData1,YData1,covvy,'update');
        mx=max(mx,y);
        rhonew=weights(covvy1); % used purely for determining x

        x2=fmincon(@(x) multi_step_negvalue(2,rhonew,x,XData1,YData1,covvy1,mx,lowrbnd,upprbnd,num_samples,varargin{:}),...
                x1+[rand,rand],[],[],[],[],lowrbnd,upprbnd,[]);

        [m,C]=gpmeancov(x2,XData1,covvy1,hs_ind);
        sigma=sqrt(C);
        ys=linspacey(m-2*sigma,m+2*sigma,num_samples)';
        width=sigma; % Any better ideas?
        w2=N(m,ys,width^2+sigma^2)*N(ys,ys,width^2);
        XData2=[XData1;x2];
        for y2_ind=1:num_samples
          y=ys(y2_ind);

          YData2=[YData1;y];
          covvy2=gpparams(XData2,YData2,covvy1,'update');
          mx=max(mx,y);
          rhonew=weights(covvy2); % used purely for determining x
          x3=fmincon(@(x) weightednegval(rhonew,x,XData2,YData2,covvy2,mx,varargin{:}),...
              x2,[],[],[],[],lowrbnd,upprbnd,[],optim_options);

          negval=negval+rho(hs_ind)*w1(y1_ind)*w2(y2_ind)*negvalue(x3,XData2,YData2,covvy2,hs_ind,mx,varargin{:});
        end
      end
    end
end
