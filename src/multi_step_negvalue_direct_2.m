function [negval,result]=multi_step_negvalue_direct_2(x1,Prob)

result = true;
%x1 = x1';
rho = Prob.user.rho;
XData = Prob.user.XData;
YData = Prob.user.YData;
covvy = Prob.user.covvy;
mx = Prob.user.mx;
num_samples = Prob.user.num_samples;
lowrbnd = Prob.user.lowrbnd;
upprbnd = Prob.user.upprbnd;

N=@(xs,ys,covariance)matrify(@(x,y)normpdf(x,y,(covariance)^0.5),xs,ys);
%two step
negval=0;
for hs_ind=1:length(rho) %hs for 'hypersample'
  [m,C]=gpmeancov(x1,XData,covvy,hs_ind);
  sigma=sqrt(diag(C));
  ys=linspacey(m-2*sigma,m+2*sigma,num_samples)';
  %width=sigma; % Any better ideas?
  width=0.5*sepn(ys);
  w=N(m,ys,width^2+sigma^2)/N(ys,ys,width^2);
  if (any(isnan(w)))
    continue;
  end
  XData1=[XData;x1];
  for y1_ind=1:num_samples
    try
      y=ys(y1_ind);
      YData1=[YData;y];
      covvy1=gpparams(XData1,YData1,covvy,'update');
      rhonew=weights(covvy1); % used purely for determining x
      [maxes, inds] = max(rhonew);
      mx=max(mx,y);

      Prob1.user.rho = rhonew;
      Prob1.user.XData = XData1;
      Prob1.user.YData = YData1;
      Prob1.user.covvy = covvy1;
      Prob1.user.mx = mx;

			Problem.f = @(XStar)negval_direct_naive(XStar, Prob1);
			bounds = [lowrbnd; upprbnd]';
			opts.maxevals = 500;
			opts.showits = 0;

			[xmin, finalatmin] = Direct(Problem, bounds, opts);
			x = finalatmin';
		
      negval = negval + ...
							 rho(hs_ind)*w(y1_ind)*...
							 negvalue(x,XData1,YData1,covvy1,hs_ind,mx);
    catch
		%	stupid = lasterror;
		%	disp(stupid.message);
		%	disp(stupid.identifier);
		%	for i=1:numel(stupid.stack)
		%		disp([stupid.stack(i).file '(' num2str(stupid.stack(i).line) ...
		%					'): ' stupid.stack(i).name]);
		%	end
    end
  end
end
