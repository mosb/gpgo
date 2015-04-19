function [f,out2,out3,out4,out5,out6,out7,out8,out9] = negvalue(XStar,XData,YData,covvy,sample,mx,personal_space,penalty,outputGP)
% Compute the expected value of taking an observation of a function at
% XStar, given that we have already observed YData at XData. covvy
% describes our covariance structure for the GP we possess over the
% function.

% the mean m and covariance C (and their gradients & Hessians) of our GP
% evaluated at XStar

NData=size(XData,1);
NStar=size(XStar,1);

if nargin<5
    sample=1;
end
if nargin<6
    mx=max(YData);
elseif isempty(mx)
    mx=max(YData);
end
if nargin<9
    outputGP=0;
    % The function will not output the computed mean and covariances of the
    % GP
end

if (~outputGP && nargout == 1) || (outputGP && nargout <= 3)
    [m,C]=gpmeancov(XStar,XData,covvy,sample);
    D=diag(C);
elseif (~outputGP && nargout == 2) || (outputGP && nargout <= 6)
    [m,C,gm,gC]=gpmeancov(XStar,XData,covvy,sample);
    
    gm=cell2mat2d(gm);
    D=diag(C);
    gD=cellfun(@diag,gC,'UniformOutput',0);
elseif (~outputGP && nargout == 3) || (outputGP && nargout <= 9)
    [m,C,gm,gC,Hm,HC]=gpmeancov(XStar,XData,covvy,sample); 
    D=diag(C);
    gD=cellfun(@diag,gC,'UniformOutput',0);
    HD=cellfun(@diag,HC,'UniformOutput',0);
end

%finds the diagonal for each derivative


D=max(D,eps^2); % effectively jitter



%have to deal specially with small covariances, otherwise causes a singularity
%smalls=(D<eps^2);

% if any(any() 
%     f=-max([mx*ones;m]);
%     if nargout>1
%         g=-zeros(size(XData,2),1);
%         if nargout>2
%             H=-0.5*cell2mat2d(Hm);
%         end
%     end
%     return
% end

% argterm=(mx-m)/(sqrt(2*D)); % scalar
% 
% f = 0.5*(mx+m)+sqrt(D/(2*pi))*exp(-argterm^2)+0.5*erf(argterm)*(mx-m);
% 
% % argterm=(mx-m)./(sqrt(2*diag(D)));
% % f = 0.5*(mx+m)+sqrt(diag(D)/(2*pi)).*exp(-argterm.^2)+0.5*erf(argterm).*(mx-m);
% f=-f;
% % Dompute the scalar value at XStar
% if nargout > 1   % fun called with two output arguments
%    g = 0.5*(2*pi*D)^(-0.5)*gD*exp(-argterm^2)+0.5*erfc(argterm)*gm; 
%    g=-g;
%    % Gradient vector of the value evaluated at XStar
%    if nargout > 2 % fun called with three output arguments
%       H = 0.5*Hm*erfc(argterm)+(2*pi*D^3)^(-0.5)*exp(-argterm^2)*(...
%           -0.25*gD*gD'+0.5*D*HD+(gm+0.5/D*gD*(mx-m))*(gm*D+0.5*gD*(mx-m))');
%       H=-H;
%       % Hessian matrix evaluated at XStar
%    end
% end



argterm=(mx-m)./(sqrt(2*D)); % scalar
f = 0.5*(mx+m)+sqrt(D/(2*pi)).*exp(-argterm.^2)+0.5*erf(argterm).*(mx-m);

% argterm=(mx-m)./(sqrt(2*diag(D)));
% f = 0.5*(mx+m)+sqrt(diag(D)/(2*pi)).*exp(-argterm.^2)+0.5*erf(argterm).*(mx-m);
f=-f;
% Compute the scalar value at XStar

if outputGP
    out2=m;
    out3=D;
end

if (~outputGP && nargout > 1) || (outputGP && nargout > 3)  % fun called with two output arguments
    
    NDims=length(gm);
    
    g = 0.5*(2*pi*repmat(D,NDims,1)).^(-0.5).*cell2mat2d(gD).*repmat(exp(-argterm.^2),NDims,1)+0.5*repmat(erfc(argterm),NDims,1).*cell2mat2d(gm); 
    g=-g;
    % Gradient vector of the value evaluated at XStar

    out2=g;
    if outputGP
        out3=m;
        out4=D;
        out5=gm;
        out6=gD;
    end
   
    if (~outputGP && nargout > 2) || (outputGP && nargout > 6) % fun called with three output arguments
                
        term1L=repmat(cell2mat2d(gD),1,NDims);
        term1R=repmat(cell2mat2d(gD'),NDims,1);
        term2L=repmat((cell2mat2d(gm)+...
                    0.5*cell2mat2d(gD).*repmat(D.^(-1).*(mx-m),NDims,1)),1,NDims);
        term2R=repmat((cell2mat2d(gm')+...
                    0.5*cell2mat2d(gD').*repmat(D.^(-1).*(mx-m),1,NDims)),NDims,1);
        
        Hm=cell2mat2d(Hm);
        HD=cell2mat2d(HD);

        H = 0.5*Hm.*repmat(erfc(argterm),NDims,NDims)+...
            repmat((2*pi*D.^3).^(-0.5).*exp(-argterm.^2),NDims,NDims).*(...
            -0.25*term1L.*term1R+repmat(D,NDims,NDims).*(0.5*HD+...
            term2L.*term2R));
        H=-H;
        % Hessian matrix evaluated at XStar

        out3=H;
        if outputGP
            out4=m;
            out5=D;
            out6=gm;
            out7=gD;
            out8=Hm;
            out9=HD;
        end
   end
end



