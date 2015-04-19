function [c,ceq,GC,GCeq] = certainabout(rho,x,XData,covvy,thresh)
% This is the constraint that the covariance of our belief about x
% is below the prespecified threshold thresh

f=0;
g=0;
for sample=1:numel(covvy.hypersamples)   
    [m,C,gm,gC]=gpmeancov(x,XData,covvy,sample,'no_mean');
    f=f+rho(sample)*C;
    g=g+rho(sample)*cell2mat2d(gC);
end


n=size(XData,2);

% Nonlinear inequalities at x
c = f-thresh;
    
% Nonlinear equalities at x
ceq = 0; 

if nargout > 2   % nonlcon called with 4 outputs
    % Gradients of the inequalities
    GC = g;    
    % Gradients of the equalities
    GCeq = zeros(n,1);
end
