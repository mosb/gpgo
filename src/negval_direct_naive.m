function [f,g,H]=negval_direct_naive(XStar, Prob)

XStar = XStar';
rho = Prob.user.rho;
XData = Prob.user.XData;
YData = Prob.user.YData;
covvy = Prob.user.covvy;
mx = Prob.user.mx;

NStar=size(XStar,1);
f=zeros(NStar,1);
g=0;
H=0;
[maxes, inds] = max(rho);
for sample=inds
    [fi,gi,Hi]=negvalue(XStar,XData,YData,covvy,sample,mx);
    f=f+rho(sample)*fi;
    g=g+rho(sample)*gi;
    H=H+rho(sample)*Hi;
end
