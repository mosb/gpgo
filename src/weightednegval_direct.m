function [f,g,H]=weightednegval_direct(XStar, Prob)

XStar = XStar;
rho = Prob.user.rho;
XData = Prob.user.XData;
YData = Prob.user.YData;
covvy = Prob.user.covvy;
mx = Prob.user.mx;

NStar=size(XStar,1);
f=zeros(NStar,1);
g=0;
H=0;
for sample=1:numel(covvy.hypersamples)   
    [fi,gi,Hi]=negvalue(XStar,XData,YData,covvy,sample,mx);
    f=f+rho(sample)*fi;
    g=g+rho(sample)*gi;
    H=H+rho(sample)*Hi;
end

%         wf=zeros(length(XPlot),1);
%         wm=zeros(length(XPlot),1);
%         wdC=zeros(length(XPlot),1);
%         for sample=1:numel(covvy.hypersamples)
%             [f,m,dC] =  negvalue(XPlot,XData,YData,covvy,sample,1);
%             wf=wf+rho(sample)*f;
%             wm=wm+rho(sample)*m;
%             wdC=wdC+rho(sample)*(dC+m.^2);
%         end
