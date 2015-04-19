function [f,g,H]=weighted_neg_gpmean(rho,XStar,XData,covvy)

NStar=size(XStar,1);
f=zeros(NStar,1);
g=0;
H=0;
for sample=1:numel(covvy.hypersamples)   
    [m,C,gm,gC,Hm]=gpmeancov(XStar,XData,covvy,sample,'no_cov');
    f=f-rho(sample)*m;
    g=g-rho(sample)*cell2mat2d(gm);
    H=H-rho(sample)*cell2mat2d(Hm);
end