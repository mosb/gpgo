NX=40;
NY=40;

X1s=linspace(lower_bound(1),upper_bound(1),NX)';
X2s=linspace(lower_bound(2),upper_bound(2),NY)';

XPlot=augment(allcombs([X1s,X2s]));
wf=zeros(length(XPlot),1);
wm=zeros(length(XPlot),1);
for sample=1:numel(covvy.hypersamples)
	[f,m,dC] = negvalue(XPlot,x_data(1:evaluation,:),y_data(1:evaluation),covvy,sample,maximum,[],[],1);
	wf=wf+rho(sample)*f;
	wm=wm+rho(sample)*m;
end

[X,Y]=meshgrid(X1s,X2s);

drawnow

figure(1)
clf(1)
title('Mean')
set(gca,'FontSize',14)
hold on

Z=reshape(m,NX,NY);

surf(X,Y,-Z); colorbar

% figure(2)
% clf(2)
% title('Standard Deviation')
% set(gca,'FontSize',14)
% hold on  
% 
% Z=reshape(sqrt(dC),NX,NY);
% surf(X,Y,Z); colorbar
% 
figure(2)
clf(2)   
title('Expected Utility')
set(gca,'FontSize',14)
hold on  
% 
Z=reshape(f,NX,NY);
surf(X,Y,Z); colorbar
% 
% realfunc = zeros(NX*NY,1);
% for i=1:numel(realfunc)
% 	realfunc(i)=fn(XPlot(i,:));
% end
% figure(4)
% clf(4)   
% title('Actual Function')
% set(gca,'FontSize',14)
% hold on  
% 
% Z=reshape(realfunc,NX,NY);
% surf(X,Y,Z); colorbar
% figure(5);
% subplot(1,2,2);
% contour(X,Y,Z); colorbar

drawnow
