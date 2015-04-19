scrsz = get(0,'ScreenSize');
figure('Position',[1 1 0.4*scrsz(3) 0.6*scrsz(4)])
% left bottom width height
hold on
xx=linspace(lowrbnd,upprbnd,500)';

for i=1:length(xx);
    fs(i)=-fn(xx(i));
end

wv=zeros(length(xx),1);
awm=zeros(length(xx),1);
awC=zeros(length(xx),1);
for sample=1:numel(covvy.hypersamples)
  [v,m,dC] =  negvalue(Augment(xx),XData,YData,covvy,sample,mx,[],[],1);
  wv=wv+rho(sample)*v;
  awm=awm-rho(sample)*m;
  awC=awC+rho(sample)*(dC+m.^2);
end
awC=awC-awm.^2;

[mv,I]=min(v);
x=xx(I);

% [awm,awC]=weighted_gpmeancov(rho,Augment(xx),XData,covvy);
% awm=-awm;

subplot(3,1,1:2);
%subplot(2,14,1:28);

[sd,mn,obs]=plot_regression(XLocs(isactualobs,:),-YData(isactualobs),xx,awm,sqrt(awC),'MarkerSize',12);

% isactualobs=abs(XDerivDirns)*ones(NDims,1)==0;
% XLocsObs=XLocs(isactualobs,:);
% YDataObs=YData(isactualobs);
% plot(XLocsObs,YDataObs,'k+');

XLocsDerivObs=XLocs(~isactualobs,:);
YDataDerivObs=YData(~isactualobs);
FDataDerivObs=FData(~isactualobs);
XFDataDerivObs=XFLocs(~isactualobs);
XDerivObsDirns=XDerivDirns(~isactualobs);

for i=1:length(YDataDerivObs)
    o=plot(XFDataDerivObs(i),-FDataDerivObs(i),'ok','MarkerSize',12);
    xi=XLocsDerivObs(i);
    zi=weighted_gpmeancov(rho,Augment(xi),XData,covvy);
    zi=-zi;
    run=0.1;
    rise=YDataDerivObs(i)*norm(run);
    l=line([xi-run,xi+run],[zi-rise,zi+rise],'Color',[0 0 0],'LineWidth',1.5);
end

fp=plot(xx,fs,'k:','LineWidth',2);

set(gca,'FontSize',14)

%axis([lowrbnd upprbnd -2.5 -1])
axis([lowrbnd upprbnd -0.3 0.1]) 
xlabel x
ylabel('y','Rotation',0.0)
%xlabel('$$\phi$$','Interpreter','latex')
%ylabel('f','Rotation',0.0)


el=plot(xx,v,'b','LineWidth',2);
mel=plot(x,mv,'bd','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','b');

legend([fp,obs,mn,sd,el,mel],'Objective function','Observation','Mean','\pm 1SD','Expected loss','Chosen posn of next observation',...
'Location','NorthEastOutside')
%

% if length(YData)==4
% legend([fp,obs,o,l,mn,sd],'Actual function','Observation','Dropped obs','Deriv obs','Mean','\pm 1SD')
% end

subplot(3,1,3);
%subplot(2,14,24:28);

ins=cat(1,covvy.hypersamples(:).hyperparameters);
% InputScales=(reshape(ins(:,1),covvy.hyperparams(2).NSamples,covvy.hyperparams(1).NSamples));
% OutputScales=(reshape(ins(:,2),covvy.hyperparams(2).NSamples,covvy.hyperparams(1).NSamples));
% weightmat=reshape(rho,covvy.hyperparams(2).NSamples,covvy.hyperparams(1).NSamples);

set(gca,'FontSize',14)

%colours=colormap(copper);L=length(colours);
colormap hot;
%colours(max(ceil(rho*L),1),:)

scatter([exp(ins(:,1));-10;-20],[exp(ins(:,2));-10;-20],[200*ones(size(rho));1;1],[max(rho,0);1;0],'s','filled');
%colorbar('YTickLabel',linspace(0,1,6),'FontSize',14)
colorbar('FontSize',14)
%axis([0.4 1.95 0.058 0.158])
axis([0.4 2 0.058 0.161])
set(gca,'ytick',linspace(0.06,0.16,3))
xlabel 'x Scale'
ylabel 'y Scale'
title 'Sample weights'
axis square

  
[ax,h3]=suplabel(['Function Evaluation #',num2str(length(YData))],'t');
set(h3,'FontSize',16)


%saveas(gcf,['xexpected_loss_',num2str(length(YData)),'.pdf'])