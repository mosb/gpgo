clear

ActualNoise=5;
ActualLS=10;
ActualTS=1;

logLS=log(ActualLS);




Xind=0;
for logTS=linspace(-7,12,40)
    Xind=Xind+1;
    Yind=0;
    for logNoise=linspace(-2,3.5,40);
        Yind=Yind+1;
        Z(Xind,Yind)=blackbox(logTS,logLS,logNoise);
        X(Xind,Yind)=logTS;
        Y(Xind,Yind)=logNoise;
    end
end

contour(X,Y,Z,linspace(-28,-24,40))
xlabel('logTS')
ylabel('logNoise')


