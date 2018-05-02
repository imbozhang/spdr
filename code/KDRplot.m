function KDRplot(f1,X,Y,B)

Z=X*B;
figure(f1);
set(gca, 'NextPlot', 'replacechildren');

cidx=find(Y(:,1)==1);
plot(Z(cidx,1),Z(cidx,2),'kd', 'MarkerFaceColor', 'r');
set(gca, 'NextPlot', 'add');
cidx=find(Y(:,1)==0);
plot(Z(cidx,1),Z(cidx,2),'ko', 'MarkerFaceColor', 'b');
drawnow;


