x=linspace(0, 5,1000);
smaller=x<1;
y=x;
y(smaller)=1-x(smaller);
y(~smaller)=1-1./x(~smaller);
plot(x,y,'Color',[0 0 1]);
axis([0 5 0 2])
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 500 200]);
saveas(gcf,'results/error_eval.png');