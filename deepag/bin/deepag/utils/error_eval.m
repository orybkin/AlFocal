x=linspace(0, 10,1000);
smaller=x<1;
y=x;
y(smaller)=1-x(smaller);
y(~smaller)=1-1./x(~smaller);
plot(x,y,'Color',[0 0 1]);
axis([0 10 0 1.25])
xlabel('f^{ob}')
ylabel('e(1,f^{ob})')
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 500 100]);
saveas(gcf,'../../../results/error_eval.eps','epsc');