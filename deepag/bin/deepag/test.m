get1=@(x) x(x<10);
get_rat_error_add=@(x,y) get1(abs(x(:,1)./x(:,2))-abs(y(:,1)./y(:,2)));

estion1=estion(:,1,1,1)./estion(:,2,1,1);
truth1=truth(:,1,1,1)./truth(:,2,1,1);
estdata=abs(estion1)-abs(truth1);

% plot ratio errors
figure();
for j=1:noises
    for i=1:axdists
        subplot(N,M,j);
        cumhist(sort(get_rat_error_add(estion(:,:,i,j),truth)),20,2,[colors{i}]);
    end
    title(['noise ' num2str(noise(j)) 'px']);
    xlabel('error magnitude');
    ylabel('frequency, %');
    xlim([-1,1]);
end
legend(arrayfun(@(x) {num2str(x)},axdist))
suptitle('additive errors in ratio w.r.t. axes distance')