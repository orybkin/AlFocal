% plots several degeneracy_hist's vs different level of noise

% settings
corr=7;
pop_size=1000;
method='|F|=0';
noise_out=0;
noises=[1 3 6 20];
% plot
figure()
for i=1:4
    subplot(2,2,i);
    imaginary_hist(corr,pop_size, method, noise_out,noises(i));
    axis([-1 1.5 0 100])
end