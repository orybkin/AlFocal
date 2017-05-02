% A script for generating scatter plot series of Bougnoux formula results
%
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

%

deepagpaths;

%% calculating
% what will be on X axis of the grpahs - either correspondences or levels of noise
key='corrs';
%key='noise';
% set parameters
if key=='corrs'
    % numbers of correspondences - x axis
    corresp=[7:17 18:3:29 30:5:49 50:10:100];
    corresp=[7:20];
    %corresp=[8:12];
    %corresp=[7:12];
    %corresp=[7 8];
    noise=1*ones(size(corresp));
    file=cell(size(corresp));
    for i=1:size(noise,2)
        file{i}=['../../data/correspondences_F_synth_1K_' num2str(noise(i)) 'noise.mat'];
    end
else
    noise=[0 0.5 1 2 3 6];
    file=cell(size(noise));
    for i=1:size(noise,2)
        file{i}=['../../data/correspondences_F_synth_1K_' num2str(noise(i)) 'noise.mat'];
    end
    corresp=8*ones(size(noise));
end


plot_data=true; % turn on if you want to see data for each focal length computation scatter plotted.
methods={'Free' '|F|=0' 'Prop' 'Prop6'};
methods={'Free' '|F|=0'};

%noise added to the first point in the first picture (in sigmas)
noise_out=0*ones(size(corresp));
pop_size=300; % number of points that will be on scatter plot
n=size(corresp,2); % x axis length
m=size(methods,2);
clear data stat_data
stat_data=NaN(2*m,n,4);
data=cell(2*m,n,2);

% size(stat_data):=(methods, n, parameters), where:
% methods are as in the legend below
% parameters are [absolute_error a_std relative_error r_std]
% size(data):=(methods, n, parameters), where:
% methods are as in the legend below
% parameters are [absolute relative]
perm32=@(x) permute(x,[1 3 2]);
resh22=@(x) perm32(reshape(x,2,2));
for i=1:n
    repS = adprintf({}, ['calculating correspondences ' num2str(i), '/', num2str(n)]);
    % beware, data are magically reshaped here
    for m=1:size(methods,2)
        [stat_temp, data_temp]=bougnoux_scatter(file{i},corresp(i),pop_size,methods{m},noise_out(i),plot_data);
        data(m*2-[1 0],i,:)=resh22(data_temp);
        stat_data(m*2-[1 0],i,:)=perm32(stat_temp);
    end
    rmprintf(repS);
end
clear data_temps stat_temp;
methods={'8pt' '8pt, real f', '8pt, imaginary f' ...
    ,'7pt', '7pt, real f', '7pt, imaginary f' ...
    ,'f-Ratio', 'f-Ratio, real f', 'f-Ratio, imaginary f' ...
    ,'6-Ratio', '6-Ratio, real f', '6-Ratio, imaginary f'...
    };
methods={'8pt, real f', '8pt, imaginary f', '7pt, real f', '7pt, imaginary f'};
methods={'8pt', '7pt'};
!notify-send hello
%% plotting
clear titles legends ylabels names

titles(1)={'Error of estimating (f1,f2)'};
legends(1,:)=methods;
ylabels(1)={'mean of f1/f1_{true}'};
names(1)={'bougnoux/error'};

titles(2)={'Std of estimating (f1,f2)'};
legends(2,:)=legends(1,:);
ylabels(2)={'std of f1/f1_{true}'};
names(2)={'bougnoux/std'};

titles(3)={'Mean error of estimating f1/f2'};
legends(3,:)=legends(1,:);
ylabels(3)={'error of f1/f2-f1_{true}/f2_{true}'};
names(3)={[names{1} '_proportion']};

titles(4)={'Std of estimating f1/f2'};
legends(4,:)=legends(1,:);
ylabels(4)={'std of f1/f2-f1_{true}/f2_{true}'};
names(4)={[names{2} '_proportion']};


for i=1:4
    colours={'b', 'k', 'r', 'g'};
    if 1
        for j=1:size(methods,2)
            if key=='corrs'
                plot(corresp,abs(stat_data(j*2-1,:,i)),[colours{j} 'x--']); %real
            else
                plot(noise,abs(stat_data(j*2-1,:,i)),[colours{j} 'x--']); %real
            end
            hold on
        end
    else
        for j=1:size(methods,2)/2
            numimag=cellfun(@(x) size(x,1),data(j*2,:,1));
            numreal=cellfun(@(x) size(x,1),data(j*2-1,:,1));
            both=numimag+numreal;
            %plot(corresp,abs([stat_data(j*2,:,i).*numimag + stat_data(j*2-1,:,i).*numreal]./both),[colours{j} 'x-']); %both
            plot(corresp,abs(stat_data(j*2-1,:,i)),[colours{j} 'x--']); %real
            hold on
            plot(corresp,abs(stat_data(j*2,:,i)),[colours{j} 'x:']); %imaginary
        end
    end
    %plot(corresp,abs(stat_data(1,:,i)), 'bx-',corresp,abs(stat_data(2,:,i)),'rx-' ...
    %    ,corresp,abs(stat_data(3,:,i)), 'mx-',corresp,abs(stat_data(4,:,i)),'kx-' ...
    %    ,corresp,abs(stat_data(5,:,i)), 'yx-',corresp,abs(stat_data(6,:,i)),'cx-');
    legend(legends(i,:));
    title(titles(i));
    if key=='corrs'
        xlabel('number of correspondences');
    else 
        xlabel('level of noise');
    end
    ylabel(ylabels(i));
    saveas(gcf,[names{i} '.fig']);
    saveas(gcf,[names{i} '.jpg']);
    %set(gcf, 'PaperUnits', 'points');
    %set(gcf, 'PaperPosition', [0 0 1024 800]);
    %set(gcf,'Position',[100, 100, 1024, 800]);
    hold off
end


%% boxplotting
clear titles legends ylabels names
titlesi=methods;
namesi={'8pt_real','8pt_imaginary','7pt_real','7pt_imaginary','fRatio_real',...
    'fRatio_imaginary','6Ratio_real','6Ratio_imaginary'};
%
titles(1)={'Estimating (f1,f2)'};
legends(1,:)={'f/f_{true}','ground truth'};
names(1)={'bougnoux/boxplot'};

titles(2)={'Estimating f1/f2'};
legends(2,:)={'f1/f2-f1_{true}/f2_{true}','ground truth'};
names(2)={'bougnoux/boxplot_proportion'};

for i=1:size(methods,2)
    for j=1:2
        %reformatting for dropbox grouping variable.
        %with preallocation
        data_temp=data(i,:,j);
        data_size=0;
        for k=1:size(data_temp,2)
            data_size=data_size+size(data_temp{k},1);
        end
        data_plot=NaN(1,data_size);
        grp=zeros(1,data_size);
        corrdisp=[];
        start=1;
        ende=0;
        for k=1:size(data_temp,2)
            d_t_size=size(data_temp{k},1);
            ende=ende+d_t_size;
            data_plot(start:ende)=data_temp{k}';
            grp(start:ende)=k*ones(1,d_t_size);
            corrdisp=[corrdisp k*ones(d_t_size>0)];
            start=start+d_t_size;
        end
        
        %plotting
        corrdisp=corrdisp+corresp(1)-1; % these are the labels of groups.
        % Here I shift it so that the first one is the first correspondence number, not '1'
        limits=[1-j 3-j];
        if all(data_plot>limits(2)) || all(data_plot<limits(1))
            limits=[-inf inf];
        end
        boxplot(data_plot,grp,'Labels',corrdisp,'Symbol','','DataLim',limits);
        hold on
        plot((2-j)*ones(1,size(corresp,2)),'gx');
        legend(legends(j,2));
        title([titles{j} '. ' titlesi{i}]);
        xlabel('number of correspondences');
        ylabel(legends(j,1));
        name=[names{j} '_' namesi{i}];
        saveas(gcf,[name '.fig']);
        saveas(gcf,[name '.jpg']);
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperPosition', [0 0 1024 800]);
        set(gcf,'Position',[100, 100, 1024, 800]);
        hold off
    end
end



