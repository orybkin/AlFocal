% A script for generating scatter plot series of Bougnoux formula results
% 
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

% 

calculating=1;
plotting=1;
boxplotting=1;
deepagpaths;

if calculating
    % set parameters
    file='../../data/paris/correspondences_F_synth_1K.mat';
    %file='../../data/paris/correspondences_F.mat';
    corr=[7:17 18:3:29 30:5:49 50:10:100];  % numbers of correspondences - x axis
    corr=[7:20];
    %corr=[8:12];
    %corr=[7:9];
    %corr=[8];
    pop_size=1000; % number of points that will be on scatter plot
    n=size(corr,2); % x axis length
    clear data stat_data
    stat_data=NaN(4,n,4);
    data=cell(4,n,2);
        
    % size(stat_data):=(methods, n, parameters), where:
    % methods are as in the legend below
    % parameters are [absolute_error a_std relative_error r_std]
    % size(data):=(methods, n, parameters), where:
    % methods are as in the legend below
    % parameters are [absolute relative]
    perm32=@(x) permute(x,[1 3 2]);
    resh22=@(x) perm32(reshape(x,2,2));
    for i=1:n
        Fparam.corr=corr(i);
        if corr(i)>7
            % beware, data are magically reshaped here
            [stat_temp, data_temp]=bougnoux_scatter(file,corr(i),pop_size,'Free');
            data(1:2,i,:)=resh22(data_temp);
            stat_data(1:2,i,:)=perm32(stat_temp);
        else
            data(1:2,i,:)={[]};
        end
        [stat_temp, data_temp]=bougnoux_scatter(file,corr(i),pop_size,'|F|=0');
        data(3:4,i,:)=resh22(data_temp);
        stat_data(3:4,i,:)=perm32(stat_temp);
    end
    clear data_temps stat_temp;
end

if plotting
    clear titles legends ylabels names
    
    titles(1)={'Error of estimating (f1,f2)'};
    legends(1,:)={'only real f', 'all f','singularized F, only real f', 'singularized F, all f'};
    ylabels(1)={'mean of f1/f1_true'};
    names(1)={'bougnoux/error'};
    
    titles(2)={'Std of estimating (f1,f2)'};
    legends(2,:)=legends(1,:);
    ylabels(2)={'std of f1/f1_true'};
    names(2)={'bougnoux/std'};
    
    titles(3)={'Mean error of estimating f1/f2'};
    legends(3,:)=legends(1,:);
    ylabels(3)={'error of f1/f2-f1_true/f2_true'};
    names(3)={[names{1} '_proportion']};
    
    titles(4)={'Std of estimating f1/f2'};
    legends(4,:)=legends(1,:);
    ylabels(4)={'std of f1/f2-f1_true/f2_true'};
    names(4)={[names{2} '_proportion']};
    
    for i=1:4
        plot(corr,abs(stat_data(1,:,i)), 'bx-',corr,abs(stat_data(2,:,i)),'rx-' ...
            ,corr,abs(stat_data(3,:,i)), 'mx-',corr,abs(stat_data(4,:,i)),'kx-');
        legend(legends(i,:));
        title(titles(i));
        xlabel('number of correspondences');
        ylabel(ylabels(i));
        saveas(gcf,[names{i} '.fig']);
        saveas(gcf,[names{i} '.jpg']);
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperPosition', [0 0 1024 800]);
        set(gcf,'Position',[100, 100, 1024, 800]);
        hold off
    end
end

if boxplotting
    clear titles legends ylabels names
    titlesi={'Only real f', 'All f','Singularized F, only real f', 'Singularized F, all f'};
    namesi={'real','all','sinreal','sinall'};
    %
    titles(1)={'Estimating (f1,f2)'};
    legends(1,:)={'f/f_true','ground truth'};
    names(1)={'bougnoux/boxplot'};
    
    titles(2)={'Estimating f1/f2'};
    legends(2,:)={'f1/f2-f1_true/f2_true','ground truth'};
    names(2)={'bougnoux/boxplot_proportion'};
    
    for i=1:4
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
            start=1;
            ende=0;
            for k=1:size(data_temp,2)
                d_t_size=size(data_temp{k},1);
                ende=ende+d_t_size;
                data_plot(start:ende)=data_temp{k}';
                grp(start:ende)=k*ones(1,d_t_size);
                start=start+d_t_size;
            end
            
            %plotting
            if i<3 && corr(1)==7
                corrdisp=corr(2:end);
            else
                corrdisp=corr;
            end
            boxplot(data_plot,grp,'Labels',corrdisp,'Symbol','','DataLim',[1-j 3-j]);
            hold on
            plot((2-j)*ones(1,size(corr,2)),'gx');
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
end


