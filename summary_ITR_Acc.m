function summary_ITR_Acc(figdir,data_dir)

load(fullfile(data_dir,'Phase.mat'),'freqs','phases');
load(fullfile(data_dir,'channel_name.mat'),'channels');
[~,freqs_I]=sort(freqs);

fs=250;
ch_no=64;
freqs=freqs(1:end);
channel_select=[48 54 55 56 57 58 61 62 63];%1:length(channels);
channels=channels(channel_select);
subject_no=35;
trial_no=40;
block_no=6;
possible_T=0.25:0.25:5;
test_block=(1:6)';

possible_T=0.25:0.25:5;
delayTime=0.14;
start_t=0.5+0.14;

method_name={};
res={};
for source_no=[4 8 12]
    source_no_i=source_no/4;
    load(['random_source_' num2str(source_no) '.mat'],'select_source')
    for source_i=1:5%1:5
        freqs_compare=setdiff(1:40,select_source(source_i,:));
        disp([num2str(source_no) '-' num2str(source_i)])
        method_i=1;
        % CCA Results
        load('cca_results.mat','res_store')
        res{method_i,source_no_i,source_i,1}=zeros(size(res_store,3),size(res_store,1));
        res{method_i,source_no_i,source_i,2}=zeros(size(res_store,3),size(res_store,1));
        for sub_no=1:size(res_store,3)
            for T_i=1:length(possible_T)
                T=possible_T(T_i);
                tmp=squeeze(res_store(T_i,:,sub_no,freqs_compare,:));
                my_acc=sum(sum(sum(tmp)))/numel(tmp);
                real_tw=length(floor(start_t*fs):floor((start_t+T)*fs-1))/fs;
                my_tw=real_tw+0.5+0.14+0.2;%my_tw=real_tw+0.5;%real_tw+delayTime+0.02+0.5;
                if my_acc==1
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc));
                elseif my_acc<1/40
                    itr=0;
                else
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc)+(1-my_acc)*log2((1-my_acc)/(40-1)));
                end
                res{method_i,source_no_i,source_i,1}(sub_no,T_i)=my_acc;
                res{method_i,source_no_i,source_i,2}(sub_no,T_i)=itr;
            end
        end
        if source_no_i==1 && source_i==1
            method_name{1}='CCA';
        end
        % ECCA Results
        method_i=method_i+1;
        load('ecca_results_all.mat','res_store')
        res{method_i,source_no_i,source_i,1}=zeros(size(res_store,3),size(res_store,1));
        res{method_i,source_no_i,source_i,2}=zeros(size(res_store,3),size(res_store,1));
        for sub_no=1:size(res_store,3)
            for T_i=1:length(possible_T)
                T=possible_T(T_i);
                tmp=squeeze(res_store(T_i,:,sub_no,freqs_compare,:));
                my_acc=sum(sum(sum(tmp)))/numel(tmp);
                real_tw=length(floor(start_t*fs):floor((start_t+T)*fs-1))/fs;
                my_tw=real_tw+0.5+0.14+0.2;%my_tw=real_tw+0.5;%real_tw+delayTime+0.02+0.5;
                if my_acc==1
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc));
                elseif my_acc<1/40
                    itr=0;
                else
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc)+(1-my_acc)*log2((1-my_acc)/(40-1)));
                end
                res{method_i,source_no_i,source_i,1}(sub_no,T_i)=my_acc;
                res{method_i,source_no_i,source_i,2}(sub_no,T_i)=itr;
            end
        end
        if source_no_i==1 && source_i==1
            method_name{end+1}='eCCA';
        end
        % TRCA only template Results
        method_i=method_i+1;
        load('trca_results.mat','res_store')
        index_i=1;%1:size(res_store,6)
        res{method_i,source_no_i,source_i,1}=zeros(size(res_store,3),size(res_store,1));
        res{method_i,source_no_i,source_i,2}=zeros(size(res_store,3),size(res_store,1));
        for sub_no=1:size(res_store,3)
            for T_i=1:length(possible_T)
                T=possible_T(T_i);
                tmp=squeeze(res_store(T_i,:,sub_no,freqs_compare,:,index_i));
                my_acc=sum(sum(sum(tmp)))/numel(tmp);
                real_tw=length(floor(start_t*fs):floor((start_t+T)*fs-1))/fs;
                my_tw=real_tw+0.5+0.14+0.2;%my_tw=real_tw+0.5;%real_tw+delayTime+0.02+0.5;
                if my_acc==1
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc));
                elseif my_acc<1/40
                    itr=0;
                else
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc)+(1-my_acc)*log2((1-my_acc)/(40-1)));
                end
                res{method_i,source_no_i,source_i,1}(sub_no,T_i)=my_acc;
                res{method_i,source_no_i,source_i,2}(sub_no,T_i)=itr;
            end
        end
        if source_no_i==1 && source_i==1
            method_name{end+1}='TRCA';
        end
        % Propose AFD based method
        method_i=method_i+1;
        load(['averaged_AFD_template_neighbor_results_N' num2str(50) '_neighber_w_random' num2str(source_i) '_of_' num2str(source_no) '_end_good1and2.mat'],'res_store');
        index_i=1;%1:size(res_store,6)
        res{method_i,source_no_i,source_i,1}=zeros(size(res_store,3),size(res_store,1));
        res{method_i,source_no_i,source_i,2}=zeros(size(res_store,3),size(res_store,1));
        for sub_no=1:size(res_store,3)
            for T_i=1:length(possible_T)
                T=possible_T(T_i);
                tmp=squeeze(res_store(T_i,:,sub_no,freqs_compare,:,index_i,4));
                my_acc=sum(sum(sum(tmp)))/numel(tmp);
                real_tw=length(floor(start_t*fs):floor((start_t+T)*fs-1))/fs;
                my_tw=real_tw+0.5+0.14+0.2;%my_tw=real_tw+0.5;%real_tw+delayTime+0.02+0.5;
                if my_acc==1
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc));
                elseif my_acc<1/40
                    itr=0;
                else
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc)+(1-my_acc)*log2((1-my_acc)/(40-1)));
                end
                res{method_i,source_no_i,source_i,1}(sub_no,T_i)=my_acc;
                res{method_i,source_no_i,source_i,2}(sub_no,T_i)=itr;
            end
        end
        if source_no_i==1 && source_i==1
            method_name{end+1}='Proposed method';
        end
        % TRCA AFD template neighbor
        method_i=method_i+1;
        load(['trca_impulse_template_neighbor_results_random' num2str(source_i) '_of_' num2str(source_no) '_end_good_2.mat'],'res_store')
        index_i=1;%1:size(res_store,6)
        res{method_i,source_no_i,source_i,1}=zeros(size(res_store,3),size(res_store,1));
        res{method_i,source_no_i,source_i,2}=zeros(size(res_store,3),size(res_store,1));
        for sub_no=1:size(res_store,3)
            for T_i=1:length(possible_T)
                T=possible_T(T_i);
                tmp=squeeze(res_store(T_i,:,sub_no,freqs_compare,:,index_i,2));
                my_acc=sum(sum(sum(tmp)))/numel(tmp);
                real_tw=length(floor(start_t*fs):floor((start_t+T)*fs-1))/fs;
                my_tw=real_tw+0.5+0.14+0.2;%my_tw=real_tw+0.5;%real_tw+delayTime+0.02+0.5;
                if my_acc==1
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc));
                elseif my_acc<1/40
                    itr=0;
                else
                    itr=60/my_tw*(log2(40)+my_acc*log2(my_acc)+(1-my_acc)*log2((1-my_acc)/(40-1)));
                end
                res{method_i,source_no_i,source_i,1}(sub_no,T_i)=my_acc;
                res{method_i,source_no_i,source_i,2}(sub_no,T_i)=itr;
            end
        end
        if source_no_i==1 && source_i==1
            method_name{end+1}='tlCCA';
        end
    end
end
for res_i=1:size(res,1)
    disp([num2str(res_i) ': ' method_name{res_i}])
end
res_store=res;
% max ITR compare
possible_T=0.25:0.25:5;
[~,t_min]=min(abs(possible_T-0.75));
possible_T=possible_T(t_min:end);
for source_no=[4 8 12]
    source_no_i=source_no/4;
    res=cell(size(res_store,1),size(res_store,4));
    for method_i=1:size(res,1)
        for type_i=1:2
            for source_i=1:5
                res{method_i,type_i}(:,:,source_i)=res_store{method_i,source_no_i,source_i,type_i};
            end
            res{method_i,type_i}=mean(res{method_i,type_i},3);
            res{method_i,type_i}=res{method_i,type_i}(:,t_min:end);
        end
    end
    tmp1=[];
    tmp2=[];
    for i=1:35
        tmp1(i,1)=max(res{1,2}(i,:));
        tmp2(i,1)=max(res{4,2}(i,:));
    end
    figure
    plot([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],[min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],'r-','linewidth',2)
    hold on
    plot(tmp1,tmp2,'b.','markersize',12)
    set(gcf,'position',[1138         538         373         288])
    axis([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5 min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5])
    set(gca,'fontsize',10)
    grid on
    saveas(gcf,fullfile(figdir,['ITR_max_compare_14_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['ITR_max_compare_14_' num2str(source_no) '.png']),'-dpng','-r300');
    close(gcf)
    %
    tmp1=[];
    tmp2=[];
    for i=1:35
        tmp1(i,1)=max(res{5,2}(i,:));
        tmp2(i,1)=max(res{4,2}(i,:));
    end
    figure
    plot([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],[min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],'r-','linewidth',2)
    hold on
    plot(tmp1,tmp2,'b.','markersize',12)
    set(gcf,'position',[1138         538         373         288])
    axis([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5 min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5])
    set(gca,'fontsize',10)
    grid on
    saveas(gcf,fullfile(figdir,['ITR_max_compare_54_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['ITR_max_compare_54_' num2str(source_no) '.png']),'-dpng','-r300');
    close(gcf)
    %
    tmp1=[];
    tmp2=[];
    for i=1:35
        tmp1(i,1)=max(res{2,2}(i,:));
        tmp2(i,1)=max(res{4,2}(i,:));
    end
    figure
    plot([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],[min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5],'r-','linewidth',2)
    hold on
    plot(tmp1,tmp2,'b.','markersize',12)
    set(gcf,'position',[1138         538         373         288])
    axis([min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5 min(min([tmp1 tmp2]))-5 max(max([tmp1 tmp2]))+5])
    set(gca,'fontsize',10)
    grid on
    saveas(gcf,fullfile(figdir,['ITR_max_compare_24_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['ITR_max_compare_24_' num2str(source_no) '.png']),'-dpng','-r300');
    close(gcf)
end
% compare different methods
possible_T=0.25:0.25:5;
[~,t_min]=min(abs(possible_T-0.75));
possible_T=possible_T(t_min:end);
for source_no=[4 8 12]
    source_no_i=source_no/4;
    res=cell(size(res_store,1),size(res_store,4));
    for method_i=1:size(res,1)
        for type_i=1:2
            for source_i=1:5
                res{method_i,type_i}(:,:,source_i)=res_store{method_i,source_no_i,source_i,type_i};
            end
            res{method_i,type_i}=mean(res{method_i,type_i},3);
            res{method_i,type_i}=res{method_i,type_i}(:,t_min:end);
        end
    end
    plot_color=colormap(jet(size(res,1)));
    plot_color=plot_color(randperm(size(plot_color,1)),:);
    plot_color=[0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.9290 0.6940 0.1250;
                0.4940 0.1840 0.5560;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330;
                0.6350 0.0780 0.1840
                plot_color];
    close(gcf)
    % Plot ITR (no calibration data)
    h_figure=figure;
    hold on
    p_val=[];
    t_val=[];
    plot_method=[1 5 4];
    method_split=[];
    method_split(1)=find(cellfun(@(x) strcmp(x,'Proposed method'),method_name)==1);
    method_split(2)=find(cellfun(@(x) strcmp(x,'tlCCA'),method_name)==1);
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_i=1:size(res,1)
            itr_val(method_i,:)=res{method_i,2}(:,T_i);
        end
        pairs=nchoosek(1:size(itr_val,1), 2);
        %
        [~,p_val(T_i,1)]=ttest(itr_val([plot_method(1)],:)',itr_val([plot_method(end)],:)');
        [~,p_val(T_i,2)]=ttest(itr_val([plot_method(2)],:)',itr_val([plot_method(end)],:)');
    end
    for method_i=plot_method%1:size(res,1)
        if method_i<method_split(1)
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'o-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        elseif method_i<method_split(2)
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'x-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        else
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'s-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        end
    end
    %
    method_i=plot_method;
    for i=1:size(p_val,2)
        for T_i=1:size(p_val,1)
            if p_val(T_i,i)<=1E-3
                p_char='***';
            elseif p_val(T_i,i)<=1E-2
                p_char='**';
            elseif p_val(T_i,i)<=0.05
                p_char='*';
            else
                p_char='';
            end
            if size(p_val,2)==1
                dd=0;
            else
                if i==1
                    dd=-1;
                else
                    dd=1;
                end
            end
            if size(p_val,2)==1
                if strcmp(p_char,'n.s.')
                    text(possible_T(T_i)+dd*0.05,160,'n.s.',...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                else
                    for j=1:length(p_char)
                        text(possible_T(T_i)+dd*0.05,165-5*j,p_char(j),...
                            'HorizontalAlignment','Center',...
                            'BackGroundColor','none',...
                            'Tag','sigstar_stars',...
                            'fontsize',12,...
                            'color',plot_color(method_i(i),:));
                    end
                end
            else
                if strcmp(p_char,'n.s.')
                    text(possible_T(T_i)+dd*0.05,160,'n.s.',...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                else
                    for j=1:length(p_char)
                        text(possible_T(T_i)+dd*0.05,165-5*j,p_char(j),...
                            'HorizontalAlignment','Center',...
                            'BackGroundColor','none',...
                            'Tag','sigstar_stars',...
                            'fontsize',12,...
                            'color',plot_color(method_i(i),:));
                    end
                end
            end
        end
    end
    %
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('ITR')
    set(gca,'xtick',possible_T(1):0.25:5)
    axis([possible_T(1)-0.25 5.25 0 170])
    legend(method_name(plot_method),'location','northeastoutside')%legend(method_name([1 4]),'location','northeastoutside')
    set(gcf,'position',[240         379        1606         420])
    saveas(h_figure,fullfile(figdir,['ITR_summary_nocalibration_random_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['ITR_summary_nocalibration_random_' num2str(source_no) '.png']),'-dpng','-r300');
    close(h_figure)
    % Plot Acc (no calibration data)
    h_figure=figure;
    hold on
    p_val=[];
    t_val=[];
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_i=1:size(res,1)
            itr_val(method_i,:)=res{method_i,1}(:,T_i);
        end
        pairs=nchoosek(1:size(itr_val,1), 2);
        [~,p_val(T_i,1)]=ttest(itr_val([plot_method(1)],:)',itr_val([plot_method(end)],:)');
        [~,p_val(T_i,2)]=ttest(itr_val([plot_method(2)],:)',itr_val([plot_method(end)],:)');
    end
    for method_i=plot_method%1:size(res,1)
        if method_i<method_split(1)
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'o-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        elseif method_i<method_split(2)
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'x-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        else
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'s-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        end
    end
    %
    method_i=plot_method;
    for i=1:size(p_val,2)
        for T_i=1:size(p_val,1)
            if p_val(T_i,i)<=1E-3
                p_char='***';
            elseif p_val(T_i,i)<=1E-2
                p_char='**';
            elseif p_val(T_i,i)<=0.05
                p_char='*';
            else
                p_char='';
            end
            if size(p_val,2)==1
                dd=0;
            else
                if i==1
                    dd=-1;
                else
                    dd=1;
                end
            end
            if size(p_val,2)==1
                if strcmp(p_char,'n.s.')
                    text(possible_T(T_i)+dd*0.05,1.15,'n.s.',...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                else
                    for j=1:length(p_char)
                        text(possible_T(T_i)+dd*0.05,1.15-0.03*j,p_char(j),...
                            'HorizontalAlignment','Center',...
                            'BackGroundColor','none',...
                            'Tag','sigstar_stars',...
                            'fontsize',12,...
                            'color',plot_color(method_i(i),:));
                    end
                end
            else
                if strcmp(p_char,'n.s.')
                    text(possible_T(T_i)+dd*0.05,1.1,'n.s.',...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                else
                    for j=1:length(p_char)
                        text(possible_T(T_i)+dd*0.05,1.15-0.03*j,p_char(j),...
                            'HorizontalAlignment','Center',...
                            'BackGroundColor','none',...
                            'Tag','sigstar_stars',...
                            'fontsize',12,...
                            'color',plot_color(method_i(i),:));
                    end
                end
            end
        end
    end
    %
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('Acc')
    set(gca,'xtick',possible_T(1):0.25:5)
    axis([possible_T(1)-0.25 5.25 0 1.15])
    legend(method_name(plot_method),'location','northeastoutside')
    set(gcf,'position',[240         379        1606         420])
    saveas(h_figure,fullfile(figdir,['Acc_summary_nocalibration_random_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['Acc_summary_nocalibration_random_' num2str(source_no) '.png']),'-dpng','-r300');
    close(h_figure)
    % Plot ITR (have calibration data)
    h_figure=figure;
    hold on
    p_val=[];
    t_val=[];
    plot_method=[2 4];
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_i=1:size(res,1)
            itr_val(method_i,:)=res{method_i,2}(:,T_i);
        end
        pairs=nchoosek(1:size(itr_val,1), 2);
        [~,p_val(T_i,1)]=ttest(itr_val([plot_method(1)],:)',itr_val([plot_method(end)],:)');
        if length(plot_method)>2
            [~,p_val(T_i,2)]=ttest(itr_val([plot_method(2)],:)',itr_val([plot_method(end)],:)');
        end
    end
    for method_i=plot_method%1:size(res,1)
        if method_i<method_split(1)
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'d-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        elseif method_i<method_split(2)
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'x-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        else
            shadedErrorBar(possible_T,res{method_i,2},{@mean,@cal_CI95},'lineProps',{'s-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        end
    end
    %
    method_i=plot_method;
    for i=1:size(p_val,2)
        for T_i=1:size(p_val,1)
            if p_val(T_i,i)<=1E-3
                p_char='***';
            elseif p_val(T_i,i)<=1E-2
                p_char='**';
            elseif p_val(T_i,i)<=0.05
                p_char='*';
            else
                p_char='';
            end
            if size(p_val,2)==1
                dd=0;
            else
                if i==1
                    dd=-1;
                else
                    dd=1;
                end
            end
            if strcmp(p_char,'n.s.')
                text(possible_T(T_i)+dd*0.05,210,'n.s.',...
                    'HorizontalAlignment','Center',...
                    'BackGroundColor','none',...
                    'Tag','sigstar_stars',...
                    'fontsize',12,...
                    'color',plot_color(method_i(i),:));
            else
                for j=1:length(p_char)
                    text(possible_T(T_i)+dd*0.05,210-5*j,p_char(j),...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                end
            end
        end
    end
    %
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('ITR')
    set(gca,'xtick',possible_T(1):0.25:5)
    axis([possible_T(1)-0.25 5.25 0 215])
    legend(method_name(plot_method),'location','northeastoutside')
    set(gcf,'position',[240         379        1606         420])
    saveas(h_figure,fullfile(figdir,['ITR_summary_calibration_random_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['ITR_summary_calibration_random_' num2str(source_no) '.png']),'-dpng','-r300');
    close(h_figure)
    % Plot Acc (have calibration data)
    h_figure=figure;
    hold on
    p_val=[];
    t_val=[];
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_i=1:size(res,1)
            itr_val(method_i,:)=res{method_i,1}(:,T_i);
        end
        pairs=nchoosek(1:size(itr_val,1), 2);
        [~,p_val(T_i,1)]=ttest(itr_val([plot_method(1)],:)',itr_val([plot_method(end)],:)');
        if length(plot_method)>2
            [~,p_val(T_i,2)]=ttest(itr_val([plot_method(2)],:)',itr_val([plot_method(end)],:)');
        end
    end
    for method_i=plot_method%1:size(res,1)
        if method_i<method_split(1)
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'d-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        elseif method_i<method_split(2)
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'x-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        else
            shadedErrorBar(possible_T,res{method_i,1},{@mean,@cal_CI95},'lineProps',{'s-','color',plot_color(method_i,:),'linewidth',2,'markersize',8});
        end
    end
    %
    method_i=plot_method;
    for i=1:size(p_val,2)
        for T_i=1:size(p_val,1)
            if p_val(T_i,i)<=1E-3
                p_char='***';
            elseif p_val(T_i,i)<=1E-2
                p_char='**';
            elseif p_val(T_i,i)<=0.05
                p_char='*';
            else
                p_char='';
            end
            if size(p_val,2)==1
                dd=0;
            else
                if i==1
                    dd=-1;
                else
                    dd=1;
                end
            end
            if strcmp(p_char,'n.s.')
                text(possible_T(T_i)+dd*0.05,1.15,'n.s.',...
                    'HorizontalAlignment','Center',...
                    'BackGroundColor','none',...
                    'Tag','sigstar_stars',...
                    'fontsize',12,...
                    'color',plot_color(method_i(i),:));
            else
                for j=1:length(p_char)
                    text(possible_T(T_i)+dd*0.05,1.15-0.03*j,p_char(j),...
                        'HorizontalAlignment','Center',...
                        'BackGroundColor','none',...
                        'Tag','sigstar_stars',...
                        'fontsize',12,...
                        'color',plot_color(method_i(i),:));
                end
            end
        end
    end
    %
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('Acc')
    set(gca,'xtick',possible_T(1):0.25:5)
    axis([possible_T(1)-0.25 5.25 0 1.15])
    legend(method_name(plot_method),'location','northeastoutside')
    set(gcf,'position',[240         379        1606         420])
    saveas(h_figure,fullfile(figdir,['Acc_summary_calibration_random_' num2str(source_no) '.fig']))
    print(fullfile(figdir,['Acc_summary_calibration_random_' num2str(source_no) '.png']),'-dpng','-r300');
    close(h_figure)
end
% compare different source number
for method_i=[4]
    res=cell(3,size(res_store,4));
    method_name={};
    for source_no=[4 8 12]
        source_no_i=source_no/4;
        method_name{1,source_no_i}=['Source No. ' num2str(source_no_i*4)];
        for type_i=1:2
            for source_i=1:5
                res{source_no/4,type_i}(:,:,source_i)=res_store{method_i,source_no_i,source_i,type_i};
            end
            res{source_no/4,type_i}=mean(res{source_no/4,type_i},3);
            res{source_no/4,type_i}=res{source_no/4,type_i}(:,t_min:end);
        end
    end
    plot_color=colormap(jet(size(res,1)));
    plot_color=plot_color(randperm(size(plot_color,1)),:);
    plot_color=[0 0.4470 0.7410;
                0.8500 0.3250 0.0980;
                0.4940 0.1840 0.5560;
                0.9290 0.6940 0.1250;
                0.4660 0.6740 0.1880;
                0.3010 0.7450 0.9330;
                0.6350 0.0780 0.1840
                plot_color];
    close(gcf)
    % Plot ITR (no calibration data)
    h_figure=figure;
    hold on
    t_val=[];
    
    p_val=[];
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_j=1:size(res,1)
            itr_val(method_j,:)=res{method_j,2}(:,T_i);
        end
        [~,p_val(1,T_i)]=ttest(itr_val(1,:)',itr_val(2,:)');
        [~,p_val(2,T_i)]=ttest(itr_val(2,:)',itr_val(3,:)');
        [~,p_val(3,T_i)]=ttest(itr_val(1,:)',itr_val(3,:)');
    end

    bar_plot=[];
    for plot_method=1:size(res,1)
        bar_plot(:,plot_method)=mean(res{plot_method,2}).';
    end

    bar(possible_T,bar_plot);
    set(gcf,'position',[1985         538        1717         271])
    tmp=get(gca,'children');
    res_tmp=res([3 2 1],:);
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('ITR')

    for T_i=1:length(possible_T)
        H=sigstar({[tmp(3).XData(T_i)+tmp(3).XOffset,tmp(2).XData(T_i)+tmp(2).XOffset],...
                 [tmp(1).XData(T_i)+tmp(1).XOffset,tmp(2).XData(T_i)+tmp(2).XOffset],...
                 [tmp(3).XData(T_i)+tmp(3).XOffset,tmp(1).XData(T_i)+tmp(1).XOffset]},...
                 [p_val(1,T_i),...
                 p_val(2,T_i),...
                 p_val(3,T_i)]);
        new_tmp=get(gca,'children');
        if method_i==4
            if strcmp(new_tmp(1).String,'')
                new_tmp(1).Position(2)=1000000;
                new_tmp(2).YData=new_tmp(2).YData+100000;
            else
                new_tmp(1).Position(2)=156+5+5+5;
                new_tmp(2).YData=new_tmp(2).YData+(155+5+5+5-0.5-new_tmp(2).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(3).Position(2)=1000000;
                new_tmp(4).YData=new_tmp(2).YData+100000;
            else
                new_tmp(3).Position(2)=148-5+5+5+5;
                new_tmp(4).YData=new_tmp(4).YData+(147+5+5+5-5-0.5-new_tmp(4).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(5).Position(2)=1000000;
                new_tmp(6).YData=new_tmp(2).YData+100000;
            else
                new_tmp(5).Position(2)=140-5*2+5+5+5;
                new_tmp(6).YData=new_tmp(6).YData+(139+5+5+5-5*2-0.5-0.5-new_tmp(6).YData(2));
            end
        else
            if strcmp(new_tmp(1).String,'')
                new_tmp(1).Position(2)=1000000;
                new_tmp(2).YData=new_tmp(2).YData+100000;
            else
                new_tmp(1).Position(2)=156+5+5+5-20-5;
                new_tmp(2).YData=new_tmp(2).YData+(155-5+5-20+5+5-0.5-new_tmp(2).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(3).Position(2)=1000000;
                new_tmp(4).YData=new_tmp(2).YData+100000;
            else
                new_tmp(3).Position(2)=148-5+5+5+5-20-5;
                new_tmp(4).YData=new_tmp(4).YData+(147-5+5-20+5+5-5-0.5-new_tmp(4).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(5).Position(2)=1000000;
                new_tmp(6).YData=new_tmp(2).YData+100000;
            else
                new_tmp(5).Position(2)=140-5*2+5+5+5-20-5;
                new_tmp(6).YData=new_tmp(6).YData+(139+5-5-20+5+5-5*2-0.5-0.5-new_tmp(6).YData(2));
            end
        end
    end
    for tmp_i=1:3
        error_bar_x=tmp(tmp_i).XData+tmp(tmp_i).XOffset;
        error_bar_y=mean(res_tmp{tmp_i,2});
        error_bar_err=cal_CI95(res_tmp{tmp_i,2});
        errorbar(error_bar_x,error_bar_y,error_bar_err,'k.');
    end
    legend(method_name,'location','northeastoutside')%legend(method_name([1 4]),'location','northeastoutside')
    xtick=unique([tmp(1).XData+tmp(1).XOffset tmp(2).XData+tmp(2).XOffset tmp(3).XData+tmp(3).XOffset tmp(1).XData tmp(2).XData tmp(3).XData]);
    xticklabel=cell(1,length(xtick));
    for x_i=1:length(xtick)
        if mod(x_i,3)==1 || mod(x_i,3)==0
            xticklabel{1,x_i}='';
        else
            xticklabel{1,x_i}=num2str(xtick(x_i));
        end
    end
    set(gca,'xtick',xtick);
    set(gca,'xticklabel',xticklabel);
    if method_i==4
        axis([0.6278    5.1222 0 182])
    else
        axis([0.6278    5.1222 0 152])
    end
    saveas(gcf,fullfile(figdir,['ITR_summary_nocalibration_random_Method' num2str(method_i) '_bar_plot.fig']))
    print(fullfile(figdir,['ITR_summary_nocalibration_random_Method' num2str(method_i) '_bar_plot.png']),'-dpng','-r300');
    close(gcf)
    
    % Plot Acc (no calibration data)
    h_figure=figure;
    hold on
    t_val=[];
    
    p_val=[];
    for T_i=1:length(possible_T)
        itr_val=[];
        for method_j=1:size(res,1)
            itr_val(method_j,:)=res{method_j,1}(:,T_i);
        end
        [~,p_val(1,T_i)]=ttest(itr_val(1,:)',itr_val(2,:)');
        [~,p_val(2,T_i)]=ttest(itr_val(2,:)',itr_val(3,:)');
        [~,p_val(3,T_i)]=ttest(itr_val(1,:)',itr_val(3,:)');
    end
    
    bar_plot=[];
    for plot_method=1:size(res,1)
        bar_plot(:,plot_method)=mean(res{plot_method,1}).';
    end

    bar(possible_T,bar_plot);
    grid on
    set(gca,'fontsize',12)
    xlabel('Sig Len (s)')
    ylabel('Acc')
    tmp=get(gca,'children');
    for T_i=1:length(possible_T)
        H=sigstar({[tmp(3).XData(T_i)+tmp(3).XOffset,tmp(2).XData(T_i)+tmp(2).XOffset],...
                 [tmp(1).XData(T_i)+tmp(1).XOffset,tmp(2).XData(T_i)+tmp(2).XOffset],...
                 [tmp(3).XData(T_i)+tmp(3).XOffset,tmp(1).XData(T_i)+tmp(1).XOffset]},...
                 [p_val(1,T_i),...
                 p_val(2,T_i),...
                 p_val(3,T_i)]);
        new_tmp=get(gca,'children');
        if method_i==4
            if strcmp(new_tmp(1).String,'')
                new_tmp(1).Position(2)=1000000;
                new_tmp(2).YData=new_tmp(2).YData+100000;
            else
                new_tmp(1).Position(2)=156/100-20/100;
                new_tmp(2).YData=new_tmp(2).YData+(155/100-20/100-0.5/100-new_tmp(2).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(3).Position(2)=1000000;
                new_tmp(4).YData=new_tmp(2).YData+100000;
            else
                new_tmp(3).Position(2)=148/100-5/100-20/100;
                new_tmp(4).YData=new_tmp(4).YData+(147/100-20/100-5/100-0.5/100-new_tmp(4).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(5).Position(2)=1000000;
                new_tmp(6).YData=new_tmp(2).YData+100000;
            else
                new_tmp(5).Position(2)=140/100-5/100*2-20/100;
                new_tmp(6).YData=new_tmp(6).YData+(139/100-20/100-5/100*2-0.5/100-new_tmp(6).YData(2));
            end
        else
            if strcmp(new_tmp(1).String,'')
                new_tmp(1).Position(2)=1000000;
                new_tmp(2).YData=new_tmp(2).YData+100000;
            else
                new_tmp(1).Position(2)=156/100-20/100;
                new_tmp(2).YData=new_tmp(2).YData+(155/100-20/100-0.5/100-new_tmp(2).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(3).Position(2)=1000000;
                new_tmp(4).YData=new_tmp(2).YData+100000;
            else
                new_tmp(3).Position(2)=148/100-5/100-20/100;
                new_tmp(4).YData=new_tmp(4).YData+(147/100-20/100-5/100-0.5/100-new_tmp(4).YData(2));
            end
            if strcmp(new_tmp(3).String,'')
                new_tmp(5).Position(2)=1000000;
                new_tmp(6).YData=new_tmp(2).YData+100000;
            else
                new_tmp(5).Position(2)=140/100-5/100*2-20/100;
                new_tmp(6).YData=new_tmp(6).YData+(139/100-20/100-5/100*2-0.5/100-new_tmp(6).YData(2));
            end
        end
    end

    for tmp_i=1:3
        error_bar_x=tmp(tmp_i).XData+tmp(tmp_i).XOffset;
        error_bar_y=mean(res_tmp{tmp_i,1});
        error_bar_err=cal_CI95(res_tmp{tmp_i,1});
        errorbar(error_bar_x,error_bar_y,error_bar_err,'k.');
    end
    legend(method_name,'location','northeastoutside')%legend(method_name([1 4]),'location','northeastoutside')
    xtick=unique([tmp(1).XData+tmp(1).XOffset tmp(2).XData+tmp(2).XOffset tmp(3).XData+tmp(3).XOffset tmp(1).XData tmp(2).XData tmp(3).XData]);
    xticklabel=cell(1,length(xtick));
    for x_i=1:length(xtick)
        if mod(x_i,3)==1 || mod(x_i,3)==0
            xticklabel{1,x_i}='';
        else
            xticklabel{1,x_i}=num2str(xtick(x_i));
        end
    end
    set(gca,'xtick',xtick);
    set(gca,'xticklabel',xticklabel);
    axis([0.6278    5.1222 0 1.45])
    set(gcf,'position',[1985         538        1717         271])
    
    saveas(gcf,fullfile(figdir,['Acc_summary_nocalibration_random_Method' num2str(method_i) '_bar_plot.fig']))
    print(fullfile(figdir,['Acc_summary_nocalibration_random_Method' num2str(method_i) '_bar_plot.png']),'-dpng','-r300');
    close(gcf)
end

