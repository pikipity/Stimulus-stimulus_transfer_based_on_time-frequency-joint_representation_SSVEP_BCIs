function study_AFD_coef_averaged_template(data_dir,rand_total_no,source_input)

load(fullfile(data_dir,'Phase.mat'),'freqs','phases');
load(fullfile(data_dir,'channel_name.mat'),'channels');
[~,freqs_I]=sort(freqs);

fs=250;
ch_no=64;
start_t=0.5+0.14;
freqs=freqs(1:end);
channel_select=[48 54 55 56 57 58 61 62 63];%1:length(channels);
channels=channels(channel_select);
subject_no=35;
trial_no=40;
block_no=6;
possible_T=0.25:0.25:5;
test_block=(1:6)';

AFD_N=50;
source_i=source_input;
load(['Averaged_template_AFD_decompose_t_move_N' num2str(AFD_N) '_random' num2str(source_i) '_of_' num2str(rand_total_no) '.mat'],...
      'x_ori_store',...
      't_store',...
      'S1_store',...
      'a_store',...
      'select_source')



re_f_abs_same_phase_same_select=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select));
re_f_1period_abs_same_phase_same_select=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select));



for T_i=1:length(possible_T)
    T=possible_T(T_i);
    for test_run=1:size(test_block,1)
        for sub_no=1:subject_no
            for ch_i=1:length(channel_select)
                disp(['Analyze AFD -> ' num2str(source_i) '/' num2str(size(select_source,1)) ' sub ' num2str(sub_no) ', T=' num2str(T) ', run ' num2str(test_run) ', ch ' num2str(ch_i)])
                
                x_diff_stim=x_ori_store{T_i,test_run,sub_no,ch_i,source_i};
                t=t_store{T_i,test_run,sub_no,ch_i,source_i};
                S1_store_tmp=S1_store{T_i,test_run,sub_no,ch_i,source_i};
                a=a_store{T_i,test_run,sub_no,ch_i,source_i};
                
                t_1period=[];
                freqs_select=1:40;
                t_tmp=t(freqs_select,:);
                t_tmp_1period=zeros(size(t_tmp,1),size(t_tmp,2)+1);
                for i=1:size(t_tmp_1period,1)
                    t_tmp_1period(i,:)=linspace(0,2*pi,size(t_tmp,2)+1);
                end
                t_tmp_1period=t_tmp_1period(:,1:end-1);
                t_1period(freqs_select,:)=t_tmp_1period;
                %
                x_diff_stim_new=[];
                x_diff_stim_new_1=[];
                re_f_recon_1period=[];
                re_f_recon_1period_1=[];
                freqs_select_target=select_source(source_i,:);
                abs_S1=abs(S1_store_tmp(:,freqs_select_target));
                keep_level=[];
                for level_i=1:size(abs_S1,1)
                    total_energy=sum(abs_S1,1);
                    percentage_energy=sum(abs_S1(1:level_i,:),1)./total_energy;
                    if ~isempty(find(percentage_energy<0.9))
                        keep_level=[keep_level level_i];
                    end
                end
                reject_level=[];
                num_outlier=[];
                IQR=[];
                span_len=[];
                for level_i=keep_level
                    [num_outlier(level_i),~,~,~,~,~,~,IQR(level_i,1)]=iqr_wiki(abs_S1(level_i,:));
                    span_len(level_i,1)=(max(abs_S1(level_i,:))-min(abs_S1(level_i,:)));
                    if length(abs_S1(level_i,:))-num_outlier(level_i)<floor(0.8*length(abs_S1(level_i,:)))
                        reject_level=[reject_level level_i];
                    end
                end
                % 
                keep_level=setdiff(keep_level,reject_level);
                for freqs_target_i=1:length(freqs_select_target)
                    freqs_select=[];
                    for freqs_i=1:40
                        [~,target_i]=min(abs(freqs(freqs_select_target)-freqs(freqs_i)));
                        if target_i==freqs_target_i
                            freqs_select=[freqs_select freqs_i];
                        end
                    end
                    x_diff_stim_tmp=x_diff_stim(freqs_select,:);
                    hilbert_x=hilbert(x_diff_stim_tmp.').';
                    t_tmp=t(freqs_select,:);
                    t_tmp_1period=t_1period(freqs_select,:);
                    S1=S1_store_tmp(:,freqs_select);
                    abs_S1=abs(S1);
                    phase_S1=angle(S1);
                    abs_S1_new=[];
                    phase_S1_new=[];
                    source_target=find(freqs_select==intersect(freqs_select_target,freqs_select));
                    for i=1:size(phase_S1,2)
                        abs_S1_new(:,i)=abs_S1(:,source_target);
                        phase_S1_new(:,i)=phase_S1(:,source_target);
                    end
                    S1_new=abs_S1_new.*exp(1j.*phase_S1_new);
                    %
                    afdcal=AFDCal(hilbert_x);
                    afdcal.setDecompMethod(3);
                    for t_re_ch=1:size(t_tmp,1)
                        afdcal.setPhase(t_re_ch,t_tmp(t_re_ch,:));
                    end
                    afdcal.set_an({a.'});
                    coef_tmp={};
                    for n=1:AFD_N
                        for ch_afd=1:size(hilbert_x,1)
                            coef_tmp{ch_afd,1}(1,n)=S1_new(n,ch_afd);
                        end
                    end
                    afdcal.set_coef(coef_tmp);
                    afdcal.init_decomp(0,0);
                    for n=1:AFD_N-1
                        afdcal.nextDecomp(0,0);
                    end
                    F=[];
                    for n=1:AFD_N
                        for ch_afd=1:size(hilbert_x,1)
                            F(n,ch_afd,:)=afdcal.deComp{ch_afd}(n,:);
                        end
                    end
                    x_diff_stim_new_1(freqs_select,:)=squeeze(real(sum(F(keep_level,:,:),1)));
                    %
                    afdcal=AFDCal(hilbert_x);
                    afdcal.setDecompMethod(3);
                    for t_re_ch=1:size(t_tmp_1period,1)
                        afdcal.setPhase(t_re_ch,t_tmp_1period(t_re_ch,:));
                    end
                    afdcal.set_an({a.'});
                    afdcal.set_coef(coef_tmp);
                    afdcal.init_decomp(0,0);
                    for n=1:AFD_N-1
                        afdcal.nextDecomp(0,0);
                    end
                    F=[];
                    for n=1:AFD_N
                        for ch_afd=1:size(hilbert_x,1)
                            F(n,ch_afd,:)=afdcal.deComp{ch_afd}(n,:);
                        end
                    end
                    re_f_recon_1period_1(freqs_select,:)=squeeze(real(sum(F(keep_level,:,:),1)));
                end

                re_f_abs_same_phase_same_select{T_i,test_run,sub_no,ch_i}=x_diff_stim_new_1;
                re_f_1period_abs_same_phase_same_select{T_i,test_run,sub_no,ch_i}=re_f_recon_1period_1;
                

                
            end
            
            

        end
    end
end

save(['study_S1_averaged_template_t_move_random' num2str(source_i) '_of_' num2str(rand_total_no) '.mat'],...
     're_f_abs_same_phase_same_select',...
     're_f_1period_abs_same_phase_same_select')
