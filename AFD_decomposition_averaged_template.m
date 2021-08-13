function AFD_decomposition_averaged_template(data_dir,rand_total_no,source_input)

load(fullfile(data_dir,'Phase.mat'),'freqs','phases');
load(fullfile(data_dir,'channel_name.mat'),'channels');
fs=250;
ch_no=64;
start_t=0.5+0.14;
freqs=freqs(1:end);
channel_select=[48 54 55 56 57 58 61 62 63];%1:length(channels);
channels=channels(channel_select);
subject_no=35;
trial_no=40;
block_no=6;


f_d=[7 90];
[filter_b, filter_a]=butter(4,f_d./(fs/2),'bandpass');

test_block=(1:6)';
train_block=zeros(size(test_block,1),6-size(test_block,2));
for test_run=1:size(test_block,1)
    train_block(test_run,:)=setdiff(1:6,test_block(test_run,:));
end

for i=1:length(freqs)
    sine_ref{i}=gen_ref_sin(freqs(i),fs,fs*6,5,phases(i));
end

if exist('subj_data.mat')
    load('subj_data.mat','subj_data')
else
    subj_data=cell(1,subject_no);
    for sub_no=1:subject_no
        sub=['S' num2str(sub_no)];
        disp(['Filter: ' sub])
        load(fullfile(data_dir,[sub '.mat']));
        y=data;
        for trial=1:trial_no
            for block=1:block_no
                for ch=1:size(y,1)
                    temp=squeeze(y(ch,:,trial,block));
                    temp=detrend(temp);
                    y(ch,:,trial,block)=filtfilt(filter_b,filter_a,temp);
                end
            end
        end
        subj_data{sub_no}=y;
    end
    save('subj_data.mat','subj_data');
end

possible_T=0.25:0.25:5;
possible_index1=1:9;



if exist('averaged_ref.mat')
    load('averaged_ref.mat','averaged_ref')
else
    averaged_ref=cell(length(possible_T),size(test_block,1),subject_no,trial_no);

    for T_i=1:length(possible_T)
        T=possible_T(T_i);
        for test_run=1:size(test_block,1)
            for sub_no=1:subject_no
                sub=['S' num2str(sub_no)];

                y=subj_data{sub_no};

                for trial=1:trial_no
                    disp(['Averaged Signal  ->  Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial)])

                    x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,train_block(test_run,:)));
                    averaged_ref{T_i,test_run,sub_no,trial}=mean(x,3);
                end

            end
        end
    end
    
    save('averaged_ref.mat','averaged_ref');
end



disp('Initial')
AFD_N=50;
index1=1;

% rand_total_no=8;
if exist(['random_source_' num2str(rand_total_no) '.mat'])
    load(['random_source_' num2str(rand_total_no) '.mat'],'select_source')
else
    s=RandStream('mlfg6331_64');
    select_source=[];
    rand_i=1;
    while rand_i<=10
        randsource=sort(randsample(s,40,rand_total_no));
        if min(diff(randsource))<2
        else
            [sort_freqs,sort_freqs_I]=sort(freqs);
            select_source(rand_i,:)=sort_freqs_I(randsource);
            rand_i=rand_i+1;
        end
    end
    save(['random_source_' num2str(rand_total_no) '.mat'],'select_source')
end
    
AFD_r=zeros(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
S1_store=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
x_recon_store=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
x_ori_store=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
a_store=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
t_store=cell(length(possible_T),size(test_block,1),subject_no,length(channel_select),size(select_source,1));
T_i_start=1;
test_run_start=1;
sub_no_start=1;


disp('Decomposition')

for source_i=source_input
    for T_i=T_i_start:length(possible_T)
            T=possible_T(T_i);
            for test_run=test_run_start:size(test_block,1)
                for sub_no=sub_no_start:subject_no
                    sub=['S' num2str(sub_no)];

                    for ch_i=1:length(channel_select)

                        x_diff_stim=[];
                        n_period=[];
                        t=[];
                        t_1period=[];
                        for trial=1:trial_no
                            x_diff_stim(trial,:,:)=averaged_ref{T_i,test_run,sub_no,trial}(ch_i,:);
                            n_period(1,trial)=T*freqs(trial);
                            total_sample=size(x_diff_stim,3);
                            sample_one_period=total_sample/n_period(1,trial);
                            t(trial,:)=0:((1/fs)*(2*pi)/(1/freqs(trial))):((1/fs)*(2*pi)/(1/freqs(trial))*(total_sample-1));%linspace(0,2*pi*n_period(1,trial),size(x_diff_stim,3));
                        end
                        x_diff_stim=squeeze(x_diff_stim);
                        x_recon=zeros(size(x_diff_stim));
                        S1_store_tmp=[];
                        t_re_store=[];
                        for freqs_i=1:2
                            if freqs_i==1
                                freqs_select=select_source(source_i,:);
                            else
                                freqs_select=setdiff(1:length(freqs),select_source(source_i,:));
                            end
                            x_diff_stim_tmp=x_diff_stim(freqs_select,:);
                            n_period_tmp=n_period(freqs_select);
                            t_tmp=t(freqs_select,:);

                            t_re=[];
                            for i=1:length(freqs_select)
                                t_re(i,:)=t(freqs_select(i),:)+phases(freqs_select(i));
                            end

                            hilbert_x=hilbert(x_diff_stim_tmp.').';
                            if freqs_i==1
                                afdcal=AFDCal(hilbert_x);
                                afdcal.setDecompMethod(3);
                                for t_re_ch=1:size(t_re,1)
                                    afdcal.setPhase(t_re_ch,t_re(t_re_ch,:));
                                end
                                afdcal.genDic(0.1,0.95);
                                afdcal.genEva();
                                afdcal.init_decomp(1);
                                for n=1:AFD_N-1
                                    afdcal.nextDecomp(1);
                                end
                                %
                                F=[];
                                S1=[];
                                for n=1:AFD_N
                                    for ch_afd=1:size(hilbert_x,1)
                                        F(n,ch_afd,:)=afdcal.deComp{ch_afd}(n,:);
                                        S1(n,ch_afd)=afdcal.coef{ch_afd}(n);
                                    end
                                end
                                a=afdcal.an{1}.';
                                level=AFD_N;
                            else
                                afdcal=AFDCal(hilbert_x);
                                afdcal.setDecompMethod(3);
                                for t_re_ch=1:size(t_re,1)
                                    afdcal.setPhase(t_re_ch,t_re(t_re_ch,:));
                                end
                                afdcal.set_an({a.'});
                                afdcal.init_decomp(0);
                                for n=1:level-1
                                    afdcal.nextDecomp(0);
                                end
                                %
                                F=[];
                                S1=[];
                                for n=1:AFD_N
                                    for ch_afd=1:size(hilbert_x,1)
                                        F(n,ch_afd,:)=afdcal.deComp{ch_afd}(n,:);
                                        S1(n,ch_afd)=afdcal.coef{ch_afd}(n);
                                    end
                                end
                                a=afdcal.an{1}.';
                            end


                            x_recon(freqs_select,:)=real(squeeze(sum(F(1:level,:,:),1)));
                            S1_store_tmp(1:level,freqs_select)=S1(1:level,:);
                            t_re_store(freqs_select,:)=t_re;
                        end
                        S1_store{T_i,test_run,sub_no,ch_i,source_i}=S1_store_tmp;
                        x_recon_store{T_i,test_run,sub_no,ch_i,source_i}=x_recon;
                        AFD_r(T_i,test_run,sub_no,ch_i,source_i)=corr2(x_diff_stim,x_recon);
                        x_ori_store{T_i,test_run,sub_no,ch_i,source_i}=x_diff_stim;
                        a_store{T_i,test_run,sub_no,ch_i,source_i}=a(1:level);
                        t_store{T_i,test_run,sub_no,ch_i,source_i}=t_re_store;
                        
                        disp([num2str(source_i) '/' num2str(size(select_source,1)) ' AFD -> ' sub ', T=' num2str(T) ', run ' num2str(test_run) ', ch ' num2str(ch_i)])
                    end




                end


            end


    end


     save(['Averaged_template_AFD_decompose_t_move_N' num2str(AFD_N) '_random' num2str(source_i) '_of_' num2str(rand_total_no) '.mat'],...
         'AFD_N','AFD_r',...
         'S1_store','x_recon_store',...
         'x_ori_store','a_store','t_store',...
         'T_i','test_run','sub_no','select_source','source_i');
 end