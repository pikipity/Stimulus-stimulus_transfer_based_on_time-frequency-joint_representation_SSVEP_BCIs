function test_ecca(data_dir)

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

u_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
v_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
r_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
R_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
res_store=zeros(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
eeg_ref=cell(length(possible_T),size(test_block,1),subject_no,trial_no);

for T_i=1:length(possible_T)
    T=possible_T(T_i);
    for test_run=1:size(test_block,1)
        for sub_no=1:subject_no
            sub=['S' num2str(sub_no)];
            y=subj_data{sub_no};
            
            
            for trial=1:trial_no
                disp(['Train ECCA  ->  Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial)])
                
                x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,train_block(test_run,:)));
                eeg_ref{T_i,test_run,sub_no,trial}=mean(x,3);
            end
            
        end
    end
end



for T_i=1:length(possible_T)
    T=possible_T(T_i);
    for test_run=1:size(test_block,1)
        for sub_no=1:subject_no
            sub=['S' num2str(sub_no)];
            y=subj_data{sub_no};

            for trial=1:trial_no
                for block=1:length(test_block(test_run,:))
                    x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,test_block(test_run,block)));
                    
                    R=[];
                    U={};
                    V={};
                    r=[];
                    for i=1:length(freqs)
                        disp(['ECCA Test -> Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial) ', block: ' num2str(test_block(test_run,block)) ', f: ' num2str(i)])

                        [u1,v1,~]=canoncorr(x.',sine_ref{i}(:,1:size(x,2)).');
                        [u2,v2,~]=canoncorr(x.',eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)).');
                        [u3,v3,~]=canoncorr(eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)).',sine_ref{i}(:,1:size(x,2)).');
                        
                        r_tmp=corrcoef((u1(:,1)'*x)',(v1(:,1)'*sine_ref{i}(:,1:size(x,2)))');
                        r(i,1)=r_tmp(2);
                        r_tmp=corrcoef((u2(:,1)'*x)',(v2(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                        r(i,2)=r_tmp(2);
                        r_tmp=corrcoef((u1(:,1)'*x)',(u1(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                        r(i,3)=r_tmp(2);
                        r_tmp=corrcoef((u3(:,1)'*x)',(u3(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                        r(i,4)=r_tmp(2);

                        R(1,i)=sign(r(i,1))*r(i,1)^2+sign(r(i,2))*r(i,2)^2+sign(r(i,3))*r(i,3)^2+sign(r(i,4))*r(i,4)^2;
                        U{1,i}=u1(:,1);
                        V{1,i}=v1(:,1);
                        U{2,i}=u2(:,1);
                        V{2,i}=v2(:,1);
                        U{3,i}=u3(:,1);
                        V{3,i}=v3(:,1);

                    end
                    u_store{T_i,test_run,sub_no,trial,block}=U;
                    v_store{T_i,test_run,sub_no,trial,block}=V;
                    r_store{T_i,test_run,sub_no,trial,block}=r;
                    R_store{T_i,test_run,sub_no,trial,block}=R;
                    [~,I_R]=max(R);
                    if I_R==trial
                        res_store(T_i,test_run,sub_no,trial,block)=1;
                    end

                end
            end

        end
    end
end


save('ecca_results_all.mat','u_store','v_store','r_store','res_store','eeg_ref')