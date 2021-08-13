function test_trca(data_dir)

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

u_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
v_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
r_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1));
res_store=zeros(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1));
eeg_ref=cell(length(possible_T),size(test_block,1),subject_no,trial_no);

for T_i=1:length(possible_T)
    T=possible_T(T_i);
    for test_run=1:size(test_block,1)
        for sub_no=1:subject_no
            sub=['S' num2str(sub_no)];

            y=subj_data{sub_no};
            
            for trial=1:trial_no
                disp(['Train TRCA  ->  Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial)])
                
                x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,train_block(test_run,:)));
                [V,~,~,~] = fun_TRCA_Matrix(x);
                u_store{T_i,test_run,sub_no,trial,1}=V;
                
                eeg_ref_tmp=[];
                for block_no=1:length(train_block(test_run,:))
                    eeg_ref_tmp(:,:,block_no)=V.'*squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,train_block(test_run,block_no)));
                end
                eeg_ref{T_i,test_run,sub_no,trial}=mean(eeg_ref_tmp,3);
            end
            
        end
    end
end


for index_i=1:length(possible_index1)
    index1=possible_index1(index_i);
    for T_i=1:length(possible_T)
        T=possible_T(T_i);
        for test_run=1:size(test_block,1)
            for sub_no=1:subject_no
                sub=['S' num2str(sub_no)];


                y=subj_data{sub_no};

                for trial=1:trial_no
                    for block=1:length(test_block(test_run,:))%1:block_no
                        x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),trial,test_block(test_run,block)));

                        R=[];
                        U=[];
                        V=[];
                        for i=1:length(freqs)
                            disp(['TRCA Test -> Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial) ', block: ' num2str(test_block(test_run,block)) ', f: ' num2str(i) ', index: ' num2str(index1)])
 
                            r=corr2(u_store{T_i,test_run,sub_no,i,1}(:,1:index1).'*x,eeg_ref{T_i,test_run,sub_no,i}(1:index1,1:size(x,2)));

                            R(1,i)=r(1);
  
                        end

                        r_store{T_i,test_run,sub_no,trial,block,index_i}=R;
                        [~,I_R]=max(R);
                        if I_R==trial
                            res_store(T_i,test_run,sub_no,trial,block,index_i)=1;
                        end

                    end
                end
                
            end
        end
    end

end

save('trca_results.mat','u_store','v_store','r_store','res_store','eeg_ref')
