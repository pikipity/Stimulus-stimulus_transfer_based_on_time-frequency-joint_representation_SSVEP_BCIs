function test_impulse_response(data_dir,source_no,source_i)

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


new_u_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
new_v_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
r_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1));
res_store=zeros(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1),2);
eeg_ref_new=cell(length(possible_T),size(test_block,1),subject_no,trial_no);


load(['study_impulse_response_2_source' num2str(source_i) '_of_' num2str(source_no) '.mat'],...
    'temlate_re')
load(['random_source_' num2str(source_no) '.mat'],'select_source')
index1=1;

AFD_r=zeros(length(possible_T),size(test_block,1),subject_no);
for T_i=1:length(possible_T)
        T=possible_T(T_i);
        for test_run=1:size(test_block,1)
            for sub_no=1:subject_no
                sub=['S' num2str(sub_no)];
                y=subj_data{sub_no};

                for freqs_i=1:40
                    eeg_ref_new{end,test_run,sub_no,freqs_i}=temlate_re{end,test_run,sub_no}(freqs_i,:);
                end
                
                freqs_select_target=select_source(source_i,:);
                for freqs_target_i=1:length(freqs_select_target)
                    freqs_select=[];
                    for freqs_i=1:40
                        [~,target_i]=min(abs(freqs(freqs_select_target)-freqs(freqs_i)));
                        if target_i==freqs_target_i
                            freqs_select=[freqs_select freqs_i];
                        end
                    end
                    
                    freqs_select=[freqs_select_target(freqs_target_i) setdiff(freqs_select,freqs_select_target(freqs_target_i))];
                    for freqs_i=1:length(freqs_select)
                        if freqs_i==1
                            x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),freqs_select(freqs_i),train_block(test_run,:)));
                            x=mean(x,3);
                            [u,v,r]=canoncorr(x.',eeg_ref_new{end,test_run,sub_no,freqs_select(freqs_i)}(1:index1,1:size(x,2)).');
                            new_u_store{T_i,test_run,sub_no,freqs_select(freqs_i),1}=u(:,1);
                            new_v_store{T_i,test_run,sub_no,freqs_select(freqs_i),1}=v(:,1);
                        else
                            new_u_store{T_i,test_run,sub_no,freqs_select(freqs_i),1}=new_u_store{T_i,test_run,sub_no,freqs_select(1),1};
                            new_v_store{T_i,test_run,sub_no,freqs_select(freqs_i),1}=new_v_store{T_i,test_run,sub_no,freqs_select(1),1};
                        end
                    end
                end


                disp(['Impulse Response Template Generation ' num2str(source_i) '-' num2str(source_no) ': ' sub ', T=' num2str(T) ', run ' num2str(test_run) ', r=' num2str(AFD_r(T_i,test_run,sub_no))])

                
                
            end
        end
end
u_store=new_u_store;
v_store=new_v_store;
eeg_ref=eeg_ref_new;


for index_i=1
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
                        for i=1:length(freqs)
                            disp(['Test impulse response ' num2str(source_i) '-' num2str(source_no) ', Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial) ', block: ' num2str(test_block(test_run,block)) ', f: ' num2str(i) ', index: ' num2str(index1)])
            
                            [u,v,r]=canoncorr((u_store{T_i,test_run,sub_no,i,1}(:,1:index1).'*x).',eeg_ref{end,test_run,sub_no,i}(1:index1,1:size(x,2)).');
                            [u2,v2,r2]=canoncorr((u_store{T_i,test_run,sub_no,i,1}(:,1:index1).'*x).',sine_ref{i}(:,1:size(x,2)).');
                            [u4,v4,r4]=canoncorr(x.',sine_ref{i}(:,1:size(x,2)).');
                            R(1,i)=sign(r(1))*r(1)^2+sign(r2(1))*r2(1)^2;
                            R(2,i)=sign(r(1))*r(1)^2+sign(r2(1))*r2(1)^2+sign(r4(1))*r4(1)^2;
                            

                        end
                        r_store{T_i,test_run,sub_no,trial,block,index_i}=R;
                        [~,I_R]=max(R(1,:));
                        if I_R==trial
                            res_store(T_i,test_run,sub_no,trial,block,index_i,1)=1;
                        end
                        [~,I_R]=max(R(2,:));
                        if I_R==trial
                            res_store(T_i,test_run,sub_no,trial,block,index_i,2)=1;
                        end
                        
                    end
                end
            end
        end
    end
end

save(['trca_impulse_template_neighbor_results_random' num2str(source_i) '_of_' num2str(source_no) '_end_good_2.mat'],...
      'u_store','v_store','r_store','res_store','eeg_ref')