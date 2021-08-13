function test_AFD_template(data_dir,rand_total_no,source_input)

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
v_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2));
r_store=cell(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1));
res_store=zeros(length(possible_T),size(test_block,1),subject_no,trial_no,size(test_block,2),length(possible_index1),4);
eeg_ref_new=cell(length(possible_T),size(test_block,1),subject_no,trial_no);
S1_store=cell(length(possible_T),size(test_block,1),subject_no);


disp('Load Averaged template')
AFD_N=50;

source_i=source_input;
source_no=rand_total_no;
load(['study_S1_averaged_template_t_move_random' num2str(source_i) '_of_' num2str(source_no) '.mat'],...
     're_f_abs_same_phase_same_select')
load(['Averaged_template_AFD_decompose_t_move_N' num2str(AFD_N) '_random' num2str(source_i) '_of_' num2str(source_no) '.mat'],...
      'select_source')
index1=1;

AFD_r=zeros(length(possible_T),size(test_block,1),subject_no);
for T_i=1:length(possible_T)
        T=possible_T(T_i);
        for test_run=1:size(test_block,1)
            for sub_no=1:subject_no
                sub=['S' num2str(sub_no)];
                y=subj_data{sub_no};
                    


                for freqs_i=1:40
                    if ~ismember(freqs_i,select_source(source_i,:))
                        for ch_i=1:length(channel_select)
                            if length(size(re_f_abs_same_phase_same_select))==4
                                eeg_ref_new{T_i,test_run,sub_no,freqs_i}(ch_i,:)=re_f_abs_same_phase_same_select{end,test_run,sub_no,ch_i}(freqs_i,:);
                            else
                                eeg_ref_new{T_i,test_run,sub_no,freqs_i}(ch_i,:)=re_f_abs_same_phase_same_select{end,test_run,sub_no,ch_i,source_i}(freqs_i,:);
                            end
                        end
                    else
                        x=squeeze(y(channel_select,floor(start_t*fs):floor((start_t+T)*fs-1),freqs_i,train_block(test_run,:)));
                        x=mean(x,3);
                        eeg_ref_new{T_i,test_run,sub_no,freqs_i}=x;
                    end
                end
                
                sig_len=length(floor(start_t*fs):floor((start_t+T)*fs-1));
                
                freqs_select_target=select_source(source_i,:);
                for freqs_target_i=1:length(freqs_select_target)
                    X_sig=[];
                    Y_ref=[];
                    for freqs_i=1:40
                        [~,target_i]=min(abs(freqs(freqs_select_target)-freqs(freqs_i)));
                        if target_i==freqs_target_i
                            X_sig=[X_sig eeg_ref_new{T_i,test_run,sub_no,freqs_i}(:,1:sig_len)];
                            Y_ref=[Y_ref sine_ref{freqs_i}(:,1:sig_len)];
                        end
                    end

                    [u,v,r]=canoncorr(X_sig.',Y_ref.');
                    new_u_store{T_i,test_run,sub_no,freqs_target_i,1}=u(:,1);
                    new_v_store{T_i,test_run,sub_no,freqs_target_i,1}=v(:,1);
                end
                
                
                disp(['AFD Template Generation -> ' num2str(source_i) '/' num2str(size(select_source,1)) ' ' sub ', T=' num2str(T) ', run ' num2str(test_run)])
                
            end
        end
end
u_store=new_u_store;
v_store=new_v_store;
eeg_ref=eeg_ref_new;


for index_i=index1
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
                        %
                        R=[];
                        for i=1:length(freqs)
                            freqs_select_target=select_source(source_i,:);
                            [~,target_i]=min(abs(freqs(freqs_select_target)-freqs(i)));
                            disp(['AFD Template Test -> ' num2str(source_i) '/' num2str(size(select_source,1)) ' Sig Len: ' num2str(T) ', Run ' num2str(test_run) ', ' sub ', trial: ' num2str(trial) ', block: ' num2str(test_block(test_run,block)) ', f: ' num2str(i) ', index: ' num2str(index1)])
                            %
                            [u,v,r]=canoncorr((u_store{T_i,test_run,sub_no,target_i,1}(:,1:index1).'*x).',sine_ref{i}(:,1:size(x,2)).');
                            [u2,v2,r2]=canoncorr((u_store{T_i,test_run,sub_no,target_i,1}(:,1:index1).'*x).',(v_store{T_i,test_run,sub_no,target_i,1}(:,1:index1).'*sine_ref{i}(:,1:size(x,2))).'); 
                            [u4,v4,r4]=canoncorr(x.',sine_ref{i}(:,1:size(x,2)).');
                            R(1,i)=sign(r(1))*r(1)^2+sign(r2(1))*r2(1)^2;
                            R(2,i)=sign(r(1))*r(1)^2+sign(r2(1))*r2(1)^2+sign(r4(1))*r4(1)^2;
                            %
                            r=[];
                            [u1,v1,~]=canoncorr(x.',sine_ref{i}(:,1:size(x,2)).');
                            [u2,v2,~]=canoncorr(x.',eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)).');
                            [u3,v3,~]=canoncorr(eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)).',sine_ref{i}(:,1:size(x,2)).');
                            r_tmp=corrcoef((u1(:,1)'*x)',(v1(:,1)'*sine_ref{i}(:,1:size(x,2)))');
                            r(1)=r_tmp(2);
                            r_tmp=corrcoef((u2(:,1)'*x)',(v2(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                            r(2)=r_tmp(2);
                            r_tmp=corrcoef((u1(:,1)'*x)',(u1(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                            r(3)=r_tmp(2);
                            r_tmp=corrcoef((u3(:,1)'*x)',(u3(:,1)'*eeg_ref{T_i,test_run,sub_no,i}(:,1:size(x,2)))');
                            r(4)=r_tmp(2);
                            %
                            R(3,i)=sign(r(1))*r(1)^2+sign(r(2))*r(2)^2+sign(r(3))*r(3)^2+sign(r(4))*r(4)^2;
                            R(4,i)=R(1,i)+R(3,i);

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
                        [~,I_R]=max(R(3,:));
                        if I_R==trial
                            res_store(T_i,test_run,sub_no,trial,block,index_i,3)=1;
                        end
                        [~,I_R]=max(R(4,:));
                        if I_R==trial
                            res_store(T_i,test_run,sub_no,trial,block,index_i,4)=1;
                        end

                    end
                end

            end
        end
    end

end

save(['averaged_AFD_template_neighbor_results_N' num2str(AFD_N) '_neighber_w_random' num2str(source_i) '_of_' num2str(source_no) '_end_good1and2.mat'],...
    'u_store','v_store','r_store','res_store','eeg_ref','AFD_r')