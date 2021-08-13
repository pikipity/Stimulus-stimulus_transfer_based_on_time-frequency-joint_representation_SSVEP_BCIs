function impulse_response_decompose(data_dir)

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

possible_T=0.25:0.25:5;
possible_index1=1:9;


load('trca_results.mat','eeg_ref');
index1=1;

disp('Initial')
AFD_N=20;
reconst_r=zeros(length(possible_T),size(test_block,1),subject_no);
impulse_reponse=cell(length(possible_T),size(test_block,1),subject_no);
impulse_store=cell(length(possible_T),size(test_block,1),subject_no);
x_ori_store=cell(length(possible_T),size(test_block,1),subject_no);
T_i_start=1;
test_run_start=1;
sub_no_start=1;


for T_i=T_i_start:length(possible_T)
        for test_run=test_run_start:size(test_block,1)
            for sub_no=sub_no_start:subject_no
                sub=['S' num2str(sub_no)];
                    
                x_diff_stim=[];
                T=possible_T(T_i);
                n_period=[];
                t=[];
                for trial=1:trial_no
                    x_diff_stim(trial,:,:)=eeg_ref{T_i,test_run,sub_no,trial}(1:index1,:);
                    n_period(1,trial)=T*freqs(trial);
                    t(trial,:)=linspace(0,2*pi*n_period(1,trial),size(x_diff_stim,3))+phases(trial);
                end
                
                x_diff_stim=squeeze(x_diff_stim);
                x_recon=zeros(size(x_diff_stim));
                impulse_res_tmp={};
                impulse_tmp={};
                
                freqs_select=[];
                freqs_count=zeros(1,8);
                for i=1:40
                    [~,target_i]=min(abs(freqs(1:8)-freqs(i)));
                    freqs_count(target_i)=freqs_count(target_i)+1;
                    freqs_select(target_i,freqs_count(target_i))=i;
                end
                x_re=[];
                for i=1:size(freqs_select,1)
                    for j=1:size(freqs_select,2)
                        if freqs_select(i,j)==0
                            continue
                        end
                        x_ori=x_diff_stim(freqs_select(i,j),:);
                        impulse_cosine=cos(t(freqs_select(i,j),:));
                        impulse=zeros(size(impulse_cosine));
                        d=1;
                        [~,start_d]=min(abs(t(freqs_select(i,j),:)-(max(phases)+0.25*pi+2*pi*(d-1))));
                        while start_d<length(impulse)
                            impulse(start_d)=impulse_cosine(start_d);
                            d=d+1;
                            [~,start_d]=min(abs(t(freqs_select(i,j),:)-(max(phases)+0.25*pi+2*pi*(d-1))));
                        end
                            
                        impulse_matrix=[];
                        total_row=ceil(((2/freqs(freqs_select(i,j)))*fs));
                        for d=1:total_row
                            impulse_matrix(d,:)=[zeros(1,d-1) impulse zeros(1,total_row-(d-1))];
                        end
                        impulse_res=inv(impulse_matrix*impulse_matrix')*impulse_matrix*[x_ori zeros(1,size(impulse_matrix,2)-length(x_ori))]';
                        x_re_tmp=impulse_res'*impulse_matrix;
                        x_re(freqs_select(i,j),:)=x_re_tmp(1:length(x_ori));
                        impulse_tmp{1,freqs_select(i,j)}=impulse_matrix;
                        impulse_res_tmp{1,freqs_select(i,j)}=impulse_res';
                    end
                end
                impulse_reponse{T_i,test_run,sub_no}=impulse_res_tmp;
                x_ori_store{T_i,test_run,sub_no}=x_diff_stim;
                impulse_store{T_i,test_run,sub_no}=impulse_tmp;
                reconst_r(T_i,test_run,sub_no)=corr2(x_diff_stim,x_re);

                disp(['impulse response decomposition -> ' sub ', T=' num2str(T) ', run ' num2str(test_run)])
                
            end
        end
end

save(['TRCA_template_impulse_decompose_2.mat'],...
      'reconst_r',...
      'impulse_reponse','impulse_store',...
      'x_ori_store',...
      'T_i','test_run','sub_no');