function study_impulse_response(data_dir,source_no,source_i)

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


load(['random_source_' num2str(source_no) '.mat'],'select_source');


load('TRCA_template_impulse_decompose_2.mat',...
      'reconst_r',...
      'impulse_reponse','impulse_store',...
      'x_ori_store',...
      'T_i','test_run','sub_no');


impulse_res_corr=cell(length(possible_T),size(test_block,1),subject_no);
temlate_re=cell(length(possible_T),size(test_block,1),subject_no);
template_re_corr=cell(length(possible_T),size(test_block,1),subject_no);

for T_i=1:length(possible_T)
    T=possible_T(T_i);
    for test_run=1:size(test_block,1)
        for sub_no=1:subject_no
            disp(['Analyze impulse reponse ' num2str(source_i) '-' num2str(source_no) ' -> sub ' num2str(sub_no) ', T=' num2str(T) ', run ' num2str(test_run)])
            % impulse_res_corr
            impulse_res_corr{T_i,test_run,sub_no}=zeros(length(freqs_I),length(freqs_I));
            for ch_i=1:length(freqs_I)
                for ch_j=1:length(freqs_I)
                    L_min=min(length(impulse_reponse{T_i,test_run,sub_no}{ch_i}),length(impulse_reponse{T_i,test_run,sub_no}{ch_j}));
                    r=corrcoef(impulse_reponse{T_i,test_run,sub_no}{ch_i}(1:L_min).',impulse_reponse{T_i,test_run,sub_no}{ch_j}(1:L_min).');
                    impulse_res_corr{T_i,test_run,sub_no}(ch_i,ch_j)=r(2);
                end
            end
            % Reconstruct 
            x_ori=x_ori_store{T_i,test_run,sub_no};
            x_re=[];
            freqs_select_target=select_source(source_i,:);
            for freqs_target_i=1:length(freqs_select_target)
                freqs_select=[];
                for freqs_i=1:40
                    [~,target_i]=min(abs(freqs(freqs_select_target)-freqs(freqs_i)));
                    if target_i==freqs_target_i
                        freqs_select=[freqs_select freqs_i];
                    end
                end
                for freqs_select_i=1:length(freqs_select)
                    template_ori=x_ori_store{T_i,test_run,sub_no}(freqs_select(freqs_select_i),:);
                    impulse_res_tmp=impulse_reponse{T_i,test_run,sub_no}{freqs_select(1)};
                    impulse_tmp=impulse_store{T_i,test_run,sub_no}{freqs_select(freqs_select_i)};
                    first_impulse=impulse_tmp(1,:);
                    impulse_index=find(first_impulse~=0);
                    tmp_more=[first_impulse(impulse_index(end)) zeros(1,length((impulse_index(1)+1):(impulse_index(2)-1)))];
                    tmp_more=tmp_more(1:length(tmp_more)-(impulse_index(1)-1));
                    new_first_impulse=[tmp_more first_impulse];
                    new_impulse_tmp=zeros(size(impulse_tmp,1),size(new_first_impulse,2));
                    new_impulse_tmp(1,:)=new_first_impulse;
                    for i=2:size(new_impulse_tmp,1)
                        new_impulse_tmp(i,:)=[zeros(1,i-1) new_first_impulse(1,1:(end-(i-1)))];
                    end
                    if length(impulse_res_tmp)>=size(new_impulse_tmp,1)
                        x_re_tmp=impulse_res_tmp(1:size(new_impulse_tmp,1))*new_impulse_tmp;
                    else
                        x_re_tmp=[impulse_res_tmp zeros(1,size(new_impulse_tmp,1)-length(impulse_res_tmp))]*new_impulse_tmp;
                    end
                    x_re_tmp=x_re_tmp((length(tmp_more)+1):(length(tmp_more)+length(template_ori)));
                    x_re(freqs_select(freqs_select_i),:)=x_re_tmp;
                end
            end
            temlate_re{T_i,test_run,sub_no}=x_re;
            x_corr=zeros(1,length(freqs_I));
            for ch=1:length(freqs_I)
                r=corrcoef(x_ori(ch,:)',x_re(ch,:)');
                x_corr(ch)=r(2);
            end
            template_re_corr{T_i,test_run,sub_no}=x_corr;
        end
    end
end

save(['study_impulse_response_2_source' num2str(source_i) '_of_' num2str(source_no) '.mat'],...
     'temlate_re')