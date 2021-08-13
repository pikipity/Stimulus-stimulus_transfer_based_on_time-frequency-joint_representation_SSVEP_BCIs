clear;clc;
data_dir='./data';
figdir='./ITR_Acc_Summary';

% CCA
test_cca(data_dir);

% ECCA
test_ecca(data_dir);

% TRCA
test_trca(data_dir);

% impulse reponse initial
impulse_response_decompose(data_dir);

for source_no=[4 8 12]
    for d=1:5
        % phase-domain based stimulus-stimulus transfer
        AFD_decomposition_averaged_template(data_dir,source_no,d);
        study_AFD_coef_averaged_template(data_dir,source_no,d);
        test_AFD_template(data_dir,source_no,d);

        % impulse response
        study_impulse_response(data_dir,source_no,d);
        test_impulse_response(data_dir,source_no,d);
    end
end

% Plot Acc and ITR
% require "shadedErrorBar" and "sigstar"
summary_ITR_Acc(figdir,data_dir);