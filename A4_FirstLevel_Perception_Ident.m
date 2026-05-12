% =========================================================================
% Script Name: A4_FirstLevel_Perception_Ident.m
% Description: 
%   SPM12 First-Level Analysis for the Perception Identification Task (Run 5).
%   This pipeline models 'Type1' conditions and computes all combinations of 
%   reverse and Unprimed contrasts.
%
% Toolboxes Required: SPM12, Parallel Computing Toolbox
% =========================================================================

clear; clc;

% Initialize parallel pool for batch acceleration
if isempty(gcp('nocreate'))
    parpool(12);
end

% 1. Define Directories & Parameters
base_path = '/Documents/Data/Resilience_fmri_pku/resilience_task_run5';
preproc_dir = fullfile(base_path, 'BIDS_preprocessed');
output_dir  = fullfile(base_path, 'First_Level_Type1');
mask_file   = '/Documents/MATLAB/Mask/MNI152_T1_3mm_brain_mask.nii';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Subject Discovery
sub_dirs = dir(fullfile(preproc_dir, 'Sub_*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);
num_subs = length(sub_dirs);

% Scanning parameters
fmri_rt    = 2;
fmri_t     = 16;
fmri_t0    = 8;
fmri_units = 'secs';

% Tracking variables
wrong_subject = zeros(1, num_subs);
sub_job = cell(1, num_subs);

%% 2. GLM Specification & Contrast Definition
fprintf('--- Beginning First-Level Batch Configuration ---\n');

for j = 1:num_subs
    try
        sub_name = sub_dirs(j).name;
        fprintf('Preparing batch for Subject: %03d - %s\n', j, sub_name);
        
        % Locate run 5 files specifically for Type1 conditions
        onset_file = dir(fullfile(preproc_dir, sub_name, '*_Type1.mat'));
        motion_file = dir(fullfile(preproc_dir, sub_name, '*rp_aSub*.txt'));
        nii_file = dir(fullfile(preproc_dir, sub_name, 'swraSub*.nii'));
        
        if ~isempty(nii_file) && (length(nii_file) == length(motion_file)) && (length(nii_file) == length(onset_file))
            
            clear matlabbatch;
            sub_out_dir = fullfile(output_dir, sub_name);
            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(sub_out_dir);
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = fmri_units;
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = fmri_rt;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = fmri_t;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;
            
            for i = 1:length(nii_file)
                % NIfTI Scans
                niifile_path = spm_select('ExtFPList', nii_file(i).folder, nii_file(i).name);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = cellstr(niifile_path);
                
                % Onsets
                onset_data = load(fullfile(onset_file(i).folder, onset_file(i).name));
                conditions = table2struct(onset_data.onset);
                
                for k = 1:length(conditions)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).name = char(conditions(k).condition);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).tmod = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).onset = round(conditions(k).onset, 10);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).duration = conditions(k).duration;
                end
                
                % Nuisance Regressors: Head Motion
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
                motion_filepath = spm_select('FPList', motion_file(i).folder, ['^', motion_file(i).name, '.*']);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = cellstr(motion_filepath);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
            end
            
            % Global Settings & Explicit Mask
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {mask_file};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
            % Model Estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            % =====================================================================
            % SPM Contrast Set
            % Assumed Order: [Baseline, No_think, Think, Unprimed]
            % Note: Contrast weights are preserved exactly as defined by the user.
            % =====================================================================
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            
            % 1. Main Effects
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Baseline';
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'No_think';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Think';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Unprimed';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
            
            % 2. Positive Contrasts vs Baseline & Think/No-Think
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'No_think-Baseline';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [-1 1 0 0];
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Think-Baseline';
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [-1 0 1 0];
            matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Think-No_think';
            matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 -1 1 0];
            matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';
            
            % 3. Reverse / Negative Contrasts
            matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Baseline-No_think';
            matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1 -1 0 0];
            matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Baseline-Think';
            matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [1 0 -1 0];
            matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'No_think-Think';
            matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 1 -1 0];
            matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'replsc';

            % 4. Unprimed Positive Contrasts
            matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'Unprimed-Baseline';
            matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [-1 0 0 1];
            matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'replsc';

            matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'Unprimed-No_think';
            matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [0 -1 0 1];
            matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'replsc';

            matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'Unprimed-Think';
            matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [0 0 -1 1];
            matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'replsc';

            % 5. Unprimed Reverse Contrasts
            matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'Baseline-Unprimed';
            matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [1 0 0 -1];
            matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'replsc';

            matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'No_think-Unprimed';
            matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [0 1 0 -1];
            matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'replsc';

            matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'Think-Unprimed';
            matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [0 0 1 -1];
            matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.delete = 0;
            sub_job{j} = matlabbatch;
        else
            wrong_subject(j) = j;
        end
    catch
        wrong_subject(j) = j;
    end
end

% Purge invalid subjects from processing queue
sub_dirs(wrong_subject ~= 0) = [];
sub_job(wrong_subject ~= 0) = [];

%% 3. Parallel Execution
fprintf('\n--- Executing First-Level GLM ---\n');
wrong_first_subject = zeros(1, length(sub_job));

parfor j = 1:length(sub_job)
    try
        spm_jobman('run', sub_job{j});
    catch
        wrong_first_subject(j) = j;
    end
end

% Cleanup failed subject directories (Robust OS-independent method)
wrong_first_subject(wrong_first_subject == 0) = [];
failed_subs = sub_dirs(wrong_first_subject);

% Remove failed subjects from the master list BEFORE copying maps
sub_dirs(wrong_first_subject) = []; 

parfor i = 1:length(failed_subs)
    fprintf('Removing failed subject folder: %s\n', failed_subs(i).name);
    failed_dir = fullfile(output_dir, failed_subs(i).name);
    if exist(failed_dir, 'dir')
        rmdir(failed_dir, 's'); % Replaces unix('rm -rf ...')
    end
end

%% 4. File Management & Archiving
fprintf('\n--- Copying and Renaming Output Files ---\n');
valid_subs = dir(fullfile(output_dir, 'Sub_*'));

parfor i = 1:length(valid_subs)
    fprintf('Processing maps for: %s\n', valid_subs(i).name);
    sub_dir = fullfile(output_dir, valid_subs(i).name);
    spm_mat_file = fullfile(sub_dir, 'SPM.mat');
    
    if exist(spm_mat_file, 'file')
        sub_onset = load(spm_mat_file);
        conditions = {sub_onset.SPM.xCon.name};
        
        for j = 1:length(conditions)
            % Sanitize condition names for file system compatibility
            safe_cond_name = strrep(conditions{j}, ' ', '_');
            
            % Handle Contrast Images
            old_con = fullfile(sub_dir, sprintf('con_%04d.nii', j));
            new_con = fullfile(sub_dir, ['con_', safe_cond_name, '.nii']);
            if exist(old_con, 'file')
                copyfile(old_con, new_con);
            end
            
            % Handle T-Statistic Images
            old_t = fullfile(sub_dir, sprintf('spmT_%04d.nii', j));
            new_t = fullfile(sub_dir, ['spmT_', safe_cond_name, '.nii']);
            if exist(old_t, 'file')
                copyfile(old_t, new_t);
            end
        end
    end
end

disp('All First-Level GLM operations and file copying are complete!');