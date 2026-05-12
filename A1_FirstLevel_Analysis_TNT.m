% =========================================================================
% First-Level fMRI Analysis for Think/No-Think (TNT) Task
% =========================================================================
% Description:
%   This script performs first-level General Linear Model (GLM) specification,
%   estimation, and contrast generation using SPM12. It is optimized for 
%   large-scale batch processing utilizing MATLAB's Parallel Computing Toolbox.
%
% Toolboxes Required:
%   - SPM12 (Statistical Parametric Mapping)
%
% Note on Design Matrix:
%   Based on the specific pipeline configuration, the constant/mean term is 
%   located at the FIRST column of the design matrix. Contrast vectors are 
%   padded with a leading 0 accordingly.
% =========================================================================

%% 1. Initialization & Environment Setup
clear; clc;

% Initialize parallel pool if not already running to accelerate processing
if isempty(gcp('nocreate'))
    parpool(12);
end

% Define directories (Ensure these match your BIDS/derivatives structure)
base_path = '/Documents/Data/Resilience_fmri_pku/resilience_task_3mm';
preproc_dir = fullfile(base_path, 'BIDS_preprocessed');
output_dir  = fullfile(base_path, 'First_Level');
mask_file   = '/Documents/MATLAB/Mask/MNI152_T1_3mm_brain_mask.nii';

% task_id = 'MID'; % Note: Consider changing to 'TNT' if this is a Think/No-Think task
task_id = 'TNT'; 

% Scanning parameters
fmri_rt    = 2;     % Repetition Time (TR) in seconds
fmri_t     = 16;    % Microtime resolution
fmri_t0    = 8;     % Microtime onset (reference slice for slice timing)
fmri_units = 'secs';

%% 2. Subject Discovery
% Locate all subject directories
sub_dirs = dir(fullfile(preproc_dir, 'Sub_*'));
sub_dirs = sub_dirs([sub_dirs.isdir]); % Filter out non-directories
num_subs = length(sub_dirs);

% Initialize tracking variables
wrong_subject = zeros(1, num_subs);
sub_job = cell(1, num_subs);

%% 3. GLM Specification & Batch Construction
fprintf('--- Beginning First-Level GLM Specification ---\n');

for j = 1:num_subs
    try
        sub_name = sub_dirs(j).name;
        fprintf('Preparing batch for Subject: %s\n', sub_name);
        
        % Locate necessary files using robust cross-platform fullfile
        onset_file = dir(fullfile(preproc_dir, sub_name, '*run*.mat'));
        motion_file = dir(fullfile(preproc_dir, sub_name, '*rp_aSub*.txt'));
        nii_file = dir(fullfile(preproc_dir, sub_name, 'swraSub*.nii'));
        
        % Validate that all files exist and match in run numbers
        if ~isempty(nii_file) && (length(nii_file) == length(motion_file)) && (length(nii_file) == length(onset_file))
            
            clear matlabbatch;
            sub_out_dir = fullfile(output_dir, sub_name);
            
            % Create output directory if it does not exist
            if ~exist(sub_out_dir, 'dir')
                mkdir(sub_out_dir);
            end
            
            % 3.1 Model Specification Parameters
            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(sub_out_dir);
            matlabbatch{1}.spm.stats.fmri_spec.timing.units   = fmri_units;
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = fmri_rt;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = fmri_t;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;
            
            % Loop through runs/sessions
            for i = 1:length(nii_file)
                % NIfTI functional images
                niifile_path = spm_select('ExtFPList', nii_file(i).folder, nii_file(i).name);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = cellstr(niifile_path);
                
                % Task Onsets
                onset_data = load(fullfile(onset_file(i).folder, onset_file(i).name));
                conditions = table2struct(onset_data.onset);
                
                for k = 1:length(conditions)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).name     = char(conditions(k).condition);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).tmod     = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).onset    = round(conditions(k).onset, 10);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).duration = conditions(k).duration;
                end
                
                % Nuisance Regressors (Head Motion)
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
                motion_filepath = spm_select('FPList', motion_file(i).folder, ['^', motion_file(i).name, '$']);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = cellstr(motion_filepath);
                
                % High-pass filter
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
            end
            
            % 3.2 Global Defaults & Masking
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {mask_file};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            
            % 3.3 Model Estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            % 3.4 Contrast Configuration
            % Note: Weights start with 0 because the first column in the design 
            % matrix represents the constant/mean term per the user's specific setup.
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name    = 'NT-T';
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 -1 1]; 
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = 'T-NT';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 -1];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
            
            matlabbatch{3}.spm.stats.con.delete = 0;
            
            % Store the job
            sub_job{j} = matlabbatch;
        else
            warning('Missing or mismatched files for Subject: %s', sub_name);
            wrong_subject(j) = j;
        end
    catch ME
        warning('Error configuring batch for Subject %s: %s', sub_dirs(j).name, ME.message);
        wrong_subject(j) = j;
    end
end

% Remove invalid subjects from execution queue
sub_dirs(wrong_subject ~= 0) = [];
sub_job(wrong_subject ~= 0) = [];

%% 4. Parallel Model Execution
fprintf('\n--- Executing First-Level Analysis (Parallel Mode) ---\n');
wrong_first_subject = zeros(1, length(sub_job));

parfor j = 1:length(sub_job)
    try
        spm_jobman('run', sub_job{j});
    catch
        wrong_first_subject(j) = j;
    end
end

% Identify subjects that failed during estimation
wrong_first_subject(wrong_first_subject == 0) = [];
failed_subs = sub_dirs(wrong_first_subject);

% Clean up directories for failed subjects using cross-platform rmdir
parfor i = 1:length(failed_subs)
    fprintf('Cleaning up failed estimation for Subject: %s\n', failed_subs(i).name);
    failed_dir = fullfile(output_dir, failed_subs(i).name);
    if exist(failed_dir, 'dir')
        rmdir(failed_dir, 's'); 
    end
end

%% 5. Output Renaming (Cross-Platform)
% Renaming con_*.nii and spmT_*.nii to explicitly include condition names
fprintf('\n--- Renaming Contrast Maps ---\n');
valid_subs = dir(fullfile(output_dir, 'Sub_*'));

parfor i = 1:length(valid_subs)
    sub_dir = fullfile(output_dir, valid_subs(i).name);
    spm_mat_file = fullfile(sub_dir, 'SPM.mat');
    
    if exist(spm_mat_file, 'file')
        sub_onset = load(spm_mat_file);
        conditions = {sub_onset.SPM.xCon.name};
        
        for j = 1:length(conditions)
            % Construct original and target file names
            orig_con = fullfile(sub_dir, sprintf('con_%04d.nii', j));
            targ_con = fullfile(sub_dir, sprintf('con_%s.nii', conditions{j}));
            
            orig_spmT = fullfile(sub_dir, sprintf('spmT_%04d.nii', j));
            targ_spmT = fullfile(sub_dir, sprintf('spmT_%s.nii', conditions{j}));
            
            % Use MATLAB native movefile for OS independence (Windows/Linux/macOS)
            if exist(orig_con, 'file')
                movefile(orig_con, targ_con);
            end
            if exist(orig_spmT, 'file')
                movefile(orig_spmT, targ_spmT);
            end
        end
        fprintf('Renamed maps for Subject: %s\n', valid_subs(i).name);
    end
end

fprintf('\n--- First-Level Analysis Completed ---\n');