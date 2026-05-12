% =========================================================================
% Script Name: A2_SecondLevel_NT_T.m
% Description: 
%   Performs a Second-Level (Group) Single-Sample T-test in SPM12.
%   This script aggregates First-Level contrast maps (NT-T) across three 
%   groups (Healthy Control, Low Resilience, High Resilience) while 
%   covarying for Age and Sex.
%
% Note on Contrast: 
%   The contrast is set to [1], explicitly targeting the first column 
%   (the Constant/Mean effect) of the design matrix, as configured in the 
%   user's specific SPM environment.
% =========================================================================

%% 1. Initialization
clear; clc;

% Initialize SPM to ensure jobman is ready
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% 2. Define Directories & Parameters
% Base directory standardized as requested
base_dir = '/Documents/Data/Resilience_fmri_pku';

% Sub-directories
dir_health     = fullfile(base_dir, 'Resilience Image data (copy 1)', 'healthy control', 'healthy control');
dir_low        = fullfile(base_dir, 'Resilience Image data (copy 1)', 'low resilience', 'low resilience');
dir_high       = fullfile(base_dir, 'Resilience Image data (copy 1)', 'high resilience', 'high resilience');
firstLevel_dir = fullfile(base_dir, 'resilience_task_3mm', 'First_Level');
behavior_dir   = fullfile(base_dir, 'resilience_task', 'behaviordata_sub');

% Output directory for Group Analysis
group_analysisDir = fullfile(base_dir, 'resilience_task_3mm', 'Second_Level', 'single_sample_test_NT_contrast_T', 'all');

if ~exist(group_analysisDir, 'dir')
    mkdir(group_analysisDir);
end

% Data problem subjects to exclude
excludeIDs = {'sub-044','sub-001','sub-015','sub-027','sub-013','sub-016','sub-017','sub-045','sub-006'};

%% 3. Retrieve Data for All Groups
fprintf('--- Extracting data for Healthy Controls ---\n');
[scans_health_t1, sex_health, age_health] = get_subject_paths_and_info(dir_health, excludeIDs, firstLevel_dir, behavior_dir);

fprintf('--- Extracting data for Low Resilience ---\n');
[scans_low_t1, sex_low, age_low] = get_subject_paths_and_info(dir_low, excludeIDs, firstLevel_dir, behavior_dir);

fprintf('--- Extracting data for High Resilience ---\n');
[scans_high_t1, sex_high, age_high] = get_subject_paths_and_info(dir_high, excludeIDs, firstLevel_dir, behavior_dir);

%% 4. Concatenate Data
% Ensure all vectors are column vectors to avoid dimension mismatch in SPM
all_scans_t1 = [scans_health_t1(:); scans_low_t1(:); scans_high_t1(:)];
sex_vector   = [sex_health(:); sex_low(:); sex_high(:)];
age_vector   = [age_health(:); age_low(:); age_high(:)];

fprintf('Total valid subjects included in analysis: %d\n', length(all_scans_t1));

%% 5. Execute SPM Batch
if ~isempty(all_scans_t1)
    fprintf('--- Building and Running SPM Second-Level Batch ---\n');
    matlabbatch = build_group_analysis_batch(group_analysisDir, all_scans_t1, sex_vector, age_vector);
    spm_jobman('run', matlabbatch);
    fprintf('--- Second-Level Analysis Completed ---\n');
else
    warning('No valid scans found. SPM batch aborted.');
end


% =========================================================================
% Local Functions
% =========================================================================

function matlabbatch = build_group_analysis_batch(groupfolder, scans, sex, age)
    % Builds the SPM matlabbatch structure for a Single-Sample T-test
    
    % Model Specification
    matlabbatch{1}.spm.stats.factorial_design.dir = {groupfolder};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans;
    
    % Covariate 1: Age
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1; % Mean centering
    
    % Covariate 2: Sex
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = sex;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1; % Mean centering
    
    % Masking and Global parameters
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Model Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % Contrast Manager
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'NT-T (Group Mean)';
    % Weight set to 1 explicitly targets the first column (Mean), ignoring covariates
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1; 
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
    matlabbatch{3}.spm.stats.con.delete = 0;
end

function [valid_scans, valid_sex, valid_age] = get_subject_paths_and_info(parentFolder, excludeIDs, firstLevel_dir, behavior_dir)
    % Extracts subject paths and behavioral data, implementing robust error
    % checking to filter out subjects with missing NIfTI or TSV files.
    
    folders = dir(parentFolder);
    folderNames = {folders([folders.isdir]).name};
    subjectFolders = folderNames(~ismember(folderNames, {'.', '..'}));
    
    % Extract subject IDs (sub-XXX)
    subjectIDs = regexp(subjectFolders, '(\d+)', 'match');
    subjectIDs = cellfun(@(x) sprintf('sub-%03d', str2double(x{1})), subjectIDs, 'UniformOutput', false);
    
    % Clean list: Unique and exclude bad IDs
    subjectIDs = unique(subjectIDs, 'stable');
    subjectIDs = subjectIDs(~ismember(subjectIDs, excludeIDs));
    
    % Preallocate output vectors (will be pruned of invalid data later)
    num_subs = length(subjectIDs);
    temp_scans = cell(num_subs, 1);
    temp_sex   = zeros(num_subs, 1);
    temp_age   = zeros(num_subs, 1);
    valid_idx  = true(num_subs, 1); % Logical index to track valid subjects
    
    for i = 1:num_subs
        try
            sub_id_bids = subjectIDs{i};
            sub_id_spm  = strrep(sub_id_bids, 'sub-', 'Sub_');
            
            % 1. Verify First-Level Contrast Map exists
            scan_file = fullfile(firstLevel_dir, sub_id_spm, 'con_0001.nii');
            if ~exist(scan_file, 'file')
                warning('Missing contrast map for %s. Skipping.', sub_id_spm);
                valid_idx(i) = false;
                continue;
            end
            temp_scans{i} = [scan_file, ',1'];
            
            % 2. Verify Behavioral TSV exists
            tsv_file = fullfile(behavior_dir, sub_id_spm, [sub_id_spm, '_run1.tsv']);
            if ~exist(tsv_file, 'file')
                warning('Missing TSV file for %s. Skipping.', sub_id_spm);
                valid_idx(i) = false;
                continue;
            end
            
            % 3. Extract Behavioral Data
            data = readtable(tsv_file, 'FileType', 'text', 'Delimiter', '\t');
            if isempty(data) || ~ismember('Sex', data.Properties.VariableNames) || ~ismember('Age', data.Properties.VariableNames)
                warning('Invalid TSV content for %s. Skipping.', sub_id_spm);
                valid_idx(i) = false;
                continue;
            end
            
            temp_sex(i) = data.Sex(1);
            temp_age(i) = data.Age(1);
            
        catch ME
            warning('Error processing %s: %s', subjectIDs{i}, ME.message);
            valid_idx(i) = false;
        end
    end
    
    % Return only subjects with complete, valid data
    valid_scans = temp_scans(valid_idx);
    valid_sex   = temp_sex(valid_idx);
    valid_age   = temp_age(valid_idx);
end