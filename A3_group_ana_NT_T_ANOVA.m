% =========================================================================
% Script Name: A3_group_ana_NT_T_ANOVA.m
% Description: 
%   Performs a Second-Level One-Way ANOVA in SPM12 to compare the NT-T 
%   contrast across three independent groups: Healthy Control, Low Resilience, 
%   and High Resilience. Age and Sex are included as mean-centered covariates.
%   An explicit ROI mask is applied to restrict the statistical search volume.
%
% Toolboxes Required:
%   - SPM12
% =========================================================================

%% 1. Initialization & Environment Setup
clear; clc;

% Initialize SPM job manager
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define base directories
base_dir       = '/public/home/yangzy/Documents/Data/Resilience_fmri_pku';
firstLevel_dir = fullfile(base_dir, 'resilience_task_3mm', 'First_Level');
behavior_dir   = fullfile(base_dir, 'resilience_task', 'behaviordata_sub');
roi_mask       = fullfile(base_dir, 'resilience_task_3mm', 'ROI', 'NT_T', 'Cluster_fwe_30', 'NT_T_cluster10_all.nii');

% Define Group Image Directories
dir_health = fullfile(base_dir, 'Resilience Image data (copy 1)', 'healthy control', 'healthy control');
dir_low    = fullfile(base_dir, 'Resilience Image data (copy 1)', 'low resilience', 'low resilience');
dir_high   = fullfile(base_dir, 'Resilience Image data (copy 1)', 'high resilience', 'high resilience');

% Output directory for Group Analysis
group_analysisDir = fullfile(base_dir, 'resilience_task_3mm', 'Second_Level', 'F_test', 'NT-T');
if ~exist(group_analysisDir, 'dir')
    mkdir(group_analysisDir);
end

% IDs to exclude (Data problem subjects & specific 3mm exclusions)
excludeIDs = {'sub-044', 'sub-001', 'sub-015', 'sub-027', 'sub-013', ...
              'sub-016', 'sub-017', 'sub-045', 'sub-006', ...
              'sub-071', 'sub-074', 'sub-075', 'sub-082'};

%% 2. Data Extraction
fprintf('--- Extracting data for Healthy Controls ---\n');
[scans_health_t1, sex_health, age_health] = get_subject_paths_and_info(dir_health, excludeIDs, firstLevel_dir, behavior_dir);

fprintf('--- Extracting data for Low Resilience ---\n');
[scans_low_t1, sex_low, age_low] = get_subject_paths_and_info(dir_low, excludeIDs, firstLevel_dir, behavior_dir);

fprintf('--- Extracting data for High Resilience ---\n');
[scans_high_t1, sex_high, age_high] = get_subject_paths_and_info(dir_high, excludeIDs, firstLevel_dir, behavior_dir);

%% 3. Prepare Covariates
% Concatenate covariates for the entire model
sex_vector = [sex_health(:); sex_low(:); sex_high(:)];
age_vector = [age_health(:); age_low(:); age_high(:)];

fprintf('Total subjects: Health (N=%d), Low (N=%d), High (N=%d)\n', ...
    length(scans_health_t1), length(scans_low_t1), length(scans_high_t1));

%% 4. Execute SPM Batch
if ~isempty(scans_health_t1) && ~isempty(scans_low_t1) && ~isempty(scans_high_t1)
    fprintf('--- Building and Running SPM ANOVA Batch ---\n');
    matlabbatch = build_anova_batch(group_analysisDir, scans_health_t1, scans_low_t1, scans_high_t1, sex_vector, age_vector, roi_mask);
    spm_jobman('run', matlabbatch);
    fprintf('--- Second-Level ANOVA Completed ---\n');
else
    error('One or more groups have no valid subjects. Aborting analysis.');
end

% =========================================================================
% Local Functions
% =========================================================================

function matlabbatch = build_anova_batch(groupfolder, scan1, scan2, scan3, sex, age, roi_mask)
    % Model Specification: One-Way ANOVA (3 levels)
    matlabbatch{1}.spm.stats.factorial_design.dir = {groupfolder};
    
    % Group Cells
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = scan1(:);
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = scan2(:);
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = scan3(:);
    
    % ANOVA Parameters
    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;     % Independent groups
    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % Unequal variance
    matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
    
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
    
    % Explicit Mask
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {[roi_mask, ',1']};
    
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    % Model Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % Contrast Manager
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    % F-Test: Main effect of group
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'NT-T F test';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 -1 0; 1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'replsc';
    
    % T-Tests: Pairwise group differences
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NT-T T_health_low';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'NT-T T_low_health';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'NT-T T_high_low';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'NT-T T_low_high';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'NT-T T_health_high';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'NT-T T_high_health';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'NT-T T_health_trauma';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1 -1 -1];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'NT-T T_trauma_health';
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-1 1 1];
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'replsc';
    
    matlabbatch{3}.spm.stats.con.delete = 0;
end

function [valid_scans, valid_sex, valid_age] = get_subject_paths_and_info(parentFolder, excludeIDs, firstLevel_dir, behavior_dir)
    % Robustly extracts subject paths and behavioral data, bypassing missing files
    
    folders = dir(parentFolder);
    folderNames = {folders([folders.isdir]).name};
    subjectFolders = folderNames(~ismember(folderNames, {'.', '..'}));
    
    % Extract subject IDs (sub-XXX)
    subjectIDs = regexp(subjectFolders, '(\d+)', 'match');
    subjectIDs = cellfun(@(x) sprintf('sub-%03d', str2double(x{1})), subjectIDs, 'UniformOutput', false);
    
    % Clean list: Unique and exclude bad IDs
    subjectIDs = unique(subjectIDs, 'stable');
    subjectIDs = subjectIDs(~ismember(subjectIDs, excludeIDs));
    
    num_subs = length(subjectIDs);
    temp_scans = cell(num_subs, 1);
    temp_sex   = zeros(num_subs, 1);
    temp_age   = zeros(num_subs, 1);
    valid_idx  = true(num_subs, 1);
    
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
    
    valid_scans = temp_scans(valid_idx);
    valid_sex   = temp_sex(valid_idx);
    valid_age   = temp_age(valid_idx);
end