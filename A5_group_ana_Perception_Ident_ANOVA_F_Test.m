% =========================================================================
% Script Name: A5_group_ana_Perception_Ident_ANOVA_F_Test.m
% Description: 
%   Second-Level One-Way ANOVA (3 independent groups) across 16 First-Level 
%   Contrasts for Run 5. It automates F-tests for main effects of group 
%   and post-hoc T-tests for pairwise/composite group comparisons.
%
% Toolboxes Required:
%   - SPM12
% =========================================================================

clear; clc;

% Initialize SPM to ensure jobman is ready
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% 1. Configuration & Parameters
% IDs to exclude (due to data quality issues, etc.)
excludeIDs = {'sub-023'};

% Define the 16 First-Level Contrasts to iterate over
contrasts_to_run = {
    % 1. Main Effects
    'Baseline'
    'No_think'
    'Think'
    'Unprimed'
    
    % 2. Positive Contrasts
    'No_think-Baseline'
    'Think-Baseline'
    'Think-No_think'
    
    % 3. Reverse Contrasts
    'Baseline-No_think'
    'Baseline-Think'    
    'No_think-Think'
    
    % 4. Unprimed Positive Contrasts
    'Unprimed-Baseline'
    'Unprimed-No_think'
    'Unprimed-Think'
    
    % 5. Unprimed Reverse Contrasts
    'Baseline-Unprimed'
    'No_think-Unprimed'
    'Think-Unprimed'
};

% Define Root Directories
second_level_root = '/public/home/yangzy/Documents/Data/Resilience_fmri_pku/resilience_task_run5/Second_Level/ANOVA_3Groups';

parentFolder_health = '/home1/yangzy/Documents/Data/Resilience_fmri_pku/Resilience Image data (copy 1)/healthy control/healthy control/';
parentFolder_low    = '/home1/yangzy/Documents/Data/Resilience_fmri_pku/Resilience Image data (copy 1)/low resilience/low resilience/';
parentFolder_high   = '/home1/yangzy/Documents/Data/Resilience_fmri_pku/Resilience Image data (copy 1)/high resilience/high resilience/';

%% 2. Batch Execution Loop Over Contrasts
for c = 1:length(contrasts_to_run)
    current_contrast = contrasts_to_run{c};
    fprintf('\n=======================================================\n');
    fprintf('Starting ANOVA Analysis for Contrast: %s\n', current_contrast);
    fprintf('=======================================================\n');
    
    % Extract valid file paths and demographic covariates for each group
    [scans_health, sex_health, age_health] = get_subject_paths_and_info(parentFolder_health, excludeIDs, current_contrast);
    [scans_low, sex_low, age_low]          = get_subject_paths_and_info(parentFolder_low, excludeIDs, current_contrast);
    [scans_high, sex_high, age_high]       = get_subject_paths_and_info(parentFolder_high, excludeIDs, current_contrast);

    % Validate data availability
    if isempty(scans_health) || isempty(scans_low) || isempty(scans_high)
        warning('Missing sufficient group data for contrast: %s. Skipping to next...', current_contrast);
        continue;
    end

    % Concatenate covariates for the SPM design matrix
    sex_vector = [sex_health(:); sex_low(:); sex_high(:)];
    age_vector = [age_health(:); age_low(:); age_high(:)];
    
    % Create output directory for the current contrast
    group_analysisDir = fullfile(second_level_root, current_contrast);
    if ~exist(group_analysisDir, 'dir')
        mkdir(group_analysisDir);
    end

    % Build and run the ANOVA batch
    try
        matlabbatch = build_anova_batch(group_analysisDir, scans_health, scans_low, scans_high, sex_vector, age_vector, current_contrast);
        spm_jobman('run', matlabbatch);
        fprintf('Successfully completed ANOVA for: %s\n', current_contrast);
    catch ME
        warning('SPM batch failed for contrast %s. Error: %s', current_contrast, ME.message);
    end
end

fprintf('\nAll ANOVA models have been successfully processed!\n');


% =========================================================================
% Local Function: Build ANOVA matlabbatch
% =========================================================================
function matlabbatch = build_anova_batch(groupfolder, scan1, scan2, scan3, sex, age, contrast_name)
    matlabbatch{1}.spm.stats.factorial_design.dir = {groupfolder};
    
    % Group Cells Configuration: 1=Health, 2=Low Resilience, 3=High Resilience
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = scan1(:);
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = scan2(:);
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = scan3(:);
    
    % Statistical Assumptions
    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;     % Independent measurements
    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % Unequal variance assumption
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
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    
    % Masking and Global parameters
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
    
    % Contrast Specification (Weights preserved strictly as defined)
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    
    % 1. Main Effect F-Test
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = sprintf('%s Main Effect of Group (F-test)', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 -1 0; 0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none'; 
    
    % 2. Pairwise T-Tests
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = sprintf('%s Health > Low', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = sprintf('%s Low > Health', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = sprintf('%s High > Low', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = sprintf('%s Low > High', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = sprintf('%s Health > High', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = sprintf('%s High > Health', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
    
    % 3. Composite T-Tests (Health vs. Combined Trauma) 
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = sprintf('%s Health > Trauma(Low+High)', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1 -0.5 -0.5];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = sprintf('%s Trauma(Low+High) > Health', contrast_name);
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-1 0.5 0.5];
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.delete = 0;
end

% =========================================================================
% Local Function: Get Subject Paths and Demographics
% =========================================================================
function [scans_t1, sex_vector, age_vector] = get_subject_paths_and_info(parentFolder, excludeIDs, contrast_name)
    % Recursively checks for valid subjects, returning their contrast images
    % and behavioral covariates securely.
    
    folders = dir(parentFolder);
    if isempty(folders)
        scans_t1 = {}; sex_vector = []; age_vector = [];
        return;
    end

    folderNames = {folders([folders.isdir]).name};
    subjectFolders = folderNames(~ismember(folderNames, {'.', '..'}));
    
    subjectIDs = regexp(subjectFolders, '(\d+)', 'match');
    empty_idx = cellfun(@isempty, subjectIDs);
    subjectIDs(empty_idx) = [];
    
    subjectIDs = cellfun(@(x) sprintf('sub-%03d', str2double(x{1})), subjectIDs, 'UniformOutput', false);
    [subjectIDs, ~] = unique(subjectIDs, 'stable');
    subjectIDs = subjectIDs(~ismember(subjectIDs, excludeIDs));

    % Target Directories
    first_level_dir = '/public/home/yangzy/Documents/Data/Resilience_fmri_pku/resilience_task_run5/First_Level_Type1';
    tsv_dir = '/public/home/yangzy/Documents/Data/Resilience_fmri_pku/resilience_task_run5/Timepoint/sub_combined';

    num_subs = length(subjectIDs);
    valid_idx = false(num_subs, 1);
    scans_t1 = cell(num_subs, 1);
    sex_vector = zeros(num_subs, 1);
    age_vector = zeros(num_subs, 1);

    for i = 1:num_subs
        subjectID_FL = strrep(subjectIDs{i}, 'sub-', 'Sub_');
        
        % Mapping explicit contrast strings generated from the first-level script
        actual_nii_name = sprintf('con_%s_-_All_Sessions.nii', contrast_name);
        con_file = fullfile(first_level_dir, subjectID_FL, actual_nii_name);

        % Locate behavior data
        tsv_num = str2double(subjectIDs{i}(5:end));
        tsv_file = fullfile(tsv_dir, sprintf('sub%03d.tsv', tsv_num));

        if exist(con_file, 'file') && exist(tsv_file, 'file')
            scans_t1{i} = [con_file, ',1'];
            
            try
                data = readtable(tsv_file, 'FileType', 'text', 'Delimiter', '\t');
                if ismember('Sex', data.Properties.VariableNames) && ismember('Age', data.Properties.VariableNames)
                    sex_vector(i) = data.Sex(1); 
                    age_vector(i) = data.Age(1);
                    valid_idx(i) = true; 
                end
            catch
                warning('Failed to read behavioral TSV for %s', subjectIDs{i});
            end
        end
    end

    % Output only complete data pairs
    scans_t1   = scans_t1(valid_idx);
    sex_vector = sex_vector(valid_idx);
    age_vector = age_vector(valid_idx);
end