%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% HAPPE Version 2.0 - Including HAPPE+ER and HAPPILEE
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of the pipeline and user options, please see 
% the following manuscripts:
%   Gabard-Durnam, et al., (2018)
%   Lopez, et al., (----) - submitted, will update ***
%   Monachino, et al., (----) - submitted, will update ***
%   ...
%
% Contributors to HAPPE:
%   Laurel Joy Gabard-Durnam (laurel.gabarddurnam@gmail.com)
%   Adriana S. Mendez Leal (asmendezleal@gmail.com)
%   Carol L. Wilkinson (carol.wilkinson@childrens.harvard.edu)
%   April R. Levin (april.levin@childrens.harvard.edu)
%   ----
%   Alexa D. Monachino (alexamonachino@gmail.com)
%   Kelsie L. Lopez (k.lopez@northeastern.edu)
%
% HAPPE includes code that is dependent on the following third-party 
% software. Please reference this third-party software in any manuscripts 
% making use of HAPPE as below:
%   EEGLAB - A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
%       analysis of single-trial EEG dynamics. Journal of Neuroscience 
%       Methods 134:9-21
%   Cleanline by Tim Mullen (as an EEGlab plug-in): Mullen, T. (2012). 
%       NITRC: CleanLine: Tool/Resource Info. Available online at: 
%       http://www.nitrc.org/projects/cleanline.
%   MARA by Irene Winkler (as an EEGlab plug-in): Winkler, et al. Automatic
%       Classification of Artifactual ICA-Components for Artifact Removal 
%       in EEG Signals. Behavioral and Brain Functions 7:30 (2011).
%   FASTER segment-level channel interpolation code: Nolan*, H., Whelan*, R.,
%       & Reilly, R.B. (2010). FASTER: Fully Automated Statistical Thresholding
%       for EEG artifact Rejection. Journal of Neuroscience Methods, 192, 
%       152-162.
%   I have modified Matlab central code for the wICA function originally 
%       posted by Jordan Sorokin (2015) 
%       https://www.mathworks.com/matlabcentral/fileexchange/55413-wica-data-varargin-
%
% Any code that is not part of the third-party dependencies is released
% under the GNU General Public License version 3.
%
% Copyright 2018, 2021 Alexa Monachino, Kelsie Lopez, Laurel Gabard-Durnam
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License (version 3) as
% published by the Free Software Foundation.
%
% This software is being distributed with the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See GNU General
% Public License for more details.
%
% In no event shall Boston Children’s Hospital (BCH), the BCH Division of 
% Developmental Medicine, the Laboratories of Cognitive Neuroscience (LCN),
% Northeastern University, the Center for Cognitive and Brain Health (CCBH),
% the Department of Psychology at Northeastern University, the Plasticity 
% in Neurodevelopment (PINE) Lab, or software contributors be liable to any
% party for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software and its 
% documentation, even if any of the above listed entities and/or contributors
% have been advised of the possibility of such damage. Software and 
% documentation is provided “as is.” The listed entities and/or software 
% contributors are under no obligation to provide maintenance, support, 
% updates, enhancements, or modifications.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear ;
%% SET FOLDERS FOR HAPPE AND EEGLAB PATHS
disp('Preparing HAPPE...') ;
% SET HAPPE AND EEGLAB PATHS
happe_directory_path = fileparts(which(mfilename('fullpath'))) ;
eeglab_path = [happe_directory_path filesep 'Packages' filesep 'eeglab2020_0'] ;

% ADD HAPPE AND REQUIRED SUBFOLDERS TO THE PATH
addpath([happe_directory_path filesep 'acquisition_layout_information'], ...
    [happe_directory_path filesep 'scripts'], ...
    [happe_directory_path filesep 'scripts' filesep 'UI_scripts'], ...
    [happe_directory_path filesep 'Packages' filesep 'MARA-master'], ...
    eeglab_path, genpath([eeglab_path filesep 'functions']));
rmpath(genpath([eeglab_path filesep 'functions' filesep 'octavefunc']));

% ADD EEGLAB AND REQUIRED SUBFOLDERS TO THE PATH  
% *** Eventually allow users to set own eeglab path. For now, assume eeglab2020_0
plugin_directories = dir([eeglab_path filesep 'plugins']) ;
plugin_directories = strcat(eeglab_path, filesep, 'plugins', filesep, ...
    {plugin_directories.name}, ';') ;
addpath([plugin_directories{:}]) ;

% ADD CLEANLINE FOLDERS TO PATH
if exist('cleanline', 'file')
    cleanline_path = which('eegplugin_cleanline.m') ;
    cleanline_path = cleanline_path(1:strfind(cleanline_path, 'eegplugin_cleanline.m')-1) ;
    addpath(genpath(cleanline_path)) ;
else; error('Please make sure cleanline is on your path') ;
end

%% DETERMINE AND SET PATH TO DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    src_folder_name = input('Enter the path to the folder containing the dataset(s):\n> ','s') ;
    if exist(src_folder_name, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd (src_folder_name) ;

%% DETERMINE IF REPROCESSING DATA
% Use command line input to determine if reprocessing data. If invalid 
% input, repeat until valid. 
disp("Select an option:") ;
disp("  raw = Run on raw data from the start") ;
disp("  reprocess = Run on HAPPE-processed data starting post-waveleting/ICA") ;
reprocessing = choose2("raw", "reprocess") ;

% IF REPROCESSING... 
if reprocessing
    currDir = dir ;
    currDir = {currDir.name} ;
    addpath([src_folder_name filesep currDir{find(contains(currDir,'quality'))}]) ;
    % LOAD PREVIOUS DATAQM TABLE: Use command line input to determine the 
    % name of the previously created dataQM file. If invalid file, repeat 
    % until valid. The none option is available in case the pipeline failed
    % before the creation of the table but after waveleting/ICA processing
    % was completed.
    while true
        loaded_dataQM_filename = input(['Name of previously created dataQC .csv file:\n' ...
            'If no file exists, enter "none" (without quotations).\n' ...
            '> '], 's') ;
        if strcmpi(loaded_dataQM_filename, 'none'); break; end
        if isfile([src_folder_name filesep currDir{find(contains(currDir, ...
                'quality'))} filesep loaded_dataQM_filename])
            loaded_dataQM = readtable(loaded_dataQM_filename, 'Delimiter', ...
                'comma', 'PreserveVariableNames', true); break;
        else; disp("Invalid input: please enter correct file name."); 
        end
    end
    
    % SAVE OR OVERWRITE FILES: Use command line input to determine whether 
    % to save a new file or overwrite potentially existing files. If new, 
    % use command line input to set a default or custom suffix for created
    % files.
    disp("Files, such as processed data and quality metrics, may already exist for this dataset.") ;
    disp("  overwrite = Overwrite existing files") ;
    disp("  new = Save new files") ;
    if choose2('overwrite', 'new')
        disp("Use default or custom suffix for processed set?") ;
        disp("  default = Default name (_rerun_dd-mm-yyyy.mat).") ;
        disp("  custom = Create your own file name.") ;
        if choose2("custom", "default"); rerun_suffix = ['_rerun_' datestr(now, 'dd-mm-yyyy')] ;
        else; rerun_suffix = input('Enter custom suffix:\n> ', 's') ;
        end
    else; rerun_suffix = '' ;
    end
    
% OTHERWISE, Set the rerun_suffix to an empty char to be used in all code
% related to saving for the remainder of the pipeline.
else; rerun_suffix = '' ;
end

%% DETERMINE IF USING PRESET PARAMETERS
% Use command line input to determine if using pre-existing parameters. If 
% invalid input, repeat until valid.
disp('Load pre-existing set of input parameters? [Y/N]') ;
pre_exist = choose2("n", "y") ;

%% LOAD PRE-EXISTING PARAMETERS (IF EXISTING)
if pre_exist
    % PARAMETER FILE FOLDER: Use command line input to collect the path of 
    % the folder containing the parameter file. If invalid path, repeat 
    % until valid. Add folder to path.
    while true
        preparam_folder_name = input('Path to the folder containing the input parameters:\n> ','s') ;
        if exist(preparam_folder_name, 'dir') == 7; break ;
        else; disp("Invalid input: please enter the complete path to the folder.") ;
        end
    end
    addpath(preparam_folder_name) ;
    cd(preparam_folder_name) ;
    
    % PARAMETER FILE NAME: Use command line input to collect file name of 
    % pre-existing parameters. If invalid file, repeat until valid. Load file.
    while true
        param_file = [preparam_folder_name filesep input('Name of file containing pre-existing parameters:\n> ', 's')] ;
        if isfile(param_file)
            disp("Loading parameters...") ;
            load(param_file) ;
            disp("Parameters loaded.") ;
            break ;
        else
            disp("Invalid input: please enter correct file name.") ;
        end
    end
    
    % COMPARE CURRENT SCRIPT VERSION TO PARAMETER VERSION: If the versions 
    % don't match, use command line input to determine whether to continue
    % or end run via error. Will also end run with an error if no version
    % is included in the parameter set.
    if isfield(params, 'HAPPE_version') && ~strcmpi(params.HAPPE_version, '2_8_2')
        fprintf(['These parameters were saved using a different version of HAPPE\n' ...
            'and may be incompatible with the running pipeline.\n' ...
            'Proceed anyway? [Y/N]\n']) ;
        if ~choose2("n", "y")
            error('Mismatch between HAPPE and input parameter versions.') ; 
        end
    elseif ~isfield(params, 'HAPPE_version')
        error('Unable to identify the HAPPE version associated with this parameter set.') ;
    end
    
    % DISPLAY PARAMETERS AND ASK IF CHANGING THEM: List out the parameters
    % in the command window. Use command line input to determine whether to
    % change them.
    listParams(params) ;
    disp("Change an existing parameter? [Y/N]") ;
    changed_params = choose2("n", "y") ;
    cd(src_folder_name) ;

% IF NOT USING EXISTING PARAMETERS...
else
    if reprocessing; disp("PLEASE ENSURE THE PARAMETERS (THROUGH SEGMENTATION) THAT YOU ENTER MATCH THE ALREADY PROCESSED DATASET."); end
    params = struct() ;
    changed_params = 0 ;
end

%% SET PARAMETERS
% Uses user input to set parameters. See setParams for further
% documentation.
params = setParams(params, pre_exist, reprocessing, changed_params) ;

%% SAVE INPUT PARAMETERS
% If created new or changed parameter set, save as a new .mat file to a new
% folder (input_parameters) added to the source folder.
if ~pre_exist || changed_params
    % CREATE "input_parameters" FOLDER AND ADD IT TO PATH, unless it
    % already exists.
    if ~isfolder([src_folder_name filesep 'input_parameters'])
        mkdir([src_folder_name filesep 'input_parameters']) ;
    end
    addpath('input_parameters') ;
    cd input_parameters ;
    
    % DETERMINE PARAMETER FILE NAME: Prompt to use a default or custom name
    % for parameter file. If file exists, ask to create new file with a 
    % different name or overwrite existing file.
    disp("Parameter file save name:") ;
    disp("  default = Default name (inputParameters_dd-mm-yyyy.mat).") ;
    disp("  custom = Create your own file name.") ;
    if choose2("custom", "default")
        param_file = paramFile_validateExist(['inputParameters_' ...
            datestr(now, 'dd-mm-yyyy') '.mat'], 2) ;
    else
        disp("File name (Do not include .mat):") ;
        param_file = paramFile_validateExist([input('> ', 's') '.mat'], 0) ;
    end

    % SAVE PARAMETERS: Save the params variable to a .mat file using the
    % name created above.
    params.HAPPE_version = '2_8_2' ;
    disp("Saving parameters...") ;
    save(param_file, 'params') ;
    disp("Parameters saved.") ;
end

%% CREATE OUTPUT FOLDERS
% Create folders within source folder to store outputs, unless they already
% exist or are not applicable to selected processing steps.
cd (src_folder_name) ;
disp('Creating output folders...') ;
folder_names = {'intermediate_processing', 'wavelet_cleaned_continuous', ...
    'ICA', 'ERP_filtered', 'segmenting', 'processed', 'quality_assessment_outputs'} ;
if ~params.ERP_analysis; folder_names(ismember(folder_names, 'ERP_filtered')) = []; end
if params.ERP_analysis || params.low_density; folder_names(ismember(folder_names, 'ICA')) = [] ; end
folderNames = {} ;
for i=1:length(folder_names)
    if ~isfolder([src_folder_name filesep folder_names{i}])
        mkdir([src_folder_name filesep num2str(i) ' - ' folder_names{i}]) ;
        folderNames = [folderNames [num2str(i) ' - ' folder_names{i}]] ;
    end
end
clear('folder_names') ;
disp('Output folders created.') ;

%% SET SENSOR LAYOUT
% For file formats where the channel locations are not imbedded in the
% file, retrieve the sensor layout from HAPPE's stored channel location
% files, if applicable.
% .mat and .raw Files:
if ismember(params.datafileformat, [0,1]) && ~params.low_density
    disp('Setting sensor layout...') ;
    chanlocs = [happe_directory_path filesep 'acquisition_layout_information' filesep] ;
    % GSN Nets:
    if params.layout_type(1,1) == 1
        if params.layout_type(1,2) == 64; chanlocs = [chanlocs 'GSN65v2_0.sfp'] ;
        else error(['Invalid sensor layout selection.', newline, 'Users ' ...
                'wishing to use an unsupported layout can run HAPPE through' ...
                'BEAPP, as described in the HAPPE manuscript.']) ;
        end
    elseif params.layout_type(1,1) == 2
        if ismember(params.layout_type(1,2), [32, 64, 128, 256])
            chanlocs = [chanlocs 'GSN-HydroCel-' num2str(params.layout_type(1,2)+1) ...
                '.sfp'] ;
        else; error(['Invalid sensor layout selection.', newline, 'Users ' ...
                'wishing to use an unsupported layout can run HAPPE through' ...
                'BEAPP, as described in the HAPPE manuscript.']) ;
        end
    % BioSemi Nets: (NOTE: only supports 32 and 128)
    elseif params.layout_type(1,1) == 3
        if params.layout_type(1,2) == 32; chanlocs = [chanlocs 'BioSemi_32.elp'] ;
        elseif params.layout_type(1,2) == 128; chanlocs = [chanlocs 'BioSemi_128.elp'] ;
        else; error(['HAPPE does not currently support this number of channels' ...
                'for BioSemi nets as a .mat file.']) ;
        end
    % Other Nets: HAPPE does not currently support other nets as .mat or .raw
    else; error(['Invalid sensor layout selection.', newline, 'Users wishing' ...
            ' to use an unsupported layout can run HAPPE through BEAPP, as' ...
            'described in the HAPPE manuscript.']) ;
    end
    disp('Sensor layout loaded.') ;
% Other Files: If the file type is not supported by HAPPE, end with an error.
elseif ~ismember(params.datafileformat, [0:4]); error(['Invalid sensor ' ...
        'layout selection.', newline, 'Users wishing to use an unsupported ' ...
        'layout can run HAPPE through BEAPP, as described in the HAPPE ' ...
        'manuscript.']) ;
end

%% DETERMINE AND COLLECT DATA TO RUN
disp('Gathering files...') ;
cd(src_folder_name);
% DETERMINE FILE EXTENSION: Uses the datafileformat to determine the
% extension for the data files.
if params.datafileformat == 0; src_file_ext = '.mat' ;
elseif params.datafileformat == 1; src_file_ext = '.raw' ;
elseif params.datafileformat == 2; src_file_ext = '.set' ;
elseif params.datafileformat == 3; src_file_ext = '.cdt' ;
elseif params.datafileformat == 4; src_file_ext = '.mff' ;
end

% COLLECT FILE NAMES: Pull the names of all files in the source folder with
% the file extension set above.
FileNames = {dir(['*' src_file_ext]).name} ;

% LOCATE STIM FILE AND NAMES: If .mat task data with a seperate stim file, 
% enter the folder indicated to hold that .txt file, and load it.
if params.datafileformat == 0 && params.task_EEG_processing == 1
    cd(params.task_event_info_location) ;
    StimNames = {dir('*.txt').name} ;
end

%% INITIALIZE QUALITY REPORT METRICS
disp('Initializing report metrics...') ;
% DATA QUALITY METRICS: Cell array to hold the data quality metrics. If 
% conducting ERP analyses, create an additional cell array to hold the 
% additional data QC associated with the ERP onset tags. Cells are filled
% in throughout the rest of the pipeline.
chan_index = [1:length(params.chan_IDs)] ;
dataQC_names = {'File_Length_in_Seconds', 'Number_User_Selected_Channels' ...
    'Number_Good_Channels_Selected', 'Percent_Good_Channels_Selected', ...
    'Interpolated_Channel_IDs', 'Number_ICs_Rejected', 'Percent_ICs_Rejected', ...
    'Percent_Variance_Kept_of_Post-Wav_Data', 'Median_Artifact_Probablity_of_Post-Wav_Data', ...
    'Mean_Artifact_Probability_of_Kept_ICs', 'Range_Artifact_Probability_of_Kept_ICs', ...
    'Min_Artifact_Probability_of_Kept_ICs', 'Max_Artifact_Probability_of_Kept_ICs', ...
    'Channels_Interpolated_for_Each_Segment', 'Number_Segments_Before_Segment_Rejection', ...
    'Number_Segments_Post_Segment_Rejection'} ;
dataQC = cell(length(FileNames), 16) ;
if params.ERP_analysis
    dataQC_erp = cell(length(FileNames), length(params.task_onset_tags)*2) ;
    dataQC_erp_names = cell(1, length(params.task_onset_tags)*2) ;
    for i=1:length(params.task_onset_tags)
        dataQC_erp_names{i*2-1} = ['Number_ ' params.task_onset_tags{i} '_Segments_Before_Segment_Rejection'] ;
        dataQC_erp_names{i*2} = ['Number_ ' params.task_onset_tags{i} '_Segments_Post_Segment_Rejection'] ;
    end     
end
 
% PIPELINE QUALITY METRICS: Variables to hold the pipeline quality metrics.
% Cells are filled in throughout the rest of the pipeline.
if ~reprocessing
    if params.line_noise == 60; lnfreqsofinterest = [50, 55, 58, 59, 60, 61, 62, 65, 70] ;
    else; lnfreqsofinterest = [40, 45, 48, 49, 50, 51, 52, 55, 60] ; end
    
    if params.ERP_analysis; freqsofinterest = [.5, 1, 2, 5, 8, 12, 20, 30, 45, 70];
    else; freqsofinterest = [2, 5, 8, 12, 20, 30, 45, 70] ;
    end
    ln_means = [] ;
    wav_means = [] ;
    if ~params.low_density % && ~params.ERP_analysis 
        ica_means = [] ;
        artifactRej_means = [] ;
    end
else
    if strcmpi(loaded_dataQM_filename, 'none')
        loaded_dataQM = cell(length(FileNames), 16) ;
    end
end

%% RUN THE PREPROCESSING PIPELINE OVER EACH DATA FILE
% For each data file, run the preprocessing steps outlined using user-selected
% parameters. Collect report metrics on data quality and pipeline performance.
for current_file = 1:length(FileNames) 
    cd(src_folder_name) ;
    try
        % IF RUNNING DATA FROM THE START...
        if ~reprocessing
            %% LOAD RAW DATA FILE AND GET SAMPLING RATE
            % Import raw data into EEGLAB. Store file length in seconds and
            % the sampling sampling rate.
            disp(['Loading ' FileNames{current_file} '...']) ;
            
            % LOAD MATLAB FORMAT: Load data from .mat file. Throw error if 
            % EEG variable cannot load correctly.
            if params.datafileformat == 0
                load(FileNames{current_file}) ;
                % LOW DENSITY MATLAB: For low-density, read the sampling
                % rate according to the user-selected method (table or
                % shared rate). Load the data in using the specified
                % channel locations file (if applicable).
                if params.low_density
                    if params.same_srate; srate = params.srate ;
                    else
                        srateTable = table2cell(readtable(params.srate)) ;
                        srate = srateTable{find(contains(list, ...
                            FileNames{current_file})), 2} ;
                    end
                    EEGloaded = pop_importdata('dataformat', 'matlab', 'data', ...
                        FileNames{current_file}, 'srate', srate, ...
                        'chanlocs', params.chanlocs) ;
                % HIGH DENSITY MATLAB: For high density matlab files, get
                % the sampling rate and EEG variable name from the
                % information in the .mat file. Try to load the file, and
                % throw an error if it fails.
                else
                    srate = intersect(who, {'samplingRate', 'EEGSamplingRate'}) ; % *** NOTE: this will not work
                    srate = double(eval(srate{1})) ;                              % if the dataset has both variables
                    file_eeg_vname = intersect(who, params.potential_eeg_var_names) ;
                    try EEGloaded = pop_importegimat(FileNames{current_file}, ...
                            srate, 0.00, file_eeg_vname{1}) ;
                    catch err_msg
                        if strcmp(err_msg.identifier, 'MATLAB:badsubscript')
                            error(['Sorry, could not read the variable name of', ...
                                ' the EEG data. Please check your file.'])
                        else; error(err_msg.message) ;
                        end
                    end
                end
                % For all task EEG processing, import the task information.
                if params.task_EEG_processing
                    cd(params.task_event_info_location) ;
                    EEGloaded = pop_importevent(EEGloaded, 'append', 'no', ...
                        'event', StimNames{current_file}, 'fields', {'type' ...
                        'latency' 'status'}, 'skipline', 1, 'timeunit', ...
                        1E-3, 'align', NaN) ;
                    orig_event_info = EEGloaded.urevent ;
                    events = EEGloaded.event ;
                end 
                
            % LOAD .RAW FORMAT: Load data from .raw file and set events.
            elseif params.task_EEG_processing && params.datafileformat == 1
                EEGloaded = pop_readegi(FileNames{current_file}, [], [], 'auto') ;
                events = EEGloaded.event ;
                orig_event_info = EEGloaded.urevent ;
                srate = double(EEGloaded.srate) ;

            % LOAD .SET FORMAT: Includes events, if task.
            elseif params.datafileformat == 2
                load('-mat', FileNames{current_file}) ;
                EEGloaded = EEG ;
                srate = double(EEGloaded.srate) ;
                if params.task_EEG_processing == 1
                    events = EEGloaded.event ;
                end
            
            % LOAD .CDT FORMAT: Includes events, if task.
            elseif params.datafileformat == 3
                EEGloaded = loadcurry([src_folder_name filesep ...
                        FileNames{current_file}], 'CurryLocations', ...
                        'false') ;
                srate = double(EEGloaded.srate) ;
                if params.task_EEG_processing == 1
                    events = EEGloaded.event ;
                end
                
            % LOAD OTHER FORMAT: No support for other file formats or baseline
            % .raw data. End run via error.
            else; error('File format type unsupported or running non-task data with .raw file');
            end

            %% RENAME AND VALIDATE LOADED FILE
            % Use eeg_checkset to confirm that the file loaded in
            % correctly.
            disp('Validating file...') ;
            EEGloaded.setname = 'rawEEG' ;
            EEG = eeg_checkset(EEGloaded) ;

            %% DETERMINE FILE LENGTH, IN SECONDS, FROM EEG STRUCT
            dataQC{current_file, 1} = EEG.xmax ;
            
            %% UPDATE CHAN-IDS
            % If this is the first file, update the chan_IDs parameter
            % using information from user-input and from the file itself,
            % when applicable.
            if current_file == 1
                % COLLECT THE FULL CHANNEL LIST
                fullchans = {} ;
                % .mat Format:
                if params.datafileformat == 0
                    % Low Density (with chanlocs): Use the chanlocs file
                    % specified by the user to collect the full list of
                    % channels.
                    if params.low_density
                        if params.has_chanlocs
                            fullchans = {EEG.chanlocs.labels} ;
                        end
                    % High Density: Get the full list using the netdata_lib
                    % provided with HAPPE.
                    else
                    	load('happe_netdata_lib.mat') ;
                        net = aquiLayout_getInfo(params.layout_type, netdata_lib) ;
                        fullchans = net.leads_all ;
                    end
                % .raw Format: Get the full list using the netdata_lib
                % provided with HAPPE.
                elseif params.datafileformat == 1
                    load('happe_netdata_lib.mat') ;
                    net = aquiLayout_getInfo(params.layout_type, netdata_lib) ;
                    fullchans = net.leads_all ;
                % For files with imbedded channel locations, use the full
                % list of channels included in the file.
                else
                    fullchans = {EEG.chanlocs.labels} ;
                end
                
                % APPLY METHOD OF CHANNEL SELECTION ON CHANNELS: Apply the
                % user-specified method of channel selection to the full
                % list of channels to create the selected list.
                if strcmpi(params.chans_all, 'coi_exclude')
                    params.chan_IDs = setdiff(fullchans, params.chan_IDs) ;
                elseif strcmpi(params.chans_all, 'coi_include')
                    params.chan_IDs = params.chan_IDs ;
                elseif strcmpi(params.chans_all, 'all')
                    params.chan_IDs = fullchans ;
                end
                
                % UPDATE CHAN-INDEX TO MATCH NEW LENGTH OF CHAN_IDS
                chan_index = [1:length(params.chan_IDs)] ;
            end
           
            %% CORRECT CHANNELS
            % For non-.set files, channel locations do not import properly from
            % netstation by default, so the locations need to be corrected.
            % This step is not applicable to low-density files.
            if ismember(params.datafileformat, [0,1]) && ~params.low_density
                disp('Correcting Channels...') ;
                EEG = pop_chanedit(EEG, 'load', {chanlocs 'filetype' 'autodetect'});
                EEG = eeg_checkset(EEG) ;
            end
            
            %% SET AND CORRECT THE REFERENCE CHANNEL
            % Collects the reference channel for later re-referencing, if
            % selected by the user. If needed, corrects 
            if params.rereference_on
                disp('Finding the reference channel...') ;
                if params.datafileformat == 0 && ~params.low_density
                    EEG.nbchan = EEG.nbchan - 1 ;
                    index = find(~cellfun(@isempty, {EEG.urchanlocs.type})) ;
                    ref_elect = EEG.chanlocs(index) ;
                    EEG.chanlocs(index) = [] ;
                    EEG.urchanlocs(index) = [] ;
                    EEG.data(index,:) = [] ;
                    ref_elect.datachan = 0 ;
                    EEG.chaninfo.nodatchans = [EEG.chaninfo.nodatchans ref_elect] ;
                else
                    ref_elect = [] ;
                    for i=1:size(EEG.chaninfo.nodatchans, 2)
                        if strcmp('REF', EEG.chaninfo.nodatchans(i).type)
                            ref_elect = EEG.chaninfo.nodatchans(i) ;
                            break ;
                        end
                    end
                end
            end

            %% SET 10-20 CHANNELS
            % Apply 10-20 electrode names to EEG data using net data library. 
            % Because ICA is not run on low-density data, skip for this
            % data type. Will also skip if the names of 10-20 channels are 
            % already included. Does not skip for ERP because the 10-20 
            % channels are included in the high-density layout channels of 
            % interest, which would affect the selecting channels step.
            tentwenty = {'FP1', 'FP2', 'F7', 'F3', 'F4', 'F8', 'C3', 'C4',  ...
                'T5', 'PZ', 'T6', 'O1', 'O2', 'T3', 'T4', 'P3', 'P4', 'Fz'} ;
            if ~params.low_density && isempty(intersect({EEG.chanlocs.labels}, tentwenty))
                disp('Setting 10-20 channels...') ;
                load('happe_netdata_lib.mat')
                net_info = aquiLayout_getInfo(params.layout_type, netdata_lib) ;
                if params.layout_type(1,1) == 3
                    for i=1:length(net_info.lead_nums_sub)
                        EEG = pop_chanedit(EEG, 'changefield', {net_info.lead_nums_sub{i} ...
                            'labels' net_info.lead_list_sub{i}}) ;
                    end
                else
                    for i=1:length(net_info.lead_nums_sub)
                        EEG = pop_chanedit(EEG, 'changefield', {net_info.lead_nums_sub(i) ...
                            'labels' net_info.lead_list_sub{i}}) ;
                    end
                end
            end
            clear('tentwenty') ;

            %% SELECT AND FILTER TO THE CHANNELS OF INTEREST
            % Removes all channels not specified in the channels of 
            % interest list designated by the user. Does not apply to
            % low-density .mat files without channel locations.
            if ~(params.low_density && params.datafileformat == 0 && ~params.has_chanlocs)
                disp('Selecting channels of interest...') ;
                EEG = pop_select(EEG,'channel', intersect({EEG.chanlocs.labels}, ...
                    params.chan_IDs)) ;
                EEG.setname = 'rawEEG_cs' ;
                full_selected_channels = EEG.chanlocs ;
            end

            %% REDUCE LINE NOISE
            % Reduces line noise in the data. This step may not completely
            % eliminate all line noise - an effect that may be mitigated by
            % re-referencing.
            disp('Reducing line noise...') ;
            EEGpreLNreduce = reshape(EEG.data, size(EEG.data,1), []) ;
            % LEGACY LINENOISE: If the user indicated to use the legacy 
            % method of linenoise reduction, the code from HAPPE v1 is 
            % applied to the data.
            if params.legacy_linenoise
                EEG = pop_cleanline(EEG, 'bandwidth', 2, 'chanlist', ...
                    chan_index, 'computepower', 1, 'linefreqs',...
                    [params.line_noise, params.line_noise*2], ...
                    'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', ...
                    0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', ...
                    100, 'verb', 0, 'winsize', 4, 'winstep', 1, ...
                    'ComputeSpectralPower', 'False');
                EEG.setname = 'rawEEG_cs_ln';
            
            % DEFAULT LINE NOISE REDUCTION:
            else
                signal = struct('data', EEG.data, 'srate', EEG.srate) ;
                lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
                    'lineNoiseChannels', 1:EEG.nbchan, 'Fs', EEG.srate, ...
                    'lineFrequencies', params.line_noise, 'p', 0.01, ...
                    'fScanBandWidth', 2, 'taperBandWidth', 2, ...
                    'taperWindowSize', 4, 'taperWindowStep', 4, ...
                    'tau', 100, 'pad', 2, 'fPassBand', [0 EEG.srate/2], ...
                    'maximumIterations', 10) ;
                [EEG, lineNoiseOut] = cleanLineNoise(EEG, lineNoiseIn) ;
            end
            EEGpostLNreduce = reshape(EEG.data, size(EEG.data,1), []) ;

            % LINE NOISE REDUCTION QM: Assesses the performance of line
            % noise reduction.
            % Cross Correlation Across Data: Evaluated between pre- and
            % post-line noise reduced EEG data.
            cc_ln = corrcoef(EEGpreLNreduce, EEGpostLNreduce);
            [cross_corr] = cc_ln(1,2);
            % Cross Correlation Across Channels: Evaluated for the 
            % frequencies of interest across all channels based on the 
            % frequency of line noise in the data.
            [cxy, f] = mscohere(EEGpostLNreduce', EEGpreLNreduce', 1000, ...
                0, lnfreqsofinterest, EEG.srate) ;
            % Stores the evaluation metrics for later output.
            ln_means = [ln_means; [(cross_corr) mean(cxy, 2, 'omitnan')']] ;

            %% CLOSE WINDOW (IF VISUALIZATIONS OFF)
            if params.visualizations == 0; close gcf; end

            %% RESAMPLE DATA
            % Resamples the data, unless conducting ERP analyses *** check
            % with Laurel to confirm that this is an option
            if params.downsample_freq ~= 0 && ~params.ERP_analysis
                disp('Resampling the data...') ;
                EEG = pop_resample(EEG, params.downsample_freq) ;
            end

            %% FILTER DATA
            % Filter the data. If ERP, just use the 100Hz filter.
            % Otherwise, filter using 1hz highpass and bandpass 1hz-100hz.
            cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
                'intermediate_processing'))}]) ;
            if params.ERP_analysis; EEG = pop_eegfiltnew(EEG, [], 100, [], ...
                    0, [], 0) ;
            else; EEG = pop_eegfiltnew(EEG, 1, 100, [], 0, [], 0) ;
            end
            EEG.setname = 'rawEEG_f_cs_ln' ;
            
            %% SAVE FILTERED, LINENOISE REDUCED FILE
            pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                src_file_ext, '_filtered_lnreduced.set')) ;

            %% BAD CHANNEL DETECTION
            % Detects bad channels in the data and removes them. Additionally
            % stores information about the good channels, and saves the names of
            % rejected channels for QC output.
            if params.badchanrej
                disp('Detecting bad channels...') ;
                % LEGACY DETECTION:
                % Conducts crude bad channel detection using spectrum criteria and
                % 3SDeviations as channel outlier threshold (twice). This
                % option is not available for low density layouts.
                if params.legacy_channels
                    EEG = pop_rejchan(EEG, 'elec', chan_index, 'threshold', [-3 3], ...
                        'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]) ;
                    EEG.setname = 'rawEEG_f_cs_ln_badc' ;
                    EEG = pop_rejchan(EEG, 'elec', [1:EEG.nbchan], 'threshold', ...
                        [-3 3], 'norm', 'on', 'measure', 'spec', 'freqrange', [1 125]);
                    EEG.setname = 'rawEEG_f_cs_ln_badc2x' ;
                
                % DEFAULT DETECTION:
                % Conducts bad channel detection optimized for HAPPE 2,
                % HAPPE+ER, and HAPPILEE.
                else
                    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, ...
                        'ChannelCriterion', .1, 'LineNoiseCriterion', ...
                        20, 'Highpass', 'off', 'BurstCriterion', 'off', ...
                        'WindowCriterion', 'off', 'BurstRejection', 'off', ...
                        'Distance', 'Euclidian') ;
                    % If low-density data, detect bad channels using the
                    % methods optimized in HAPPILEE.
                    if params.low_density
                         EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', ...
                             'off', 'ChannelCriterion', .7, 'LineNoiseCriterion', ...
                             2.5, 'Highpass', 'off', 'BurstCriterion', 'off', ...
                             'WindowCriterion', 'off', 'BurstRejection', ...
                             'off', 'Distance', 'Euclidian') ;
                        EEG = pop_rejchan(EEG, 'elec', [1:EEG.nbchan], ...
                            'threshold', [-2.75 2.75], 'norm', 'on', 'measure', ...
                            'spec', 'freqrange', [1 100]) ;
                    % Otherwise, detect bad channels using methods
                    % optimized in HAPPE 2 and HAPPE+ER.
                    else
                        EEG = pop_rejchan(EEG, 'elec', [1:EEG.nbchan], 'threshold', [-5 3.5], ...
                            'norm', 'on', 'measure', 'spec', 'freqrange', [1 100]) ;
                        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', ...
                            .8, 'LineNoiseCriterion', 6, 'Highpass', 'off', 'BurstCriterion', ...
                            'off', 'WindowCriterion', 'off', 'BurstRejection', 'off', 'Distance', ...
                            'Euclidian') ;
                    end
                end
                
                % SAVE PRE-WAVLET, POST-BAD CHAN REJ EEG
                pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, '_badchanrej.set')) ;
            end

            % BAD CHANNEL DETECTION QC: Assesses the performance of bad
            % channel rejection. Measures include the names of the channels
            % labeled as bad vs good, and the channel locations of the good
            % channels.
            if ~(params.low_density && params.datafileformat == 0 && ~params.has_chanlocs)
                disp('Evaluating bad channel detection...') ;
                selected_good_channel_locations = EEG.chanlocs ;
                selected_channel_labels = {selected_good_channel_locations.labels} ;
                bad_channels_removed = setdiff(params.chan_IDs, selected_channel_labels, ...
                    'stable') ;
                if exist('params.ROI_channels', 'var')
                    [~,ROI_indices_in_selected_chanlocs] = intersect(selected_channel_labels, ...
                        params.ROI_channels,'stable') ;
                end
            end
            
            % Hold current EEG for later artifact rejection QC - only
            % applicable to data running through ICA.
            if ~params.low_density && ~params.ERP_analysis
                EEG_preAR = EEG ;
            end

            %% WAVELET THRESHOLDING
            % Applies wavelet thresholding to remove artifact from the data.
            disp('Wavelet thresholding...') ;
            cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
                'wavelet_cleaned_continuous'))}]) ;
            
            % Set wavelet level depending on the task paradigm *** for low
            % density
            if params.ERP_analysis
                if EEG.srate > 500; wavLevel = 13 ;
                elseif EEG.srate > 250 && EEG.srate <= 500; wavLevel= 12 ;
                elseif EEG.srate <= 250; wavLevel = 11 ;
                else; wavLevel = 7 ; 
                end
            else; wavLevel = 9 ;
            end

            % LEGACY WAVELET:
            % Kept from HAPPEv1 so comparison/unfinished analyses can be
            % run without needing to switch versions.
            % ICA for clustering data. Uses a soft, global threshold for 
            % the wavelets. The wavelet family is coiflet (level 5). 
            % Threshold multiplier is used to remove more high frequency 
            % noise or for ERP analyses. For details, refer to wICA.m.
            if params.legacy_wavelet
                if params.ERP_analysis == 1; threshmultiplier = 3;
                elseif params.ERP_analysis == 0; threshmultiplier = 1;
                end
                try
                    if params.visualizations == 0
                        [wIC, A, W, IC] = wICA(EEG, 'runica', threshmultiplier, ...
                            0, [], 5, 'coif5');
                    elseif params.visualizations == 1
                        [wIC, A, W, IC] = wICA(EEG, 'runica', threshmultiplier, ...
                            1, srate, 5, 'coif5');
                    end
                catch wica_err
                    if strcmp (['Output argument "wIC" (and maybe others) ' ...
                            'not assigned during call to "wICA".'], wica_err.message)
                        error('Error during wICA, most likely due to memory settings. Please confirm your EEGLAB memory settings are set according to the description in the HAPPE ReadMe')
                    else; rethrow(wica_err)
                    end
                end
                artifacts = A * wIC ;
                
            % DEFAULT WAVELET:
            % Uses a global threshold for the wavelets. Wavelet family is
            % coiflet (level depending). Threshold the wavelet coefficients
            % to generate artifact signals, reconstructing signal as 
            % channels x samples format.
            else
                % Low-Density Data: Uses hard threshold - not used for ERPs
                if params.low_density && ~params.ERP_analysis
                    artifacts = wdenoise(reshape(EEG.data, size(EEG.data, 1), ...
                        [])', 10, 'Wavelet', 'coif4', 'DenoisingMethod', ...
                        'Bayes', 'ThresholdRule', 'Hard', 'NoiseEstimate', ...
                        'LevelDependent')' ;
                % Otherwise: Uses soft threshold
                else
                    artifacts = wdenoise(reshape(EEG.data, size(EEG.data, 1), ...
                        [])', wavLevel, 'Wavelet', 'coif4', 'DenoisingMethod', ...
                        'Bayes', 'ThresholdRule', 'Soft', 'NoiseEstimate', ...
                        'LevelDependent')' ;
                end
            end

            
            % REMOVE ARTIFACT FROM DATA: Subtract out the wavelet artifact
            % signal from the EEG signal and save the wavcleaned data into
            % EEGLAB structure. If conducting ERP analyses, filter the data
            % to the user-specified frequency range for analyses purposes
            % only.
            if params.ERP_analysis   
                preWavEEG = reshape(pop_eegfiltnew(EEG, params.ERP_highpass_cutoff, ...
                    params.ERP_lowpass_cutoff, [], 0, [], 0).data, ...
                    size(EEG.data,1), []) ;
                EEG2D = reshape(EEG.data, size(EEG.data,1), []) ;
                wavcleanEEG = EEG2D - artifacts ;
                EEG.data = wavcleanEEG ;
                EEG.setname = 'wavcleanedEEG' ;
                postWavEEG = reshape(pop_eegfiltnew(EEG, params.ERP_highpass_cutoff, ...
                    params.ERP_lowpass_cutoff, [], 0, [], 0).data, ...
                    size(EEG.data,1), []) ;
            else
                preWavEEG = reshape(EEG.data, size(EEG.data,1), []) ;
                postWavEEG = preWavEEG - artifacts ;
                EEG.data = postWavEEG ;
                EEG.setname = 'wavcleanedEEG' ;
            end

            % SAVE WAVELET-CLEANED EEG
            pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                src_file_ext, '_wavclean.set')) ;
            
            % WAVELETING QC METRICS
            disp('Evaluating wavelet thresholding...') ;
            % Root Mean Squared Error (RMSE):
            RMSE_alldata = calcRMSE(preWavEEG, postWavEEG, 2) ;
            % Mean Absolute Error (MAE):
            MAE_alldata = calcMAE(preWavEEG, postWavEEG, 2) ;
            % Signal to Noise Ratio (SNR) AND Peak Signal to Noise Ratio (PSNR):
            [SNR_alldata, PeakSNR_alldata] = calcSNR_PSNR(preWavEEG, postWavEEG, 2) ;
            % Cross Correlation Across Data and by Channel:
            cc = corrcoef(postWavEEG, preWavEEG) ;
            [cross_corr] = cc(1, 2) ;
            [cxy, f] = mscohere(postWavEEG', preWavEEG', 1000, 0, ...
                freqsofinterest, EEG.srate) ;
            % Store waveleting pipeline QC
            wav_means = [wav_means; [(RMSE_alldata) (MAE_alldata) ...
                (PeakSNR_alldata) (cross_corr) mean(cxy, 2, 'omitnan')']] ;

            %% INDEPENDENT COMPONENT ANALYSIS (ICA)
            % Performs ICA and uses the components to remove artifact from
            % the data. Uses MARA with an artifact rejection threshold of 
            % 50%. Not conducted on low-density data or ERP analysis. If
            % ICA is not run, set all relevent QC to NA.
            if ~params.low_density && ~params.ERP_analysis
                disp('Independent component analysis...') ;
                cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
                'ICA'))}]) ;
                % RUN ICA: Save the weights and spheres.
                EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on') ;
                ICAweightstotransfer = EEG.icaweights ;
                ICAspheretotransfer = EEG.icasphere ;
                
                % SAVE: Includes the ICA decomposition.
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, '_ICA.set')) ;

                % FLAG ARTIFACT ICs USING MARA: Automatically flag artifact
                % ICs if probability > 0.5.
                [~, EEG, ~] = processMARA(EEG,EEG,EEG, [0, 0, params.visualizations, ...
                    params.visualizations, params.visualizations]) ;
                EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject)) ;
                EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob >= .5) = 1 ;
                ICs_with_reject_flagged = EEG.reject.gcompreject ;
                EEG.setname = 'wavcleanedEEG_ICA_MARA' ;

                % STORE DATA/PIPELINE QM: Use MARA variables as measures of
                % data quality, and store IC variables and variance 
                % retained post-IC rejection.
                index_ICs_kept = find(EEG.reject.gcompreject == 0) ;
                dataQC{current_file, 9} = num2str(median(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept))) ;
                dataQC{current_file, 10} = num2str(mean(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept))) ;
                dataQC{current_file, 11} = num2str(range(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept))) ;
                dataQC{current_file, 12} = num2str(min(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept))) ;
                dataQC{current_file, 13} = num2str(max(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept))) ;
                [projWav, varianceWav] = compvar(EEG.data, EEG.icaact, ...
                    EEG.icawinv, index_ICs_kept) ;
                dataQC{current_file, 8} = varianceWav ;

                % REJECT ARTIFACT ICS
                artifact_ICs = find(EEG.reject.gcompreject == 1) ;
                EEG = pop_subcomp(EEG, artifact_ICs, 0) ;
                EEG.setname = 'wavcleanedEEG_ICA_MARA_rej' ;

                % SAVE POST-MARA CLEANED INTERMEDIATE FILE
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, '_ICAcleanedwithMARA.set')) ;
                
                % Hold post-artifact rejection variable for QC
                EEG_postAR = EEG ;

                % ICA QUALITY MEASURES
                disp('Evaluating ICA...') ;
                % Number and Percent of ICs Rejected:
                dataQC{current_file, 6} = num2str(length(artifact_ICs)) ;
                dataQC{current_file, 7} = str2double(dataQC(current_file,6))/str2double(dataQC(current_file,3)) * 100 ;
                % Reshape EEG to channels x samples format
                EEG2DICA = reshape(EEG.data, size(EEG.data,1), []) ;
                % RMSE:
                RMSE_alldata = calcRMSE(EEG2DICA, postWavEEG, 1) ;
                % MAE:
                MAE_alldata = calcMAE(EEG2DICA, postWavEEG, 1) ;
                % SNR and PSNR:
                [SNR_alldata, PeakSNR_alldata] = calcSNR_PSNR(EEG2DICA, postWavEEG, 1) ;
                % Cross Correlation Across Data and by Channel:
                cc = corrcoef(EEG2DICA, postWavEEG) ;
                [cross_corr] = cc(1, 2) ;
                [cxy, f] = mscohere(EEG2DICA', postWavEEG', 1000, 0, ...
                    freqsofinterest, EEG.srate) ;
                % Store ICA QC
                ica_means = [ica_means; [(RMSE_alldata) (MAE_alldata) (SNR_alldata) ...
                    (PeakSNR_alldata) (cross_corr) mean(cxy, 2, 'omitnan')']] ;
                
                %% PIPELINE QUALITY MEASURES ON ALL ARTIFACT REJECTION
                disp('Evaluating all artifact rejection...') ;
                EEG_preAR = reshape(EEG_preAR.data, size(EEG_preAR.data,1), []) ;
                EEG_postAR = reshape(EEG_postAR.data, size(EEG_postAR.data, 1), []) ;
                % RMSE:
                RMSE_alldata = calcRMSE(EEG_postAR, EEG_preAR, 1) ;
                % MAE:
                MAE_alldata = calcMAE(EEG_postAR, EEG_postAR, 1) ;
                % SNR AND PSNR:
                [SNR_alldata, PeakSNR_alldata] = calcSNR_PSNR(EEG_postAR, EEG_preAR, 1) ;
                % Cross Correlation across Data and by Channel:
                cc = corrcoef(EEG_postAR, EEG_preAR);
                [cross_corr] = cc(1,2);
                [cxy, f] = mscohere(EEG_postAR', EEG_preAR', 1000, 0, ...
                    freqsofinterest, EEG.srate);
                % Store QC before and after all artifact rejection
                artifactRej_means = [artifactRej_means; [(RMSE_alldata) ...
                    (MAE_alldata) (SNR_alldata) (PeakSNR_alldata) (cross_corr) ...
                    mean(cxy, 2, 'omitnan')']] ;
            else
                for i=6:13
                    dataQC{current_file, i} = 'NA' ;
                end
            end
            
            %% FILTERING FOR ERP ANALYSES
            % Apply the user-specified high- and low-pass filters to the
            % data.
            if params.ERP_analysis == 1
                disp('Filtering using ERP cutoffs...') ;
                cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
                'ERP_filtered'))}]) ;
                EEG = pop_eegfiltnew(EEG, params.ERP_highpass_cutoff, ...
                    params.ERP_lowpass_cutoff, [], 0, [], 0) ;
                EEG.setname = 'wavcleanedEEGforERP_rej_filt' ;

                % SAVE ERP-PREPPED FILE
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, '_wavcleaned_filteredERP.set'));
            end
        end
        
        %% LOAD DATA IF REPROCESSING (REPROCESSING ONLY)
        % If data is being reprocessed, load in the files with
        % already-processed data needed to continue.
        if reprocessing
            % LOAD FILTERED, LN-REDUCED FILE: From this, get information
            % about the full set of selected channels in the data.
            disp('Loading filtered, line-noise reduced file...') ;
            EEG = pop_loadset('filename', strrep(FileNames{current_file}, ...
                src_file_ext, '_filtered_lnreduced.set'), 'filepath', ...
                [src_folder_name filesep folderNames{find(contains(folderNames, ...
                'intermediate_processing'))}]) ;
            full_selected_channels = EEG.chanlocs ;
            
            % LOAD POST-ARTIFACT REJECTION FILE
            disp('Loading artifact-cleaned file...') ;
            % For ERPs, this is the file post-waveleting and post-filtering
            % using user-selected filters.
            if params.ERP_analysis
                EEG = pop_loadset('filename', strrep(FileNames{current_file}, ...
                    '.raw', '_wavcleaned_filteredERP.set'), 'filepath', ...
                    [src_folder_name filesep folderNames{find(contains(folderNames, ...
                    'ERP_filtered'))}]) ;
            else
                % Since ICA is not run on low-density data, load the
                % wavelet-cleaned file.
                if params.low_density
                    EEG = pop_loadset('filename', strrep(FileNames{current_file}, ...
                        src_file_ext, '_wavclean.mat'), 'filepath', ...
                        [src_folder_name filesep folderNames{find(contains(folderNames, 'wavelet_cleaned_continuous'))}]) ;
                % Load the waveleted and ICA cleaned data.
                else
                    EEG = pop_loadset('filename', strrep(FileNames{current_file}, ...
                        src_file_ext, '_ICAcleanedwithMARA.set'), 'filepath', ...
                        [src_folder_name filesep folderNames{find(contains(folderNames, 'ICA'))}]) ;
                end
            end

            % STORE INFORMATION ABOUT LOADED EEGS
            srate = double(EEG.srate) ;
            selected_good_channel_locations = EEG.chanlocs ;
            selected_channel_labels = {selected_good_channel_locations.labels} ;
            bad_channels_removed = setdiff(params.chan_IDs, selected_channel_labels, 'stable') ;
        end
        
        %% SEGMENT DATA
        % If option selected by user, segment the data based on whether the
        % data is time-frequency or ERP.
        cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
            'segmenting'))}]) ;
        if params.segment_data
            disp('Segmenting...') ;
            % IF NOT ERP PROCESSING, SEGMENT USING USER-SELECTED LENGTH
            if ~params.task_EEG_processing
                EEG = eeg_regepochs(EEG, 'recurrence', params.segment_length, ...
                    'limits', [0 params.segment_length], 'rmbase', [NaN]);
            
            % IF ERP PROCESSING... 
            else
                % Transform task offset from milliseconds to samples
                samples_offset = srate * params.task_offset/1000 ;

                % Correct for timing offset (in samples) between event 
                % initiation and presentation.
                for i = 1:size(EEG.event, 2)
                    EEG.event(i).latency = EEG.event(i).latency + samples_offset ;
                end

                % Generate segments around the corrected stimulus presentation
                % timing.        
                EEG = pop_epoch(EEG, params.task_onset_tags, ...
                    [params.task_segment_start params.task_segment_end], ...
                    'verbose', 'no', 'epochinfo', 'yes');
            end

            % SAVE SEGMENTED DATA
            EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                src_file_ext, ['_segmented' rerun_suffix '.set'])) ;
        end
        
        % SEGMENTATION QC METRICS
        dataQC{current_file, 15} = EEG.trials ;
        % ERP-specific segmentation qc
        if params.task_EEG_processing
            for i=1:length(params.task_onset_tags)
                try
                    dataQC_erp{current_file, i*2-1} = length(pop_selectevent(EEG, ...
                        'type', params.task_onset_tags{i}).epoch) ;
                catch
                    dataQC_erp{current_file, i*2-1} = 'ERROR' ;
                    fprintf("No instances of tag '%s' for this file.\n", ...
                        params.task_onset_tags{i}) ;
                end
            end
        end

        %% BASELINE CORRECT TASK (Only if ERP processing & user-selected)
        if params.task_EEG_processing && isfield(params, 'baseline_correction') ...
                && params.baseline_correction
            disp('Performing baseline correction...') ;
            EEG = pop_rmbase(EEG, [params.baseline_corr_start params.baseline_corr_end]) ;
            
            % SAVE BASELINE CORRECTED DATA
            EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                src_file_ext, ['_segmented_BLcorrected' rerun_suffix '.set'])) ;
        end

        %% INTERPOLATE BAD DATA USING CHANNELS
        % Evaluate the channels for each segment and interpolate channels 
        % with bad data for each segment using the FASTER program, 
        % interpolating channels scoring above/below z threshold of 3 for 
        % a segment: Code adapted from FASTER program (Nolan et al., 2010).
        % Not available if bad channel detection was not performed.
        if params.segment_interpolation
            disp("Interpolating bad data...") ;
            eeg_chans = [1:length(selected_good_channel_locations)] ;
            ext_chans = [] ;
            o.epoch_interp_options.rejection_options.measure = [1 1 1 1] ;
            o.epoch_interp_options.rejection_options.z = [3 3 3 3] ;
            if length(size(EEG.data)) > 2
                status = '' ;
                lengths_ep = cell(1, size(EEG.data,3)) ;
                for v=1:size(EEG.data, 3)
                    list_properties = single_epoch_channel_properties(EEG, ...
                        v, eeg_chans);
                    lengths_ep{v} = eeg_chans(logical(min_z(list_properties, ...
                        o.epoch_interp_options.rejection_options)));
                    status = [status sprintf('[%d:',v) sprintf(' %d', ...
                        lengths_ep{v}) ']'] ;
                end
                EEG = h_epoch_interp_spl(EEG, lengths_ep, ext_chans) ;
                EEG.saved = 'no' ;

                % Add info about which channels were interpolated for each 
                % segment to the dataQM output csv.
                EEG.etc.epoch_interp_info = [status] ;
                dataQC{current_file, 14} = cellstr(status) ;
            end

            % SAVE INTERPOLATED DATA
            EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                src_file_ext, ['_segments_interp' rerun_suffix '.set'])) ;
        
        % If not interpolating, note that no channels were interpolated in
        % the QM
        else; dataQC{current_file, 14} = 'NA' ;
        end

        %% SEGMENT REJECTION
        % Using amplitude-based and joint probability artifact detection,
        % reject segments. If indicated by user, will only use the
        % specified region of interest. Otherwise, will examine all
        % channels.
        if params.segment_rejection
            disp('Rejecting segments...')
            % APPLY AMPLITUDE CRITERIA:
            if strcmpi(params.seg_rej_method, 'amplitude') || ... 
                    strcmpi(params.seg_rej_method, 'both')
                if params.ROI_channels_only
                    EEG = pop_eegthresh(EEG, 1, [ROI_indices_in_selected_chanlocs]', ...
                        [params.reject_min_amp], [params.reject_max_amp], [EEG.xmin], ...
                        [EEG.xmax], 2, 0) ;
                else
                    EEG = pop_eegthresh(EEG, 1, [1:EEG.nbchan], [params.reject_min_amp], ...
                        [params.reject_max_amp], [EEG.xmin], [EEG.xmax], 2, 0) ;
                end
            end
            
            % APPLY SIMILARITY CRITERIA:
            if strcmpi(params.seg_rej_method, 'similarity') || ... 
                    strcmpi(params.seg_rej_method, 'both')
                if params.low_density; num = 2; else; num = 3 ; end
                if ~params.ROI_channels_only
                    EEG = pop_jointprob(EEG, 1, [1:EEG.nbchan], num, num, ...
                        params.visualizations, 0, params.visualizations, [], ...
                        params.visualizations) ;
                else
                    EEG = pop_jointprob(EEG, 1, [ROI_indices_in_selected_chanlocs]', ...
                        num, num, params.visualizations, 0, params.visualizations, ...
                        [], params.visualizations) ;
                end
            end

            % REJECT TRIALS: Use the trials that the user has flagged as
            % bad, if the option was selected.
            if params.task_EEG_processing == 1 && params.datafileformat == 0 ...
                    && params.user_selected_trials == 1
                EEG = pop_selectevent(EEG, 'status', 'good', 'deleteevents', ...
                    'on', 'deleteepochs', 'on', 'invertepochs', 'off') ;
            end
            EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1) ;

            % CONFIRM TRIALS AND SAVE: Check to see if all trials have been
            % rejected. If so, print a warning to the command line and move
            % onto the next file. Otherwise, save the file.
            if (isfield(EEG, 'reject') && all(EEG.reject.rejglobal)) || ...
                    (isfield(EEG, 'rej') && all(EEG.rej.rejglobal))
                EEG.trials = 0 ;
                error('MATLAB:happeAllTrialsRej', ['All trials rejected.' ...
                    ' No further processing on this file is possible.']) ;
            else
                EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal], 0) ;
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, ['_segments_postreject' rerun_suffix '.set'])) ;
            end  
        end

        %% INTERPOLATE BAD CHANNELS
        if ~(params.low_density && params.datafileformat == 0 && ~params.has_chanlocs)
            EEG = pop_interp(EEG, full_selected_channels, 'spherical') ;
            EEG.setname = 'wavcleanedEEG_ICA_MARA_rej_chan_int' ;
        end

        %% REREFERENCE DATA
        if params.rereference_on
            disp('Re-referencing...') ;
            % AVERAGE REREFERENCE:
            if params.average_rereference == 1 
                EEG = pop_reref(EEG, [], 'keepref', 'on', 'refloc', ref_elect) ;
                EEG.setname = 'wavcleanedEEG_ICA_MARA_rej_chan_int_avgreref' ;
            
            % REREFERENCE DATA BY SUBSET:
            else
                [~, ref_chan_indices_in_full_selected_chanlocs] = intersect({full_selected_channels.labels}, ...
                    params.no_av_reref_channel_subset, 'stable') ;
                EEG = pop_reref(EEG, ref_chan_indices_in_full_selected_chanlocs) ;
                EEG.setname = 'wavcleanedEEG_ICA_MARA_rej_chan_int_chansubsetreref' ;
            end
        end

        %% SPLIT BY TAGS (IF ERP)
        % Create seperate EEG files for each user-specified ERP stimulus.
        % If the tag of interest doesn't exist, mark an error in the QM and
        % alert the user through the command line.
        if params.ERP_analysis && size(params.task_onset_tags, 2) > 1
            disp('Selecting by tags...') ;
            eeg_byTags = [] ;
            for i=1:length(params.task_onset_tags)
                try
                    eeg_byTags = [eeg_byTags pop_selectevent(EEG, 'type', ...
                        params.task_onset_tags{i})] ;
                    dataQC_erp{current_file, i*2} = length(eeg_byTags(i).epoch) ;
                catch
                    dataQC_erp{current_file, i*2} = 'ERROR' ;
                    fprintf('No instances of tag %s appear in this file', ...
                        params.task_onset_tags{i}) ;
                end
            end
        end

        %% STORE REMAINING OUTPUTS AND REPORT METRICS
        disp('Evaluating remaining QC metrics...') ;
        if ~reprocessing
            % If the file has channel names, include the number of channels
            % selected by the user and the number of good channels
            % remaining after bad channel rejection. If no channel names
            % are included, instead list out the number of channels present
            % in the dataset (this is only applicable to chanloc-less .mat
            % low-density data).
            if (params.low_density && params.datafileformat == 0 && ~params.has_chanlocs)
                dataQC{current_file, 2} = num2str(size(EEG.data,1)) ;
                dataQC{current_file, 3} = num2str(size(EEG.data,1)) ;
            else
                dataQC{current_file, 2} = num2str(size(params.chan_IDs, 2)) ;
                dataQC{current_file, 3} = num2str(size(selected_good_channel_locations, 2)) ;
            end
            % List the percent of good channels kept out of the number of
            % channels selected by the user.
            dataQC{current_file, 4} = str2double(dataQC(current_file, 3))/str2double(dataQC(current_file, 2)) * 100 ;
        end
        % Fill in the bad channels removed from the data. If the file has
        % no channel locations, instead use NA.
        if (params.low_density && params.datafileformat == 0 && ~params.has_chanlocs)
            dataQC{current_file, 5} = 'NA' ;
        elseif ~isempty(bad_channels_removed)
            dataQC{current_file, 5} = [sprintf('%s ', bad_channels_removed{1:end-1}), ...
                bad_channels_removed{end}] ;
        end
        % The number of segments retained post-processing. This should only
        % be different if segments were rejected.
        dataQC{current_file, 16} = num2str(EEG.trials) ;

        %% SAVE PREPROCESSED DATASET
        % Save the completely pre-processed dataset in a format specified
        % by the user. Can be in .txt, .mat, or .set format.
        disp('Saving preprocessed dataset(s)...') ;
        cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
                'processed'))}]) ;
        switch params.save_as_format
            % SAVE AS TEXT FILES: Usually used for ERP data. Prints out a
            % txt file with the Individual Trials, a .txt file with the
            % Average over Trials, a .set file for the EEG including all
            % tags, and a .set file for each ERP tag.
            case 1
                % Text file containing individual trials:
                pop_export(EEG, strrep(FileNames{current_file}, src_file_ext, ...
                    ['_processed_IndivTrial' rerun_suffix '.txt']), 'transpose', ...
                    'on', 'precision', 8);
                % Text file containing averages over trials:
                pop_export(EEG, strrep(FileNames{current_file}, src_file_ext, ...
                    ['_processed_AveOverTrials' rerun_suffix '.txt']), 'transpose', ...
                    'on', 'erp', 'on', 'precision', 8);
                % Text files for each EEG by tag (Individual and Average):
                if size(params.task_onset_tags, 2) > 1
                    for i=1:length(eeg_byTags)
                        pop_export(eeg_byTags(i), strrep(FileNames{current_file}, ...
                            src_file_ext, ['_processed_IndivTrial' ...
                            params.task_onset_tags{i} rerun_suffix '.txt']), ...
                            'transpose', 'on', 'precision', 8);
                        pop_export(eeg_byTags(i), strrep(FileNames{current_file}, ...
                            src_file_ext, ['_processed_AveOverTrials_' ...
                            params.task_onset_tags{i} rerun_suffix '.txt']), ...
                            'transpose', 'on', 'erp', 'on', 'precision', 8);
                    end
                end
                % Set file for the complete EEG
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext, ['_processed' rerun_suffix '.set'])) ;
                % Set files for each EEG by tag.
                if size(params.task_onset_tags,2) > 1
                    for i=1:length(eeg_byTags)
                        eeg_byTags(i) = pop_saveset(eeg_byTags(i), 'filename', ...
                            strrep(FileNames{current_file}, src_file_ext, ...
                            ['_' params.task_onset_tags{i} '_processed' rerun_suffix ...
                            '.set']));
                    end
                end
            
            % SAVE AS MAT FILE:    
            case 2
                save(strrep(FileNames{current_file}, src_file_ext, ...
                    ['_processed' rerun_suffix '.mat']), 'EEG') ;
            
            % SAVE AS SET FILE:    
            case 3
                EEG = pop_saveset(EEG, 'filename', strrep(FileNames{current_file}, ...
                    src_file_ext,['_processed' rerun_suffix '.set'])) ;
        end

        %% GENERATE POWER SPECTRUM AND TOPOPLOT VISUALIZATION (if enabled)
        % Plot the spectrum across channels to evaluate pipeline performance
        if params.visualizations == 1
            if params.ERP_analysis == 0
                figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', ...
                    [[params.freq_to_plot]], 'freqrange', [[params.vis_freq_min] ...
                    [params.vis_freq_max]], 'electrodes', 'off') ;
                saveas(gcf, strrep(FileNames{current_file}, src_file_ext, ...
                    ['_processedspectrum.' rerun_suffix 'jpg'])) ;
            
            elseif params.ERP_analysis == 1
                figure; pop_timtopo(EEG, [params.vis_time_start params.vis_time_end], ...
                    [params.time_to_plot], '');
                saveas(gcf, strrep(FileNames{current_file}, src_file_ext, ...
                    ['_processedERPspectrum' rerun_suffix '.jpg']));
            end
        end
        
    %% IF AN ERROR OCCURS DURING THE PROCESSING OF ANY FILE...
    % Record the errors in all QM outputs.
    catch ME
        if strcmp(ME.identifier, 'MATLAB:happeAllTrialsRej')
            fill = 'ALL_SEG_REJ' ;
        else; fill = 'ERROR' ;
        end
        
        % ERRORS IN DATA QM:
        for i=1:size(dataQC, 2)
            if isempty(dataQC{current_file, i})
                dataQC{current_file, i} = fill ;
            end
        end
        % ERRORS IN DATA QM FOR ERP:
        if params.ERP_analysis
            for i=1:size(dataQC_erp, 2)
                if isempty(dataQC_erp{current_file, i})
                    dataQC_erp{current_file, i} = fill ;
                end
            end
        end
        
        if ~reprocessing
            % ERRORS IN LINENOISE:
            ln_means = qm_ERROR(ln_means, 1 + length(lnfreqsofinterest), ...
                current_file) ;
            % ERRORS IN WAVELET THRESHOLDING:
            wav_means = qm_ERROR(wav_means, 5 + length(freqsofinterest), ...
                current_file) ;
            if ~params.low_density && ~params.ERP_analysis
                % ERRORS IN ICA:
                ica_means = qm_ERROR(ica_means, 5 + length(freqsofinterest), ...
                    current_file) ;
                % ERRORS IN ARTIFACT REJ:
                artifactRej_means = qm_ERROR(artifactRej_means, ...
                    5 + length(freqsofinterest), current_file) ;
            end
        end
    end
end

%% GENERATE OUTPUT TABLES
disp('Generating quality assessment outputs...') ;
cd([src_folder_name filesep folderNames{find(contains(folderNames, ...
    'quality_assessment_outputs'))}]) ;
rmpath(genpath(cleanline_path));
if ~reprocessing
    % CREATE VARIABLE NAMES
    % Line Noise Metrics:
    ln_names = {'r all freqs pre/post linenoise removal'};
    for i = 2:size(lnfreqsofinterest,2)+1
        ln_names{i} = ['r ' num2str(lnfreqsofinterest(i-1)) 'hz pre/post linenoise removal'] ;
    end
    % Wavelet Thresholding Metrics:
    wav_names = {'RMSE post/pre wav-threshold', 'MAE post/pre wav-threshold', ...
        'SNR post/pre wav-rejection', 'PeakSNR post/pre wav-threshold', ...
        'r alldata pre/post wav-threshold'} ;
    for i = 1:size(freqsofinterest,2)
        wav_names{i+4} = ['r ' num2str(freqsofinterest(i)) 'hz post/pre wav-threshold'] ;
    end
    if ~params.low_density && ~params.ERP_analysis
        % ICA Metrics:
        ica_names = {'RMSE post/pre ICA-rejection', 'MAE post/pre ICA-rejection', ...
            'SNR post/pre ICA-rejection', 'PeakSNR post/pre ICA-rejection', ...
            'r alldata pre/post ICA-rejection'} ;
        for i = 1:size(freqsofinterest, 2)
            ica_names{i+5} = ['r ' num2str(freqsofinterest(i)) 'hz post/pre ICA-rejection'] ;
        end
        % All Artifact Rejection Metrics:
        artifactRej_names = {'RMSE post/pre artifact rejection', ...
            'MAE post/pre artifact rejection', 'SNR post/pre artifact rejection', ...
            'PeakSNR post/pre artifact rejection', ...
            'r alldata post/pre artifact rejection'} ;
        for i = 1:size(freqsofinterest, 2)
            artifactRej_names{i+5} = ['r ' num2str(freqsofinterest(i)) ...
                'hz post/pre artifact rejection'] ;
        end
    end
        
    % CONCATENATE ARRAYS
    if ~params.low_density && ~params.ERP_analysis
        pipelineQM_names = [(ln_names) (wav_names) (ica_names) (artifactRej_names)] ;
        pipelineQM = [(ln_means) (wav_means) (ica_means) (artifactRej_means)] ;
    else
        pipelineQM_names = [(ln_names) (wav_names)] ;
        pipelineQM = [(ln_means) (wav_means)] ;
    end
    if params.ERP_analysis && size(params.task_onset_tags,2) > 1
        dataQC = [dataQC dataQC_erp] ;
        dataQC_names = [dataQC_names dataQC_erp_names] ;
    end

    % CONVERT TO TABLE; NAME VARIABLES AND ROWS
    pipelineQM = array2table(pipelineQM, 'VariableNames', pipelineQM_names, ...
        'RowNames', FileNames) ;
    dataQC = cell2table(dataQC, 'VariableNames', dataQC_names, 'RowNames', ...
        FileNames) ;
    
    % CREATE SAVE NAMES FOR QUALITY METRIC FILES
    % Data Quality Measures:
    dataQM_saveName = ['HAPPE_dataQC_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
    indx = 2 ;
    while isfile(dataQM_saveName)
        dataQM_saveName = ['HAPPE_dataQC_' datestr(now, 'dd-mm-yyyy') '_' ...
            num2str(indx) '.csv'] ;
        indx = indx + 1 ;
    end
    
    % Pipeline Quality Measures:
    indx = 2 ;
    pipelineQM_saveName = ['HAPPE_pipelineQC_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
    while isfile(pipelineQM_saveName)
        pipelineQM_saveName = ['HAPPE_pipelineQC_' datestr(now, 'dd-mm-yyyy') ...
            '_' num2str(indx) '.csv'] ;
        indx = indx + 1 ;
    end
    
    % WRITE TABLE TO .CSV
    writetable(dataQC, dataQM_saveName, 'WriteRowNames', true, 'QuoteStrings', true) ;
    writetable(pipelineQM, pipelineQM_saveName, 'WriteRowNames', true, ...
        'QuoteStrings', true);

% IF REPROCESSING...
else
    % REPLACE OLD DATAQM METRICS WITH THOSE OF CURRENT RUN
    loaded_dataQM.Number_Segments_Before_Segment_Rejection = dataQC(:, 15) ;
    loaded_dataQM.Number_Segments_Post_Segment_Rejection = dataQC(:, 16) ;
    if params.ERP_analysis && size(params.task_onset_tags,2) > 1
        dataQC_erp = cell2table(dataQC_erp) ;
        loaded_dataQM(:, (width(loaded_dataQM)+1) - width(dataQC_erp):width(loaded_dataQM)) = dataQC_erp ;
    end
    
    % CREATE DATA SAVE NAME
    dataQM_saveName = ['HAPPE_dataQC' rerun_suffix '.csv'] ;
    indx = 2 ;
    while isfile(dataQM_saveName)
        dataQM_saveName = ['HAPPE_dataQC' rerun_suffix '_' num2str(indx) '.csv'] ;
        indx = indx + 1 ;
    end
    
    % WRITE TABLE TO .CSV
    writetable(loaded_dataQM, dataQM_saveName, 'WriteRowNames', true, ...
        'QuoteStrings', true) ;
end
