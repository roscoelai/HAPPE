% setParams() - A helper function for HAPPE that prints
%                out the user-entered parameters to the command window in a
%                way that is easy to read. Allows the user to check their
%                parameters prior to confirming them.
%
% Usage: 
%   >> params = setParams(params, pre_exist, reprocessing, changed_params)
%
% Inputs:
%   params         - A struct containing all of the parameters needed to 
%                    run HAPPE.
%   pre_exist      - A flag indicating whether or not these are pre-existing
%                    paramaters (0 for no, 1 for yes).
%   reprocessing   - A flag indicating whether or not the data is being
%                    reprocessed (0 for no, 1 for yes).
%   changed_params - A flag indicating whether or not the parameters are
%                    being changed.
%
% Outputs:
%   params         - A struct containing all the parameters needed to run
%                    HAPPE
%
% NOTE: Changes to the way any parameters are encoded, or changes in the
% number of parameters (adding/subtracting) will result in this code
% needing to be updated to reflect those changes. There is no easy way to
% make it reliant on other scripts... Sorry :(
%
% Author: Alexa D. Monachino, PINE Lab at Northeastern University, 2021
%
% This file is part of HAPPE.
%
% HAPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% HAPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with HAPPE. If not, see <https://www.gnu.org/licenses/>.

function params = setParams(params, pre_exist, reprocessing, changed_params)
if reprocessing
    changeMessage = ['segmentation, interpolation, segment rejection, re-referencing, save format.\n' ...
        '---------------------------------------------------------------------\n' ...
        'If you do not see an option, quit (Ctrl+C or Cmd+.) and re-run HAPPE.\n' ...
        'Enter "done" (without quotations) when finished changing parameters.\n'] ;
else
    changeMessage = ['data file format, acquisition layout, channels of interest,\n' ...
        'bad channel detection, line noise frequency, line noise reduction, visualizations,\n' ...
        'resampling, wavelet thresholding, segmentation, interpolation, segment rejection,\n' ...
        're-referencing, save format.\n' ...
        '--------------------------------------------------------------------------------\n' ...
        'If you do not see an option, quit (Ctrl+C or Cmd+.) and re-run HAPPE.\n' ...
        'Enter "done" (without quotations) when finished changing parameters.\n'] ;
end
param_choice = 'na' ;
while true
    %% BREAK IF NOT CHANGING ANY PRE-EXISTING PARAMETERS
    if pre_exist && ~changed_params; break ; end
    
    %% LIST CHANGEABLE PARAMETERS:
    if changed_params
        fprintf('Parameter to change: ') ;
        fprintf(changeMessage) ;
        param_choice = input('> ', 's') ;
    end
    
    if ~pre_exist
        %% DENSITY:
        disp('Low-Density data? [Y/N]') ;
        disp('For HAPPE, low-density data contains 30 or less channels.') ;
        params.low_density = choose2("N", "Y") ;

        %% REST VS TASK:
        [params.task_EEG_processing, params.ERP_analysis] = determ_rVt() ;
        if params.ERP_analysis
            fprintf(['Enter the task onset tags, one at a time, pressing ' ...
                'enter/return between each entry.\nWhen you have entered ' ...
                'all tags, input "done" (without quotations).\n']) ;
            params.task_onset_tags = userInput_cellArray(1,{}) ;
            params.ERP_lowpass_cutoff = input(['Enter low-pass filter, in Hz:\n' ...
                'Common low-pass filter is 30 - 45 Hz\n> ']) ;
            params.ERP_highpass_cutoff = input(['Enter high-pass filter, in Hz:\n' ...
                'Common high-pass filter is 0.1 - 0.3 Hz.\n> ']) ;
        end
    end

    %% DATA FILE FORMAT AND AQUISITION LAYOUT
    if ~pre_exist || (strcmpi(param_choice, 'data file format') && ~reprocessing) ...
            || (strcmpi(param_choice, 'acquisition layout') && ~reprocessing)
        params.datafileformat = determ_fileFormat() ;
        if params.low_density
            params.layout_type = [0, 0] ;
            %% LOW DENSITY .MAT
            if params.datafileformat == 0
                % Specify Sampling Rate:
                fprintf(['Do all your files share the same sampling rate? ' ...
                    '[Y/N]\n']) ;
                params.same_srate = choose2('n','y') ;
                if params.same_srate
                    fprintf('Sampling rate:\n') ;
                    while true
                        params.srate = input(['> ']) ;
                        if isnumeric(params.srate); break;
                        else; disp('Invalid input: please enter a real number.') ;
                        end
                    end
                else
                    fprintf(['Enter the name of the file containing the ' ...
                        'sampling rates for each file, including the path' ...
                        ' and file extension:\nSee the HAPPE user guide for' ...
                        ' how this file should be formatted.\n']) ;
                    while true
                        params.srate = input('> ', 's') ;
                        if isfile(params.srate); break;
                        else; disp('Invalid input: please enter an existing file.');
                        end
                    end
                end

                % Get Chanlocs, if available:
                fprintf(['Do you have a channel locations file for your ' ...
                    'data? [Y/N]\nNOTE: A list of supported files can be ' ...
                    'found in the HAPPE user guide.\n']) ;
                params.has_chanlocs = choose2('n', 'y') ;
                if params.has_chanlocs
                    params.chanlocs = input(['Enter the name of the file ' ...
                        'containing the chanlocs, including the full path and ' ...
                        'file extension.\n> '], 's') ;
                else; params.chanlocs = [] ;
                end

                % If task, get task event info
                if params.task_EEG_processing
                    params.task_event_info_location = matFormat_Task() ;
                end

            elseif params.datafileformat == 2
            else; error(['HAPPE currently does not support low-density data in ' ...
                    'this file format.']) ;
            end

        %% HIGH DENSITY
        else
            fprintf(['Acquisition layout net type:\n  1 = EGI Geodesic Sensor ' ...
                'Net\n  2 = EGI HydroCel Geodesic Sensor Net\n  3 = BioSemi\n' ...
                '  4 = Brain Products Standard BrainCap (BC)\n  5 = Brain ' ...
                'Products Wet-Sponge R-Net for actiCHamp Plus (RNP-AC)\n  6' ...
                ' = Neuroscan Quik-Cap\n  7 = Other\n']) ;
            while true
                params.layout_type(1,1) = input('> ') ;
                if ismember(params.layout_type(1,1), [1:7]); break;
                else; fprintf(['Invalid input: please enter 1, 2, 3, 4, 5, 6 or 7']) ;
                end
            end

            % .MAT FILES:
            if params.datafileformat == 0
                % List files that are supported - need to have chanlocs files
                % imbedded in HAPPE
                if params.layout_type(1,1) == 1; fprintf(['For EGI GSN nets,' ...
                        ' HAPPE supports 64 channels in .mat format.\n']) ;
                elseif params.layout_type(1,1) == 2; fprintf(['For EGI HydroCel ' ...
                        'GSN nets, HAPPE supports 32, 64, 128, and 256 ' ...
                        'channels in .mat format.\n']) ;
                elseif params.layout_type(1,1) == 3; fprintf(['For BioSemi ' ...
                        'nets, HAPPE supports 32, and 128 channels in .mat ' ...
                        'format.\n']) ;
                % For any other net, must enter a chanlocs file
                elseif ismember(params.layout_type(1,1), [4:7])
                    fprintf(['Do you have a channel locations file for your ' ...
                        'data? [Y/N]\nNOTE: A list of supported files can be ' ...
                        'found in the HAPPE user guide.\n']) ;
                    params.has_chanlocs = choose2('n', 'y') ;
                    if params.has_chanlocs
                        while true
                            params.chanlocs = input(['Enter the name of the file ' ...
                                'containing the chanlocs, including the full path and ' ...
                                'file extension:\n> '], 's') ;
                            if isfile(params.chanlocs); break;
                            else; fprintf(['Invalid input: please enter the file ' ...
                                    'containing the channel locations.\n']) ;
                            end
                        end
                    else; error(['HAPPE requires channel locations for high-density' ...
                            ' .mat file formats.']) ;
                    end
                else; error(['HAPPE does not support this layout for .mat files. For ' ...
                        'other nets, run through BEAPP, as described in the HAPPE ' ...
                        'manual.']) ;
                end
                % Collect the number of channels
                params.layout_type(1,2) = input('Number of channels:\n> ') ;
                if (params.layout_type(1,1) == 1 && params.layout_type(1,2) ~= 64) ...
                        || (params.layout_type(1,1) == 2 && ...
                        ~ismember(params.layout_type(1,2), [32,64,128,256])) ...
                        || (params.layout_type(1,1) == 3 && ...
                        ismember(params.layout_type(1,2), [32,128]))
                    error(['The entered number of channels is not supported for this' ...
                        ' net in a .mat file.']) ;
                else
                    fprintf(['Enter the potential EEG variable names, one at a ' ...
                        'time.\nPress enter/return between each entry.\nNOTE: ' ...
                        'variable names containing "segment" may cause issues.\n']) ;
                    params.potential_eeg_var_names = userInput_cellArray(1, {}) ;
                end

            % .RAW files
            elseif params.datafileformat == 1
                if params.layout_type(1,1) == 1; fprintf(['For EGI GSN nets,' ...
                        'HAPPE supports 64 channels in .raw format.\n']) ;
                elseif params.layout_type(1,1) == 2; fprintf(['For EGI ' ...
                        'HydroCel GSN nets, HAPPE supports 32, 64, 128, ' ...
                        'and 256 channels in .raw format.\n']) ;
                else; error(['HAPPE does not support this layout for .raw files. ' ...
                        'For other nets, run through BEAPP, as described in the ' ...
                        'HAPPE manual.']) ;
                end
                params.layout_type(1,2) = input('Number of channels: \n> ') ;
                if (params.layout_type(1,1) ~= 1 && params.layout_type(1,2) ~= 64) && ...
                        (params.layout_type(1,1) ~= 1 && ~ismember(params.layout_type(1,2), ...
                        [32, 64, 128, 256]))
                    error(['The entered number of channels is not supported for this' ...
                        ' net as a .raw file.']) ;
                end

            % .SET files    
            elseif params.datafileformat == 2
                % probably needs to contain 10-20 channels
                fprintf(['Does your file have the 10-20 channels labeled? [Y/N]\n']) ;
                if ~choose2('n','y'); error(['For .set files, the 10-20 channels ' ...
                        'must be labeled.']) ;
                end
                fprintf('Number of channels: \n') ;
                params.layout_type(1,2) = input('> ') ;

            % .CDT files
            elseif params.datafileformat == 3
                if params.layout_type(1,1) ~= 6; error(['To run a .cdt file, the net ' ...
                        'must be a NeuroScan layout.']) ;
                else
                    fprintf('Number of channels: \n') ;
                    params.layout_type(1,2) = input(['For Neuroscan Quik-Caps, ' ...
                        'HAPPE supports 32, 64, and 128 channels.\n> ']) ;
                end

            % .MFF files    
            elseif params.datafileformat == 4

            end
        end
    end
    
    noChans = params.low_density && params.datafileformat == 0 && ...
        ~params.has_chanlocs ;

    %% CHANNELS OF INTEREST
    if ~pre_exist || (strcmpi(param_choice, 'channels of interest') && ~reprocessing)
        if noChans
            params.chan_IDs = {} ;
            params.chans_all = 'all' ;
        else; [params.chans_all, params.chan_IDs] = determ_chanIDs(params.low_density) ;
        end
    end

    %% BAD CHANNEL DETECTION
    if (~pre_exist || strcmpi(param_choice, 'bad channel detection')) && ~reprocessing
        if noChans; params.badchanrej = 0;
        else
            fprintf(['Perform bad channel detection? [Y/N]\n']) ;
            params.badchanrej = choose2('n','y') ;
            if params.badchanrej && ~params.low_density
                disp("Bad channel detection method:") ;
                disp("  default = Default method optimized in HAPPE v2.") ;
                disp("  legacy = Method from HAPPE v1 (NOT RECOMMENDED).") ;
                params.legacy_channels = choose2("default", "legacy") ;
            else
                params.legacy_channels = 0 ;
            end
        end
    end

    %% LINE NOISE FREQUENCY AND METHOD
    if (~pre_exist || strcmpi(param_choice, 'line noise frequency')) && ~reprocessing
        params.line_noise = input(['Frequency of electrical (line) noise in Hz:\n' ...
            'USA data probably = 60; Otherwise, probably = 50\n' ...
            '> ']) ; 
    end

    if (~pre_exist || strcmpi(param_choice, 'line noise reduction')) && ~reprocessing
        disp("Line noise reduction method:") ;
        disp("  default = Default method optimized in HAPPE v2.") ;
        disp("  legacy = Method from HAPPE v1 (NOT RECOMMENDED).") ;
        params.legacy_linenoise = choose2("default", "legacy") ;
    end

    %% VISUALIZATIONS
    if ~pre_exist || (strcmpi(param_choice, 'visualizations') && ~reprocessing)
        disp("Run HAPPE with visualizations? [Y/N]") ;
        params.visualizations = choose2("N", "Y") ;
        if params.visualizations
            % POWER SPECTRUM:
            % Min and Max
            params.vis_freq_min = input("Minimum value for power spectrum figure:\n> ") ;
            params.vis_freq_max = input("Maximum value for power spectrum figure:\n> ") ;

            % Frequencies for spatial topoplots
            disp("Enter the frequencies, one at a time, to generate spatial topoplots for:") ;
            disp("When you have entered all frequencies, input 'done' (without quotations).") ;
            indx = 1 ;
            while true
                user_input = input('> ', 's') ;
                if strcmpi(user_input, 'done')
                    params.freq_to_plot = unique(params.freq_to_plot, 'stable') ;
                    break ;
                else
                    params.freq_to_plot(indx) = str2num(user_input) ;
                    indx = indx + 1 ;
                end
            end
            
            if params.ERP_analysis
                % DETERMINE TIME RANGE FOR THE TIMESERIES FIGURE        
                params.vis_time_start = input('Start time, in MILLISECONDS, for the ERP timeseries figure:\n> ') ;
                params.vis_time_end = input(['End time, in MILLISECONDS, for the ERP timeseries figure:\n' ...
                    'NOTE: This should end 1 millisecond before your segmentation parameter ends. (e.g. 299 for 300)\n' ...
                    '> ']) ;
                
                % Frequencies for spatial topoplots
                disp("Enter the latencies, one at a time, to generate spatial topoplots for:") ;
                disp("When you have entered all latencies, input 'done' (without quotations).") ;
                indx = 1 ;
                while true
                    user_input = input('> ', 's') ;
                    if strcmpi(user_input, 'done')
                        params.time_to_plot = unique(params.time_to_plot, 'stable') ;
                        break ;
                    else
                        params.time_to_plot(indx) = str2num(user_input) ;
                        indx = indx + 1 ;
                    end
                end
            end
        end
    end

    %% RESAMPLE
    if (~pre_exist || strcmpi(param_choice, 'resampling')) && ~reprocessing
        params.downsample_freq = determ_downsample() ;
    end

    %% WAVELET METHODOLOGY
    if (~pre_exist || strcmpi(param_choice, 'legacy wavelet')) && ~reprocessing
        disp("Method of wavelet thresholding:") ;
        disp("  default = Default method optimized in HAPPE v2.") ;
        disp("  legacy = Method from HAPPE v1 (NOT RECOMMENDED).") ;
        params.legacy_wavelet = choose2("default", "legacy") ;
    end
    
    %% SEGMENTATION
    if ~pre_exist || strcmpi(param_choice, 'segmentation')
        disp("Segment data? [Y/N]") ;
        params.segment_data = choose2("N", "Y") ;
        if params.segment_data
            if params.task_EEG_processing
                % SET SEGMENT START AND END
                params.task_segment_start = input(['Segment start, in MILLISECONDS, ' ...
                    'relative to stimulus onset:\nExample: -500\n> '])/1000 ;
                params.task_segment_end = input(['Segment end, in MILLISECONDS, ' ...
                    'relative to stimulus onset:\n> '])/1000 ;
                if params.ERP_analysis
                    % DETERMINE TASK OFFSET
                    % *** For this, maybe make it possible to upload a list
                    % of offset delays?
                    params.task_offset = input(['Offset delay, in MILLISECONDS, ' ...
                        'between stimulus initiation and presentation:\n' ...
                        'NOTE: Please enter the total offset (combined system' ...
                        ' and task-specific offsets).\n' ...
                        '> ']) ;
                    % DETERMINE IF WANT BASELINE CORRECTION
                    disp("Perform baseline correction (by subtraction)? [Y/N]") ;
                    params.baseline_correction = choose2("n", "y") ;
                    if params.baseline_correction
                        % DETERMINE BASELINE START AND END
                        params.baseline_corr_start = input(['Enter, in MILLISECONDS,' ...
                            ' where the baseline segment begins:\nExample: -100\n> ']) ;
                        params.baseline_corr_end = input(['Enter, in MILLISECONDS,' ...
                            ' where the baseline segment ends:\n' ...
                            'NOTE: 0 indicates stimulus onset.\n> ']) ;
                    end
                end
            % DETERMINE SEGMENT LENGTH
            elseif ~params.task_EEG_processing; params.segment_length = ...
                    input("Segment length, in SECONDS:\n> ") ;
            end
        end
    end
    
    %% INTERPOLATION
    if ~pre_exist || strcmpi(param_choice, 'interpolation')
        if noChans || ~params.badchanrej; params.segment_interpolation = 0 ;
        else
            disp("Interpolate the specific channels' data determined to be artifact/bad within each segment? [Y/N]") ;
            params.segment_interpolation = choose2("n", "y") ;
        end
    end
    
    %% SEGMENT REJECTION
    if ~pre_exist || strcmpi(param_choice, 'segment rejection')
        disp("Perform segment rejection? [Y/N]") ;
        params.segment_rejection = choose2("n", "y") ;
        if params.segment_rejection
            disp("Choose a method of segment rejection: ") ;
            disp("  amplitude = Amplitude criteria only") ;
            disp("  similarity = Segment similarity only") ;
            disp("  both = Both amplitude criteria and segment similarity") ;
            while true
                params.seg_rej_method = input('> ','s') ;
                if strcmpi(params.seg_rej_method, 'amplitude') || ...
                        strcmpi(params.seg_rej_method, 'both')
                    params.reject_min_amp = input(['Minimum signal amplitude' ...
                        ' to use as the artifact threshold:\n> ']) ;
                    params.reject_max_amp = input(['Maximum signal amplitude' ...
                        'to use as the artifact threshold:\n> ']) ;
                    break ;
                elseif strcmpi(params.seg_rej_method, 'similarity'); break ;
                else
                    disp("Invalid input: please enter 'amplitude', 'similarity', or 'both' (without quotations)") ;
                end
            end
            disp("Use all channels or a region of interest for segment rejection?") ;
            disp("  all = all channels") ;
            disp("  roi = region of interest") ;
            params.ROI_channels_only = choose2("all", "roi") ;
            if params.ROI_channels_only
                disp("Enter the channels in the ROI, one at a time.") ;
                disp("When you have finished entering all channels, enter 'done' (without quotations).") ;
                params.ROI_channels = userInput_cellArray(1,{}) ;
            end
            if params.task_EEG_processing && params.datafileformat == 0
                disp("Use pre-selected 'usable' trials to restrict analysis? [Y/N]") ;
                params.user_selected_trials = choose2("n", "y") ;
            end
        end
    end
    
    %% RE-REFERENCING
    if ~pre_exist || strcmpi(param_choice, 're-referencing')
        if noChans; params.rereference_on = 0;
        else
            disp("Re-reference data? [Y/N]") ;
            params.rereference_on = choose2("n", "y") ;
            if params.rereference_on
                disp("Re-Referencing Type:") ;
                disp("  subset = Re-referencing to another channel/subset of channels") ;
                disp("  average = Average re-referencing") ;
                params.average_rereference = choose2("subset", "average") ;
                if ~params.average_rereference
                    disp("Enter channel/subset of channels to re-reference to, one at a time.") ;
                    disp("When you have entered all channels, input 'done' (without quotations).") ;
                    params.no_av_reref_channel_subset = userInput_cellArray(1,{}) ;
                end
            end
        end
    end
    
    %% SAVE FORMAT
    if ~pre_exist || strcmpi(param_choice, 'save format')
        params.save_as_format = determ_saveFormat(params.datafileformat, ...
            params.low_density) ;
    end
    
   %% DONE
   if ~pre_exist || strcmpi(param_choice, 'done')
        disp("Please check your parameters before continuing.") ;
        listParams(params) ;
        disp("Are the above parameters correct? [Y/N]") ;
        if choose2("n","y"); break ;
        elseif ~pre_exist
            changed_params = 1 ;
            pre_exist = 1 ;
        end
   end   
end
end