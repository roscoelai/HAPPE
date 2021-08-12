% listParams() - A helper function for setParams.m and HAPPE that prints
%                out the user-entered parameters to the command window in a
%                way that is easy to read. Allows the user to check their
%                parameters prior to confirming them.
%
% Usage: 
%   >> listParams(params)
%
% Inputs:
%   params - A struct containing all of the parameters needed to run HAPPE.
%
% Outputs:
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

function listParams(params)
disp("---------------------------------------------") ;
disp("PARAMETER SETTINGS:") ;
%% Density
fprintf('Density: ') ;
if params.low_density; fprintf('Low (<= 30 channels)\n') ;
else; fprintf('High (> 30 channels)\n') ;
end

%% Rest VS Task
fprintf('Resting State or Task: ') ;
if params.task_EEG_processing
    fprintf('Task\nERP Analysis: ') ;
    if params.ERP_analysis
        fprintf('Yes\n') ;
        fprintf(' - Lowpass Cutoff: %g\n', params.ERP_lowpass_cutoff) ;
        fprintf(' - Highpass Cutoff: %g\n', params.ERP_highpass_cutoff) ;
    else; fprintf('No\n') ;
    end
else; fprintf('Resting State\n') ;
end

%% Data File Format
fprintf('Data File Format: ') ;
if params.datafileformat == 0
    fprintf('.mat (Matlab array)\n') ;
    if params.low_density
        fprintf(' - Same Sampling Rate Across Files: ');
        if params.same_srate; fprintf('Yes\n') ;
        else
            fprintf('No\n    - List of Sampling Rates: %s\n', params.srate) ;
        end
        fprintf(' - Channel Locations Provided: ') ;
        if params.has_chanlocs
            fprintf('Yes\n    - Channel Locations File: %s\n', params.chanlocs) ;
        else; fprintf('No\n') ;
        end
    else
        fprintf(' - Potential EEG Variable Names: ') ;
        listCells(params.potential_eeg_var_names) ;
        if params.low_density
            fprintf(' - Sampling Rate Variable Name: %s\n', params.sampling_rate_varname{1}) ;
        end
        if params.task_EEG_processing == 1
           fprintf('Task Event Information Location: %s\n', params.task_event_info_location) ;
        end
    end
elseif params.datafileformat == 1
    fprintf('.raw (Netstation simple binary)\n - Task Onset Tags: ') ;
    listCells(params.task_onset_tags) ;
elseif params.datafileformat == 2; fprintf('.set (EEGLab)\n') ;
elseif params.datafileformat == 3; fprintf('.cdt (Neuroscan)\n') ;
elseif params.datafileformat == 4; fprintf('.mff (EGI)\n') ;
end

%% Aquisition Layout
fprintf('Acquisition Layout: ') ;
if params.layout_type(1,1) == 1; fprintf(['%i channel EGI Geodesic Sensor ' ...
        'Net\n'], params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 2; fprintf(['%i channel EGI HydroCel ' ...
        'Geodesic Sensor Net\n'], params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 3; fprintf('%i channel BioSemi Net\n', ...
        params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 4; fprintf(['%i channel Brain Products ' ...
        'Standard BrainCap Net\n'], params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 5; fprintf(['%i channel Brain Products ' ...
        'Wet Sponge R-Net for actiCHamp Plus\n'], params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 5; fprintf(['%i channel Neuroscan SynAmps ' ...
        'Net\n'], params.layout_type(1,2)) ;
elseif params.layout_type(1,1) == 7; fprintf('Unspecified\n') ;
end

%% Channels of Interest
fprintf('Channels: ') ;
if strcmpi(params.chans_all, 'all'); fprintf('All\n') ;
elseif strcmpi(params.chans_all, 'coi_include')
    listCells(params.chan_IDs) ;
elseif strcmpi(params.chans_all, 'coi_exclude')
    fprintf(' All Except ') ;
    listCells(params.chan_IDs) ;    
end

%% Line Noise
% FREQUENCY
if isfield(params, 'line_noise')
    fprintf('Line Noise Frequency: %i Hz\n', params.line_noise) ;
end
% LEGACY
if isfield(params, 'legacy_linenoise')
    fprintf('Line Noise Reduction: ') ;
    if params.legacy_linenoise; fprintf('Legacy\n') ;
    else; fprintf('Default\n') ;
    end
end

%% Resampling
if isfield(params, 'downsample_freq')
    fprintf('Resample: ')
    if params.downsample_freq == 0; fprintf('Off\n') ;
    else; fprintf('On\n') ;
    end
end

%% Legacy Bad Channel Selection
if isfield(params, 'legacy_channels')
    fprintf('Bad Channel Detection: ') ;
    if params.badchanrej
        fprintf('On\n - Bad Channel Detection Method: ') ;
        if params.legacy_channels; fprintf('Legacy\n') ;
        else; fprintf('Default\n') ;
        end
    else; fprintf('Off\n') ;
    end  
end

%% Legacy Wavelet
if isfield(params, 'legacy_wavelet')
    fprintf('Wavelet Thresholding: ') ;
    if params.legacy_wavelet; fprintf('Legacy\n') ;
    else; fprintf('Default\n') ;
    end
end

%% Segmentation
fprintf('Segmentation: ') ;
if params.segment_data
    fprintf('On\n') ;
    if params.task_EEG_processing
        fprintf(' - Starting Parameter for Stimulus: %g seconds\n', params.task_segment_start) ;
        fprintf(' - Ending Parameter for Stimulus: %g seconds\n', params.task_segment_end) ;
        if params.ERP_analysis
            fprintf(' - Task Offset: %g milliseconds\n - Baseline Correction: ', params.task_offset) ;
            if params.baseline_correction
                fprintf('On\n    - Baseline Correction Start: %g milliseconds\n', params.baseline_corr_start) ;
                fprintf('    - Baseline Correction End: %g milliseconds\n', params.baseline_corr_end) ;
            else
                fprintf('Off\n') ;
            end
        end
    else
        fprintf(' - Segment Length: %g seconds\n', params.segment_length) ;
    end
end

%% Interpolation
fprintf('Interpolation: ') ;
if params.segment_interpolation; fprintf('On\n') ;
else; fprintf('Off\n') ;
end

%% Segment Rejection
fprintf('Segment Rejection: ') ;
if params.segment_rejection
    fprintf('On\n - Segment Rejection Method: ') ;
    if strcmpi(params.seg_rej_method, 'both'); fprintf('Both amplitude and similarity criteria\n') ;
    elseif strcmpi(params.seg_rej_method, 'amplitude'); fprintf('Amplitude criteria only\n');
    elseif strcmpi(params.seg_rej_method, 'similarity'); fprintf('Similarity criteria only\n') ;
    end
    if strcmpi(params.seg_rej_method, 'both') || ...
            strcmpi(params.seg_rej_method, 'amplitude')
        fprintf('    - Minimum Segment Rejection Threshold: %i\n', ...
            params.reject_min_amp) ;
        fprintf('    - Maximum Segment Rejection Threshold: %i\n', ...
            params.reject_max_amp) ;
        fprintf('    - Segment Rejection based on All Channels or ROI: ') ;
        if params.ROI_channels_only
            fprintf('ROI\n     - ROI Channels: ') ;
            listCells(params.ROI_channels) ;
        else
            fprintf('All\n') ;
        end
        if params.task_EEG_processing && params.datafileformat == 0
            fprintf('Restrict Analysis Using Pre-Selected Trials: ') ;
            if params.user_selected_trials; fprintf('On\n') ;
            else; fprintf('Off\n') ;
            end
        end
    end
else; fprintf('Off\n') ;
end

%% Re-Reference
fprintf('Re-Referencing: ') ;
if params.rereference_on
    fprintf('On\n - Re-Reference Method: ') ;
    if params.average_rereference
        fprintf('Average\n') ;
    else
        fprintf('To Subset - ') ;
        listCells(params.no_av_reref_channel_subset) ;
    end
end

%% Visualizations
fprintf('Visualizations: ') ;
if params.visualizations
    fprintf('On\n') ;
    if params.ERP_analysis
        fprintf(' - Start Time: %i milliseconds\n', params.vis_time_start) ;
        fprintf(' - End Time: %i milliseconds\n', params.vis_time_end) ;
        fprintf(' - Times to Plot: ');
        for i=1:length(params.time_to_plot)
            if i == length(params.time_to_plot)
                fprintf('%i\n', params.time_to_plot(i)) ;
            else
                fprintf('%i, ', params.time_to_plot(i)) ;
            end
        end
    else
        fprintf(' - Power Spectrum Minimum: %i\n', params.vis_freq_min) ;
        fprintf(' - Power Spectrum Maximum: %i\n', params.vis_freq_max) ;
        fprintf(' - Frequencies to Plot: ');
        for i=1:length(params.freq_to_plot)
            if i == length(params.freq_to_plot)
                fprintf('%i\n', params.freq_to_plot(i)) ;
            else
                fprintf('%i, ', params.freq_to_plot(i)) ;
            end
        end
    end
else
    fprintf('Off\n') ;
end

%% Save Format
fprintf('Save Format: ') ;
if params.save_as_format == 1
    fprintf('.txt file\n') ;
elseif params.save_as_format == 2
    fprintf('.mat file\n') ;
elseif params.save_as_format == 3
    fprintf('.set file\n') ;
end
disp("---------------------------------------------") ;    
end