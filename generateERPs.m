%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate ERPs - a post-processing script in association with HAPPE+ER to
% create ERP waveforms and calculate common ERP measures.
%
% Developed at Northeastern University's PINE Lab
%
% For a detailed description of this script and user options, please see 
% the following manuscript(s):
%   Monachino, et al., (----) - submitted, will update ***
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
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% SET FOLDERS FOR HAPPE AND EEGLAB PATHS
clear ;
disp('Preparing HAPPE - ERP ADD-ON...') ;

%% DETERMINE AND SET PATH TO DATA
% Use input from the command line to set the path to the data. If an 
% invalid path is entered, repeat until a valid path is entered.
while true
    src_folder_name = input('Enter the path to the folder containing the processed dataset(s):\n> ','s') ;
    if exist(src_folder_name, 'dir') == 7; break ;
    else; disp("Invalid input: please enter the complete path to the folder containing the dataset(s).") ;
    end
end
cd (src_folder_name) ;

%% COLLECT FILE NAMES
fprintf(['Enter the suffix used for this dataset, including stimulus tag' ...
    ' (if applicable).\nIf no extension beyond "AveOverTrials", press ' ...
    'enter/return.\n']) ;
suffix = input('> ', 's') ;
FileNames = {dir(['*' '_AveOverTrials' suffix '.txt']).name} ;

%% DETERMINE CHANNELS OF INTEREST
% Ask the user if they want to include entered channels or
% exclude entered channels. Adapted from the setParams function.
fprintf('Examine all channels (all) or only channels of interest (coi)?\n') ;
while true
    chans_all = input('> ', 's') ;
    if strcmpi(chans_all, 'all'); chan_IDs = {} ; break ;
   % If the user wants a subset, request user input to collect coi.
    elseif strcmpi(chans_all, "coi")
        disp("Choose an option for entering channels:") ;
        disp("  include - Include ONLY the entered channel names.") ;
        disp("  exclude - Include every channel EXCEPT the entered channel names.") ;
        while true
            chans_all = ['coi_' input('> ', 's')] ;
            if strcmpi(chans_all, 'coi_include') || strcmpi(chans_all, 'coi_exclude')
                 break ;
            else
                disp("Invalid input: please enter 'include' or 'exclude' (without quotations)") ;
            end
        end
        chan_IDs = {} ;
        indx = 1 ;
        fprintf(['Enter channels, including the preceding letter, one at a' ...
            ' time.\nPress enter/return between each entry.\nExamples: E17' ...
            '\n          M1\nWhen you have entered all channels, input ' ...
            '"done" (without quotations).\n']) ;
        while true
            u_input = input('> ', 's') ;
            if strcmpi(u_input, 'done')
                chan_IDs = unique(chan_IDs, 'stable') ;
                break ;
            else
                chan_IDs{indx} = u_input ;
                indx = indx + 1 ;
            end
        end
        break ;
    else; fprintf('Invalid input: please enter "all" or "coi" (without quotations)\n') ;    
    end
end

%% DETERMINE IF EXCLUDING BAD CHANNELS
% Ask the user whether or not to include bad/interpolated channels in the
% analyses.
fprintf(['Include bad channels in calculating ERP?\n  include = keep bad ' ...
    'channels\n  exclude = remove bad channels\n']) ;
ignoreBadChans = choose2('include', 'exclude') ;
if ignoreBadChans
    while true
        fprintf(['Enter the file containing the bad channels, including ' ...
            'the complete path.\nRefer to the HAPPE User Guide for ' ...
            'instructions on creating this file and an example.\n']) ;
        badChan_file = input('> ', 's') ;
        if isfile(badChan_file); break ;
        else; fprintf('Invalid input: please enter an existing file.') ;
        end
    end
end

%% DETERMINE IF CALCULATING VALUES
% Determine, using user input from the command line, whether calculating
% ERP values. If so, collect the latency windows, and the method(s) for
% calculating area under the curve/50% area under the curve.
fprintf(['Calculate ERP values? [Y/N]\n']) ;
calcVals = choose2('n', 'y') ;
if calcVals
    values = {'peaks', 'area under the curve', '50% area under the curve'} ;
    % COLLECT LATENCY WINDOWS
    windows = [] ;
    fprintf(['Enter latency windows of interest with anticipated peak:\n' ...
        'Enter each window as two consecutive numbers followed by "max" or ' ...
        '"min" (without quotations).\nPress Enter/Return between entries.\n' ...
        'Use latencies present in your dataset or your windows will be ' ...
        'corrected to the nearest existing values.\nExample: -100 100 max\n']) ;
    while true
        temp = split(input('> ', 's')) ;
        if length(temp) == 1 && strcmpi(temp, 'done'); break ;
        elseif size(temp, 1) == 3; temp = reshape(temp, 1, 3) ;
        else
            fprintf(['Invalid input: please enter two numbers and "max"/"min"' ...
                ' or "done" (without quotations).\n']) ;
            continue ;
        end
        if str2num(temp{1}) > str2num(temp{2})
            fprintf(['Invalid input: please make sure that the entries are consecutive.\n']) ;
        elseif str2num(temp{1}) == str2num(temp{2})
            fprintf(['Invalid input: you cannot have a latency window of 0.\n']) ;
        elseif ~strcmpi(temp{3}, 'max') && ~strcmpi(temp{3}, 'min')
            fprintf(['Invalid input: please specify "max" or "min" ' ...
                '(without quotations).\n']) ;
        else
            windows = [windows; temp] ;
        end
    end
    
    % DETERMINE METHOD FOR AUC/50% AUC
    fprintf(['Choose a method for calculating area under the curve:\n' ...
        '  windows = restrict calculations to the specified latency window\n' ...
        '  zeros = calculate area under the curve using points where the amplitude = 0\n' ...
        '  both = calculate both by windows and by zeros\n']) ;
    aucMethod = [0,0] ;
    while true
        user_input = input('> ', 's') ;
        if strcmpi(user_input, 'windows'); aucMethod(1) = 1; break;
        elseif strcmpi(user_input, 'zeros'); aucMethod(2) = 1; break ;
        elseif strcmpi(user_input, 'both'); aucMethod = [1,1] ; break ;
        else; fprintf(['Invalid input: please enter "windows", "zeros", ' ...
                'or "both" (without quotations).\n']) ;
        end
    end
    
    % CREATE TABLE TO HOLD VALUES FOR OUTPUT
    peakVals = zeros(size(FileNames,2)+1, 2*size(windows,1) + 4) ;
    allPeaks = cell(size(FileNames,2)+1, 2) ;
    if aucMethod(1)
        aucValsWindows = zeros(size(FileNames,2)+1, size(windows,1) + 1) ;
        fiftyAL = zeros(size(FileNames,2)+1, 2*size(windows,1) + 2) ;
    end
    if aucMethod(2)
        aucValsZeros = cell(size(FileNames,2)+1,1) ;
        fiftyALZeros = cell(size(FileNames,2)+1,1) ;
    end
end

%% SETTING UP THE TABLE FOR BAD CHANNELS:
if ignoreBadChans
    % Load the File
    loadedBadChans = readtable(badChan_file) ;

    % Next, convert table to cell
    badChans = table2cell(loadedBadChans) ;

    % Split the string into a cell array
    for i=1:size(badChans, 1)
        badChans{i,2} = split(badChans{i,2})' ;
    end
end

%% CREATE ERP ARRAY FOR EACH FILE
allSubsAve = [] ;
for currfile = 1:size(FileNames, 2)
     try
        % LOAD THE FILE
        currsub = importdata(FileNames{currfile}) ;
        if currfile == 1 && calcVals
            lats = currsub.data(:,1) ;
            % CORRECT LATENCIES:
            % If the requested latency boundary does not exist as a
            % timepoint within the dataset, correct to the closest latency
            % value for each user-specified window.
            for i=1:size(windows,1)
                for j=1:2
                    if ~any(ismember(lats, str2num(windows{i,j})))
                        [minVal, closestIndx] = min(abs(lats-str2num(windows{i,j}))) ;
                        windows{i,j} = num2str(lats(closestIndx)) ;
                    end
                end
            end
        end
        
        % COMPILE LIST OF BAD CHANNELS
        if ignoreBadChans
            subBadChans = badChans{find(contains({badChans{:,1}}, ... 
                strrep(FileNames{currfile}, ['_processed_AveOverTrials' ...
                suffix '.txt'], ''))), 2} ;
        end
        
        % COLLECT CHANNEL NAMES:
        % Get the channel names based on coi method and whether to exclude
        % bad channels.
        chanCols = [] ;
        if strcmpi(chans_all, 'coi_exclude')
            subChanIDs = setdiff(currsub.colheaders, chan_IDs) ;
            if ignoreBadChans; subChanIDs = setdiff(subChanIDs, subBadChans) ; end
        elseif strcmpi(chans_all, 'coi_include')
            subChanIDs = chan_IDs ;
            if ignoreBadChans; subChanIDs = setdiff(subChanIDs, intersect(subBadChans, chan_IDs)) ; end
        elseif strcmpi(chans_all, 'all')
            subChanIDs = currsub.colheaders ;
            if ignoreBadChans; subChanIDs = setdiff(subChanIDs, subBadChans) ; end
        end
        
        % GET COLUMN NUMBERS FOR COI:
        for int=1:length(subChanIDs)
            chanCols = [chanCols find(contains(currsub.colheaders, subChanIDs{int}))] ;
        end
        
        % CREATE ARRAY TO HOLD COI
        subjSet = [] ;
        for i=1:size(chanCols, 2)
            subjSet = [subjSet currsub.data(:, chanCols(i))] ;
        end
        
        % COMPILE MEANS
        currSubAve = mean(subjSet,2) ;
        allSubsAve = [allSubsAve currSubAve] ;
        currSubAve_noBL = currSubAve(find(lats==0):size(currSubAve,1),:) ;
        
        % CALCULATE VALUES
        if calcVals
            % CREATE WINDOWS
            [currSubWindows, currSubGlobal] = createWindows(lats, ...
                currSubAve_noBL, windows) ;
            
            % CALCULATE PEAKS:
            peakVals = calcPeaks(peakVals, windows, currSubWindows, ...
                currSubGlobal, currfile) ;
            [listedmaxes, listedmins] = getAllPeaks(currSubGlobal) ;
            allPeaks(currfile,1) = {listedmaxes} ;
            allPeaks(currfile,2) = {listedmins} ;
            
            % CALCULATE AREA UNDER CURVE & 50% AREA UNDER CURVE
            if aucMethod(1)
                [aucValsWindows, fiftyAL] = calcAUC(windows, currSubWindows, ...
                    aucValsWindows, fiftyAL, currSubGlobal, currfile) ;
            end
            if aucMethod(2)
                [aucValsZeros, fiftyALZeros] = calcAUCZeros(lats, currSubGlobal, ...
                    aucValsZeros, fiftyALZeros, currfile) ;
            end
        end
    catch
        fprintf("Error in file %s.\n", FileNames{currfile}) ;
    end
end
currfile = currfile + 1;
%% PLOT ERP WAVEFORMS
disp("Figure 1: All subjects' ERP without mean") ;
plot(allSubsAve) ;
figure() ;
disp("Figure 2: Mean ERP across all subjects") ;
totalMean = mean(allSubsAve, 2) ;
totalMean_noBL = totalMean(find(lats==0):size(totalMean,1),:) ;
if calcVals
    % Create the windows
    [currSubWindows, currSubGlobal] = createWindows(lats, totalMean_noBL, ...
        windows) ;
    
    % PEAKS:
    peakVals = calcPeaks(peakVals, windows, currSubWindows, currSubGlobal, ...
        currfile) ;
    [listedmaxes, listedmins] = getAllPeaks(currSubGlobal) ;
    allPeaks(currfile,1) = {listedmaxes} ;
    allPeaks(currfile,2) = {listedmins} ;

    % AREA UNDER CURVE & 50% AREA UNDER CURVE
    if aucMethod(1)
        [aucValsWindows, fiftyAL] = calcAUC(windows, currSubWindows, aucValsWindows, fiftyAL, ...
            currSubGlobal, currfile) ;
    end
    if aucMethod(2)
        [aucValsZeros, fiftyALZeros] = calcAUCZeros(lats, currSubGlobal, ...
            aucValsZeros, fiftyALZeros, currfile) ;
    end
end
plot(totalMean) ;
figure() ;
disp("Figure 3: All subjects' ERP with mean") ;
allSubsAve = [allSubsAve totalMean] ;
plot(allSubsAve) ;
disp("To view details, use the 'View > Plot Browser' function in the figure window.") ;

%% COMPILE STAT TABLES
if calcVals
    peakNames = cell(1,size(peakVals,2)) ;
    for i=1:size(windows,1)
        peakNames{i*2-1} = [windows{i,3} ' Value for Window ' windows{i,1} '-' ...
            windows{i,2}] ;
        peakNames{i*2} = ['Latency at ' windows{i,3} ' for Window ' ...
            windows{i,1} '-' windows{i,2}] ;
    end
    peakNames(size(peakNames,2)-3:size(peakNames,2)) = {'Global Max Value', ...
        'Latency at Global Max', 'Global Min Value', 'Latency at Global Min'} ;
    allPeakNames = {'All Maxes (with Values)' 'All Mins (with Values)'} ;
    
    if aucMethod(1)
        aucNames = cell(1,size(aucValsWindows,2)) ;
        for i=1:size(windows,1)
            aucNames{i} = ['Area Under the Curve for Window ' windows{i,1} '-' ...
                windows{i,2}] ;
        end
        aucNames{size(aucValsWindows,2)} = 'Global Area Under the Curve' ;

        fiftyALNames = cell(1,size(fiftyAL,2)) ;
        for i=1:size(windows,1)
            fiftyALNames{i*2-1} = ['50% Area Under Curve for Window ' ...
                windows{i,1} '-' windows{i,2}] ;
            fiftyALNames{i*2} = ['Latency at 50% Area Under Curve for Window ' ...
                windows{i,1} '-' windows{i,2}] ;
        end
        fiftyALNames{size(fiftyAL,2)-1} = 'Global 50% Area Under the Curve' ;
        fiftyALNames{size(fiftyAL,2)} = 'Latency at 50% Area Under the Curve' ;
    end
    if aucMethod(2)
        aucZerosNames = {'AUC for 0-Bound Windows'} ;
        fiftyALZerosNames = {'50% AUC for 0-Bound Windows and Associated Latencies'} ;
    end
    
    if aucMethod(1) && ~aucMethod(2)
        erpVals = [array2table(peakVals, 'VariableNames', peakNames, 'RowNames', ...
            [FileNames'; 'Average']), cell2table(allPeaks, 'VariableNames', ...
            allPeakNames, 'RowNames', [FileNames'; 'Average']), array2table(aucValsWindows, ...
            'VariableNames', aucNames, 'RowNames', [FileNames'; 'Average']), ...
            array2table(fiftyAL, 'VariableNames', fiftyALNames, 'RowNames', ...
            [FileNames'; 'Average'])] ;
    elseif aucMethod(2) && ~aucMethod(1)
        erpVals = [array2table(peakVals, 'VariableNames', peakNames, 'RowNames', ...
            [FileNames'; 'Average']), cell2table(allPeaks, 'VariableNames', ...
            allPeakNames, 'RowNames', [FileNames'; 'Average']), array2table(aucValsZeros, ...
            'VariableNames', aucZerosNames, 'RowNames', [FileNames'; 'Average']), ...
            array2table(fiftyALZeros, 'VariableNames', fiftyALZerosNames, 'RowNames', ...
            [FileNames'; 'Average'])] ;
    elseif aucMethod(1) && aucMethod(2)
        erpVals = [array2table(peakVals, 'VariableNames', peakNames, 'RowNames', ...
            [FileNames'; 'Average']), cell2table(allPeaks, 'VariableNames', ...
            allPeakNames, 'RowNames', [FileNames'; 'Average']), array2table(aucValsWindows, ...
            'VariableNames', aucNames, 'RowNames', [FileNames'; 'Average']), ...
            array2table(aucValsZeros, 'VariableNames', aucZerosNames, ...
            'RowNames', [FileNames'; 'Average']), array2table(fiftyAL, ...
            'VariableNames', fiftyALNames, 'RowNames', [FileNames'; 'Average']), ...
            array2table(fiftyALZeros, 'VariableNames', fiftyALZerosNames, ...
            'RowNames', [FileNames'; 'Average'])] ;
    end
    
    erpVals_save = ['allSubjectsERPvals_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
    indx = 2 ;
    while isfile(erpVals_save)
        erpVals_save = ['allSubjectsERPvals_' datestr(now, 'dd-mm-yyyy') '_' ...
            num2str(indx) '.csv'] ;
        indx = indx + 1 ;
    end
    
    writetable(erpVals, erpVals_save, 'WriteRowNames', true, 'QuoteStrings', ...
        true);
end
    
%% SAVE OUT AS CSV
erps_save = ['allSubjects_generatedERPs' suffix '_' datestr(now, 'dd-mm-yyyy') '.csv'] ;
indx = 2 ;
while isfile(erps_save)
    erps_save = ['allSubjectsERPvals_' suffix '_' datestr(now, 'dd-mm-yyyy') ...
        '_' num2str(indx) '.csv'] ;
    indx = indx + 1 ;
end

allSubsAve = array2table(allSubsAve, 'VariableNames', [FileNames 'Average']) ;
writetable(allSubsAve, erps_save) ;

%-------------------------------------------------------------------------%
%% FUNCTIONS USED TO CALCULATE ER
% CREATE WINDOWS OF INTEREST:
% Creates arrays containing the data points for the windows of interest as
% indicated previously by the user.
function [currSubWindows, currSubGlobal] = createWindows(lats, currSubAve_noBL, windows)
    currSubGlobal = [lats(find(lats==0):size(lats,1),:) currSubAve_noBL] ;
    currSubWindows =  {};
    for i=1:size(windows, 1)
        currSubWindows{i} = currSubGlobal(find(currSubGlobal==str2num(windows{i,1})):find(currSubGlobal==str2num(windows{i,2})),:) ;
    end
end

% CALCULATE PEAKS
% Find the peak as indicated (max or min) in each of the windows, as well
% as for the ERP waveform as a whole.
function peakVals = calcPeaks(peakVals, windows, currSubWindows, currSubGlobal, currfile)
    % Find Peak (Max or Min) in Windows
    for i=1:size(windows, 1)
        currWin = currSubWindows{i} ;
        if strcmpi(windows{i, 3}, "max")
            [currPeak, currPeakIndx] = max(currWin(:,2)) ;
        elseif strcmpi(windows{i,3}, "min")
            [currPeak, currPeakIndx] = min(currWin(:,2)) ;
        end
        peakVals(currfile, i*2-1) = currWin(currPeakIndx,2) ;
        peakVals(currfile,i*2) = currWin(currPeakIndx,1) ;
    end
    % Find Global Max and Min
    [globalMax, globalMaxIndx] = max(currSubGlobal(:,2)) ;
    [globalMin, globalMinIndx] = min(currSubGlobal(:,2)) ;
    % Store Info
    peakVals(currfile, size(peakVals,2)-3) = currSubGlobal(globalMaxIndx,2) ;
    peakVals(currfile, size(peakVals,2)-2) = currSubGlobal(globalMaxIndx,1) ;
    peakVals(currfile, size(peakVals,2)-1) = currSubGlobal(globalMinIndx,2) ;
    peakVals(currfile, size(peakVals,2)) = currSubGlobal(globalMinIndx,1) ;
end

% CALCULATE AUC AND 50% AUC USING WINDOWS
% Calculates the area under the curve and 50% area under the curve using
% the user-specified latency values as the latency boundaries. Additionally
% calculates the total area under the curve and 50% area under the curve.
function [aucVals, fiftyAL] = calcAUC(windows, currSubWindows, aucVals, fiftyAL, currSubGlobal, currfile)
    % Calculate AUC & 50% AUC
    for i=1:size(windows,1)
        currWin = currSubWindows{i} ;
        auc = cumtrapz(currWin(:,1), abs(currWin(:,2))) ;
        aucVals(currfile, i) = auc(length(auc)) ;
        [minVal, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
        fiftyAL(currfile, i*2) = currWin(closestIndx,1) ;
        fiftyAL(currfile, i*2-1) = auc(closestIndx) ;
    end
    % Calculate Global AUC and 50% AUC
    auc = cumtrapz(currSubGlobal(:,1), abs(currSubGlobal(:,2))) ;
    aucVals(currfile, size(aucVals,2)) = auc(length(auc)) ;
    [minVal, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
    fiftyAL(currfile, size(fiftyAL,2)) = currSubGlobal(closestIndx,1) ;
    fiftyAL(currfile, size(fiftyAL,2)-1) = auc(closestIndx) ;
end

% FIND ALL PEAKS IN THE DATA
function [listedmaxes, listedmins] = getAllPeaks(currSubGlobal)
    maxes = [currSubGlobal(:,1) islocalmax(currSubGlobal(:,2))] ;
    mins = [currSubGlobal(:,1) islocalmin(currSubGlobal(:,2))] ;

    listedmaxes = '' ;
    listedmins = '' ;
    for i=1:size(maxes,1)
        if maxes(i,2)
            listedmaxes = [listedmaxes num2str(maxes(i,1)) '(' ...
                num2str(round(currSubGlobal(i,2),2)) ') '] ;
        end
        if mins(i,2)
            listedmins = [listedmins num2str(mins(i,1)) '( ' ...
                num2str(round(currSubGlobal(i,2),2)) ') '] ;
        end
    end
end

% CALCULATE AUC & 50% AUC USING ZEROS
% Calculates both values by first finding zeros crossings in the dataset
% and using them as the limits for the calculations. Includes the starting
% and ending latencies as boundary-points.
function [aucValsZeros, fiftyALZeros] = calcAUCZeros(lats, currSubGlobal, aucValsZeros, fiftyALZeros, currfile)
    % Find the Zero Crossings
    zeroCrossings = [0] ;
    for i=2:size(currSubGlobal,1)
        prev = currSubGlobal(i-1,2) ;
        curr = currSubGlobal(i,2) ;
        if curr == 0
            zeroCrossings = [zeroCrossings currSubGlobal(i,1)] ;
        elseif (prev < 0 && curr > 0) || (prev > 0 && curr < 0)
            zeroCrossings = [zeroCrossings currSubGlobal(i-1,1)] ;
        end
    end
    zeroCrossings = [zeroCrossings currSubGlobal(length(currSubGlobal),1)] ;

    % Create Windows from the Zero Crossings
    zeroCrossWinds =  {} ;
    for i=2:size(zeroCrossings, 2)
        zeroCrossWinds{i-1} = currSubGlobal(find(currSubGlobal==zeroCrossings(i-1)):find(currSubGlobal==zeroCrossings(i)),:) ;
    end

    % Calculate AUC & 50% AUC for each Window
    currAUC = '' ;
    currFiftyAL = '' ;
    for i=1:size(zeroCrossWinds,2)
        currWin = zeroCrossWinds{i} ;
        if size(currWin, 1) == 1; auc = currWin(1,2) ;
        else; auc = cumtrapz(currWin(:,1), abs(currWin(:,2))) ;
        end
        [minVal, closestIndx] = min(abs(auc-auc(length(auc))/2)) ;
        currAUC = [currAUC '{' num2str(currWin(1,1)) '-' num2str(currWin(size(currWin,1), ...
            1)) ': ' num2str(auc(length(auc))) '} '] ;
        currFiftyAL = [currFiftyAL '{' num2str(currWin(1,1)) '-' ...
            num2str(currWin(size(currWin,1),1)) ': ' num2str(auc(closestIndx)) ...
            ' at ' num2str(currWin(closestIndx,1)) '} '] ;
    end
    
    % Save Outputs
    aucValsZeros{currfile} = {currAUC} ;
    fiftyALZeros{currfile} = {currFiftyAL} ;
end