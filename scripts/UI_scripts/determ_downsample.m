% determ_downsample() - A helper function for setParams.m from HAPPE.
%                       Determines whether the user intends to downsample
%                       the data, and if so, to what frequency, through user
%                       input via the command line.
%
% Usage: 
%   >> downsample_freq = determ_downsample()
%
% Inputs:
%
% Outputs:
%   downsample_freq - An integer representing the frequency to downsample
%                     the data to. If the frequency is 0, that indicates
%                     that downsampling should not be performed.
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

function downsample_freq = determ_downsample()
    disp("Resample data? [Y/N]") ;
    disp("NOTE: Resampling is recommended for files <= 60 seconds long.") ;
    if choose2("n", "y")
%         downsample_freq = 250 ;
        disp("HAPPE supports resampling to 250, 500, and 1000.") ;
        disp("Resample frequency:") ;
        while true
            downsample_freq = input('> ') ;
            if ismember(downsample_freq, [250, 500, 1000]); break ;
            else; disp("Invalid input: please enter 250, 500, or 1000.") ;
            end
        end
    else; downsample_freq = 0 ;
    end
end