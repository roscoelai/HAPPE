% determ_fileFormat() - A helper function for setParams.m from HAPPE.
%                       Determines the file format of the data to be loaded
%                       in HAPPE through user input via the command line.
%                       Does not throw errors for or accept invalid input.
%
% Usage: 
%   >> datafileformat = determ_fileFormat()
%
% Inputs:
%
% Outputs:
%   datafileformat - An integer (0-4) representing the file format of the 
%                    data, with each number symbolizing a different file
%                    format.
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

function datafileformat = determ_fileFormat()
    fprintf(['File Format:\n  0 = .mat (MATLAB array)\n  1 = .raw' ...
        ' (Netstation simple binary)\n  2 = .set (EEGLAB format)\n' ...
        '  3 = .cdt (Neuroscan)\n  4 = .mff (EGI)\n']) ;
    while true
        datafileformat = input('> ') ;
        if ismember(datafileformat, [0:4]); break ;
        else; disp("Invalid input: please enter 0, 1, 2, 3, or 4.") ;
        end
    end
end