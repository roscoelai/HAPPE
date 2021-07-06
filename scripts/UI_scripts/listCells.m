% listCells() - A helper function for listParams.m from HAPPE. Prints out a
%               cell array so there is a comma and a space between each
%               element in the array, excluding following the final
%               element. Allows for neat, compact printing of a cell array.
%
% Usage: 
%   >> listCells(cellArray)
%
% Inputs:
%   cellArray - The cell array whose elements the user wishes to print out
%               to the command window. Each element should be a string or a
%               character array.
%
% Outputs:
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

function listCells(cellArray)
    for i=1:length(cellArray)
        if i == length(cellArray); fprintf('%s\n', cellArray{i}) ;
        else; fprintf('%s, ', cellArray{i}) ;
        end
    end
end