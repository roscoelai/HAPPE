% calcRMSE() - Calculate Root Mean Squared Error in a dataset pre/post- any
%              particular processing step.
%
% Usage: 
%   >> RMSE = calcRMSE(preEEG, postEEG, order)
%
% Inputs:
%   preEEG  - EEG signal pre-processing step in channels x samples format
%   postEEG - EEG signal post-processing step in channels x samples format
%
% Outputs:
%   RMSE - Root Mean Squared Error over all channels and timepoints
%
% Author: Alexa D. Monachino, PINE Lab at Northeastern University, 2021
%         Laurel J. Gabard-Durnam, PINE Lab at Northeastern University, 2021
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

function RMSE = calcRMSE(preEEG, postEEG)
    RMSE = mean(realsqrt(mean((preEEG - postEEG) .^ 2, 2)));
    
%    % Step-by-step
%    differences = preEEG - postEEG;
%    squared_errors = differences .^ 2;
%    mean_squared_errors = mean(squared_errors, 2);
%    root_mean_squared_errors = realsqrt(mean_squared_errors);
%    RMSE = mean(root_mean_squared_errors);
end
