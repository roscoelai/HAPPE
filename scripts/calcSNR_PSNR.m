% calcSNR_PSNR() - Calculate the Signal to Noise Ratio and Peak Signal to 
%                  Noise Ratio in a dataset pre/post- any particular 
%                  processing step.
%
% Usage: 
%   >> [SNR_alldata, PeakSNR_alldata] = calcSNR_PSNR(preEEG, postEEG, order)
%
% Inputs:
%   preEEG  - EEG signal pre-processing step in channels x samples format
%   postEEG - EEG signal post-processing step in channels x samples format
%   order   - The order in which to subtract the pre and post EEGs.
%                 {For 1, preEEG - postEEG; For 2, postEEG - preEEG} 
%
% Outputs:
%   SNR     - Signal to noise ratio over all channels and timepoints
%   PeakSNR - Peak signal to noise ratio over all channels and
%                     timepoints
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

function [SNR, PeakSNR] = calcSNR_PSNR(preEEG, postEEG, order)
    % Set up needed variables:
    NUM = zeros(size(preEEG,1),1) ;
    DEN = zeros(size(preEEG,1),1) ;
    SNR = zeros(size(preEEG,1),1) ;
    PeakSNR = zeros(size(preEEG, 1), 1) ;

    % Calculate mean squared error per channel as intermediate variable:
    MSE = zeros(size(preEEG, 1), 1);
    for j=1:size(preEEG, 1)
        mse = 0;
        for i=1:length(preEEG)
            if order == 1; mse = mse + (preEEG(j,i) - postEEG(j,i))^2 ;
            elseif order == 2; mse = mse + (postEEG(j,i) - preEEG(j,i))^2 ;
            end
        end
        [mse] = mse/length(preEEG) ;
        MSE(j) = [mse] ;
    end

    % Complete SNR and PSNR calculations:
    for j = 1:size(preEEG, 1)
        num = 0 ;
        den = 0 ;
        for i=1:length(preEEG)
            if order == 1; den = den + (preEEG(j,i) - postEEG(j,i))^2 ;
            elseif order == 2; den = den + (postEEG(j,i) - preEEG(j,i))^2 ;
            end
        end
        DEN(j) = [den] ;
        for i = 1:length(preEEG)
            if order == 1; num = num + preEEG(j,i)^2;
            elseif order == 2; num = num + postEEG(j,i)^2 ;
            end
        end
        NUM(j) = [num];
        SNR(j) = 20*log10(sqrt(NUM(j))/sqrt(DEN(j)));
        if order == 1; PeakSNR(j) = 20*log10(max(preEEG(j,:))/sqrt(MSE(j))) ;
        elseif order == 2; PeakSNR(j) = 20*log10(max(postEEG(j,:))/sqrt(MSE(j)));
        end
    end

    % SNR over all channels and timepoints:
    SNR = mean(SNR);

    % PSNR over all channels and timepoints:
    PeakSNR = mean(PeakSNR);
end