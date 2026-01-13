function demo_roi_quantification_mip()
% DEMO_ROI_QUANTIFICATION_MIP
%
% Public-release demo for ROI-based quantification on 3D abundance maps.
% The workflow is designed for quantifying signal changes in regions such as
% the spinal cord and cisterna magna from spectrally unmixed 3D volumes.
%
% This demo:
%   1) Loads 3D "abundance_map" volumes (.mat) for each time point
%   2) Computes an XY maximum-intensity projection (MIP) along Z
%   3) Lets the user draw polygon ROIs (roipoly) on the MIP (optional)
%   4) Quantifies the mean intensity within a percentile band of ROI pixels
%      (e.g., mean of pixels between top 10% and top 30%)
%   5) Saves ROI masks and quantification results for reuse
%   6) Plots normalized time courses
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 Seongwook Choi et al.
% SPDX-License-Identifier: MIT
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% -------------------------------------------------------------------------
%
% Requirements:
%   - Image Processing Toolbox (roipoly)
%
% Expected input:
%   - For each time point, a .mat file containing variable:
%       abundance_map : 3D numeric array [Nx x Ny x Nz]
%
% Output:
%   - Saved ROI masks per lesion (per mouse) in outputDir
%   - Saved quantification table (Qtf) and normalized values (Qtf_norm)
%
% Notes:
%   - Intensity is reported in arbitrary units (arb. units), because it is an
%     instrument-scaled quantity after reconstruction/unmixing.
%
% File naming convention (editable below):
%   e.g., CSF_<mouse>_unmixed_<time>_1.mat  containing abundance_map

close all;
clc;

%% -------------------- User configuration --------------------
inputDir  = fullfile(pwd, "unmixed_maps");          % <-- change
outputDir = fullfile(inputDir, "Quantification");   % <-- change if desired
if ~exist(outputDir, "dir"), mkdir(outputDir); end

mouseIDs  = {"m3"};                 % biological replicates (e.g., {"m1","m2","m3"})
timePoints = {"30m","24h"};         % must match your filenames (ordered)
lesions   = {"spinal_cord","cisterna_magna"};

% Visualization range for the MIP display (edit based on your data scale)
displayCLim = [0 10000];

% Quantification: use ROI pixels sorted descending; take mean of pixels
% between topPercentLow and topPercentHigh (exclusive/inclusive via indices)
topPercentLow  = 10;    % e.g., 10
topPercentHigh = 30;    % e.g., 30
assert(topPercentLow < topPercentHigh && topPercentHigh <= 100, "Invalid percentile settings.");

% ROI mode:
%   "draw"  : draw and save ROI masks
%   "reuse" : load previously saved ROI masks and quantify
roiMode = "draw";

% Filename pattern function
makeFileName = @(mouseID, tstr) sprintf("CSF_%s_unmixed_%s_1.mat", mouseID, tstr);

%% -------------------- Preallocate outputs --------------------
nM = numel(mouseIDs);
nT = numel(timePoints);
nL = numel(lesions);

% Qtf dimensions: [time x lesion x mouse]
Qtf = nan(nT, nL, nM);

% Store ROI masks: [time x lesion] per mouse (logical image in MIP coordinates)
ROI_masks = cell(nT, nL, nM);

%% -------------------- Main loop --------------------
for m = 1:nM
    mouseID = mouseIDs{m};

    % Load first volume to infer image sizes and coordinate convention
    firstFile = fullfile(inputDir, makeFileName(mouseID, timePoints{1}));
    S0 = load_required(firstFile, "abundance_map");
    vol0 = S0.abundance_map;
    validateattributes(vol0, {"numeric"}, {"3d","nonempty"});
    [Nx, Ny, Nz] = size(vol0); %#ok<ASGLU>

    % (Optional) If you want physical axes (mm), define step sizes here.
    % They are not required for quantification itself.
    % stepX = 0.1; stepY = 0.1;  % mm per pixel (example)

    % If reuse mode, load ROI masks saved earlier
    roiSaveFile = fullfile(outputDir, sprintf("ROI_masks_%s.mat", mouseID));
    if roiMode == "reuse"
        if ~isfile(roiSaveFile)
            error("ROI mask file not found for reuse mode: %s", roiSaveFile);
        end
        tmp = load_required(roiSaveFile, "ROI_masks");
        ROI_masks(:,:,m) = tmp.ROI_masks(:,:,m); %#ok<NASGU>
    end

    for t = 1:nT
        tstr = timePoints{t};
        inFile = fullfile(inputDir, makeFileName(mouseID, tstr));
        S = load_required(inFile, "abundance_map");
        vol = double(S.abundance_map);
        if ~isequal(size(vol), size(vol0))
            error("Volume size mismatch for %s at %s.", mouseID, tstr);
        end

        % Compute XY MIP along Z
        mipXY = max(vol, [], 3);

        % Match your original display orientation (optional but consistent)
        mipDisp = format_for_display(mipXY);

        % Display for ROI drawing
        if roiMode == "draw"
            figure("Name", sprintf("MIP - %s - %s", mouseID, tstr));
            imagesc(mipDisp);
            axis image off;
            colormap("hot");
            caxis(displayCLim);
        end

        for l = 1:nL
            lesionName = lesions{l};

            % Obtain ROI mask
            if roiMode == "draw"
                title(sprintf("%s | %s | %s (draw ROI)", mouseID, tstr, strrep(lesionName,"_"," ")));
                BW = roipoly();  % mask in displayed coordinates
                if isempty(BW)
                    error("ROI drawing cancelled or empty for %s | %s | %s.", mouseID, tstr, lesionName);
                end
                BW(isnan(BW)) = 0;
                BW = logical(BW);

                % Store mask
                ROI_masks{t,l,m} = BW;

            elseif roiMode == "reuse"
                BW = ROI_masks{t,l,m};
                if isempty(BW)
                    error("Missing ROI mask for %s | %s | %s in reuse mode.", mouseID, tstr, lesionName);
                end
            else
                error("Unknown roiMode: %s", roiMode);
            end

            % Quantify signal in ROI using percentile-band mean
            Qtf(t,l,m) = quantify_percentile_band(mipDisp, BW, topPercentLow, topPercentHigh);

            % (Optional) save ROI mask per lesion/time
            roiMaskOut = fullfile(outputDir, sprintf("ROI_%s_%s_%s.mat", mouseID, lesionName, tstr));
            save(roiMaskOut, "BW");
        end
    end

    % Save ROI masks for this mouse (to reuse later)
    if roiMode == "draw"
        save(roiSaveFile, "ROI_masks");
    end
end

%% -------------------- Save quantification outputs --------------------
outQtfFile = fullfile(outputDir, "Qtf_results.mat");
save(outQtfFile, "Qtf", "mouseIDs", "timePoints", "lesions", ...
    "topPercentLow", "topPercentHigh", "roiMode");

%% -------------------- Normalize and plot --------------------
% Normalize to (time=1, lesion=1, same mouse) by default
% You can change normalization reference as needed.
Qtf_norm = nan(size(Qtf));
for m = 1:nM
    ref = Qtf(1,1,m);
    if ~isfinite(ref) || ref == 0
        error("Invalid normalization reference for mouse %s: %g", mouseIDs{m}, ref);
    end
    Qtf_norm(:,:,m) = Qtf(:,:,m) ./ ref;
end

save(fullfile(outputDir, "Qtf_normalized.mat"), "Qtf_norm", "Qtf", "mouseIDs", "timePoints", "lesions");

% Plot (each lesion as a separate line; show individual mice)
figure("Name", "ROI quantification (normalized)");
hold on;
x = 1:nT;

for l = 1:nL
    y = squeeze(Qtf_norm(:,l,:)); % [time x mouse]
    plot(x, y, "-o", "LineWidth", 1.5, "MarkerSize", 4); %#ok<*PLOT>
end

xticks(x);
xticklabels(timePoints);
xlabel("Time");
ylabel("Normalized PA amplitude (arb. units)");
legend(strrep(lesions,"_"," "), "Location", "best");
set(gca, "FontSize", 12, "FontWeight", "bold");
grid on;
title(sprintf("ROI quantification (mean of pixels between top %d%% and %d%%)", topPercentLow, topPercentHigh));

end % demo_roi_quantification_mip


% =====================================================================
% Helper functions
% =====================================================================

function S = load_required(filePath, varName)
% LOAD_REQUIRED loads a .mat file and verifies a required variable exists.
if ~isfile(filePath)
    error("Missing file: %s", filePath);
end
S = load(filePath);
if ~isfield(S, varName)
    error("File %s does not contain required variable '%s'.", filePath, varName);
end
end


function v = quantify_percentile_band(img, mask, pLow, pHigh)
% QUANTIFY_PERCENTILE_BAND
% Sort ROI pixel intensities in descending order, then average a percentile band:
%   mean( pixels in (pLow, pHigh] of sorted intensities )
%
% Inputs:
%   img  : 2D image (display-oriented MIP)
%   mask : 2D logical mask in same coordinates as img
%   pLow : lower percentile (e.g., 10)
%   pHigh: upper percentile (e.g., 30)

mask = logical(mask);
roiPixels = img(mask);
roiPixels = roiPixels(isfinite(roiPixels));

if isempty(roiPixels)
    error("ROI contains no valid pixels.");
end

roiPixels = sort(roiPixels, "descend");
n = numel(roiPixels);

% Define indices corresponding to percentile band
i1 = floor(n * (pLow/100)) + 1;   % first included index
i2 = ceil(n * (pHigh/100));      % last included index
i1 = max(1, min(i1, n));
i2 = max(1, min(i2, n));

if i2 < i1
    error("Percentile index range invalid: i1=%d, i2=%d, n=%d", i1, i2, n);
end

v = mean(roiPixels(i1:i2));
end


function imgOut = format_for_display(imgIn)
% FORMAT_FOR_DISPLAY
% Keeps compatibility with the orientation used in the original script:
% imagesc(fliplr(rot90(max(vol,[],3),-1)))
imgOut = fliplr(rot90(imgIn, -1));
end
