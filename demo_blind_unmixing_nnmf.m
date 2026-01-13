function demo_blind_unmixing_nnmf()
% DEMO_BLIND_UNMIXING_NNMF
%
% Public-release demo for NNMF-based blind spectral unmixing on multispectral
% photoacoustic (PA) volumes.
%
% This demo:
%   1) Loads Hb spectra and laser energy calibration from .xlsx files
%   2) Loads multispectral 3D PA volumes from .mat files (one file per wavelength)
%   3) Energy-normalizes each wavelength
%   4) Runs NNMF to estimate endmember spectra and abundance maps
%   5) Computes abundance maps via nonnegative least squares (NNLS)
%   6) Visualizes maximum intensity projections (MIPs) of abundance volumes
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
%   - Statistics and Machine Learning Toolbox: nnmf
%   - Optimization Toolbox: lsqnonneg
%
% Notes:
%   - HbSpectra.xlsx is loaded for completeness/documentation. In this blind
%     unmixing demo, Hb spectra are not directly used (endmembers are learned
%     from data). You may use HbSpectra to identify which estimated endmember
%     corresponds to HbO/HbR/ICG in a post-hoc step if desired.
%
% File convention:
%   Each wavelength volume is stored as:
%     <baseName><wavelength>.mat
%   and must contain a 3D numeric array variable named "image3D_tot".
%
% Example:
%   baseName = "WT_m1_1h_";
%   files: WT_m1_1h_756.mat, WT_m1_1h_780.mat, WT_m1_1h_796.mat, WT_m1_1h_900.mat

%% -------------------- User configuration --------------------
wavelengths = [756 780 796 900];   % nm (order must match your filenames)
numEndmembers = 3;

% Input folder containing .mat volumes
inputDir  = fullfile(pwd, "example_data");   % <-- change to your folder
baseName  = "WT_m1_1h_";                     % <-- change to your prefix

% Calibration spreadsheets (ship these with the repository)
hbSpectraXlsx   = fullfile(pwd, "HbSpectra.xlsx");
laserEnergyXlsx = fullfile(pwd, "LaserEnergy.xlsx");

% NNMF options
optsNnmf.replicates = 5;     % increase for stability (e.g., 10)
optsNnmf.maxIter    = 500;   % adjust if needed
optsNnmf.seed       = 1;     % for reproducibility

% Abundance quantification options
optsAbundance.normalizeEndmembers = true; % matches your original normalization

%% -------------------- Load calibration tables --------------------
hbSpectra = readmatrix(hbSpectraXlsx); %#ok<NASGU>
% hbSpectra is currently not used in the blind unmixing computation.
% Keep it for transparency and optional endmember identification.

energyTable = readmatrix(laserEnergyXlsx);
% Expected format: [wavelength_nm, energy_value]
assert(size(energyTable,2) >= 2, "LaserEnergy.xlsx must have at least 2 columns: [wavelength, energy].");

%% -------------------- Compute energy at requested wavelengths --------------------
energyAtWvl = interp1(energyTable(:,1), energyTable(:,2), wavelengths, "linear", "extrap");
energyAtWvl = energyAtWvl(:)';  % row vector
if any(~isfinite(energyAtWvl)) || any(energyAtWvl <= 0)
    error("Invalid energy interpolation result. Check LaserEnergy.xlsx coverage and values.");
end

%% -------------------- Load multispectral data cube --------------------
% data: [numVoxels x numWavelengths]
[data, volSize] = load_multispectral_volumes(inputDir, baseName, wavelengths);

% Energy normalization per wavelength
for i = 1:numel(wavelengths)
    data(:,i) = double(data(:,i)) ./ energyAtWvl(i);
end

%% -------------------- Blind unmixing via NNMF --------------------
% NNMF factorization: data ≈ W * H
% where:
%   W: [numVoxels x numEndmembers]  -> abundance-like (not final)
%   H: [numEndmembers x numWavelengths] -> endmember spectra (rows)
%
% We expose endmembers as: E = H' [numWavelengths x numEndmembers]
rng(optsNnmf.seed, "twister");
nnmfOpts = statset("MaxIter", optsNnmf.maxIter, "Display", "final");

[W, H] = nnmf(max(data,0), numEndmembers, ...
    "replicates", optsNnmf.replicates, ...
    "options", nnmfOpts);

E = H';  % [numWavelengths x numEndmembers]

%% -------------------- Abundance map via NNLS --------------------
A = calculate_abundance_nnls(data, E, optsAbundance); % [numVoxels x numEndmembers]

%% -------------------- Visualize abundance MIPs --------------------
Nx = volSize(1); Ny = volSize(2); Nz = volSize(3);

for k = 1:numEndmembers
    abundanceVol = reshape(A(:,k), [Nx, Ny, Nz]);

    mip = max(abundanceVol, [], 3); % XY MIP
    figure("Name", sprintf("Endmember %d abundance (XY MIP)", k));
    imagesc(format_for_display(mip));
    axis image off; colormap(gray);
    title(sprintf("Abundance of Endmember %d (XY MIP)", k));
end

disp("Done. Endmember spectra matrix E is [numWavelengths x numEndmembers].");
disp("You can compare columns of E with HbSpectra/ICG spectra to label endmembers.");

end % demo_blind_unmixing_nnmf


% =====================================================================
% Helper functions
% =====================================================================

function [data, volSize] = load_multispectral_volumes(inputDir, baseName, wavelengths)
% LOAD_MULTISPECTRAL_VOLUMES
% Loads one 3D volume per wavelength and stacks into a 2D data matrix.
%
% Inputs:
%   inputDir   : folder containing .mat files
%   baseName   : filename prefix (string/char)
%   wavelengths: array of wavelengths (nm)
%
% Output:
%   data   : [numVoxels x numWavelengths]
%   volSize: [Nx Ny Nz]

numWvl = numel(wavelengths);
vol = [];

for i = 1:numWvl
    matName = sprintf("%s%d.mat", baseName, wavelengths(i));
    matPath = fullfile(inputDir, matName);

    if ~isfile(matPath)
        error("Missing file: %s", matPath);
    end

    S = load(matPath);
    if ~isfield(S, "image3D_tot")
        error("File %s does not contain variable 'image3D_tot'.", matPath);
    end

    img = S.image3D_tot;
    if i == 1
        vol = img;
        volSize = size(img);
        if numel(volSize) ~= 3
            error("image3D_tot must be a 3D array. Got size: %s", mat2str(volSize));
        end
        numVox = prod(volSize);
        data = zeros(numVox, numWvl);
    else
        if ~isequal(size(img), volSize)
            error("Volume size mismatch at %dnm. Expected %s, got %s.", ...
                wavelengths(i), mat2str(volSize), mat2str(size(img)));
        end
    end

    data(:,i) = reshape(img, [], 1);
end

end


function A = calculate_abundance_nnls(D, E, opts)
% CALCULATE_ABUNDANCE_NNLS
% Computes nonnegative abundances A for each voxel:
%   minimize || E * a - d ||_2  subject to a >= 0
%
% Inputs:
%   D: [numVoxels x numWavelengths] (each row is a spectrum)
%   E: [numWavelengths x numEndmembers] (endmember spectra)
%   opts.normalizeEndmembers: if true, normalizes each endmember column
%
% Output:
%   A: [numVoxels x numEndmembers]

arguments
    D double
    E double
    opts.normalizeEndmembers (1,1) logical = true
end

E_use = E;
if opts.normalizeEndmembers
    % Normalize columns to unit norm (matches your original intent)
    colNorm = vecnorm(E_use, 2, 1);
    colNorm(colNorm == 0) = 1;
    E_use = E_use ./ colNorm;
end

numVoxels = size(D,1);
numEnd = size(E_use,2);
A = zeros(numVoxels, numEnd);

% Loop per voxel (simple and transparent for public release)
for i = 1:numVoxels
    A(i,:) = lsqnonneg(E_use, D(i,:)')';
end

end


function imgOut = format_for_display(imgIn)
% FORMAT_FOR_DISPLAY
% Optional rotation/flip to match your original visualization style.
% Original code used: imagesc(fliplr(rot90(max(vol,[],3),-1)))
imgOut = fliplr(rot90(imgIn, -1));
end