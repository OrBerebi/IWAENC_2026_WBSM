% Batch Audio Resampler
% Reads .wav files from a base directory, resamples them to 48kHz,
% and saves them in a parallel directory while preserving folder structure.

% 1. Define your base directory here
baseDir = '/Users/orberebi/Documents/GitHub/IWAENC_2026_WBSM/results/04_03_24_test_signals/'; % <-- REPLACE THIS PATH

% Clean up trailing slashes for safer path manipulation
if endsWith(baseDir, filesep)
    baseDir = baseDir(1:end-1);
end

% Define the target sample rate
targetFs = 48000;

% 2. Automatically generate the parallel output folder path
[parentDir, baseName, ~] = fileparts(baseDir);
outDir = fullfile(parentDir, [baseName, sprintf('_resampled_%dkHz', targetFs/1000)]);

% 3. Find all .wav files recursively
fprintf('Searching for .wav files in: %s\n', baseDir);
wavFiles = dir(fullfile(baseDir, '**', '*.wav'));
numFiles = length(wavFiles);

if numFiles == 0
    disp('No .wav files found in the specified directory.');
    return;
end

fprintf('Found %d files. Starting resampling process...\n', numFiles);

% 4. Loop through and process each file
for i = 1:numFiles
    % Get full path of the current input file
    inFile = fullfile(wavFiles(i).folder, wavFiles(i).name);
    
    % Extract the relative path to maintain subdirectory structure
    relativePath = strrep(wavFiles(i).folder, baseDir, '');
    if startsWith(relativePath, filesep)
        relativePath = relativePath(2:end);
    end
    
    % Define the corresponding output subdirectory and create it if needed
    outSubDir = fullfile(outDir, relativePath);
    if ~exist(outSubDir, 'dir')
        mkdir(outSubDir);
    end
    
    % Define full output file path
    outFile = fullfile(outSubDir, wavFiles(i).name);
    
    try
        % Read the audio file
        [y, fs] = audioread(inFile);
        
        % Check if the file is already at the target sample rate
        if fs == targetFs
            fprintf('[%d/%d] Skipping resampling (already %dHz): %s\n', i, numFiles, targetFs, wavFiles(i).name);
            copyfile(inFile, outFile); % Copy to keep the new dataset complete
            continue;
        end
        
        % Resample the audio array
        % MATLAB's resample applies an automatic polyphase anti-aliasing filter
        y_resampled = resample(y, targetFs, fs);
        
        % Write the new file
        audiowrite(outFile, y_resampled, targetFs);
        fprintf('[%d/%d] Resampled and saved: %s\n', i, numFiles, wavFiles(i).name);
        
    catch ME
        % Error handling so one corrupted wav file doesn't kill the whole batch
        warning('Failed to process %s.\nReason: %s', wavFiles(i).name, ME.message);
    end
end

disp('Batch resampling complete!');
fprintf('All files saved to: %s\n', outDir);