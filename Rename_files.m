%% Renaming file names from 488_YYYY-MM-DD-HH_MM_SS.tiff --> 1.tiff, 2.tiff and so on

% Get the list of .tiff files in the current folder
files = dir('*.tiff');

% Preallocate a struct array to store filenames and their corresponding time
fileData = struct('name', {}, 'time', {});

% Loop through each file and extract the time portion
for i = 1:length(files)
    % Get the current file name
    fileName = files(i).name;
    
    % Extract the time part (HH_MM_SS) from the filename
    % Assuming the format is: 488_YYYY-MM-DD-HH_MM_SS.tiff
    tokens = regexp(fileName, '(\d{2}_\d{2}_\d{2})', 'match');
    timeStr = tokens{1}; % Extract the matched time string
    
    % Convert the time string to a serial date number for easy sorting
    fileTime = datenum(timeStr, 'HH_MM_SS');
    
    % Store the filename and its time
    fileData(i).name = fileName;
    fileData(i).time = fileTime;
end

% Sort the files based on the extracted time
[~, sortedIdx] = sort([fileData.time]);

% Rename the files in ascending order as 1.tiff, 2.tiff, etc.
for i = 1:length(sortedIdx)
    oldName = fileData(sortedIdx(i)).name;
    newName = sprintf('%d.tiff', i);
    movefile(oldName, newName);  % Rename the file
end

disp('Files renamed successfully.');