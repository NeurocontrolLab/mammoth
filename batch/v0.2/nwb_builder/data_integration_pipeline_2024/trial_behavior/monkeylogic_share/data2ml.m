% data_path = '/AMAX/cuihe_lab/share_rw/Neucyber-NC-2023-A-01/Nezha/Data_recording/20240315_centerOut_002/bhv/240315_Nezhat_nezha.bhv2';
% data_path=argv(1);
warning('off','all');

disp(data_path);

currentFolder = '/AMAX/cuihe_lab/cuilab_share/monkeylogic';
 
addpath(currentFolder);

subfolders = dir(currentFolder);
subfolders = {subfolders(~[subfolders.isdir]).name};
for i = 1:length(subfolders)
    addpath(fullfile(currentFolder, subfolders{i}));
end

data = mlread(data_path);

[path, name, ext] = fileparts(data_path);

disp(path)

% disp(data)

% disp([path, name,'.mat'])

save([path, '/', name,'.mat'], 'data');
