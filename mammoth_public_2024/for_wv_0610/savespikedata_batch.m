addpath(genpath('/path/to/NPMK'));

root_dir = '/your/root/folder';  % 改成你的根目录
ext = '.nev';                    % 你要查找的文件后缀

% 获取所有子文件夹路径
all_folders = genpath(root_dir);
folder_list = strsplit(all_folders, pathsep);  % 用 pathsep 分割成 cell 数组

for i = 1:length(folder_list)
    Brfilepath = folder_list{i};
    if isempty(Brfilepath)
        continue;
    end

    % 获取该文件夹下所有指定后缀文件
    files = dir(fullfile(Brfilepath, ['*' ext]));
    for j = 1:length(files)
        Brfilename = files(j).name;

        % 打印当前处理文件
        fprintf('Processing: %s\n', fullfile(Brfilepath, Brfilename));

        % === 调用你已有的代码块 ===
        openNEV(fullfile(Brfilepath, Brfilename));
        timestamp = NEV.Data.Spikes.TimeStamp;
        electrode = NEV.Data.Spikes.Electrode;
        spikedata = cell(max(electrode)+1,1);
        marker = NEV.Data.SerialDigitalIO.UnparsedData;
        markertime = NEV.Data.SerialDigitalIO.TimeStampSec * 1000;  % ms
        spikedata{max(electrode)+1} = cat(2, marker, markertime);
        
        % 你可以在这里保存结果，例如：
        % save(fullfile(Brfilepath, [Brfilename, '_processed.mat']), 'spikedata');

    end
end


Brfilepath='/home/cuihe_lab/chenyun/share_rw/CuiLab-Database/double_reach/Caesar/data_recording/20201010_DoubleReach_001_TestUDP/';
cd(Brfilepath)
Brfilename='caeser20201010001.nev';
matfilename=fullfile(Brfilepath,Brfilename);
matfilename=strrep(matfilename,'.nev','_BRspikeinfo.mat');
%%
openNEV(fullfile(Brfilepath,Brfilename));
% matfilename=strrep(fullfile(Brfilepath,Brfilename),'nev','mat');
% load(filename);
timestamp=NEV.Data.Spikes.TimeStamp;
electrode=NEV.Data.Spikes.Electrode;
spikedata=cell(max(electrode)+1,1);
marker=NEV.Data.SerialDigitalIO.UnparsedData;
markertime=NEV.Data.SerialDigitalIO.TimeStampSec*1000; %ms
spikedata{max(electrode)+1}=cat(2,marker,markertime');
%%
%spiketime
for i=1:max(electrode)
    a=timestamp(electrode==i);
    a=double(a)/30;   % ms 
    spikedata{i}=a;
end
