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
