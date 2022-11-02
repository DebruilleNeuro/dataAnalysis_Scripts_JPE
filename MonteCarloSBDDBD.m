%Script to create files needed for the Monte-Carlo procedure 
%Exctracts the trials that were accepted after artifact rejection and
%homogenize the conditions as to have the same number of trials in all
%conditions to do the Monte-Carlo. 
%
%This script is for a 5 conditions JPE experiment, but can be modified to
%fit any number of conditions. 
%
%
% Add eeglab path to matlab path
%addpath(genpath('/Users\jeula\Documents\Debruille2022\CurrentFolderVersion\Enora Matlab Current subjects\eeglab14_1_1b'))
addpath(genpath('/Users\jeula\Documents\current subjects\eeglab2021.1'))

% define the filename
name='ICs_ICA_0.1HZ_C2J12H1.set';%load manually if does not work and run from line 19 redraw
currentPath = pwd;

% eeglab command to open the GUI
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%load the data set in variable called 'name'
EEG = pop_loadset('filename',name,'filepath',currentPath); 

%eegcheckset, check the consistency of the fields of an EEG dataset
EEG = eeg_checkset( EEG ); 
eeglab redraw;

%% ARTIFACT REJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF YOUR DATA DIDNT GO THROUGH ARTIFACT REJECTION UNCOMMENT THE FOLLOWING
%SECTION
%
%electrodes=[1 2 7:28]; %
%       frontals = [3:6];%
% 
% %    Mark epochs containing activity above an upper threshold and below a
% %    lower threshold [ -75 75] uV for channels defined aboved in
% %    'electrodes', in time window  [ -204 1200]. Will flag the trials  in the EEG structure. 
%     EEG  = pop_artextval( EEGSET , 'Channel', electrodes, 'Flag',  1, 'Threshold', [ -75 75], 'Twindow',...
%         [ -204 1200] );
%     
%     %Same but for for fp1 fp2 f7 f8, threshold is [ -100 100] instead of [ -75 75]
%     EEG  = pop_artextval( EEG , 'Channel', frontals, 'Flag',  1, 'Threshold', [ -100 100], 'Twindow',...
%         [ -204 1200] );
%         
%     %removes data where there is a flatline at any electrode for more
%     %than 100ms
%     EEG  = pop_artflatline( EEG , 'Channel', electrodes, 'Duration',  100, 'Flag',  1, 'Threshold', [ -1e-07 1e-07], 'Twindow',...
%         [ -204 1200] );
%     %same for frontal electrodes
%     EEG  = pop_artflatline( EEG , 'Channel', frontals, 'Duration',  100, 'Flag',  1, 'Threshold', [ -1e-07 1e-07], 'Twindow',...
%         [ -204 1200] );
% 
%     EEG = pop_saveset(EEG, [name], [currentPath]); 
%     save EEG
% EEG structure named eeg and saved 
%%
 eeg=EEG;
 save eeg eeg

close (eeglab)%close 
%clear ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY LASTCOM STUDY eeg name subjectPath

%remove eeglab path from matlab path
%rmpath(genpath('/Users\jeula\Documents\current subjects\eeglab2021.1'))
rmpath(genpath('/Users\jeula\Documents\current subjects\eeglab2021.1'))

%clear all
%% 
addpath(genpath('/Users\jeula\Documents\current subjects\fieldtrip\fieldtrip-20190410'))    

%load('eeg.mat');%load structure

filename = eeg.filename(1:6)

% create a variable "code" including all S and D events and rejected trials
for i=1:numel(eeg.event)
code(i,1)=eeg.EVENTLIST.eventinfo(i).bini;
code(i,2)=eeg.EVENTLIST.eventinfo(i).flag; 
end

%create list of events per condition (number 1:5 in bini) with only
%accepted trials (flag=0 not 1)
for i=1:numel(eeg.event)
    if code(i,1)==1 && code(i,2)==0
        ar_code(i,1)=1; %bin1 condition SBD
    elseif code(i,1)==2 && code(i,2)==0
        ar_code(i,1)=2; %bin2 condition DBD
%     elseif code(i,1)==3 && code(i,2)==0
%         ar_code(i,1)=3; %bin3 condition INI
%     elseif code(i,1)==4 && code(i,2)==0
%         ar_code(i,1)=4; %bin4 condition NII
%     elseif code(i,1)==5 && code(i,2)==0
%         ar_code(i,1)=5; %bin5 condition NN
%     else
        ar_code(i,1)=0;
    end
end

% finding the row numbers corresponding to artifact-free trials for each
% condition

% SBD=find(ar_code==1);
% DBD=find(ar_code==2);
SBD=find(code==1);
DBD=find(code==2);
% INI=find(ar_code==3);
% NII=find(ar_code==4);
% NN=find(ar_code==5); 

%Rejects trials so that each condition has the same amount of trials 
    
%BETWEEN ALL CONDITIONS
trials_array = [numel(DBD),numel(SBD)]
min_val = min(min(trials_array))

if min_val<35;
    msg = 'condition(s) may be missing or not enough trials may be present.';
    error(msg)
end

conditions = {SBD DBD}

for i = 1:numel(conditions)
    %bin = conditions(i)
    if numel(conditions{1,i})~=min_val
    gap= numel(conditions{1,i}) - min_val
    for n=1:abs(gap)
        Index = randi(length(conditions{1,i}), 1)
        conditions{1,i}(Index,:) = []
    end
    end 
end

% extracting the trials
for i = 1:numel(conditions);
    for n=1:length(conditions{1,i});
     currentdata = eeg.data(:,:,conditions{1,i});
     currentBin(:,:,n)= currentdata(:,:,n);
    end
    dataBins{1,i}=currentBin;
end

eegperbin=struct('SBD', eeg, 'DBD', eeg)
fieldNames = fieldnames(eegperbin); 
cfg=[];

for loopBins= 1:numel(fieldNames) 
    eegperbin.(fieldNames{loopBins}).data = dataBins{1:loopBins};
    eegperbin.(fieldNames{loopBins}).trials=size(eegperbin.(fieldNames{loopBins}).data,3);
    eegperbin.(fieldNames{loopBins}).data = eeglab2fieldtrip( eegperbin.(fieldNames{loopBins}), 'preprocessing' );
    [monteCarloFiles.(fieldNames{loopBins})] = eegperbin.(fieldNames{loopBins}).data
    cfg.keeptrials='yes';
    %computes the timelocked average ERP and computes the covariance matrix
    [monteCarloFiles.(fieldNames{loopBins})] = ft_timelockanalysis(cfg, monteCarloFiles.(fieldNames{loopBins}))
    %adds filename to monteCarloFiles struct
    [monteCarloFiles.(fieldNames{loopBins}).filename] = eegperbin.(fieldNames{loopBins}).filename
    [monteCarloFiles.(fieldNames{loopBins}).filebin] = (fieldNames{loopBins});
    nameAVG = [eegperbin.(fieldNames{loopBins}).filename(end-10:end-4), '_',(fieldNames{loopBins}), '_AVG'];
    currentfile = monteCarloFiles.(fieldNames{loopBins});
   save(nameAVG, 'currentfile');
end 


%clearvars -except monteCarloFiles
fprintf(':) done :)\n');