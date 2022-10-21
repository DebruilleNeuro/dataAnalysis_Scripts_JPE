% Add eeglab path
addpath(genpath('/Users\jeula\Documents\Debruille2022\CurrentFolderVersion\Enora Matlab Current subjects\eeglab14_1_1b'))
%addpath(genpath('C:\Users\jeula\Documents\current
%subjects\eeglab2021.1\plugins\ERPLAB9.00'))

% define the files' labels with respect to cycles (subj nb)
name='AR1MH1.set';%load manually if does not work and run from line 19 redraw
currentPath = pwd;

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;                        % eeglab command to open the GUI
EEG = pop_loadset('filename',name,'filepath',currentPath); % load the data set, the dataset's name is stored in the variable we called it 'name'
EEG = eeg_checkset( EEG );
eeglab redraw;
% This is how to convert eeglab data to fieldtrip structure
% first need to open a set file in eeglab

% at this point, you can apply your own scripts to reject artefacts... you
% can copy paste in the space below the lines of codes


% once this is done,  you save the structure EEG that you can see in
% Matlab workspace, and give it a new name (e.g. eeg)
eeg=EEG;
save eeg eeg
close (eeglab)
clear ALLEEG EEG CURRENTSET ALLCOM CURRENTSTUDY LASTCOM STUDY eeg name subjectPath
%cd ..

% remove eeglab path
%rmpath(genpath('/Users\jeula\Documents\current subjects\eeglab2021.1'))
rmpath(genpath('/Users\jeula\Documents\Debruille2022\CurrentFolderVersion\Enora Matlab Current subjects\eeglab14_1_1b'))
clear all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% converting eeg matfiles into set files in eeglab and saving  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After closing eeglab,  you need to add fieldtrip in Matlab path,
addpath(genpath('/Users\jeula\Documents\current subjects\fieldtrip\fieldtrip-20190410'))
ft_defaults

load('eeg.mat');

%%

filename = eeg.filename(1:6)
%j=1;

% create a variable "code" including all S and D events and rejected trials
for i=1:numel(eeg.event)
    code(i,1)=eeg.EVENTLIST.eventinfo(i).bini;
    code(i,2)=eeg.EVENTLIST.eventinfo(i).flag;
end

% recode the events so as to exclude rejected trials ==> 0
% recode the events so as to exclude rejected trials ==> 0
for i=1:numel(eeg.event)
    if code(i,1)==1 && code(i,2)==0
        ar_code(i,1)=1;
    elseif code(i,1)==2 && code(i,2)==0
        ar_code(i,1)=2;
    else
        ar_code(i,1)=0;
    end
end

% finding the row numbers corresponding to artefact-free S and D trials
SBD=find(ar_code==1);
DBD=find(ar_code==2);

%for Discordant/Concordant conditions only
gap = numel(DBD) - numel(SBD)
if gap < 0
    for n=1:abs(gap)
        randomIndex = randi(length(SBD), 1)
        %selected_DBD_value = DBD(randomIndex)
        SBD(randomIndex,:) = []
    end
    
else if gap > 0
        for n=1:abs(gap)
            randomIndex = randi(length(DBD), 1)
            %selected_DBD_value = DBD(randomIndex)
            DBD(randomIndex,:) = []
        end
    end
end

% extracting the SBD data trials
for i=1:length(SBD)
    data_SBD(:,:,i)=eeg.data(:,:,SBD(i));
end
% extracting the DBD data trials
for i=1:length(DBD)
    data_DBD(:,:,i)=eeg.data(:,:,DBD(i));
end

eeg_SBD=eeg;
eeg_DBD=eeg;

eeg_SBD.data=data_SBD;
eeg_DBD.data=data_DBD;

eeg_SBD.trials=size(data_SBD,3);
eeg_DBD.trials=size(data_DBD,3);


%clear eeg* code data_SBD data_DBD data_INI data_NII data_NN
% the following function will convert the eeg structure into a fieldtrip
% compatible format, and give it a new name (e.g. data).

data_SBD = eeglab2fieldtrip( eeg_SBD, 'preprocessing' );
data_DBD = eeglab2fieldtrip( eeg_DBD, 'preprocessing' );

% Since the file I have access to contains the data of two participants, I
% will create in the following 2 files, containing the data of partipant 1
%clear eeg* code SBD DBD INI NII NN
%clear eeg* code SBD DBD INI NII NN data_SBD data_DBD data_INI data_NII data_NN
% remove the field elec
data_S=rmfield(data_SBD,'elec');
data_D=rmfield(data_DBD,'elec');

%clear eeg* code SBD DBD data_SBD data_DBD data_INI data_NII data_NN


%% Here once we created the two fieldtrip files, we are going to seperate the s1 and s2
load('label');

%str2double(eeg.filename(end-4))
% organize the data of S1
% s1_s
if str2double(eeg.filename(end-4))==1;
    % participant 1 data,
    H_SBD=data_S;
    H_SBD.label=label;
    H_DBD=data_D;
    H_DBD.label=label;
    
    for kk=1:size(H_SBD.trial,2)
        aa=H_SBD.trial{kk}; % extract the first half of the data matrix from trial kk=1 till the end
        a=aa(1:35,:);
        b=a(3:8,:);
        c=a(10:17,:);
        d=a(19:26,:);
        e=a(28:35,:);
        f=[b ; c ; d ; e];
        H_SBD.trial{kk}=f;
    end
    % organize the data of S1
    % s1_d
    for kk=1:size(H_DBD.trial,2)
        aa=H_DBD.trial{kk}; % extract the first half of the data matrix from trial kk=1 till the end
        a=aa(1:35,:);
        b=a(3:8,:);
        c=a(10:17,:);
        d=a(19:26,:);
        e=a(28:35,:);
        f=[b ; c ; d ; e];
        H_DBD.trial{kk}=f;
    end
    
else if str2double(eeg.filename(end-4))==2;
        % participant 1 data,
        H_SBD=data_S;
        H_SBD.label=label;
        H_DBD=data_D;
        H_DBD.label=label;
        
        for kk=1:size(H_SBD.trial,2)
            aa=H_SBD.trial{kk}; % extract the first half of the data matrix from trial kk=1 till the end
            a=aa(37:71,:); % extract the second half of the data matrix from trial kk=1 till the end
            b=a(3:8,:);
            c=a(10:17,:);
            d=a(19:26,:);
            e=a(28:35,:);
            f=[b ; c ; d ; e];
            H_SBD.trial{kk}=f;
        end
        % organize the data of S1
        % s1_d
        for kk=1:size(H_DBD.trial,2)
            aa=H_DBD.trial{kk}; % extract the first half of the data matrix from trial kk=1 till the end
            a=aa(37:71,:); % extract the second half of the data matrix from trial kk=1 till the end
            b=a(3:8,:);
            c=a(10:17,:);
            d=a(19:26,:);
            e=a(28:35,:);
            f=[b ; c ; d ; e];
            H_DBD.trial{kk}=f;
        end
    end
end
        
        % averaging trials and keeping all trials
        cfg=[];
        cfg.keeptrials='yes';
        [H_SBD_avg] = ft_timelockanalysis(cfg, H_SBD)
        
        cfg=[];
        cfg.keeptrials='yes';
        [H_DBD_avg] = ft_timelockanalysis(cfg, H_DBD)
        
        
        name1=[filename, '_SBD_AVG'];
        name2=[filename, '_DBD_AVG'];
        
        save(name1,'H_SBD_avg');
        save(name2, 'H_DBD_avg');
        
        clear
        fprintf(':) done :) \n');