%Script to preprocess Tarlan's data.
%Version 4
%Input: Iwave files (ic.asc.eeglab.data.txt and ic.asc.eeglab.events.txt).
%!!Input file name: Pair number needs to be characters 1:2 (ie paire 1= 01...txt)
%Output: Set &a erp files, preprocessed (Pruned with ICA, 0.1hz/50hz filters applied, artifact rejection done...)

clear all
%close all
clc

txtfilename = 'BS110.ic.asc.eeglab.events.txt';
datafilename = 'BS110.ic.asc.eeglab.data.txt';

path2eeglab ='/Users\jeula\Documents\current subjects\eeglab2021.1' %Here, put your OWN path to EEGLAB
path2fieldtrip ='/Users\jeula\Documents\current subjects\fieldtrip\fieldtrip-20190410' %Here, put your OWN path to FIELDTRIP

addpath(genpath(path2eeglab))
%open eeglab
eeglab

%initializing variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ERPSET =[];
ERP =[];
nameset=[];
nameerp = [];
name_temp=[];
events=[];
EEGSET=[];
h=0;
EEG=[];
exportname=[];
electrodes=[];
placingelectrode=[];
nameset = datafilename;
EEG=EEGSET

EEGSET = pop_importdata('dataformat','ascii','nbchan',28,'data',datafilename,'srate',248,'pnts',0,'xmin',0,'chanlocs',[pwd '/Standard-10-20-Cap8128.ced']);
EEGSET = pop_importevent( EEGSET, 'event', txtfilename,'fields',{'latency' 'type'},'skipline',1,'timeunit',1, 'append','no','optimalign','on');
events = [pwd '/EVENTLIST-SHEYLA.txt'];
EEGSET  = pop_editeventlist( EEGSET , 'BoundaryNumeric', { -99}, 'BoundaryString', { 'boundary' }, 'List', events, 'SendEL2', 'EEG', 'UpdateEEG', 'codelabel', 'Warning', 'on' ); % GUI: 28-Jan-2020 12:15:23
EEGSET = pop_epochbin( EEGSET , [-204.0  1000.0],  [ -204 -4]);

% savefile
EEGSET = pop_saveset( EEGSET, nameset,[pwd]);

name_temp = datafilename;
% Selecting the subject

% human 1 0,1hz 50hz filters
EEG=EEGSET

nameset = ['0.1HZ_' name_temp(3:5) '_S1' '.'];
nameerp = [name_temp(3:5) '_S1' '.erp'];
electrodes =[5:28];
frontals = [1:4];
placingelectrode = {'nch1 = ch1 label Fp2',  'nch2 = ch2 label Fp1',  'nch3 = ch3 label F8',  'nch4 = ch4 label F7',...
    'nch5 = ch5 label Fz',  'nch6 = ch6 label Cz',  'nch7 = ch7 label Pz',  'nch8 = ch8 label P4',  'nch9 = ch9 label P3',  'nch10 = ch10 label T6',...
    'nch11 = ch11 label T5',  'nch12 = ch12 label T4',  'nch13 = ch13 label T3',  'nch14 = ch14 label F4',  'nch15 = ch15 label F3',...
    'nch16 = ch16 label Ft8',  'nch17 = ch17 label Ft7',  'nch18 = ch18 label Fc4',  'nch19 = ch19 label Fc3',  'nch20 = ch20 label Fcz',  'nch21 = ch21 label C4',...
    'nch22 = ch22 label C3',  'nch23 = ch23 label Tp8',  'nch24 = ch24 label Tp7',  'nch25 = ch25 label Cp4',  'nch26 = ch26 label Cp3',...
    'nch27 = ch27 label O2',  'nch28 = ch28 label O1'};
EEG = pop_eegchanoperator(EEG, placingelectrode);%placing electrodes
EEG = pop_editset(EEG, 'run', [], 'chanlocs', [pwd '/Standard-10-20-Cap8128.ced']);
% filter data
EEGSET = pop_eegfiltnew( EEGSET, [], 50, [], false, [], 0); %50hz
EEGSET = pop_eegfiltnew( EEGSET, 0.1, [], [], false, [], 0); %0.1hz
EEG = pop_saveset(EEG,[nameset],[pwd]) %save
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

EEG=EEGSET
nameset = ['1HZ_' name_temp(3:5) '_S1' '.'];
electrodes=[5:28];
frontals = [1:4];
placingelectrode = {'nch1 = ch1 label Fp2',  'nch2 = ch2 label Fp1',  'nch3 = ch3 label F8',  'nch4 = ch4 label F7',...
    'nch5 = ch5 label Fz',  'nch6 = ch6 label Cz',  'nch7 = ch7 label Pz',  'nch8 = ch8 label P4',  'nch9 = ch9 label P3',  'nch10 = ch10 label T6',...
    'nch11 = ch11 label T5',  'nch12 = ch12 label T4',  'nch13 = ch13 label T3',  'nch14 = ch14 label F4',  'nch15 = ch15 label F3',...
    'nch16 = ch16 label Ft8',  'nch17 = ch17 label Ft7',  'nch18 = ch18 label Fc4',  'nch19 = ch19 label Fc3',  'nch20 = ch20 label Fcz',  'nch21 = ch21 label C4',...
    'nch22 = ch22 label C3',  'nch23 = ch23 label Tp8',  'nch24 = ch24 label Tp7',  'nch25 = ch25 label Cp4',  'nch26 = ch26 label Cp3',...
    'nch27 = ch27 label O2',  'nch28 = ch28 label O1'};
EEG = pop_eegchanoperator(EEG, placingelectrode); %placing electrodes
% filter data
EEGSET = pop_eegfiltnew( EEGSET, 1, [], [], false, [], 0);
EEG = pop_editset(EEG, 'run', [], 'chanlocs', [pwd '/Standard-10-20-Cap8128.ced']);%load channel location info
EEG = pop_saveset(EEG,[nameset],[pwd]) %save set
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%% Human 1 (ICA Artifact rejection, Creation of ERP) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    j=1; % participant 1
    
    EEG = pop_loadset('filename',['1HZ_' name_temp(3:5) '_S' int2str(j) '.set'],'filepath',[pwd]); %load 1hz dataset for ICA
    nameerp = [name_temp(3:5) '_S' int2str(j) '.erp'];
    nameset = [name_temp(3:5) '_S' int2str(j) '.set'];
    EEG = pop_editset(EEG, 'run', [], 'chanlocs', [pwd '/Standard-10-20-Cap8128.ced']);%load channel location info
    
    %automatic channel rejection
    %pop_rejchan(EEG)
    EEG = pop_rejchan(EEG, 'elec',[1:28],'measure','prob','norm','on','threshold',5); %automatic rejection parameters
    %  fprintf('If *bad* channels exist, remove them from brackets in pop_runica')
    %% run ICA
    %!!! Chanind remove *bad* electrode in Pop up (or if ICA automated, remove from brackets to run e.g.[1:23 25:28] %24)
    EEG = pop_runica( EEG )
    EEG = pop_saveset( EEG, 'filename',['ICA_1HZ_' name_temp(3:5) '_S' int2str(j) '.set'],'filepath',[pwd ]); %save set
    
    %ICA activation matrix
    TMP.icawinv = EEG.icawinv;
    TMP.icasphere = EEG.icasphere;
    TMP.icaweights = EEG.icaweights;
    TMP.icachansind = EEG.icachansind;
    
    % apply matrix to 0.1hz dataset
    clear EEG;
    
    EEG = pop_loadset('filename', ['0.1HZ_' name_temp(3:5) '_S' int2str(j) '.set'], 'filepath', [pwd]); %load 0.1hz .set
    EEG.icawinv = TMP.icawinv;
    EEG.icasphere = TMP.icasphere;
    EEG.icaweights = TMP.icaweights;
    EEG.icachansind = TMP.icachansind;
    clear TMP;
    EEG = pop_saveset(EEG, 'filename',['ICA_0.1HZ_' name_temp(3:5) '_S' int2str(j) '.set'], 'filepath', [pwd]); %save 0.1hz+ICA matrix .set
    
    %% !!! when 'reject component' window pops up, before rejecting need to label components manually (precaution)
    EEG = pop_loadset('filename', ['ICA_0.1HZ_' name_temp(3:5) '_S' int2str(j) '.set'], 'filepath', [pwd]);
    %IC component rejection
    EEG=iclabel(EEG);
    noisethreshold = [0 0;0.9 1; 0.9 1; 0 0; 0 0; 0 0; 0 0]; %IC label parameters: 90% Muscle and Eye probability;
    EEG = pop_icflag(EEG, noisethreshold);
    pop_viewprops (EEG, 0)
    pause
    % remove bad component(s)
    EEG = pop_subcomp( EEG ); %
    % save
    EEG = pop_saveset(EEG, 'filename',['ICs_ICA_0.1HZ_' name_temp(3:5) '_S' int2str(j) '.set'], 'filepath', [pwd]); %set 0.1hz filter + ICA + bad ICs removed
    
    % check bad channels again
    EEG = pop_rejchan(EEG)
   
%% artifact detection loop

EEG=EEGSET
    EEG = pop_loadset('filename',['ICs_ICA_0.1HZ_' name_temp(3:5) '_S' int2str(j) '.set'],'filepath',[pwd]); %load 1hz dataset for ICA
    
    nameerp = [name_temp(3:5) '_S' int2str(j) '.erp'];
    namesum =[name_temp(3:5) '_S' int2str(j) ];
    nameset = [name_temp(3:5) '_S' int2str(j) '.set'];
    %%%!! exclude *bad* electrodes, comment which electrode(s) and restore
    %%%after participant is done
    frontals = [1:4]; %1 2
    electrodes=[5:28];  %28 17take note of which electrode is removed
    
    %peak to peak (frontal elec and other elec)
    EEG  = pop_artextval( EEG , 'Channel', electrodes, 'Flag',  1, 'Threshold', [ -75 75], 'Twindow',[ -204 1000] );
    EEG  = pop_artextval( EEG , 'Channel', frontals, 'Flag',  1, 'Threshold', [ -100 100], 'Twindow',[ -204 1000] );
    
    %flat line (frontal elec and other elec)
    EEG  = pop_artflatline( EEG , 'Channel', electrodes, 'Duration',  100, 'Flag',  1, 'Threshold', [ -1e-07 1e-07], 'Twindow', [ -204 1000] );
    EEG  = pop_artflatline( EEG , 'Channel', frontals, 'Duration',  100, 'Flag',  1, 'Threshold', [ -1e-07 1e-07], 'Twindow', [ -204 1000] );
    
    %close;
    EEG = pop_saveset( EEG, [nameset] ,[pwd]);
    
    %% compute erp
    ERP = pop_averager( EEG , 'Criterion', 'good', 'DSindex',1, 'ExcludeBoundary', 'on', 'SEM', 'on' );
    
    % load channel location information
    ERP = pop_erpchanedit( ERP, [pwd '/Standard-10-20-Cap8128.ced']);
    
    % Save the erp
    ERP = pop_savemyerp(ERP, 'erpname', nameerp, 'filename', nameerp, 'filepath', [pwd], 'Warning', 'on');
    ERP = pop_summary_AR_erp_detection(ERP, [pwd '\' namesum '.txt'])
    
    ERP = pop_summary_rejectfields(EEG);
    
    clear
    
    fprintf(':) Participant done :)');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%