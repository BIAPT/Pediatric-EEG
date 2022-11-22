%% Spectrogram and Topographic Map for PD EEG
% Run them first prior to running the rest of the analysis
set(groot, 'DefaultFigureVisible', 'off');
%% Load clean EEG data set
% select the folder with the data set
waitfor(msgbox('Select the folder which contains the data in BIDS format.'));
datafolder=uigetdir(path);

waitfor(msgbox('Select the Saving directory'));
resultsfolder = uigetdir(path);

cd(datafolder)
cd("eeg")
files = dir('*.set*');

for f = 1:numel(files)
    %% Load data and get information about the state
    recording = load_set(files(f).name,pwd);
    sampling_rate = recording.sampling_rate;
    info = split(files(f).name,'_');
    ID = info{1}(5:end);
    task = info{2}(6:end);
    hemispheres = "Whole";
    disp("load complete: " + ID + '_' + task)
    
    %% Spectrogram 
    spectopo_prp = spectopo_prp_struct_PD;
    disp('spectopo_prp load complete')
    outdir_spectrogram = fullfile(resultsfolder, ID,"Whole",'Spectrogram');
    spectrogram = spectrogram_function(recording, spectopo_prp, ID, task, outdir_spectrogram); %figure is not shown because set to 'off'
    
    %% Topographic Maps of Alpha and Theta Power
    outdir_topographicmap = fullfile(resultsfolder, ID, "Whole",'Topographic Maps');
    topographic_map = topographic_map_function_PD(recording, spectopo_prp, ID, task, outdir_topographicmap); %figure is not shown because set to 'off'
end