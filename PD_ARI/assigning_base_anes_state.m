%% Prior to running step2_code of ARI pipeline [2]

% Miriam Han July 1st, 2022
% Use this code to automatically select all the necessary functional

% connectivity data (.mat files) to run step2 of ARI code

%% Assign ID, baseline and anesthesia state

% select the folder with the 
waitfor(msgbox('Select the folder which contains the calculated values from STEP 1 of ARI code.'));
datafolder=uigetdir(path);

waitfor(msgbox('Select the Saving directory'));
resultsfolder = uigetdir(path);

cd(datafolder)
cd("fc_data")
%cd("fc_data_DR") or cd("fc_data_ARI")
files = dir('*.mat*');
all_state_num = numel(files)/18;  % number of all the states

ID = "";
all_states = string(zeros(1,all_state_num));
anes_states = string(zeros(1,all_state_num-1));

for f = 1: all_state_num
    info = split(files(f).name,'_');
    ID = info{4}(1:5);
    task = info{5}(1:6);
    all_states(1,f) = task;
    disp (string(f)+'. '+task);
end

baseline_num = input ('\nSelect the number to indicate the baseline: \n\n');
baseline_state = all_states(baseline_num);

% create array of all pairs of 1 baseline to 1 anesthesia state
ID_pair = string(zeros(all_state_num-1,3)); % ex: 003PD_1
ID_pair(:,2) = baseline_state;

anes_states = setdiff(all_states, baseline_state);

for f = 1:length(anes_states)
        ID_pair(f,1) = strcat(ID,'_',string(f));
        ID_pair(f,3) = anes_states(1,f);
end

%% Save files in pairs of 1 baseline to 1 anesthesia state (CHANGE OUTDIR1!)

frequencies = ["alpha","theta","delta"];

%outdir1 = fullfile(resultsfolder,ID); % For the separation of : 2state_comparison by ID
outdir1 = fullfile(resultsfolder) % For the all 2states together in one folder / not separated by ID

for i = 1: length(ID_pair(:,1))
    mkdir(fullfile(outdir1,ID_pair(i)));
    textfile_directory = strcat(fullfile(outdir1,ID));
    states_ID_pair = [ID_pair(i,2),ID_pair(i,3)];
    
    for fre = 1:(length(frequencies))
        outdir2= fullfile(outdir1,ID_pair(i));
        
        for ii = 2:3
            
            if ii == 2
                state = "base";
            else
                state = "anes";
            end
            % wPLI (non-overlapping window of 10) (not time-resolved)
            %wpli_file_save = strcat('wpli_ARI_',frequencies (fre),'_',ID_pair(:,1),'_',state);
            % nov23: wpli_file_save = strcat('wpli_ARI_',frequencies (fre),'_',ID,'_',state);
            wpli_file_save = strcat('wpli_',frequencies (fre),'_',ID,'_',state);
            %wpli_file_directory = strcat(fullfile(outdir2,filesep, wpli_file_save, '.mat'));
            wpli_file_directory = strcat(outdir2,filesep, wpli_file_save, '.mat');

            disp ("Saving wpli of "+ID_pair(:,1)+' :'+ frequencies (fre) + '_'+state);
            % nov23:load_wpli = strcat('wpli_ARI_',frequencies (fre),'_',ID,'_',ID_pair(i,ii),'.mat');
            load_wpli = strcat('wpli_',frequencies (fre),'_',ID,'_',ID_pair(i,ii),'.mat');
            load (load_wpli);

            result_wpli_ARI = result_wpli.data.wpli;
            save (wpli_file_directory, 'result_wpli_ARI');
            
            % dPLI (non-overlapping window of 10) (not time-resolved)
            % nov23:dpli_file_save = strcat('dpli_ARI_',frequencies (fre),'_',ID,'_',state);
            dpli_file_save = strcat('dpli_',frequencies (fre),'_',ID,'_',state);
            %dpli_file_save = strcat('dpli_ARI_',frequencies (fre),'_',ID,'_',state)
            dpli_file_directory = strcat(outdir2,filesep, dpli_file_save, '.mat');

            disp ("Saving dpli of "+ID_pair(:,1)+' :'+ frequencies (fre) + '_'+state);
            % nov23:load_dpli = strcat('dpli_ARI_',frequencies (fre),'_',ID,'_',ID_pair(i,ii),'.mat');
            load_dpli = strcat('dpli_',frequencies (fre),'_',ID,'_',ID_pair(i,ii),'.mat');
            load (load_dpli);

            result_dpli_ARI = result_dpli.data.dpli;
            save (dpli_file_directory, 'result_dpli_ARI');
            
        end
    end
end
%% save the table of ID_pairs

RESULTS = table(ID_pair(:,1),ID_pair(:,2),ID_pair(:,3), 'VariableNames', {'pair_number', 'baseline','anesthesia'});
writetable(RESULTS, strcat(textfile_directory,' ARI_pair_info.txt'));

disp("Pipeline done :D!")


