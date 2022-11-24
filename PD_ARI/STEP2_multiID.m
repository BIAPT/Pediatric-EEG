%% Miriam Han June 30th, 2022 [3]
% this script is mainly based on the previous material from Yacin Mahdid
% and Charlotte Maschke
% Uses 2 state comparison to calculate the ARI

% need to loop over alpha and theta frequencies
% use two states comparison

%% File directories

%IN_DIR = '/Users/biapt/Desktop/PD_Analysis/ARI_results/step1_fc';
IN_DIR = '/Users/biapt/Desktop/PD_Analysis/Results';
MAP_FILE = '/Users/biapt/Desktop/PD_28_Whole_label.csv';
OUT_DIR = '/Users/biapt/Desktop/PD_Analysis/ARI_results/step2_ari/delta/new'; % saved in 1 state from each folder

%%  NEEDS TO BE DONE
% 1. create a new table with interest
% 2. table with different GCS, GOSE scales
% 3. Code to convert GOSE to recover vs non-recovered (setting limit)
%% create a a table from the excel file
%'2state_pair' excel file

pair_file = '/Users/biapt/Desktop/PD_Analysis/ARI_results/2statepair_eligible.csv';
state_pairs = table2array(readtable(pair_file)); % should be a file with (# of pairs) x 4 entries (ID, pair_num, base, anes)

P_ID = string(state_pairs(:,2)'); % for all states

P_LABEL = zeros(1,13); %for all states
%P_LABEL(1,1) = 1; % for all states

P_LABEL = [7,7,7,7,5,5,5,5,1,1,1,1,3];

threshold_range = 0.70:-0.01:0.01; % More connected to less connected

COLOR = 'jet';

outcome = "GOSE";

%interested_PairID = {'003PD_1','005PD_1','005PD_2','006PD_1','006PD_2','006PD_3','006PD_4','006PD_5','006PD_6','006PD_7',...
 %   '006PD_8','006PD_9','006PD_10','007PD_1','007PD_2','007PD_3','007PD_4',...
 %   '009PD_1','011PD_1','011PD_2'};
 
%%

%P_ID_interest = {'003PD_1','005PD_1','005PD_2','006PD_8','007PD_3','009PD_1','011PD_1'};

%state_pairs_interest = state_pairs(string(state_pairs(:,2)'))== P_ID_interest);

%state_pairs_interest = zeros(length(P_ID_interest),4);
%for i = 1:length(P_ID_interest)
   % state_pairs_interest(i,:) = ((state_pairs(:,2)==P_ID_interest(i)),:));
%end
%state_pairs_interest = state_pairs(string(state_pairs(:,2)')==string(P_ID_interest),:);


%% Loop through frequency 

%FREQUENCY = ["alpha","theta"]; 
FREQUENCY = ["delta"]; 


for fre = 1: length(FREQUENCY)
    % Here we iterate over each participant and each epochs to create the 2 subplots per figure
    dpli_dris = zeros(1, length(P_ID));
    hub_dris = zeros(1, length(P_ID));

    for p = 1:length(P_ID)
        participant = P_ID{p}(1:5); % ex: 003PD
        participant_pair = P_ID{p};
        disp(strcat("Participant: ", participant , "_dPLI"));
        disp(strcat("Participant pair: ", participant_pair));
    
        %% Calculate the dpli-dris
        %filename for baseline and anesthesia_file
        baseline_state = string(state_pairs(p,3)); %from the pair-info array
        anesthesia_state = string(state_pairs(p,4)); %from the pair-info array
        
        %nov23: baseline_file_dpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data_ARI',filesep, 'dpli_ARI_',FREQUENCY(fre),'_',participant,'_',baseline_state,'.mat');
        baseline_file_dpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data',filesep, 'dpli_',FREQUENCY(fre),'_',participant,'_',baseline_state,'.mat');
        %nov23: anesthesia_file_dpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data_ARI',filesep,'dpli_ARI_',FREQUENCY(fre),'_',participant,'_',anesthesia_state,'.mat');
        anesthesia_file_dpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data',filesep,'dpli_',FREQUENCY(fre),'_',participant,'_',anesthesia_state,'.mat');

        [baseline_r_dpli,baseline_r_location,baseline_r_regions] = process_dpli_original(baseline_file_dpli);
        [anesthesia_r_dpli,anesthesia_r_location,anesthesia_r_regions] = process_dpli_original(anesthesia_file_dpli);
        
        % Get the common location
        [common_labels, common_region] = get_subset_two(baseline_r_location, anesthesia_r_location, baseline_r_regions, anesthesia_r_regions);

        % Filter the matrices to have the same size
        baseline_f_dpli = filter_matrix(baseline_r_dpli, baseline_r_location, common_labels);
        anesthesia_f_dpli = filter_matrix(anesthesia_r_dpli, anesthesia_r_location, common_labels);

        % Calculate a contrast matrix, this is different than before.
        % now we want to have high score for large differences and low
        % score for similarities
        baseline_vs_anesthesia = abs(baseline_f_dpli - anesthesia_f_dpli);

        % new feature: remove the diagonal from the dPLI (wd = without diagonal)
        % temporary value
        tmp = baseline_vs_anesthesia';
        baseline_vs_anesthesia_wd = reshape(tmp(~eye(size(tmp))), size(baseline_vs_anesthesia, 2)-1, [])';

        dpli_dris(p) = sum(baseline_vs_anesthesia_wd(:));
        % normalize it so that headset with more channels aren't artificially
        % inflated
        dpli_dris(p) = dpli_dris(p)/length(baseline_vs_anesthesia_wd(:));
        % At this point we have the dpli_dri
        
        %% Calculate the Hub-DRI
        disp(strcat("Participant: ", participant , "_HUB"));
        disp(strcat("Participant pair: ", participant_pair));
        
        baseline_file_wpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data',filesep, 'wpli_',FREQUENCY(fre),'_',participant,'_',baseline_state,'.mat');
        %nov23: baseline_file_wpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data_ARI',filesep, 'wpli_ARI',FREQUENCY(fre),'_',participant,'_',baseline_state,'.mat');
        %baseline_file_dpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data',filesep, 'dpli_',FREQUENCY(fre),'_',participant,'_',baseline_state,'.mat');
        anesthesia_file_wpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data',filesep,'wpli_',FREQUENCY(fre),'_',participant,'_',anesthesia_state,'.mat');
        %nov23:anesthesia_file_wpli = strcat(IN_DIR,filesep,participant, filesep,'fc_data_ARI',filesep,'wpli_ARI',FREQUENCY(fre),'_',participant,'_',anesthesia_state,'.mat');
        
     
        [baseline_r_wpli, baseline_r_labels, baseline_r_regions, baseline_r_location] = process_wpli(baseline_file_wpli);
        [anesthesia_r_wpli, anesthesia_r_labels, anesthesia_r_regions, anesthesia_r_location] = process_wpli(anesthesia_file_wpli);
    
        % Filter the matrix to have the same size
        % Get the common location
        [common_labels, common_region] = get_subset_two(baseline_r_labels, anesthesia_r_labels, baseline_r_regions, anesthesia_r_regions);           

        % Filter the matrices to have the same size
        baseline_f_wpli = filter_matrix(baseline_r_wpli, baseline_r_labels, common_labels);
        anesthesia_f_wpli = filter_matrix(anesthesia_r_wpli, anesthesia_r_labels, common_labels);

        % Regenerate the common location from the common labels that was
        % lost in the get_subset process
        common_location = generate_common_location(common_labels, baseline_r_labels, baseline_r_location);

        %% Binarize the three states using the minimally spanning tree and calculate the hub location
        % thereby we calculate the threshold only for baseline data and keep it
        % for the other conditions
        disp("Baseline Threshold: ")
        [threshold] = find_smallest_connected_threshold(baseline_f_wpli, threshold_range);
        [baseline_b_wpli] = binarize_matrix(threshold_matrix(baseline_f_wpli, threshold));

        % here we are using only the degree and not the betweeness centrality
        [~, baseline_weights] = binary_hub_location(baseline_b_wpli, common_location, 1.0, 0.0);
        baseline_norm_weights = (baseline_weights - mean(baseline_weights))  / std(baseline_weights);

        disp("Anesthesia Threshold: ")    
        %[threshold] = find_smallest_connected_threshold(anesthesia_f_wpli, threshold_range);
        [anesthesia_b_wpli] = binarize_matrix(threshold_matrix(anesthesia_f_wpli, threshold));

        % here we are using only the degree and not the betweeness centrality
        [~, anesthesia_weights] = binary_hub_location(anesthesia_b_wpli, common_location,  1.0, 0.0);
        anesthesia_norm_weights = (anesthesia_weights - mean(anesthesia_weights))  / std(anesthesia_weights);

        % Calculate a contrast matrix, this is same thing as for dpli-dri
        % now we want to have high score for large differences and low
        % score for similarities
        bva = abs(baseline_norm_weights - anesthesia_norm_weights);
        disp(strcat("Cosine BvA: ", string(bva)))

        disp("-----")
        % Calculate the hub dri
        hub_dris(p) = sum(bva(:))/length(bva(:));

        plot_dpli_hub_two(baseline_f_dpli, anesthesia_f_dpli, baseline_norm_weights, anesthesia_norm_weights, common_location, participant_pair, OUT_DIR);

    end
 
    %% Plot figure for the dpli-dri
    % outcome labels
    % 0 recovered
    % 1 non-recovered
    
    
    if outcome == "recovery" % when it is recovered vs non-recovered
        handle = figure;
        scatter(dpli_dris(P_LABEL ==0),hub_dris(P_LABEL ==0),[],'red','filled')
        hold on;
        scatter(dpli_dris(P_LABEL ==1),hub_dris(P_LABEL ==1),[],'blue','filled')
    elseif outcome == "GOSE" % when it is scores from 1-8 of GOSE scale
         handle = figure;
        scatter(dpli_dris(P_LABEL ==1),hub_dris(P_LABEL ==1),[],'red','filled')
        hold on;
        scatter(dpli_dris(P_LABEL ==2),hub_dris(P_LABEL ==2),[],'white','filled') 
        scatter(dpli_dris(P_LABEL ==3),hub_dris(P_LABEL ==3),[],'blue','filled')
        scatter(dpli_dris(P_LABEL ==4),hub_dris(P_LABEL ==4),[],'cyan','filled')
        scatter(dpli_dris(P_LABEL ==5),hub_dris(P_LABEL ==5),[],'magenta','filled')
        scatter(dpli_dris(P_LABEL ==6),hub_dris(P_LABEL ==6),[],'yellow','filled')
        scatter(dpli_dris(P_LABEL ==7),hub_dris(P_LABEL ==7),[],'black','filled')
        scatter(dpli_dris(P_LABEL ==8),hub_dris(P_LABEL ==8),[],'green','filled') 
    end

    xlabel("dPLI-DRI");
    ylabel("Hub-DRI");
    %axis([0.005 0.08 0.4 1.4])
    axis([0.005 0.16 0.4 1.4])
    leg = legend('1','2', '3','4','5','6','7','8','Location','northeastoutside');
    title(leg,'GOSE scale')
    
    %legend('Recovered','Non-Recovered');
    
    % Save it to disk
    filename = strcat(OUT_DIR,"/DPLI_DRI_PDEEG_",FREQUENCY(fre),".png");
    saveas(handle,filename);
    close all; 

    % Print out the result in a good format
    for p = 1:length(P_ID)
        participant = P_ID{p};
        label = P_LABEL(p);
        dpli_dri = dpli_dris(p);
        hub_dri = hub_dris(p);

        msg = sprintf("Participant: %s\nLabel: %d\ndPLI-DRI: %f\nHub-DRI: %f\n------\n",participant,label,dpli_dri,hub_dri);
        disp(msg);
    end

    RESULTS = table(P_ID(:), dpli_dris(:) , hub_dris(:), 'VariableNames', { 'ID', 'dPLI','Hub'});
    % Write data to text file
    writetable(RESULTS, strcat(OUT_DIR,'/ARI_',FREQUENCY(fre),'.txt'))
end


