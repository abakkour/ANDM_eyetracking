function [trialsPerRun] = cat_form_probe_pairs(subjectID, order, numRuns)

% function [trialsPerRun] = cat_form_probe_pairs(subjectID, order, numRuns)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik December 2014 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function organizes the matrices for each block of the probe session of the boost
% (cue-approach) task, divided to number of runs as requested (1 or 2 or 4 would
% work. Pay attention that number of runs should be a divisor of number of
% comparisons.

% This function is for the version where only 48 items are being trained,


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order*.txt'' --> created by sort_ratings


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stimuliForProbe_order%d_block_%d_run%d.txt'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'test999';
% order = 1;
% test_comp = 4;
% mainPath = '/Users/schonberglabimac1/Documents/Boost_Israel_New_Rotem_mac';
% numRunsPerBlock = 1;
% block = 1;

rng shuffle

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================
mainPath=pwd;
outputPath = [mainPath '/Output'];
numRunsPerBlock=1;
block=1;

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order' num2str(order) '*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
% bidIndex = data{3};
% bidValue = data{4};
fclose(fid);

%==============================================
%%   'DATA ORGANIZATION'
%==============================================

% determine stimuli to use based on order number
%-----------------------------------------------------------------
switch order
    case 1
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [ 9 12 13  16]; % HV_beep
        HV_nobeep = [10 11 14  15]; % HV_nobeep
        
        LV_beep =   [93 96 97 100]; % LV_beep
        LV_nobeep = [94 95 98  99]; % LV_nobeep
        
        
        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - - 
        sanityHVGO = [17 20]; % sanity_HVGO
        sanityLVGO = [89 92]; % sanity_LVGO
        sanityHVNOGO = [18 19]; % sanity_HVNOGO
        sanityLVNOGO = [90 91]; % sanity_LVNOGO
        
    case 2
        
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [10 11 14  15]; % HV_beep
        HV_nobeep = [ 9 12 13  16]; % HV_nobeep
        
        
        LV_beep =   [94 95 98  99]; % LV_beep
        LV_nobeep = [93 96 97 100]; % LV_nobeep
        
                
        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - - -
        sanityHVGO = [18 19]; % sanity_HVGO
        sanityLVGO = [90 91]; % sanity_LVGO
        sanityHVNOGO = [17 20]; % sanity_HVNOGO
        sanityLVNOGO = [89 92]; % sanity_LVNOGO
        
end % end switch order


%   add multiple iterations of each item presentation
%-----------------------------------------------------


%   TRIAL TYPE 1: HighValue Go vs. HighValue NoGo(Stop)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHVbeepItems = length(HV_beep);
numHVnobeepItems = length(HV_nobeep);

HV_beep_new = repmat(HV_beep,numHVbeepItems,1);
HV_beep_new = HV_beep_new(:)';
HV_nobeep_new = repmat(HV_nobeep,1,numHVnobeepItems);
HV_nobeep_new = HV_nobeep_new(:)';

[shuffle_HV_beep_new,shuff_HV_beep_new_ind] = Shuffle(HV_beep_new);
shuffle_HV_nobeep_new = HV_nobeep_new(shuff_HV_beep_new_ind);



%   TRIAL TYPE 2: LowValue Go vs. LowValue NoGo(Stop)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numLVbeepItems = length(LV_beep);
numLVnobeepItems = length(LV_nobeep);

LV_beep_new = repmat(LV_beep,numLVbeepItems,1);
LV_beep_new = LV_beep_new(:)';
LV_nobeep_new = repmat(LV_nobeep,1,numLVnobeepItems);
LV_nobeep_new = LV_nobeep_new(:)';

[shuffle_LV_beep_new,shuff_LV_beep_new_ind] = Shuffle(LV_beep_new);
shuffle_LV_nobeep_new = LV_nobeep_new(shuff_LV_beep_new_ind);


%   TRIAL TYPE 3+4: HighValue vs. LowValue GO + NOGO
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numSanityHVGOItems = length(sanityHVGO);
numSanityLVGOItems = length(sanityLVGO);

numSanityHVNOGOItems = length(sanityHVNOGO);
numSanityLVNOGOItems = length(sanityLVNOGO);

sanityHVGO_new = repmat(sanityHVGO,numSanityHVGOItems,1);
sanityHVGO_new = sanityHVGO_new(:)';
sanityLVGO_new = repmat(sanityLVGO,1,numSanityLVGOItems);
sanityLVGO_new = sanityLVGO_new(:)';
sanityHVNOGO_new = repmat(sanityHVNOGO,numSanityHVNOGOItems,1);
sanityHVNOGO_new = sanityHVNOGO_new(:)';
sanityLVNOGO_new = repmat(sanityLVNOGO,1,numSanityLVNOGOItems);
sanityLVNOGO_new = sanityLVNOGO_new(:)';


[shuffle_sanityHVGO_new,shuff_sanityHVGO_new_ind] = Shuffle(sanityHVGO_new);
shuffle_sanityLVGO_new = sanityLVGO_new(shuff_sanityHVGO_new_ind);
[shuffle_sanityHVNOGO_new,shuff_sanityHVNOGO_new_ind] = Shuffle(sanityHVNOGO_new);
shuffle_sanityLVNOGO_new = sanityLVNOGO_new(shuff_sanityHVNOGO_new_ind);


%   randomize all possible comparisons for all trial types
%-----------------------------------------------------------------
numComparisonsHV = numHVbeepItems^2;
numComparisonsLV = numLVbeepItems^2;
numComparisons = numComparisonsHV + numComparisonsLV;
numSanity = numSanityHVGOItems^2 + numSanityHVNOGOItems^2;
total_num_trials = numComparisons + numSanity;
trialsPerRun = total_num_trials/numRunsPerBlock;

stimnum1 = zeros(numRuns,trialsPerRun);
stimnum2 = zeros(numRuns,trialsPerRun);
leftname = cell(numRuns,trialsPerRun);
rightname = cell(numRuns,trialsPerRun);
pairType = zeros(numRuns,trialsPerRun);


numComparisonsPerRun = numComparisons/numRunsPerBlock;
numSanityPerRun = numSanity/numRunsPerBlock;
pairType(1:numRuns,1:numComparisonsPerRun/2) = 1;
pairType(1:numRuns,numComparisonsPerRun/2+1:numComparisonsPerRun) = 2;
pairType(1:numRuns,numComparisonsPerRun+1:numComparisonsPerRun+numSanityPerRun/2) = 3;
pairType(1:numRuns,numComparisonsPerRun+numSanityPerRun/2+1:numComparisonsPerRun+numSanityPerRun) = 4;


leftGo = ones(numRuns,total_num_trials./numRunsPerBlock);
leftGo(:,[1:numComparisonsPerRun/4 numComparisonsPerRun/2+1:numComparisonsPerRun*3/4 1+numComparisonsPerRun:numComparisonsPerRun+numSanityPerRun/2]) = 0;
    
for numRun = 1:numRuns
    pairType(numRun,:) = Shuffle(pairType(numRun,:));
    leftGo(numRun,:) = Shuffle(leftGo(numRun,:));
end % end for numRun = 1:numRunsPerBlock

HV_beep = shuffle_HV_beep_new;
HV_nobeep = shuffle_HV_nobeep_new;
LV_beep = shuffle_LV_beep_new;
LV_nobeep = shuffle_LV_nobeep_new;

sanityHVGO = shuffle_sanityHVGO_new;
sanityLVGO = shuffle_sanityLVGO_new;
sanityHVNOGO = shuffle_sanityHVNOGO_new;
sanityLVNOGO = shuffle_sanityLVNOGO_new;


% % Divide the matrices of each comparison to the number of trials
% HV_beep_allRuns = reshape(HV_beep,2,length(HV_beep)/2);
% HV_nobeep_allRuns = reshape(HV_nobeep,2,length(HV_nobeep)/2);
% LV_beep_allRuns = reshape(LV_beep,2,length(LV_beep)/2);
% LV_nobeep_allRuns = reshape(LV_nobeep,2,length(LV_nobeep)/2);
% sanityHV_nobeep_allRuns = reshape(sanityHV_nobeep,2,length(sanityHV_nobeep)/2);
% sanityLV_nobeep_allRuns = reshape(sanityLV_nobeep,2,length(sanityLV_nobeep)/2);



for numRun = 1:numRuns
    HH = 1;
    LL = 1;
    HL_NG = 1;
    HL_G = 1;
    % Create stimuliForProbe.txt for this run
    fid1 = fopen([outputPath '/' sprintf('%s_stimuliForProbe_order%d_block_%d_run%d.txt',subjectID,order,block,numRun)], 'w');
    
    for trial = 1:trialsPerRun % trial num within block      
        switch pairType(numRun,trial)
            case 1

                % HighValue Go vs. HighValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = HV_beep(HH);
                stimnum2(numRun,trial) = HV_nobeep(HH);
                HH = HH+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 2
                
                % LowValue Go vs. LowValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = LV_beep(LL);
                stimnum2(numRun,trial) = LV_nobeep(LL);
                LL = LL+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 3
                
                % HighValue Go vs. LowValue Go
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = sanityHVGO(HL_G);
                stimnum2(numRun,trial) = sanityLVGO(HL_G);
                HL_G = HL_G+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 4
                
                % HighValue NoGo(Stop) vs. LowValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = sanityHVNOGO(HL_NG);
                stimnum2(numRun,trial) = sanityLVNOGO(HL_NG);
                HL_NG = HL_NG+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
        end % end switch pairtype

        fprintf(fid1, '%d\t %d\t %d\t %d\t %s\t %s\t \n', stimnum1(numRun,trial),stimnum2(numRun,trial),leftGo(numRun,trial),pairType(numRun,trial),leftname{numRun,trial},rightname{numRun,trial});
    end % end for trial = 1:total_num_trials
    
    fprintf(fid1, '\n');
    fclose(fid1);
end % end for numRun = 1:numRunsPerBlocks


%---------------------------------------------------------------------
% create a data structure with info about the run and all the matrices
%---------------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_stimuliForProbe_order%d_block_%d_%d_trials_%s.mat',subjectID,order,block,total_num_trials,date));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.script_name = mfilename;

save(outfile);


end % end function

