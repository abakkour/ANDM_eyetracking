function form_food_choice_pairs(subjectID, order,numRuns)

% function [trialsPerRun] = form_choice_pairs(subjectID, numRuns)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ======================= by Akram Bakkour June 2016 ======================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function is for the version where 60 items are used for food choice

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order*.txt'' --> created by sort_ratings


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stimuliFor_food_choice_run%d.txt'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'test';
% order = 1;
% mainPath = '~/Dropbox/ANDM';
% numRunsPerBlock = 1;
% block = 1;

rng shuffle

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================
mainPath=pwd;
outputPath = [mainPath '/Output'];
trialsPerRun=70;

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order' num2str(order) '*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sort_ratings
stimName = data{1};
rankOrder = data{3};
rating = data{4};
fclose(fid);

file = [mainPath '/FoodPairs.txt'];
fid = fopen(file);
foodpairs = textscan(fid, '%d %d %d %d') ;% these contain everything from the sort_ratings
pairnumber=foodpairs{1};
stim1rankorder=foodpairs{2};
stim2rankorder=foodpairs{3};


%==============================================
%%   'DATA ORGANIZATION'
%==============================================

for i=1:length(pairnumber)
    stim1name(i)=stimName(rankOrder==stim1rankorder(i));
    stim1rating(i)=rating(rankOrder==stim1rankorder(i));
    stim2name(i)=stimName(rankOrder==stim2rankorder(i));
    stim2rating(i)=rating(rankOrder==stim2rankorder(i));
end

[shuff_pairnumber,shuff_index]=Shuffle(pairnumber);
shuff_stim1rankorder=stim1rankorder(shuff_index);
shuff_stim2rankorder=stim2rankorder(shuff_index);
shuff_stim1name=stim1name(shuff_index);
shuff_stim2name=stim2name(shuff_index);
shuff_stim1rating=stim1rating(shuff_index);
shuff_stim2rating=stim2rating(shuff_index);

leftHVs=[ones(1,35) zeros(1,35)];
leftHV = [Shuffle(leftHVs) Shuffle(leftHVs) Shuffle(leftHVs)];

count=0;

for numRun = 1:numRuns
    
    % Create stimuliForProbe.txt for this run
    fid1 = fopen([outputPath '/' sprintf('%s_stimuliForFoodChoice_run%d.txt',subjectID,numRun)], 'w');
    
    for trial = 1:trialsPerRun % trial num within block      
        count=count+1;        
        if leftHV(count) == 1
            leftname(numRun,trial) = shuff_stim1name(count);
            rightname(numRun,trial) = shuff_stim2name(count);
            leftrankOrder(numRun,trial) = shuff_stim1rankorder(count);
            rightrankOrder(numRun,trial) = shuff_stim2rankorder(count);
            leftrating(numRun,trial) = shuff_stim1rating(count);
            rightrating(numRun,trial) = shuff_stim2rating(count);
        else
            leftname(numRun,trial) = shuff_stim2name(count);
            rightname(numRun,trial) = shuff_stim1name(count);
            leftrankOrder(numRun,trial) = shuff_stim2rankorder(count);
            rightrankOrder(numRun,trial) = shuff_stim1rankorder(count);
            leftrating(numRun,trial) = shuff_stim2rating(count);
            rightrating(numRun,trial) = shuff_stim1rating(count);
        end
        pairnum(numRun,trial)=shuff_pairnumber(count);
                
        fprintf(fid1, '%d\t %d\t %s\t %s\t %d\t %d\t %.2f\t %.2f \n', leftHV(count),pairnum(numRun,trial),leftname{numRun,trial},rightname{numRun,trial},leftrankOrder(numRun,trial),rightrankOrder(numRun,trial),leftrating(numRun,trial),rightrating(numRun,trial));
    end % end for trial = 1:total_num_trials
    
    fprintf(fid1, '\n');
    fclose(fid1);
end % end for numRun = 1:numRunsPerBlocks


%---------------------------------------------------------------------
% create a data structure with info about the run and all the matrices
%---------------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_stimuliForFoodChoice_%s.mat',subjectID,date));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.script_name = mfilename;

save(outfile);


end % end function

