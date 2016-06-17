function [] = sort_ratings(subjectID,order)

% function [] = sort_ratings(subjectID,order)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik May 2015 =========================
% =============== modified by Akram Bakkour June 2016 =====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function sorts the stimuli according to the BDM results.
% This function is a version in which only 40 of the items are included
% in the training


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   [mainPath '\Output\' subjectID '_food_rating_*.txt']


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order
%   'stopGoList_trainingstim.txt' ---> The file for training 48 items


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% outputPath = '~/Dropbox/ANDM/Output';
% subjectID = 'test';
% order = 1;

%=========================================================================
%%  read in info from ratings.txt
%=========================================================================

outputPath='Output/';
files=dir([outputPath '/' subjectID '_food_rating_run2*.txt']);
fid = fopen([outputPath files(length(files)).name]); %in case there are multiple files, take the last one
rating_data = textscan(fid, '%s%d%s%.2f%d' , 'HeaderLines', 1); %read in data as new matrix   
fclose(fid);


%=========================================================================
%%  Create matrix sorted by descending ratings value
%========================================================================

[ratings_sorted,trialnum_sort_byrating] = sort(rating_data{4},'descend');

ratings_sortedM(:,1) = trialnum_sort_byrating; % trialnums organized by descending rating
ratings_sortedM(:,2) = ratings_sorted; % ratings sorted large to small
ratings_sortedM(:,3) = 1:108; % stimrank

stimnames_sorted_by_rating = rating_data{3}(trialnum_sort_byrating);


%=========================================================================
%%   The ranking of the stimuli determine the stimtype
%=========================================================================

if order == 1
    
    ratings_sortedM([              9 12 13  16     17 20                                                                    ], 4) = 11; % HV_beep
    ratings_sortedM([     3  6    10 11 14  15     18 19     25 26 31 32 37 38 43 44 49 50                                  ], 4) = 12; % HV_nobeep
    ratings_sortedM([    89 92    93 96 97 100                                                                              ], 4) = 22; % LV_beep
    ratings_sortedM([    59 60 65 66 71 72 77 78 83 84     90 91     94 95 98 99   103 106                                  ], 4) = 24; % LV_nobeep
    ratings_sortedM([ 1:2 4:5 7:8 21:24 27:30 33:36 39:42 45:48 51:58 61:64 67:70 73:76 79:82 85:88 101:102 104:105 107:108 ], 4) = 0; % notTrained
    
    
elseif order == 2
    
    ratings_sortedM([             10 11 14 15     18 19                                                                     ], 4) = 11; % HV_beep
    ratings_sortedM([     3  6     9 12 13 16     17 20     25 26 31 32 37 38 43 44 49 50                                   ], 4) = 12; % HV_nobeep
    ratings_sortedM([    90 91    94 95 98 99                                                                               ], 4) = 22; % LV_beep
    ratings_sortedM([    59 60 65 66 71 72 77 78 83 84     89 92    93 96 97 100    103 106                                 ], 4) = 24; % LV_nobeep
    ratings_sortedM([ 1:2 4:5 7:8 21:24 27:30 33:36 39:42 45:48 51:58 61:64 67:70 73:76 79:82 85:88 101:102 104:105 107:108 ], 4) = 0; % notTrained
    
else
    print('\n\n order number must be 1 or 2 \n\n');
end % end if order == 1

itemsForTraining = ratings_sortedM([3 6 9:20 25:26 31:32 37:38 43:44 49:50 59:60 65:66 71:72 77:78 83:84 89:100 103 106],:);
itemsNamesForTraining = stimnames_sorted_by_rating([3 6 9:20 25:26 31:32 37:38 43:44 49:50 59:60 65:66 71:72 77:78 83:84 89:100 103 106]);

%=========================================================================
%%  create stopGoList_allstim.txt
%   this file is used during probe
%=========================================================================

fid2 = fopen([outputPath '/' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)], 'w');    

for i = 1:length(ratings_sortedM)
    fprintf(fid2, '%s\t%d\t%d\t%d\t%d\t\n', stimnames_sorted_by_rating{i,1},ratings_sortedM(i,4),ratings_sortedM(i,3),ratings_sortedM(i,2),ratings_sortedM(i,1)); 
end
fprintf(fid2, '\n');
fclose(fid2);

fid3 = fopen([outputPath '/' subjectID '_stopGoList_trainingstim.txt'], 'w');    

for i = 1:length(itemsForTraining)
    fprintf(fid3, '%s\t%d\t%d\t%d\t%d\t\n', itemsNamesForTraining{i,1},itemsForTraining(i,4),itemsForTraining(i,3),itemsForTraining(i,2),itemsForTraining(i,1)); 
end
fprintf(fid3, '\n');
fclose(fid3);

end % end function