function food_rating(subjectID,run,use_eyetracker)
%food_rating(subjectID,use_eyetracker) Run food rating task.
%   food_rating(subjectID,use_eyetracker) runs the food tating task,
%   food_rating outputs a text file named 'subjectID_food_rating.txt'
%   that contains the ratings for each item for subject subjectID.
%

% Code written by Akram Bakour 06/2016 based on bits of existing code.
% Contributors include Tom Schonberg, Rotem Botvinik, Tom Salomon

%=============================================================================

%---------------------------------------------------------------
%%   'GLOBAL VARIABLES'
%---------------------------------------------------------------

outputPath = 'Output/';

% essential for randomization
rng('shuffle');

% about timing
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

% -----------------------------------------------
%% 'INITIALIZE SCREEN'
% -----------------------------------------------

Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
Screen('Preference', 'SkipSyncTests', 1); %ONLY FOR TESTING
screennum = min(Screen('Screens')); %select external screen

pixelSize=32;
[w, windowRect] = Screen('OpenWindow',screennum,[],[],pixelSize);
[xCenter, yCenter] = RectCenter(windowRect);


%   colors
% - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
blue = [0 0 255];

Screen('FillRect', w, black);
Screen('Flip', w);

%   text
% - - - - - -
theFont='Arial';
Screen('TextSize',w,36);
Screen('TextFont',w,theFont);
Screen('TextColor',w,white);

WaitSecs(1);

%---------------------------------------------------------------
%%   Load in food images
%---------------------------------------------------------------

shuff_names=Shuffle(dir('stim/*.jpg'));
Images = cell(1, length(shuff_names));
imx=0;
for i=1:length(shuff_names)
    imx=imx+2;
    Images{i}=imread(['stim/',shuff_names(i).name]);
    CenterText(w,'Loading images...',white,0,-100);
    Screen('FillRect', w, white, [xCenter-108 yCenter-20 xCenter-108+imx yCenter+20]); 
    Screen(w,'Flip');
end

baseRect = [0 0 400 20];
pointerRect = [0 0 5 30];

%-----------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
% use_eyetracker=1; % set to 1/0 to turn on/off eyetracker functions
if use_eyetracker
    dummymode=0;
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(w);
    % Disable key output to Matlab window:
    
    el.backgroundcolour = black;
    el.backgroundcolour = black;
    el.foregroundcolour = white;
    el.msgfontcolour    = white;
    el.imgtitlecolour   = white;
    el.calibrationtargetcolour = el.foregroundcolour;
    EyelinkUpdateDefaults(el);
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end;
    
    [~,vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile=[subjID '_food_rating_run' num2str(run) '.edf'];
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    WaitSecs(2);
    %     % STEP 5
    %     % start recording eye position
    %     Eyelink('StartRecording');
    %     % record a few samples before we actually start displaying
    %     WaitSecs(0.1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish Initialization %
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------
KbQueueCreate;
Screen('TextSize',w, 40);

CenterText(w,'Rate how much you like these items from 0 to 10', white, 0,-100);
CenterText(w,'using the mouse.', white, 0,0);
CenterText(w,'Click anywhere to continue', white, 0,150);
Screen(w,'Flip');
WaitSecs(0.01);

[~,~,~,~]=GetClicks(w,0);
WaitSecs(0.01);

Screen('TextSize',w, 60);
CenterText(w,'+', white,0,0);
runStartTime=Screen(w,'Flip');

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    WaitSecs(.05);

    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message', ['SYNCTIME start: ', num2str(runStartTime)]); % mark start time in file
end

WaitSecs(2); % Wait 2 sec before first stimulus


%   'Write output file header'
%---------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID '_food_rating_run' num2str(run) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\t onsetTime\t itemName\t rating\t RT\n'); %write the header line

%	pre-allocating matrices and setting defaults
%---------------------------
respTime(1:length(shuff_names),1) = 999;
rating(1:length(shuff_names),1) = 999;
image_start_time(1:length(shuff_names),1) = 999;
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter+350); %main scale

%Run trial loop
for trialNum = 1:length(shuff_names)   % To cover all the items in one run.
    xpos=Shuffle(-100:100); %want the pointer to start in a random position every trial
    centeredPointer=CenterRectOnPointd(pointerRect, xCenter+xpos(1), yCenter+350);
    SetMouse(xCenter, yCenter,w); %place cursor in middle of screen
    Screen('TextSize',w, 40);
    Screen('PutImage',w,Images{trialNum});
    Screen('FillRect', w, white, centeredRect);
    Screen('FillRect', w, blue, centeredPointer);
    CenterText(w,'0',white,-200,360);
    CenterText(w,'5',white,0,360);
    CenterText(w,'10',white,200,360);
    image_start_time(trialNum) = Screen(w,'Flip'); % display images according to Onset times    
    vbl=image_start_time;
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message', ['trial: ' num2str(trialNum) ' stim: ',shuff_names(trialNum).name,' start_time: ',num2str(image_start_time)]); % mark start time in file
    end
    
    
    %---------------------------------------------------
    %% Move pointer with cursor
    %---------------------------------------------------
    noclick=1;
    while noclick
        % Get the current position of the mouse
        [mx, my, buttons] = GetMouse(w);
        
        % See if the mouse cursor is inside the square
        inside = IsInRect(mx, my, centeredRect);
        
        if inside==1
            centeredPointer=CenterRectOnPointd(pointerRect, mx, yCenter+350);
            Screen('PutImage',w,Images{trialNum});
            Screen('FillRect', w, white, centeredRect);
            Screen('FillRect', w, blue, centeredPointer);
            CenterText(w,'0',white,-200,360);
            CenterText(w,'5',white,0,360);
            CenterText(w,'10',white,200,360);
            vbl=Screen(w,'Flip');
            if sum(buttons) > 0
                rating(trialNum)=(mx-xCenter+200)/40;
                respTime(trialNum)=vbl-image_start_time(trialNum);
                noclick=0;
            end
        end
        
        
    
    end %%% End big while waiting for response 
    
    
    %   Show fixation
    %---------------------------
    Screen('TextSize',w, 60);
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip');
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message', ['trial: ' num2str(trialNum) ' Fixation_ITI_start: ',num2str(GetSecs)]); % mark start time in file
    end
 
    %   'Save data'
    %---------------------------
    fprintf(fid1,'%s\t %d\t %s\t %.2f\t %d\n', subjectID, image_start_time(trialNum)-runStartTime, shuff_names(trialNum).name, rating(trialNum), respTime(trialNum)*1000);
    
    WaitSecs(1); %1 sec ITI
    
end; %	End the big trialNum loop showing all the images in one run.


%---------------------------------------------------------------
%%   save data to a .mat file & close out
%---------------------------------------------------------------
outfile = strcat(outputPath, '/', subjectID,'_food_rating_run_', num2str(run), '_', timestamp,'.mat');
% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.revision_date = revision_date;
run_info.script_name = mfilename;
clear Images;

save(outfile);


if use_eyetracker
    %---------------------------------------------------------------
    %%   Finishing eye tracking system %
    %---------------------------------------------------------------
    
    % STEP 7
    %---------------------------
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    WaitSecs(.1);
    Eyelink('CloseFile');
    
    % download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch rdf
        fprintf('Problem receiving data file ''%s''\n', edfFile );
        rdf;
    end
    
    
    if dummymode==0
        movefile(edfFile,['./Output/', subjectID,'_food_rating_run_',num2str(run),'_',timestamp,'.edf']);
    end;
end

%   outgoing msg & closing
% ------------------------------
CenterText(w,'Great Job. Thank you!',white, 0,-270);
CenterText(w,'Now we will continue to the next part', white, 0, -180);
Screen('Flip',w);

WaitSecs(4);

Screen('CloseAll');

clear all

