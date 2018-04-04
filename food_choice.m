function food_choice(subjectID, numRun, use_eyetracker)

% function cat_probe(subjectID, numRun, use_eyetracker)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ====================== by Rotem Botvinik May 2015 =======================
% =============== Modified by Tom Salomon, September 2015 =================
% ================= Modified by Akram Bakkour, June 2016 ==================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the probe session of the boost (cue-approach) task. In
% this session, the subject is presented with two items in each trial, and
% should choose which one he prefers within 1.5 seconds. At the end of the
% experiment the function 'probeResolve_Israel' chooses a random trial from
% this probe session and the subject is given the item he chose on that
% chosen trial.
% This function runs one run each time. Each block is devided to
% 'numRunsPerBlock' runs. The stimuli for each run are shuffled, chosen and
% organized in the 'organizeProbe_Israel' function.

% This function is a version in which only the 40 of the items are included
% in the training


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''stopGoList_allstim_order*.txt'' --> created by sortBDM_Israel
%   ''stimuliForProbe_order%d_block_%d_run%d.txt'' --> Created by
%   organizeProbe_Israel
%   onset lists


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''probe_block_' block '_' timestamp '.txt''
%


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'test999';
% order = 1;
% test_comp = 4;
% sessionNum = 1;
% mainPath = '/Users/schonberglabimac1/Documents/Boost_Israel_New_Rotem_mac';
% block = 1;
% numRun = 1;
% trialsPerRun = 8; % for debugging

rng('Default')

% =========================================================================
% Get input args and check if input is ok
% =========================================================================

sessionNum=1;

trialsPerRun = 70; %70 trials per run x 3 runs = 210 trials total
numRunsPerBlock = 1;

% =========================================================================
% set the computer and path
% =========================================================================
test_comp=2; % 1-MRI experiment

% Set main path
mainPath=pwd;

outputPath = [mainPath '/Output'];

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================

% 'block' indicates how many times the subject has completed all
% trials of the probe experiment. Therefore, the first time the subject completes
% a probe block on his/her 2nd session, it is actually probe block 3 for
% that person:

block = 1;

%   'timestamp'
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

%   'set phase times'
% - - - - - - - - - - - - - - - - -
maxtime = 2.5;      % 2.5 second limit on each selection
baseline_fixation_dur = 2; % Need to modify based on if first few volumes are saved or not
afterrunfixation = 2;

fixationTime = zeros(trialsPerRun,1); % for later saving fixation times for each trial

% define the size factor by which the image size and the distance between the images will be reduced
sizeFactor = 1;

trial_time_fixated_left = 999;
trial_time_fixated_right = 999;
trial_time_fixated_mid = 999;
trial_num_left_fixations = 999;
trial_num_right_fixations = 999;
trial_num_mid_fixations = 999;
first_fixation_duration = 999;
first_fixation_area = 'x';


%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
%Screen('Preference', 'SkipSyncTests', 1); %ONLY FOR TESTING
screennum = min(Screen('Screens'));

pixelSize = 32;
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

HideCursor;


% Define Colors
% - - - - - - - - - - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);


% setting the text
% - - - - - - - - - - - - - - -
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

% stack locations
% - - - - - - - - - - - - - - -
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

penWidth = 10;

HideCursor;

%==============================================
%% 'ASSIGN response keys'
%==============================================
KbName('UnifyKeyNames');

if test_comp == 1
    leftstack = 'b';
    rightstack = 'y';
    badresp = 'x';
else
    leftstack = 'u';
    rightstack = 'i';
    badresp = 'x';
end

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_rated_food_allstim_order*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
% bidIndex = data{3};
%bidValue = data{4};
fclose(fid);

%   'read in organized list of stimuli for this run'
% - - - - - - - - - - - - - - - - - - - - - - - - - -

fid = fopen([outputPath '/' sprintf('%s_stimuliForFoodChoice_run%d.txt',subjectID,numRun)]);
stimuli = textscan(fid, '%d %d %s %s %d %d %.2f %.2f') ;% predetermined choice trials in form_food_choice_pairs
stimnum1 = stimuli{5};
stimnum2 = stimuli{6};
stim1rating = stimuli{7};
stim2rating = stimuli{8};
leftHV = stimuli{1};
pairNumber = stimuli{2};
leftname = stimuli{3};
rightname = stimuli{4};
fclose(fid);


%   'load image arrays'
% - - - - - - - - - - - - - - -
Images = cell(1,length(stimName));
for i = 1:length(stimName)
    Images{i} = imread([mainPath sprintf('/stim/%s',stimName{i})]);
end

% Define image scale - Change according to your stimuli
% - - - - - - - - - - - - - - -
stackH= sizeFactor*size(Images{1},1);
stackW= sizeFactor*size(Images{1},2);
%
% stackW = 576*sizeFactor;
% stackH = 432*sizeFactor;
distcent = 300*sizeFactor; % half of the distance between the images
leftRect = [xcenter-stackW-distcent ycenter-stackH/2 xcenter-distcent ycenter+stackH/2];
rightRect = [xcenter+distcent ycenter-stackH/2 xcenter+stackW+distcent ycenter+stackH/2];
midRect = [xcenter-150 ycenter-150 xcenter+150 ycenter+150];


%   'load onsets'
% - - - - - - - - - - - - - - -
onsetlist = 0:3.5:245; %70 trials at 3.5 sec per trial

%-----------------------------------------------------------------
%% 'Write output file header'
%-----------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID sprintf('_food_choice_run%d_', numRun) timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\trun\ttrial\tonsettime\tImageLeft\tImageRight\tratingOrderLeft\tratingOrderRight\tIsleftHV\tResponse\tPairNumber\tOutcome\tRT\tratingLeft\tratingRight\t\fixationTime\ttimeFixMid\ttimeFixLeft\ttimeFixRight\tnumMidFix\tnumLeftFix\tnumRightFix\tfirstFix\tfirstFixTime\n'); %write the header line    
    
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
    
    [~, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % open file to record data to
    edfFile=['foodch' num2str(numRun) '.edf'];
    i=Eyelink('Openfile', edfFile);
    
    if i~=0
        fprint('Cannot create EDF file ''%s'' ', edfFile);
        Eyelink('Shutdown');
        sca;
        return;
    end
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERS,STATUS,INPUT');
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERS,AREA,HTARGET,STATUS,INPUT');
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
    %     % STEP 5
    %     % start recording eye position
    %     Eyelink('StartRecording');
    %     % record a few samples before we actually start displaying
    %     WaitSecs(0.1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish Initialization %
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

%==============================================
%% 'Display Main Instructions'
%==============================================
ListenChar(2); %suppress keyboard output to command line
KbQueueCreate;
Screen('TextSize',w, 40);

%%% While they are waiting for the trigger
if test_comp == 1
    
    CenterText(w,'We will continue with more of the food chioce task.',white, 0, -100);
    
    CenterText(w,'Please press any button to continue ...',white, 0, 100);
    Screen(w,'Flip');
    
    % wait for the subject to press the button
    KbPressWait(-1);
    
    CenterText(w,'Waiting for trigger...GET READY....', white, 0, 50);
    Screen(w,'Flip');
    
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    
    DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected
    
else
    if numRun == 1
        if  block == 1 % If this is the first run of the first block, show instructions
            
            CenterText(w,'We will continue with more of the food chioce task.',white, 0, -100);
            CenterText(w,'Please press any button to continue ...',white, 0, 100);
            Screen(w,'Flip');
            % wait for the subject to press the button
            KbPressWait(-1);
            
        else % this is the first run but not the first block
            CenterText(w,'We will continue with more of the food chioce task.',white, 0, -100);
            CenterText(w,'Please press any button to continue ...', white,0,100);
            Screen('Flip',w);
            
            % wait for the subject to press the button
            KbPressWait(-1);
            
        end % end if block == 1
    else % if this is not the first run of the block
        Screen('TextSize',w, 40);
        CenterText(w,'Another run begins now.', white,0,-200);
        Screen('Flip',w);
        WaitSecs(2);
    end % end if numRun == 1
    
end % end if test_comp==1


%   baseline fixation cross
% - - - - - - - - - - - - -
prebaseline = GetSecs;
% baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
while GetSecs < prebaseline+baseline_fixation_dur
    %    Screen(w,'Flip', anchor);
    Screen('TextSize',w, 60);
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip');
    
end
postbaseline = GetSecs;
baseline_fixation = postbaseline - prebaseline;


%==============================================
%% 'Run Trials'
%==============================================

runtrial = 1;
runStart = GetSecs;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    WaitSecs(.1);
    
    %   Eyelink MSG
    % ---------------------------
    % messages to save on each trial ( trial number, onset and RT)
    Eyelink('Message',['SYNCTIME at run start:',num2str(GetSecs)]); % mark start time in file
    if ~dummymode
        eye_used = Eyelink('EyeAvailable');
        if eye_used == -1
            fprintf('Eyelink aborted - could not find which eye being used.\n');
            cleanup;
        end
    end
end

% for trial = 1:5 % for debugging
for trial = 1:trialsPerRun
    
    
    % initial box outline colors
    % - - - - - - -
    colorLeft = black;
    colorRight = black;
    out = 999;
    
    
    %-----------------------------------------------------------------
    % display images
    %-----------------------------------------------------------------
    Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
    Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
    %eyelink_message=['block: ' num2str(block), ' run: ' num2str(numRun),' trial: ',num2str(trial),', StimLeft: ',stimName{stimnum1(trial)},' StimRight: ',stimName{stimnum2(trial)},' time: ',num2str(GetSecs)];
    eyelink_message=['run: ' num2str(numRun),' trial: ',num2str(trial)];
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',eyelink_message);
        trial_time_fixated_left = 0;
        trial_time_fixated_right = 0;
        trial_time_fixated_mid = 0;
        trial_num_left_fixations = 0;
        trial_num_right_fixations = 0;
        trial_num_mid_fixations = 0;
        
        % current_area determines which area eye is in (left, right, neither)
        % xpos and ypos are used for eyepos_debug
        [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,midRect,leftRect,rightRect);
        
        % last_area will track what area the eye was in on the previous loop
        % iteration, so we can determine when a change occurs
        % fixation_onset_time stores the time a "fixation" into an area began
        last_area=current_area;
        fixation_onset_time = GetSecs;
        
        % tracking first fixation
        first_fixation_duration = 0;
        first_fixation_area = current_area; % this will report 'n' in output if they never looked at an object
        first_fixation_flag = (first_fixation_area=='m' || first_fixation_area=='l' || first_fixation_area=='r'); % flags 1 once the first fixation has occurred, 2 once the first fixation has been processed
        first_fixation_onset = fixation_onset_time;
    end
    
    CenterText(w,'+', white,0,0);
    StimOnset = Screen(w,'Flip', runStart+onsetlist(runtrial)+baseline_fixation);
    
    KbQueueFlush;
    KbQueueStart;
    
    
    %-----------------------------------------------------------------
    % get response
    %-----------------------------------------------------------------
    
    noresp = 1;
    goodresp = 0;
    while noresp
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck;
        
        if keyIsDown && noresp
            if use_eyetracker
                %   Eyelink MSG
                % ---------------------------
                Eyelink('Message',['block: ',num2str(block),' run: ',num2str(numRun),' trial: ',num2str(trial),' Press_time: ',num2str(GetSecs)]);
            end
            
            keyPressed = KbName(firstPress);
            if ischar(keyPressed) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed = char(keyPressed);
                keyPressed = keyPressed(1);
            end
            switch keyPressed
                case leftstack
                    respTime = firstPress(KbName(leftstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
                case rightstack
                    respTime = firstPress(KbName(rightstack))-StimOnset;
                    noresp = 0;
                    goodresp = 1;
            end
        end % end if keyIsDown && noresp
        
        if use_eyetracker==1
            % get eye position
            [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,midRect,leftRect,rightRect);
            
            % they are looking in a new area
            % Currently has initial fixation problems? (color, count, etc.)
            if current_area~=last_area
                % update timings
                switch last_area
                    case 'm'
                        trial_time_fixated_mid = trial_time_fixated_mid + (GetSecs-fixation_onset_time);
                        trial_num_mid_fixations = trial_num_mid_fixations + 1;
                    case 'l'
                        trial_time_fixated_left = trial_time_fixated_left + (GetSecs-fixation_onset_time);
                        trial_num_left_fixations = trial_num_left_fixations + 1;
                    case 'r'
                        trial_time_fixated_right = trial_time_fixated_right + (GetSecs-fixation_onset_time);
                        trial_num_right_fixations = trial_num_right_fixations + 1;
                end
                
                fixation_onset_time=GetSecs;
                
                % they have looked away from their first fixation: record its
                % duration and the target (left/right)
                if(first_fixation_flag==1)
                    %outstr=['first fixation lasted ' GetSecs-first_fixation_onset ' seconds'];
                    %Eyelink('Message',outstr);
                    first_fixation_duration = GetSecs-first_fixation_onset;
                    first_fixation_flag = 2;
                end
                
                % this is their first time fixating on an object this trial
                if(first_fixation_flag==0 && (current_area=='m' || current_area=='l' || current_area=='r'))
                    %outstr=['first fixation on ' last_area];
                    %Eyelink('Message',outstr);
                    first_fixation_flag = 1;
                    first_fixation_onset = fixation_onset_time;
                    first_fixation_area = current_area;
                end
            end
            
            last_area = current_area;
            fixation_duration = GetSecs-fixation_onset_time;
        end
        
        % check for reaching time limit
        if noresp && GetSecs-runStart >= onsetlist(runtrial)+baseline_fixation+maxtime
            noresp = 0;
            keyPressed = badresp;
            respTime = maxtime;
        end
    end % end while noresp
    
    %-----------------------------------------------------------------
    % determine what bid to highlight
    %-----------------------------------------------------------------
    if use_eyetracker==1
        % Eyelink msg
        % - - - - - - -
        rtmsg = strcat('RT = ',num2str(respTime));
        Eyelink('Message',rtmsg);
        switch last_area
            case 'm'
                trial_time_fixated_mid = trial_time_fixated_mid + fixation_duration;
                trial_num_mid_fixations = trial_num_mid_fixations + 1;
            case 'l'
                trial_time_fixated_left = trial_time_fixated_left + fixation_duration;
                trial_num_left_fixations = trial_num_left_fixations + 1;
            case 'r'
                trial_time_fixated_right = trial_time_fixated_right + fixation_duration;
                trial_num_right_fixations = trial_num_right_fixations + 1;
        end
        % time limit reached while fixating on first fixated object
        if(first_fixation_flag==1)
            %outstr=['first fixation lasted ' GetSecs-first_fixation_onset ' seconds'];
            %Eyelink('Message',outstr);
            first_fixation_duration = GetSecs-first_fixation_onset;
            first_fixation_flag = 2;
        end
    end
    
    switch keyPressed
        case leftstack
            colorLeft = green;
            if leftHV(trial) == 0
                out = 0;
            else
                out = 1;
            end
        case rightstack
            colorRight = green;
            if leftHV(trial) == 1
                out = 0;
            else
                out = 1;
            end
    end
    
    if goodresp==1
        
        Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
        Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
        
        switch keyPressed
            case leftstack
                Screen('FrameRect', w, green, leftRect, penWidth);
            case rightstack
                Screen('FrameRect', w, green, rightRect, penWidth);
        end
        
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(trial)+respTime+baseline_fixation);
        
    else
        %         Screen('DrawText', w, 'You must respond faster!', xcenter-400, ycenter, white);
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
        if use_eyetracker
            %   Eyelink MSG
            % ---------------------------
            Eyelink('Message',['block: ',num2str(block),' run: ',num2str(numRun),' trial: ',num2str(trial),' Respond_faster_time: ',num2str(GetSecs)]);
        end
        
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+baseline_fixation);
    end % end if goodresp==1
    
    
    %-----------------------------------------------------------------
    % show fixation ITI
    %-----------------------------------------------------------------
    
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5+baseline_fixation);
    fixationTime(trial) = GetSecs - runStart;
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['block: ',num2str(block),' run: ',num2str(numRun),' trial: ',num2str(trial),' Fixation_ITI_time: ',num2str(GetSecs)]);
    end
    
    if goodresp ~= 1
        respTime = 999;
    end
    
    %-----------------------------------------------------------------
    % 'Save data'
    %-----------------------------------------------------------------
    
    fprintf(fid1,'%s\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %f\t %.2f\t %.2f\t %d\t %f\t %f\t %f\t %d\t %d\t %d\t %c\t %f\n', ...
        subjectID, numRun, runtrial, StimOnset-runStart, char(leftname(trial)), char(rightname(trial)), stimnum1(trial), stimnum2(trial), ...
        leftHV(trial), keyPressed, pairNumber(trial), out, respTime, stim1rating(trial), stim2rating(trial), fixationTime(trial), ... 
        trial_time_fixated_mid, trial_time_fixated_left, trial_time_fixated_right, ...
        trial_num_mid_fixations, trial_num_left_fixations, trial_num_right_fixations,...
        first_fixation_area, first_fixation_duration);
    
    runtrial = runtrial+1;
    %     KbQueueFlush;
    
end % loop through trials
fclose(fid1);


Postexperiment = GetSecs;

if use_eyetracker
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message',['block: ',num2str(block),' run: ',num2str(numRun),' trial: ',num2str(trial),' Part_end_time: ',num2str(GetSecs)]);
end


while GetSecs < Postexperiment + afterrunfixation;
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end

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
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
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
        movefile(edfFile,['./Output/', subjectID,'_Food_Choice_EyeTracking_run_',num2str(numRun),'_' timestamp,'.edf']);
    end
end

%-----------------------------------------------------------------
%	display outgoing message
%-----------------------------------------------------------------
WaitSecs(2);
Screen('FillRect', w, black);
Screen('TextSize',w, 40);

CenterText(w,sprintf('Great Job! We will continue shortly.') ,white,0,-100);
CenterText(w,sprintf('Take a short break then get the experimenter.') ,white,0,100);
Screen('Flip',w);

WaitSecs(3);


%---------------------------------------------------------------
% create a data structure with info about the run
%---------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_probe_block_%2d_run_%2d_%s.mat',subjectID,block,numRun,timestamp));

% create a data structure with info about the run
run_info.subject=subjectID;
run_info.date=date;
run_info.outfile=outfile;

run_info.script_name=mfilename;
clear Images;
save(outfile);

KbQueueFlush;
ShowCursor;
ListenChar(0); %listen to keyboard again.
Screen('CloseAll');
end

function [current_area,  xpos, ypos] = get_current_fixation_area(dummymode,el,eye_used,midRect,leftRect,rightRect)
xpos = 0;
ypos = 0;
if ~dummymode
    evt=Eyelink('NewestFloatSample');
    x=evt.gx(eye_used+1);
    y=evt.gy(eye_used+1);
    if(x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0)
        xpos=x;
        ypos=y;
    end
else % in dummy mode use mousecoordinates
    [xpos,ypos] = GetMouse;
end

% check what area the eye is in
if IsInRect(xpos,ypos,midRect)
    current_area='m';
elseif IsInRect(xpos,ypos,leftRect)
    current_area='l';
elseif IsInRect(xpos,ypos,rightRect)
    current_area='r';
else
    current_area='n';
end
return
end


