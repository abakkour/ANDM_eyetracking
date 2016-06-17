function cat_probe(subjectID, order, run, use_eyetracker)

% function cat_probe(subjectID, order, run, use_eyetracker)
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

rng shuffle

% =========================================================================
% Get input args and check if input is ok
% =========================================================================

trialsPerRun = 40; % 4x4 HGOvHNOGO + 4x4 LGOvLNOGO + 2x2 HGOvLGO + 2x2 HNOGOvLNOGO

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

%   'timestamp'
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

%   'set phase times'
% - - - - - - - - - - - - - - - - -
maxtime = 2.5;      % 1.5 second limit on each selection

fixationTime = zeros(trialsPerRun,1); % for later saving fixation times for each trial

% define the size factor by which the image size and the distance between the images will be reduced
sizeFactor = 1;

%==============================================
%% 'INITIALIZE Screen variables'
%==============================================
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
Screen('Preference', 'SkipSyncTests', 1); %ONLY FOR TESTING
screennum = min(Screen('Screens'));

pixelSize = 32;
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

HideCursor;
ListenChar(2);

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
bidValue = data{4};
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

%   'load onsets'
% - - - - - - - - - - - - - - -
%r = Shuffle(1:4);
onsetlist = 2:3.5:170; % 48 trials x 3.5sec trial length = 168 + 2sec = 170

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
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile='cat_probe.edf';
    Eyelink('Openfile', edfFile);
    
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
KbQueueCreate;
Screen('TextSize',w, 40);

CenterText(w,'Choose one item using buttons `u` and `i`',white, 0, -100);
CenterText(w,'for left and right respectively.',white, 0, -50);
CenterText(w,'You have 2 seconds to make a choice.',white, 0, 0);
CenterText(w,'Please press any button to continue ...',white, 0, 100);
Screen(w,'Flip');
% wait for the subject to press the button
KbPressWait(-1);


for numRun=run:2
    
    %   'read in organized list of stimuli for this run'
    % - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    fid = fopen([outputPath '/' sprintf('%s_stimuliForcatProbe_order%d_block_1_run%d.txt',subjectID,order,numRun)]);
    stimuli = textscan(fid, '%d %d %d %d %s %s') ;% these contain everything from the organizedProbe
    stimnum1 = stimuli{1};
    stimnum2 = stimuli{2};
    leftGo = stimuli{3};
    pairType = stimuli{4};
    leftname = stimuli{5};
    rightname = stimuli{6};
    fclose(fid);
    
    %-----------------------------------------------------------------
    %% 'Write output file header'
    %-----------------------------------------------------------------
    
    fid1 = fopen([outputPath '/' subjectID sprintf('_cat_probe_run%d_', numRun) timestamp '.txt'], 'a');
    fprintf(fid1,'subjectID\torder\trun\ttrial\tonsettime\tImageLeft\tImageRight\tratingOrderLeft\tratingOrderRight\tIsleftGo\tResponse\tPairType\tOutcome\tRT\tratingLeft\tratingRight\t\fixationTime\t\n'); %write the header line    
    
    %   baseline fixation cross
    % - - - - - - - - - - - - -

    Screen('TextSize',w, 60);
    CenterText(w,'+', white,0,0);
    Screen(w,'Flip');

    
    
    %==============================================
    %% 'Run Trials'
    %==============================================
    
    runtrial = 1;
    runStart = GetSecs;
    
    if use_eyetracker
        % start recording eye position
        %---------------------------
        Eyelink('StartRecording');
        WaitSecs(.05);
        
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message',['SYNCTIME at run start:',num2str(GetSecs)]); % mark start time in file
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
        if leftGo(trial) == 1
            Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
            Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
            eyelink_message=['run: ' num2str(numRun),' trial: ',num2str(trial),', StimLeft: ',stimName{stimnum1(trial)},' StimRight: ',stimName{stimnum2(trial)},' time: ',num2str(GetSecs)];
        else
            Screen('PutImage',w,Images{stimnum2(trial)}, leftRect);
            Screen('PutImage',w,Images{stimnum1(trial)}, rightRect);
            eyelink_message=['run: ' num2str(numRun),' trial: ',num2str(trial),', StimLeft: ',stimName{stimnum2(trial)},' StimRight: ',stimName{stimnum1(trial)},' time: ',num2str(GetSecs)];
        end
        
        if use_eyetracker
            %   Eyelink MSG
            % ---------------------------
            Eyelink('Message',eyelink_message);
        end
        
        CenterText(w,'+', white,0,0);
        StimOnset = Screen(w,'Flip', runStart+onsetlist(runtrial));
        
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
                    Eyelink('Message',['run: ',num2str(numRun),' trial: ',num2str(trial),' Press_time: ',num2str(GetSecs)]);
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
            
            
            % check for reaching time limit
            if noresp && GetSecs-runStart >= onsetlist(runtrial)+maxtime
                noresp = 0;
                keyPressed = badresp;
                respTime = maxtime;
            end
        end % end while noresp
        
        %-----------------------------------------------------------------
        % determine what bid to highlight
        %-----------------------------------------------------------------
        
        switch keyPressed
            case leftstack
                colorLeft = green;
                if leftGo(trial) == 0
                    out = 0;
                else
                    out = 1;
                end
            case rightstack
                colorRight = green;
                if leftGo(trial) == 1
                    out = 0;
                else
                    out = 1;
                end
        end
        
        if goodresp==1
            if leftGo(trial)==1
                Screen('PutImage',w,Images{stimnum1(trial)}, leftRect);
                Screen('PutImage',w,Images{stimnum2(trial)}, rightRect);
            else
                Screen('PutImage',w,Images{stimnum2(trial)}, leftRect);
                Screen('PutImage',w,Images{stimnum1(trial)}, rightRect);
            end
            
            switch keyPressed
                case leftstack
                    Screen('FrameRect', w, green, leftRect, penWidth);
                case rightstack
                    Screen('FrameRect', w, green, rightRect, penWidth);
            end
            
            CenterText(w,'+', white,0,0);
            Screen(w,'Flip',runStart+onsetlist(trial)+respTime);
            
        else
            %         Screen('DrawText', w, 'You must respond faster!', xcenter-400, ycenter, white);
            CenterText(w,sprintf('You must respond faster!') ,white,0,0);
            if use_eyetracker
                %   Eyelink MSG
                % ---------------------------
                Eyelink('Message',['run: ',num2str(numRun),' trial: ',num2str(trial),' Respond_faster_time: ',num2str(GetSecs)]);
            end
            
            Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime);
        end % end if goodresp==1
        
        
        %-----------------------------------------------------------------
        % show fixation ITI
        %-----------------------------------------------------------------
        
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5);
        fixationTime(trial) = GetSecs - runStart;
        
        if use_eyetracker
            %   Eyelink MSG
            % ---------------------------
            Eyelink('Message',['run: ',num2str(numRun),' trial: ',num2str(trial),' Fixation_ITI_time: ',num2str(GetSecs)]);
        end
        
        if goodresp ~= 1
            respTime = 999;
        end
        
        %-----------------------------------------------------------------
        % 'Save data'
        %-----------------------------------------------------------------
        if leftGo(trial)==1
            fprintf(fid1,'%s\t %d\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %.2f\t %.2f\t %.2f\t %d\t \n', subjectID, order, numRun, runtrial, StimOnset-runStart, char(leftname(trial)), char(rightname(trial)), stimnum1(trial), stimnum2(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum1(trial)), bidValue(stimnum2(trial)), fixationTime(trial) );
        else
            fprintf(fid1,'%s\t %d\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %.2f\t %.2f\t %.2f\t %d\t \n', subjectID, order, numRun, runtrial, StimOnset-runStart, char(leftname(trial)), char(rightname(trial)), stimnum2(trial), stimnum1(trial), leftGo(trial), keyPressed, pairType(trial), out, respTime*1000, bidValue(stimnum2(trial)), bidValue(stimnum1(trial)), fixationTime(trial) );
        end
        
        runtrial = runtrial+1;
        %     KbQueueFlush;
        
    end % loop through trials
    fclose(fid1);

    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['run: ',num2str(numRun),' trial: ',num2str(trial),' Part_end_time: ',num2str(GetSecs)]);
    end
    
 
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');

    
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
            movefile(edfFile,['./Output/', subjectID,'_cat_probe_EyeTracking_run_',num2str(numRun),'_' timestamp,'.edf']);
        end;
    end
    
    %-----------------------------------------------------------------
    %	display outgoing message
    %-----------------------------------------------------------------
    WaitSecs(2);
    Screen('FillRect', w, black);
    Screen('TextSize',w, 40);
    
    %---------------------------------------------------------------
    % create a data structure with info about the run
    %---------------------------------------------------------------
    outfile = strcat(outputPath,'/', sprintf('%s_cat_probe_run_%2d_%s.mat',subjectID,numRun,timestamp));
    
    % create a data structure with info about the run
    run_info.subject=subjectID;
    run_info.date=date;
    run_info.outfile=outfile;
    
    run_info.script_name=mfilename;
    clear Images;
    save(outfile);
    Images = cell(1,length(stimName));
    for i = 1:length(stimName)
        Images{i} = imread([mainPath sprintf('/stim/%s',stimName{i})]);
    end
    %   'timestamp'
    % - - - - - - - - - - - - - - - - -
    c = clock;
    hr = sprintf('%02d', c(4));
    minutes = sprintf('%02d', c(5));
    timestamp = [date,'_',hr,'h',minutes,'m'];
end %end run loop

CenterText(w,sprintf('Great Job! We will continue with another run.') ,white,0,-100);
CenterText(w,sprintf('Please get the experimenter') ,white,0,100);
Screen('Flip',w);

WaitSecs(3);

KbQueueFlush;
ListenChar(0);
ShowCursor;
sca;

% end % end function

