function ColorDots_test(subjid,test_comp,exp_init,eye,scan,run,button_order)

% function ColorDots_test(subjid,test_comp,exp_init,eye,scan,run,task_order, button_order)
% This demo shows color dots three times and returns their information.
%
% info{trial} has the following fields:
% - color_coh
% : corresponds to the difficulty and the answer.
%   Negative means yellow, positive blue.
%   The larger absolute value, the easier it is.
%
% - prop
% : the probability of a dot being blue.
%   Note that color_coh == logit(prop).
%
% - xy_pix{fr}(xy, dot) has the dot position on that frame in pixel.
%     The first row is x, the second y.
%
% - col2{fr}(dot) = 1 means that the dot on that frame was blue.
%
% Dots contains more information, but perhaps they wouldn't matter
% in most cases.

% 2016 YK wrote the initial version. hk2699 at columbia dot edu.
% Feb 2016 modified by AB. ab4096 at columbia dot edu.

Screen('Preference', 'VisualDebugLevel', 0);
Screen('Preference', 'SkipSyncTests', 1); %FOR TESTING ONLY

c=clock;
hr=num2str(c(4));
min=num2str(c(5));
timestamp=[date,'_',hr,'h',min,'m'];

ColorDots_init_path;

% Initialization
scr = 0;
background_color = 0;
win = Screen('OpenWindow', scr, background_color);

green=[0 255 0];
white=[255 255 255];
black=[0 0 0];

KbQueueCreate;
HideCursor;
ListenChar(2);

if eye==1
    %==============================================
    %% 'INITIALIZE Eyetracker'
    %==============================================
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initializing eye tracking system %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ListenChar(2);
    dummymode=0;
    eyepos_debug=0;
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(win);
    % Disable key output to Matlab window:
    %%%%%%%%%%%%%ListenChar(2);
    
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
t        cleanup;  % cleanup function
        return;
    end;
    
    [~, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile='recdata.edf';
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    
    
    Screen('TextSize',win, 40);
    CenterText(win,'Continue using eyetracker?', white,0,0);
    CenterText(win,'(y)es', white,-75,100);
    CenterText(win,'/', white,0,100);
    CenterText(win,'(n)o', white,75,100);
    Screen(win,'Flip');
    
    noresp=1;
    while noresp
        [keyIsDown, firstPress] = KbQueueCheck;
        if keyIsDown && noresp
            keyPressed=KbName(firstPress);
            if ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed=char(keyPressed);
                keyPressed=keyPressed(1);
            end
            switch keyPressed
                case 'y'
                    noresp=0;
                    eye=1;
                    ycol=green;
                    ncol=white;
                    CenterText(win,'Continue using eyetracker?', white,0,0);
                    CenterText(win,'(y)es', ycol,-75,100);
                    CenterText(win,'/', white,0,100);
                    CenterText(win,'(n)o', ncol,75,100);
                    Screen(win,'Flip');
                    WaitSecs(.5);
                    % do a final check of calibration using driftcorrection
                    EyelinkDoDriftCorrection(el);
                case 'n'
                    noresp=0;
                    eye=0;
                    ycol=white;
                    ncol=green;
                    CenterText(win,'Continue using eyetracker?', white,0,0);
                    CenterText(win,'(y)es', ycol,-75,100);
                    CenterText(win,'/', white,0,100);
                    CenterText(win,'(n)o', ncol,75,100);
                    Screen(win,'Flip');
                    WaitSecs(.5);
            end
        end
    end
    
    
    %ListenChar(0);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish Initialization %
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

switch scan
    case 1
        switch button_order
            case 1
                blue='3#';
                bluebutton='1';
                yellow='4$';
                yellowbutton='2';
            case 2
                blue='4$';
                bluebutton='2';
                yellow='3#';
                yellowbutton='1';
        end
    case 0
        switch button_order
            case 1
                blue='u';
                bluebutton='u';
                yellow='i';
                yellowbutton='i';
            case 2
                blue='i';
                bluebutton='i';
                yellow='u';
                yellowbutton='u';
        end
end

% You may want to change the following three parameters
% with the ones measured from the experimental setup.
Dots = ColorDots( ...
    'scr', scr, ...
    'win', win, ...
    'dist_cm', 55, ...
    'monitor_width_cm', 30);

n_trial = 70;
info = cell(1, n_trial);
outcomes=zeros(1,n_trial);
%load('onsets/dots_onset.mat');
iti=ones(1,70); %iti fixed 1 sec

% I recommend the pool of color coherences in the code.
% You might omit one of the zeros (i.e., leave only one zero)
% - that would slightly reduce the power for reverse correlation
% later.
% You might also omit the -2 and 2, but that might reduce the
% range of RTs.
color_coh_pool = [-2, -1, -0.5, -0.25, -0.125, 0, 0, 0.125, 0.25, 0.5, 1, 2];
n_fr = 150;


%ListenChar(2);

%INTRUCTIONS
Screen('TextSize',win,40);
CenterText(win,'You will see a cloud of flickering dots that are either yellow or blue.',white,0,-300);
CenterText(win,['If you think the cloud contains more yellow dots than blue on average, press key `' yellowbutton '`.'],white,0,-250);
CenterText(win,['If you think the cloud contains more blue dots, press key `' bluebutton '`.'],white,0,-200);
CenterText(win,'Do not try to count the exact number of dots in each color,',white,0,-150);
CenterText(win,'because the number fluctuates rapidly over time,',white,0,-100);
CenterText(win,'and because each dot appears only briefly.',white,0,-50);
CenterText(win,'Rather, try to estimate the rough average.',white,0,0);
CenterText(win,'Please respond as soon as you have an answer.',white,0,50);
CenterText(win,'Press any button to continue...',white,0,200);
Screen('Flip',win);
KbQueueWait;
KbQueueFlush;

if scan==1
    CenterText(win,'GET READY!', white, 0, 0);    %this is for the MRI scanner, it waits for a 't' trigger signal from the scanner
    Screen('Flip',win);
    KbTriggerWait(KbName('t'),triggerkbid);
end



if eye==1
    % STEP 5
    % start recording eye position
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    Eyelink('Message', 'SYNCTIME after fixations'); % mark start time in file
    if ~dummymode
        eye_used = Eyelink('EyeAvailable');
        if eye_used == -1
            fprintf('Eyelink aborted - could not find which eye being used.\n');
            cleanup;
        end
    end
end

KbQueueCreate;
KbQueueStart;

CenterText(win,'+',white,0,0);
runStart=Screen('Flip',win);
WaitSecs(2);

fid1=fopen(['Output/' subjid '_dots_test_run_' num2str(run) '_' timestamp '.txt'], 'a');

%write the header line
fprintf(fid1,'subjid\t scanner\t test_comp\t experimenter\t runtrial\t onsettime\t color_coh\t prop\t response\t outcome\t disptime\t RT\t button_order\n');

% Iterating trials
for trial=1:n_trial;
    tstime=0;
    fprintf('-----\n');
    fprintf('Trial %d:\n', trial);
    
    t_fr = zeros(1, n_fr + 1);
    
    color_coh = randsample(color_coh_pool, 1);
    prop = invLogit(color_coh);
    
    % Dots.init_trial must be called before Dots.draw.
    Dots.init_trial(prop);
    KbQueueFlush;
    
    for fr = 1:n_fr
        % Since the dots should update every frame,
        % draw other components (e.g., fixation point)
        % before each flip, around Dots.draw.
        
        % Draw components here to have dots draw over them.
        
        Dots.draw;
        
        % Draw components here to draw over the dots.
        
        t_fr(fr) = Screen('Flip', win);
        
        if eye==1 && fr==1
            % Eyelink msg
            % - - - - - - -
            onsetmessage=strcat('Trial ',num2str(trial),' Onset = ',num2str(t_fr(fr)-runStart));
            Eyelink('Message',onsetmessage);
        end
        
        [keyIsDown, firstPress] = KbQueueCheck;
        keyPressed=KbName(firstPress);
        if keyIsDown && ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
            keyPressed=char(keyPressed);
            keyPressed=keyPressed(1);
        end
        if keyIsDown && (strcmp(keyPressed,blue) || strcmp(keyPressed,yellow))
            break;
        end
    end
    
    % One more flip is necessary to erase the dots,
    % e.g., after the button press.
    t_fr(n_fr + 1) = Screen('Flip', win);
    disp(t_fr(end) - t_fr(1)); % Should be ~2.5s on a 60Hz monitor.
    
    if eye==1
        % Eyelink msg
        % - - - - - - -
        rtmsg = strcat('RT = ',num2str(t_fr(end) - t_fr(1)));
        Eyelink('Message',rtmsg);
    end
    
    info{trial} = Dots.finish_trial;
    info{trial}.t_fr = t_fr;
    
    if isempty(keyPressed)
        keyPressed='x';
    end
    
    info{trial}.keypressed = keyPressed;
    if keyPressed ~= 'x'
        info{trial}.rt=firstPress(KbName(keyPressed))-t_fr(1);
    else
        info{trial}.rt=NaN;
    end
    info{trial}.disptime=t_fr(end) - t_fr(1);
    
    if color_coh > 0 && strcmp(keyPressed,blue)
        outcome = 1;
    elseif color_coh < 0 && strcmp(keyPressed,yellow)
        outcome = 1;
    elseif color_coh == 0
        outcome = NaN;
    else
        outcome = 0;
    end
    
    info{trial}.outcome=outcome;
    outcomes(trial)=outcome;
    
    fprintf(fid1,'%s\t %d\t %s\t %s\t %d\t %f\t %f\t %f\t %s\t %d\t %f\t %f\t %d\n', subjid, scan, test_comp, exp_init, trial, t_fr(1)-runStart, color_coh, prop, keyPressed, outcome, info{trial}.disptime, info{trial}.rt, button_order);
    
    if keyPressed=='x'
        CenterText(win,'TOO SLOW!',white,0,0);
        tstime=.5;
        Screen('Flip', win);
        WaitSecs(.5);
    end
    
    CenterText(win,'+',white,0,0);
    fixtime=Screen('Flip', win);
    
    if eye==1
        % Eyelink msg
        % - - - - - - -
        fixcrosstime = strcat('fixcrosstime = ',num2str(fixtime-runStart));
        Eyelink('Message',fixcrosstime);
    end
    
    disp(info{trial});
    fprintf('\n\n');
    
    intertrial_interval = 2.5 - info{trial}.disptime - tstime + iti(trial);
    WaitSecs(intertrial_interval);
    
end

save(['Output/' subjid '_dots_test_run_' num2str(run) '_' timestamp '.mat'],'Dots','info')

fclose(fid1);

%==============================================
%% 'BLOCK over, close out and save data'
%==============================================

%---------------------------------------------------------------
%   close out eyetracker
%---------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%
% finishing eye tracking %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 7
% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
if eye==1
    Eyelink('StopRecording');
    WaitSecs(.1);
    Eyelink('CloseFile');
    
    % download data file
    % - - - - - - - - - - - -
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
        movefile('recdata.edf',strcat('Output/', subjid,'_dots_test_run_', num2str(run), '_', timestamp,'.edf'));
    end;
end

CenterText(win,'Great job!', white,0,-100);
if scan==1
    CenterText(win,'We will continue to the next run shortly.', white,0,0);
elseif run<3
    CenterText(win,'We will continue to the next run shortly.', white,0,0);
    CenterText(win,'Please get the experimenter.', white,0,100);
else
    CenterText(win,'Great job! Please get the experimenter.', white,0,0);
end
Screen('Flip', win);
WaitSecs(5);

% Finishing up
ListenChar(0);
ShowCursor;
Screen('Close', win);
end

% Cleanup routine:
function cleanup

% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
Eyelink('Stoprecording');
Eyelink('CloseFile');
Eyelink('Shutdown');

% Close window:
Screen('CloseAll');

% Restore keyboard output to Matlab:
%ListenChar(0);
ShowCursor;
end