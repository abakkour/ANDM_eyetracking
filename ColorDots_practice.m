function ColorDots_practice(subjid,test_comp,exp_init,eye,scan,button_order)
% This demo shows color dots three times and returns their information.
%
% ColorDots_practice(subjid,test_comp,exp_init,eye,scan,button_order)
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
%Screen('Preference', 'SkipSyncTests', 1); % FOR TESTING PURPOSES ONLY!

c=clock;
hr=num2str(c(4));
min=num2str(c(5));
timestamp=[date,'_',hr,'h',min,'m'];

ColorDots_init_path;

% Initialization
scr = 0;
background_color = 0;
[win, windowRect] = Screen('OpenWindow', scr, background_color);
[xCenter, yCenter] = RectCenter(windowRect);

buffer=20; %20 pixel buffer zone around rects of interest
dotsRect=CenterRectOnPointd([0 0 400+buffer 400+buffer], xCenter, yCenter);

green=[0 255 0];
red=[255 0 0];
white=[255 255 255];
black=[0 0 0];

KbQueueCreate;
KbQueueStart;
HideCursor;


trial_time_fixated_dots=999;
trial_num_dots_fixations=999;

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
        cleanup;  % cleanup function
        return;
    end;
    
    [~, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
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

ListenChar(2);

switch scan
    case 1
        switch button_order
            case 1
                blue='1!';
                bluekey='1';
                yellow='2@';
                yellowkey='2';
            case 2
                blue='2@';
                bleukey='2';
                yellow='1!';
                yellowkey='1';
        end
    case 0
        switch button_order
            case 1
                blue='u';
                bluekey='u';
                yellow='i';
                yellowkey='i';
            case 2
                blue='i';
                bluekey='i';
                yellow='u';
                yellowkey='u';
        end
end

% You may want to change the following three parameters
% with the ones measured from the experimental setup.
Dots = ColorDots( ...
    'scr', scr, ...
    'win', win, ...
    'dist_cm', 55, ...
    'monitor_width_cm', 30);

n_trial = 200;
info = cell(1, n_trial);
outcomes=zeros(1,n_trial);

% I recommend the pool of color coherences in the code.
% You might omit one of the zeros (i.e., leave only one zero)
% - that would slightly reduce the power for reverse correlation
% later.
% You might also omit the -2 and 2, but that might reduce the
% range of RTs.
%color_coh_pool = [-2, -1, -0.5, -0.25, -0.125, 0, 0, 0.125, 0.25, 0.5, 1, 2];
coh_pool{1} = [-2, 2];
coh_pool{2} = [-2, -1, 1, 2];
coh_pool{3} = [-2, -1, -0.5, 0.5, 1, 2];
coh_pool{4} = [-2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2];
coh_pool{5} = [-2, -1, -0.5, -0.25, -0.125, 0, 0, 0.125, 0.25, 0.5, 1, 2];

n_fr = 150;


%ListenChar(2);
KbQueueFlush;

%INTRUCTIONS
Screen('TextSize',win,40);
CenterText(win,'You will see a cloud of flickering dots that are either yellow or blue.',white,0,-300);
CenterText(win,['If you think the cloud contains more yellow dots than blue on average, press key `' yellowkey '`.'],white,0,-250);
CenterText(win,['If you think the cloud contains more blue dots, press key `' bluekey '`.'],white,0,-200);
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
    Screen('Flip',w);
    escapeKey = KbName('t');
    while 1
        [keyIsDown, firstPress] = KbQueueCheck;
        if keyIsDown && firstPress(escapeKey)
            break;
        end
    end
end

if eye==1
    % STEP 5
    % start recording eye position
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    Eyelink('Message', 'SYNCTIME after fixation run1'); % mark start time in file
    if ~dummymode
        eye_used = Eyelink('EyeAvailable');
        if eye_used == -1
            fprintf('Eyelink aborted - could not find which eye being used.\n');
            cleanup;
        end
    end
end


CenterText(win,'+',white,0,0);
runStart=Screen('Flip',win);
WaitSecs(1);

fid1=fopen(['Output/' subjid '_dots_practice_' timestamp '.txt'], 'a');

%write the header line
fprintf(fid1,'subjid\t scanner\t test_comp\t experimenter\t runtrial\t onsettime\t color_coh\t prop\t response\t outcome\t disptime\t RT\t button_order\t time_fix_dots\t num_dots_fix\n');

% Iterating trials
trial=0;
for c=1:5
    color_coh_pool=coh_pool{c};
    a=1;
    while a
        
        trial=trial+1;
        
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
                
                trial_time_fixated_dots = 0;
                trial_num_dots_fixations = 0;
                
                % current_area determines which area eye is in (left, right, neither)
                % xpos and ypos are used for eyepos_debug
                [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,dotsRect);
                
                % last_area will track what area the eye was in on the previous loop
                % iteration, so we can determine when a change occurs
                % fixation_onset_time stores the time a "fixation" into an area began
                first_fixation_duration = 0;
                first_fixation_area = current_area; % this will report 'n' in output if they never looked at an object
                first_fixation_flag = (first_fixation_area=='f'); % flags 1 once the first fixation has occurred, 2 once the first fixation has been processed
                last_area=current_area;
                fixation_onset_time = GetSecs;
            elseif eye==1 && fr > 1
                % get eye position
                [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,dotsRect);
                
                % they are looking in a new area
                % Currently has initial fixation problems? (color, count, etc.)
                if current_area~=last_area
                    % update timings
                    switch last_area
                        case 'f'
                            trial_time_fixated_dots = trial_time_fixated_dots + (GetSecs-fixation_onset_time);
                            trial_num_dots_fixations = trial_num_dots_fixations + 1;
                    end
                    
                    fixation_onset_time=GetSecs;
                    
                    % they have looked away from their first fixation: record its
                    % duration and the target (food/scale)
                    if(first_fixation_flag==1)
                        %outstr=['first fixation lasted ' GetSecs-first_fixation_onset ' seconds'];
                        %Eyelink('Message',outstr);
                        first_fixation_duration = GetSecs-first_fixation_onset;
                        first_fixation_flag = 2;
                    end
                    
                    % this is their first time fixating on an object this trial
                    if(first_fixation_flag==0 && current_area=='f')
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
            
            [keyIsDown, firstPress] = KbQueueCheck;
            keyPressed=KbName(firstPress);
            if keyIsDown && ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                keyPressed=char(keyPressed);
                keyPressed=keyPressed(1);
            end
            if keyIsDown && (keyPressed == blue || keyPressed == yellow)
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
            
            switch last_area
                case 'd'
                    trial_time_fixated_dots = trial_time_fixated_dots + fixation_duration;
                    trial_num_dots_fixations = trial_num_dots_fixations + 1;
            end
            % time limit reached while fixating on first fixated object
            if(first_fixation_flag==1)
                %outstr=['first fixation lasted ' GetSecs-first_fixation_onset ' seconds'];
                %Eyelink('Message',outstr);
                first_fixation_duration = GetSecs-first_fixation_onset;
                first_fixation_flag = 2;
            end
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
        
        if color_coh > 0 && keyPressed == blue
            outcome = 1;
        elseif color_coh < 0 && keyPressed == yellow
            outcome = 1;
        elseif color_coh == 0 && ~keyPressed=='x'
            outcome = NaN;
        else
            outcome = 0;
        end
        
        outcomes(trial)=outcome;
        
        if keyPressed=='x'
            CenterText(win,'TOO SLOW!',white,0,0);
        else
            if color_coh == 0
                outcome=randsample([1 0],1);
            end
            if outcome==1
                CenterText(win,'CORRECT!',green,0,0);
            else
                CenterText(win,'INCORRECT',red,0,0);
            end
        end
        
        Screen('Flip', win);
        
        info{trial}.outcome=outcome;
        
        fprintf(fid1,'%s\t %d\t %s\t %s\t %d\t %f\t %f\t %f\t %s\t %d\t %f\t %f\t %d\t %.4f\t %d\n', ...
            subjid, scan, test_comp, exp_init, trial, t_fr(1)-runStart, ...
            color_coh, prop, keyPressed, outcome, info{trial}.disptime, info{trial}.rt, button_order, ...
            trial_time_fixated_dots, trial_num_dots_fixations);
        
        WaitSecs(.5);
        
        CenterText(win,'+',white,0,0);
        fixtime=Screen('Flip', win);
        
        if eye==1
            % Eyelink msg
            % - - - - - - -
            fixcrosstime = strcat('fixcrosstime = ',num2str(fixtime-runStart));
            Eyelink('Message',fixcrosstime);
        end
        
        intertrial_interval = 1;
        WaitSecs(intertrial_interval);
        
        disp(info{trial});
        fprintf('\n\n');
        
        if c==1 && trial >= 10 && nanmean(outcomes(trial-9:trial))==1
            a=0;
        elseif c==2 && trial >= 20 && nanmean(outcomes(trial-9:trial))==1
            a=0;
        elseif c==3 && trial >= 30 && nanmean(outcomes(trial-9:trial))>=.9
            a=0;
        elseif c==4 && trial >= 40 && nanmean(outcomes(trial-9:trial))>=.8
            a=0;
        elseif trial== 200
            a=0;
        end
        
        if trial==100
            CenterText(win,'Great job! Take a break.',white,0,-100);
            CenterText(win,'Take as long as you need.',white,0,0);
            CenterText(win,'When you are ready, please get the experimenter...',white,0,100);
            Screen('Flip', win);
            switch trial
                case 100
                    run=1;
                case 200
                    run=2;
                case 300
                    run=3;
            end
            if eye==1
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
                        movefile('recdata.edf',strcat('Output/', subjid,'_dots_practice_run',num2str(run),'_',timestamp,'.edf'));
                    end;
                end
            end
            KbQueueFlush;
            KbQueueWait;
            if eye==1
                ListenChar(0);
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
                    cleanup;  % cleanup function
                    return;
                end;
                
                [~, vs]=Eyelink('GetTrackerVersion');
                fprintf('Running experiment on a ''%s'' tracker.\n', vs );
                
                % make sure that we get gaze data from the Eyelink
                Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
                Eyelink('command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
                Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERS,STATUS,INPUT');
                Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERS,AREA,HTARGET,STATUS,INPUT');
                
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
                
                
                ListenChar(2);
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % Finish Initialization %
                %%%%%%%%%%%%%%%%%%%%%%%%%
            end
            CenterText(win,'Please press any key to continue...',white,0,0);
            Screen('Flip', win);
            KbQueueFlush;
            KbQueueWait;
            CenterText(win,'+',white,0,0);
            Screen('Flip', win);
            if eye==1
                % STEP 5
                % start recording eye position
                Eyelink('StartRecording');
                % record a few samples before we actually start displaying
                WaitSecs(0.1);
                Eyelink('Message', 'SYNCTIME after fixation run2'); % mark start time in file
                if ~dummymode
                    eye_used = Eyelink('EyeAvailable');
                    if eye_used == -1
                        fprintf('Eyelink aborted - could not find which eye being used.\n');
                        cleanup;
                    end
                end
            end
            WaitSecs(1);
        end
        
    end
end

save(['Output/' subjid '_dots_practice_' timestamp '.mat'],'Dots','info')

fclose(fid1);
if trial == 200
    run=2;
else
    run=4;
end
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
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
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
        movefile('recdata.edf',strcat('Output/', subjid,'_dots_practice_run',num2str(run),'_',timestamp,'.edf'));
    end
end

CenterText(win,'Thank you! Great job!', white,0,-100);
CenterText(win,'Please get the experimenter.', white,0,0);
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
WaitSecs(0.1);
Eyelink('Command', 'set_idle_mode');
WaitSecs(0.5);

Eyelink('CloseFile');
Eyelink('Shutdown');

% Close window:
Screen('CloseAll');

% Restore keyboard output to Matlab:
ListenChar(0);
ShowCursor;
end

function [current_area,  xpos, ypos] = get_current_fixation_area(dummymode,el,eye_used,dotsRect)
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
if IsInRect(xpos,ypos,dotsRect)
    current_area='d';
else
    current_area='n';
end
return
end