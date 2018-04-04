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
rng('default');

% about timing
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];


baseRect = [0 0 400 20];
pointerRect = [0 0 5 30];

trial_time_fixated_food=999;
trial_time_fixated_scale=999;
trial_num_food_fixations=999;
trial_num_scale_fixations=999;
first_fixation_area='x';
first_fixation_duration=999;

% -----------------------------------------------
%% 'INITIALIZE SCREEN'
% -----------------------------------------------

Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
%Screen('Preference', 'SkipSyncTests', 1); %ONLY FOR TESTING
screennum = min(Screen('Screens')); %select external screen

pixelSize=32;
[w, windowRect] = Screen('OpenWindow',screennum,[],[],pixelSize);
[xCenter, yCenter] = RectCenter(windowRect);

buffer=20; %20 pixel buffer zone around rects of interest
foodRect=CenterRectOnPointd([0 0 400+buffer 400+buffer], xCenter, yCenter);
scaleRect=CenterRectOnPointd([0 0 400+buffer 20+buffer], xCenter, yCenter+350);

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
    Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERS,STATUS,INPUT');
    Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERS,AREA,HTARGET,STATUS,INPUT');
    
    % open file to record data to
    edfFile=['rating' num2str(run) '.edf'];
    i=Eyelink('Openfile', edfFile);
        if i~=0
        fprint('Cannot create EDF file ''%s'' ', edfFile);
        Eyelink('Shutdown');
        sca;
        return;
    end
    
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

if run==1
    CenterText(w,'You will see a series of pictures of food.', white, 0,-300);
    CenterText(w,'Imagine you had to eat one of these foods today.', white, 0,-250);
    CenterText(w,'For each picture, please rate how much you would prefer to eat that food.', white, 0,-200);
    CenterText(w,'You will rate each picture on a scale from 0 to 10,', white, 0,-150);
    CenterText(w,'with 0 being that you would not want to eat that food at all', white, 0,-100);
    CenterText(w,'and 10 being that you most strongly prefer to eat that food.', white, 0,-50);
    CenterText(w,'Use the mouse to move the blue rating indicator bar', white, 0,0);
    CenterText(w,'along the scale to indicate your preference.', white, 0,50);
    CenterText(w,'There are no right answers. Please rate only according to your own preference.', white, 0,100);
    CenterText(w,'Take as much time as you would like.', white, 0,150);
    CenterText(w,'Click anywhere to continue', white, 0,250);
else
    CenterText(w,'Now that you have seen all the possible food items,', white, 0,-25);
    CenterText(w,'please rate again how much you would prefer to eat the food in each picture.', white, 0,25);
    CenterText(w,'Click anywhere to continue', white, 0,125);
end
    
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
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.05);
    Eyelink('StartRecording');
    WaitSecs(.05);

    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message', ['SYNCTIME start: ', num2str(runStartTime)]); % mark start time in file
    if ~dummymode
        eye_used = Eyelink('EyeAvailable');
        if eye_used == -1
            fprintf('Eyelink aborted - could not find which eye being used.\n');
            cleanup;
        end
    end
end

WaitSecs(2); % Wait 2 sec before first stimulus


%   'Write output file header'
%---------------------------------------------------------------

fid1 = fopen([outputPath '/' subjectID '_food_rating_run' num2str(run) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\t onsetTime\t itemName\t rating\t RT\t time_onfood\t time_onscale\t num_food_fix\t num_scale_fix\t first_fix_area\t first_fix_duration\n'); %write the header line

%	pre-allocating matrices and setting defaults
%---------------------------
respTime(1:length(shuff_names),1) = 999;
rating(1:length(shuff_names),1) = 999;
image_start_time(1:length(shuff_names),1) = 999;
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter+350); %main scale

%Run trial loop
for trialNum = 1:length(shuff_names)   % To cover all the items in one run.
    xpos=Shuffle(-200:200); %want the pointer to start in a random position every trial
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
        %Eyelink('Message', ['trial: ' num2str(trialNum) ' stim: ',shuff_names(trialNum).name,' start_time: ',num2str(image_start_time(trialNum))]); % mark start time in file
        Eyelink('Message', ['trial: ' num2str(trialNum)]); % mark start time in file
        trial_time_fixated_food = 0;
        trial_time_fixated_scale = 0;
        trial_num_food_fixations = 0;
        trial_num_scale_fixations = 0;
        
        % current_area determines which area eye is in (left, right, neither)
        % xpos and ypos are used for eyepos_debug
        [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,foodRect,scaleRect);
        
        % last_area will track what area the eye was in on the previous loop
        % iteration, so we can determine when a change occurs
        % fixation_onset_time stores the time a "fixation" into an area began
        first_fixation_duration = 0;
        first_fixation_area = current_area; % this will report 'n' in output if they never looked at an object
        first_fixation_flag = (first_fixation_area=='f' || first_fixation_area=='s'); % flags 1 once the first fixation has occurred, 2 once the first fixation has been processed
        last_area=current_area;
        fixation_onset_time = GetSecs;
        first_fixation_onset=fixation_onset_time;
    end
    
    
    %---------------------------------------------------
    %% Move pointer with cursor
    %---------------------------------------------------
    noclick=1;
    while noclick
        
        if use_eyetracker==1
            % get eye position
            [current_area, ~, ~] = get_current_fixation_area(dummymode,el,eye_used,foodRect,scaleRect);
            
            % they are looking in a new area
            % Currently has initial fixation problems? (color, count, etc.)
            if current_area~=last_area
                % update timings
                switch last_area
                    case 'f'
                        trial_time_fixated_food = trial_time_fixated_food + (GetSecs-fixation_onset_time);
                        trial_num_food_fixations = trial_num_food_fixations + 1;
                    case 's'
                        trial_time_fixated_scale = trial_time_fixated_scale + (GetSecs-fixation_onset_time);
                        trial_num_scale_fixations = trial_num_scale_fixations + 1;
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
                if(first_fixation_flag==0 && (current_area=='f' || current_area=='s'))
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
            while any(buttons)
                [mx, my, buttons] = GetMouse(w);
                centeredPointer=CenterRectOnPointd(pointerRect, mx, yCenter+350);
                Screen('PutImage',w,Images{trialNum});
                Screen('FillRect', w, white, centeredRect);
                Screen('FillRect', w, blue, centeredPointer);
                CenterText(w,'0',white,-200,360);
                CenterText(w,'5',white,0,360);
                CenterText(w,'10',white,200,360);
                vbl=Screen(w,'Flip');
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
        switch last_area
            case 'f'
                trial_time_fixated_food = trial_time_fixated_food + fixation_duration;
                trial_num_food_fixations = trial_num_food_fixations + 1;
            case 's'
                trial_time_fixated_scale = trial_time_fixated_scale + fixation_duration;
                trial_num_scale_fixations = trial_num_scale_fixations + 1;
        end
        % time limit reached while fixating on first fixated object
        if(first_fixation_flag==1)
            %outstr=['first fixation lasted ' GetSecs-first_fixation_onset ' seconds'];
            %Eyelink('Message',outstr);
            first_fixation_duration = GetSecs-first_fixation_onset;
            first_fixation_flag = 2;
        end
    end
 
    %   'Save data'
    %---------------------------
    fprintf(fid1,'%s\t %d\t %s\t %.2f\t %.4f\t %.4f\t %.4f\t %d\t %d\t %c\t %.4f\n',...
        subjectID, image_start_time(trialNum)-runStartTime, shuff_names(trialNum).name, rating(trialNum), respTime(trialNum), ...
        trial_time_fixated_food, trial_time_fixated_scale, trial_num_food_fixations, trial_num_scale_fixations,...
        first_fixation_area, first_fixation_duration);
    
    WaitSecs(1); %1 sec ITI
    
end %	End the big trialNum loop showing all the images in one run.


%---------------------------------------------------------------
%%   save data to a .mat file & close out
%---------------------------------------------------------------
outfile = strcat(outputPath, '/', subjectID,'_food_rating_run_', num2str(run), '_', timestamp,'.mat');
% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
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
        WaitSecs(1);
        movefile(edfFile,['./Output/', subjectID,'_food_rating_run_',num2str(run),'_',timestamp,'.edf']);
    end
end

%   outgoing msg & closing
% ------------------------------
CenterText(w,'Great Job. We will continue shortly.',white, 0,-100);
CenterText(w,'Please get the experimenter when you are ready.', white, 0, 100);
Screen('Flip',w);

WaitSecs(4);
ListenChar(0);
Screen('CloseAll');

end

function [current_area,  xpos, ypos] = get_current_fixation_area(dummymode,el,eye_used,foodRect,scaleRect)
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
if IsInRect(xpos,ypos,foodRect)
    current_area='f';
elseif IsInRect(xpos,ypos,scaleRect)
    current_area='s';
else
    current_area='n';
end
return
end
