classdef ColorDots < DeepCopyable
    % Color dots packaged for the value/perceptual decision project.
    %
    % Run either one of the demos:
    % >> info = ColorDots.demo_PTB
    % >> info = ColorDots.demo_OOP
    %
    % The latter code might be easier to blend with your existing code.
    
    % 2016 YK wrote the initial version. hk2699 at columbia dot edu.
properties
    Scr = PsyScr;
    RDK = PsyRDKConst;
    
    scr = []; % leave empty to set to the last screen.
    win = []; % This should be provided if you do not use Dots.open().
    dist_cm = 55;
    monitor_width_cm = 30;
end
%% Demo
methods (Static)
    function [info, Dots] = demo_PTB
        % This demo might blend better with the existing code.
        
        % Initialization
        scr = 0;
        background_color = 0;
        win = Screen('OpenWindow', scr, background_color);
        
        % You may want to change the following three parameters
        % with the ones measured from the experimental setup.        
        Dots = ColorDots( ...
            'scr', scr, ... 
            'win', win, ...
            'dist_cm', 55, ...
            'monitor_width_cm', 30);

        n_trial = 3;
        info = cell(1, n_trial);

        % Iterating trials
        for trial = 1:n_trial
            n_fr = 150;
            t_fr = zeros(1, n_fr + 1);
            
            color_coh_pool = ...
                [-2, -1, -0.5, -0.25, -0.125, 0, 0.125, 0.25, 0.5, 1, 2];
            color_coh = randsample(color_coh_pool, 1);
            prop = invLogit(color_coh);

            Dots.init_trial(prop);
            
            for fr = 1:n_fr
                Dots.draw;
                
                t_fr(fr) = Screen('Flip', win);
            end
            t_fr(n_fr + 1) = Screen('Flip', win); % One more flip is necessary to erase dots.
            disp(t_fr(end) - t_fr(1)); % Should be ~2.5s on a 60Hz monitor.

            info{trial} = Dots.finish_trial;
            info{trial}.t_fr = t_fr;
            
            intertrial_interval = 1;
            WaitSecs(intertrial_interval);
        end
        
        % Finishing up
        Screen('Close', win);
    end
    function [info, Dots] = demo_OOP
        % This is just another way of running the code.
        % It returns identical info as demo_PTB.
        % It might retain extra information inside the objects, 
        % but most likely they are irrelevant..
        
        % You may want to change the following three parameters
        % with the ones measured from the experimental setup.        
        Dots = ColorDots( ...
            'scr', 0, ... 
            'dist_cm', 55, ...
            'monitor_width_cm', 30);
        Dots.open; % OOP-specific
        
        n_trial = 3;
        info = cell(1, n_trial);
        
        for trial = 1:n_trial
            n_fr = 150;
            color_coh_pool = ...
                [-2, -1, -0.5, -0.25, -0.125, 0, 0.125, 0.25, 0.5, 1, 2];
            color_coh = randsample(color_coh_pool, 1);
            prop = invLogit(color_coh);

            Dots.init_trial(prop);

            for fr = 1:n_fr
                Dots.draw;
                Dots.flip; % OOP-specific
            end
            Dots.RDK.hide; % OOP-specific
            Dots.flip; % OOP-specific

            info{trial} = Dots.finish_trial;
            disp(info{trial}.t_fr(end) - info{trial}.t_fr(1)); % Should be ~2.5s on a 60Hz monitor.
            
            intertrial_interval = 1;
            WaitSecs(intertrial_interval);
        end
        Dots.close; % OOP-specific
    end
end
%% Script interface
methods
    function Dots = ColorDots(varargin)
        varargin2props(Dots, varargin);
        
        Dots.add_deep_copy({'Scr', 'RDK'});
        Dots.Scr = PsyScr;
        Dots.Scr.init( ...
            'scr', Dots.scr, ...
            'win', Dots.win, ...
            'distCm', Dots.dist_cm, ...
            'widthCm', Dots.monitor_width_cm);
        Dots.RDK = PsyRDKConst(Dots.Scr);
        Dots.RDK.win = Dots.win;
        
        Dots.Scr.addObj('Vis', Dots.RDK);
        
        Dots.scr = Dots.Scr.info.scr;
    end
    function init_trial(Dots, prop)
        % Ignore warning from DrawDots.
        if Dots.is_Scr_open
            Dots.Scr.initLogTrial;
        else
            Dots.RDK.initLogTrial;
        end
        Dots.RDK.init(0, prop, {'shuffle', 'shuffle', 'shuffle'});
        Dots.RDK.show;
    end
    function draw(Dots)
        Dots.RDK.update('befDraw');
        Dots.RDK.draw;
    end
    function info = finish_trial(Dots)
        Dots.Scr.closeLog;
        info = Dots.get_info;
    end
end
%% Internal
methods
    function open(Dots)
        Dots.Scr.open;
    end
    function tf = is_Scr_open(Dots)
        tf = Dots.Scr.opened;
    end
    function update_only(Dots)
        Dots.RDK.update('befDraw');
    end
    function flip(Dots)
        Dots.Scr.flip;
    end
    function info = get_info(Dots)
        info = struct;
        info.prop = Dots.RDK.prop;
        info.color_coh = logit(info.prop);
        info.xy_pix = Dots.RDK.vTrim('xyPix');
        info.col2 = Dots.RDK.vTrim('col2');
        
        if Dots.is_Scr_open
            info.t_fr = Dots.RDK.absSec('xyPix');
            info.t_fr(end + 1) = Dots.RDK.absSec('off');
        else
            info.t_fr = [];
        end
    end
    function close(Dots)
        Dots.Scr.close;
    end
end
end