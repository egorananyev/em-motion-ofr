function em

%% Preparing the variables.

% Current test variables:
subj = 1001;
domEye = 0; % 0=right, 1=left
condFileName = 'cond01';

% Keyboard:
KbName('UnifyKeyNames');
quitkey = 'c';
space = 'space';
% targetUsageName = 'Keyboard'; % change accordingly
% targetProduct = 'Dell USB Keyboard'; % change accordingly
% targetProduct = 'Apple Keyboard'; % temp
% dev = PsychHID('Devices');
% deviceIndex = find(strcmpi(targetUsageName, {dev.usageName}) & ...
%     strcmpi(targetProduct, {dev.product}));
deviceIndex = -3;
KbQueueCreate(deviceIndex);
KbQueueStart(deviceIndex);

% Output file name:
dateNtime = datestr(now,'yyyy-mm-dd_HHMMSS');
sessionName = strcat('em_', condFileName', '_s', mat2str(subj), ...
    '_d', mat2str(domEye)', '_', dateNtime);
outFileName = strcat('bcfs/', sessionName, '.csv');

% Instructions text:
textInstr = 'Press any button to continue';
textNextTrial = 'Press spacebar to continue';

%% Preparing PsychToolBox and screen.

% Prepping PsychToolBox:
% Screen('Preference', 'SkipSyncTests', 1); % a necessary evil

AssertOpenGL; % for 3D rendering 
screenid = max(Screen('Screens'));
% InitializeMatlabOpenGL(1); % for 3D rendering
% PsychImaging('PrepareConfiguration'); 

% Some verification for colour names (for the gratings):
white = WhiteIndex(screenid);
black = BlackIndex(screenid);
gray = round((white+black)/2);
if gray == white
    gray = white / 2;
end
inc = white - gray; % increment
backgroundCol = black;

try
[wPtr, rect] = Screen('OpenWindow', screenid, backgroundCol);
% wPtr = 10 % the number that designates the created window;
% rect = [0 0 1920 1080] % RectLeft=1, RectTop=2, RectRight=3, RectBottom=4
% Screen('BlendFunction', wPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Screen('CloseAll'); %temp

%% Resolution and display locations.

% screen resolution, in pixels:
disp.resX = rect(3); % 1920
disp.resY = rect(4); % 1080

% display center locations, in pixels:
disp.boxColour = [255 255 255]; % white
disp.boxSize = 90;
disp.distX = 170; % 150; % display center distance from the vertical midline (left/right)
disp.distY = -100; % display center distance from the horizontal midline (up)
disp.centX(1:2) = [disp.resX/2-disp.distX disp.resX/2+disp.distX]; % 230(R) & 730(L)
disp.centY = disp.resY/2 + disp.distY; % 540-200=340

%Screen('CloseAll'); %temp

%% Experimental conditions.

% Reading the conditions file with the settings
[~,~,condTable] = xlsread(strcat(condFileName, '.xlsx'));
% condVars = condCells(1,:);
numofConds = size(condTable,1)-1; % number of conditions
% condTable = condCells(2:numofConds+1,:);

% Also reading screen dimensions for px2cm and cm2px conversions:
[~,~,sdims] = xlsread('screenDims.xlsx');
% sdimsVars = sdims(1,:);
% sdimsTable = sdims(2,:);

% Basic settings:
numofBlocks = 2;
numofTrials = numofBlocks * numofConds;

% randomizing presentation of conditions within each block of trials
trialCond = [];
for i=1:numofBlocks
    trialCond = [trialCond randperm(numofConds)]; %#ok<AGROW>
end

%% Presenting the instructins window.
Screen('TextFont', wPtr, 'Cambria');
Screen('TextSize', wPtr, 32);
% Defining the edges of the left and right box:
boxL = [disp.centX(1)-disp.boxSize disp.centY-disp.boxSize...
    disp.centX(1)+disp.boxSize disp.centY+disp.boxSize];
boxR = [disp.centX(2)-disp.boxSize disp.centY-disp.boxSize...
    disp.centX(2)+disp.boxSize disp.centY+disp.boxSize];
% Drawing the fixation box on the left:
drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, disp.boxColour);
% Drawing the text box centered in the left box:
DrawFormattedText(wPtr, textInstr, 'center', 'center', [255 255 50], 15, ...
    [], [], [], [], boxL)
% Drawing the fixation and text boxes on the right:
drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, disp.boxColour);
DrawFormattedText(wPtr, textInstr, 'center', 'center', [255 255 50], 15, ...
    [], [], [], [], boxR)
Screen(wPtr, 'Flip');

%% Monitoring for the key presses during the instruction:
while 1,
%     [keyIsDownPrev, secsPrev, keyCodePrev] = KbQueueCheck(deviceIndex);
    [keyIsDown, ~, keyCode] = KbCheck(deviceIndex); % not sure if this is the
        % optimal command in terms of waiting.
    if keyIsDown,
        % If the quit key is pressed, quit:
        if keyCode(KbName(quitkey)),
            Screen('CloseAll');
            ShowCursor;
            return;
        else
        % If any other key is pressed, proceed after .5 seconds:
            break;
        end
    end
end
WaitSecs(0.5);

%% Preparing the gratings and screen.

% Run the movie animation for a fixed period.
% curTrial = 5; %TEMP
frameRate=Screen('FrameRate',screenid);
% If MacOSX does not know the frame rate the 'FrameRate' will return 0.
% That usually means we run on a flat panel with 60 Hz fixed refresh
% rate:
if frameRate == 0
    frameRate = 60;
end

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(wPtr);
Priority(priorityLevel);

%% Going through the trials.

for curTrial=1:numofTrials
    %% Drawing the gratings for this trial.
    curCond = trialCond(curTrial);

    % Compute each frame of the movie and convert the those frames, stored in
    % MATLAB matices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    numFrames = round(c2n(condTable,'tL',curCond) / 1000 * frameRate);
%     round( ((1/c2n(condTable,'sfL',curTrial))/ ...
%         c2n(condTable,'speedL',curTrial)) * frameRate );

    for i = 1:numFrames
        %% Drawing the left grating.
        phase = (i/numFrames)*2*pi; % this is the same for both eyes
        stimsz = round(cm2px(c2n(condTable,'sizeL',curCond), sdims));
        [x,y] = meshgrid(-stimsz:stimsz, -stimsz:stimsz);
        angle = c2n(condTable,'dirL',curCond)*pi/180;
        f = px2cm (c2n(condTable,'sfL',curCond), sdims)*2*pi; % 0.0142 cy/px = .5 cy/cm
        a = cos(angle) * f;
        b = sin(angle) * f;
        m = exp(-((x/90).^2)-((y/90).^2)) .* sin(a*x+b*y+phase);
        % multiplying by a flat top window (for black surround & artefact
        % reduction):
        [r, c] = size(m);
        wc = window(@flattopwin, c); % TODO: nuttallwin or flattopwin?
        wr = window(@flattopwin, r);
        [maskr,maskc]=meshgrid(wr,wc);
        ftwindow = maskr.*maskc;
        ftwindow(ftwindow<0)=0; % not really needed with nuttallwin
        gratL(i) = Screen('MakeTexture', wPtr, ftwindow.*(gray+inc*m) ); %#ok<AGROW>
        %% Preparing the right grating.
        stimsz = round(cm2px(c2n(condTable,'sizeR',curCond), sdims));
        [x,y] = meshgrid(-stimsz:stimsz, -stimsz:stimsz);
        angle = c2n(condTable,'dirR',curCond)*pi/180;
        f = px2cm (c2n(condTable,'sfR',curCond), sdims)*2*pi; % 0.0142 cy/px = .5 cy/cm
        a = cos(angle) * f;
        b = sin(angle) * f;
        m = exp(-((x/90).^2)-((y/90).^2)) .* sin(-a*x+b*y+phase);
        % multiplying by a function window for black surround & artefact reduction):
        [r, c] = size(m);
        wc = window(@flattopwin, c); % TODO: nuttallwin or flattopwin?
        wr = window(@flattopwin, r);
        [maskr,maskc]=meshgrid(wr,wc);
        ftwindow = maskr.*maskc;
        ftwindow(ftwindow<0)=0; % not really needed with nuttallwin
        gratR(i) = Screen('MakeTexture', wPtr, ftwindow.*(gray+inc*m) ); %#ok<AGROW>
    end

    %% Animation loop.

    % Convert movieDuration in seconds to duration in frames to draw:
    %TEMP; should be div/1000
    % stimDur_frames=round(condTable.tL(curTrial) / 100 * frameRate);
    stimDur_frames=numFrames;
    stimFrameIndcs=mod(0:(stimDur_frames-1), numFrames) + 1;
    
    % Random starting phase:
    randStartPhase = randi(stimDur_frames);
    frameSet = [randStartPhase:stimDur_frames 1:randStartPhase-1];

    for i=frameSet
        % Draw the left image:
        boxMulti = 2.2;
        Screen('DrawTexture', wPtr, gratL(stimFrameIndcs(i)), [], ...
            [disp.centX(1)-disp.boxSize*boxMulti disp.centY-disp.boxSize*boxMulti ...
            disp.centX(1)+disp.boxSize*boxMulti disp.centY+disp.boxSize*boxMulti]);
        drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, ...
            disp.boxColour);
        % Draw the right image:
        Screen('DrawTexture', wPtr, gratR(stimFrameIndcs(i)), [], ...
            [disp.centX(2)-disp.boxSize*boxMulti disp.centY-disp.boxSize*boxMulti ...
            disp.centX(2)+disp.boxSize*boxMulti disp.centY+disp.boxSize*boxMulti]);
        drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, ...
            disp.boxColour);
        Screen('Flip', wPtr);

        %% Monitoring for keypresses.
        [keyIsDown, ~, keyCode] = KbCheck(deviceIndex);
        if keyIsDown,
            if keyCode(KbName(quitkey)),
                Screen('CloseAll');
                ShowCursor;
                return;
            end
        end
    end
    
    %% Continue screen.
    % Drawing the fixation box on the left:
    drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, disp.boxColour);
    % Drawing the text box centered in the left box:
    DrawFormattedText(wPtr, textNextTrial, 'center', 'center', [255 255 50], 15, ...
        [], [], [], [], boxL)
    % Drawing the fixation and text boxes on the right:
    drawFixationBox(wPtr, disp.centX(2), disp.centY, disp.boxSize, disp.boxColour);
    DrawFormattedText(wPtr, textNextTrial, 'center', 'center', [255 255 50], 15, ...
        [], [], [], [], boxR)
    Screen(wPtr, 'Flip');
    % Monitoring for a keypress:
    while 1,
        [keyIsDown, ~, keyCode] = KbCheck(deviceIndex); % not sure if this is the
            % ... optimal command in terms of waiting.
        if keyIsDown,
            % If the 'quit' key is pressed, quit:
            if keyCode(KbName(quitkey)),
                Screen('CloseAll');
                ShowCursor;
                return;
            elseif keyCode(KbName(space)),
            % If the 'space' key is pressed, proceed to the next trial:
                break;
            end
        end
    end

end

Screen('CloseAll');

catch %#ok<CTCH>
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end