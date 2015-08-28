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
Screen('Preference', 'SkipSyncTests', 1); % a necessary evil

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

% Screen('CloseAll'); %temp

for curTrial=1:numofTrials
    %% Drawing the gratings for this trial.
    curCond = trialCond(curTrial);
    display('=====');
    display(sprintf('current trial: %3i', curTrial));
    display(sprintf('current condition: %i', curCond));
    
    % The number of frames for rendering and presentation:
    numFrames = round(c2n(condTable,'tL',curCond) / 1000 * frameRate);
    
    % Random starting phase:
    randStartPhase = rand(1);
    
    % Preparing the stimulus configurations based on the trial:
    sizeL = c2n(condTable,'sizeL',curCond);
    dirL = c2n(condTable,'dirL',curCond);
    sfL = c2n(condTable,'sfL',curCond);
%     speedL = c2n(condTable,'speedL',curCond) * (numFrames/frameRate);
    speedL = c2n(condTable,'speedL',curCond);
    display(sprintf('left stim size (px): %i', sizeL));
    display(sprintf('left orientation/direction (deg) = %i', dirL));
    display(sprintf('left spatial frequency (cy/px) = %i', sfL));
    sizeR = c2n(condTable,'sizeL',curCond);
    dirR = c2n(condTable,'dirR',curCond);
    sfR = c2n(condTable,'sfR',curCond);
    speedR = c2n(condTable,'speedR',curCond);
    display(sprintf('right stim size (cm): %i', sizeR));
    display(sprintf('right orientation/direction (deg) = %i', dirR));
    display(sprintf('right spatial frequency (cpd)= %i', sfR));
    % Modifying configs into the right dimensions:
    dirL = dirL*pi/180;
    dirR = dirR*pi/180;
    sfL = px2cm(sfL, sdims)*2*pi; % 0.0142 cy/px = .5 cy/cm
    sfR = px2cm(sfR, sdims)*2*pi;
    sizeL = round(cm2px(sizeL, sdims));
    sizeR = round(cm2px(sizeR, sdims));
    
    %% Rendering the gratings:
    for i = 1:numFrames
        % Drawing the left grating.
        phaseL = ((i/numFrames)+randStartPhase)*2*pi*speedL*sfL;
        m = renderGrating(sizeL, dirL, sfL, phaseL);
        ftwindow = renderWindow(m);
        gratL(i) = Screen('MakeTexture', wPtr, ftwindow.*(gray+inc*m) ); %#ok<AGROW>
        % Preparing the right grating.
        phaseR = ((i/numFrames)+randStartPhase)*2*pi*speedR*sfR;
        m = renderGrating(sizeR, dirR, sfR, phaseR);
        ftwindow = renderWindow(m);
        gratR(i) = Screen('MakeTexture', wPtr, ftwindow.*(gray+inc*m) ); %#ok<AGROW>
    end

    %% Animation loop.

    for i=1:numFrames %frameSet
        % Draw the left image:
        boxMulti = 1.8;
        Screen('DrawTexture', wPtr, gratL(i), [], ...
            [disp.centX(1)-disp.boxSize*boxMulti disp.centY-disp.boxSize*boxMulti ...
            disp.centX(1)+disp.boxSize*boxMulti disp.centY+disp.boxSize*boxMulti]);
        drawFixationBox(wPtr, disp.centX(1), disp.centY, disp.boxSize, ...
            disp.boxColour);
        % Draw the right image:
        Screen('DrawTexture', wPtr, gratR(i), [], ...
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