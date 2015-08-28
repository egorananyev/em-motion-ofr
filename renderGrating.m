function [ gratM ] = renderGrating( stimSize, stimDir, stimsf, stimPhase )
% Creates a grating based on stimulus size (in pixels), stimulus direction

[x,y] = meshgrid(-stimSize:stimSize, -stimSize:stimSize);
a = cos(stimDir) * stimsf;
b = sin(stimDir) * stimsf;
gratM = exp(-((x/90).^2)-((y/90).^2)) .* sin(a*x+b*y+stimPhase);

end