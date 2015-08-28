function [ imageWindow ] = renderWindow ( m )
% Renders a nuttallwin / flattopwin window based on matrix dimensions of
% gratM.

[r, c] = size(m);
wc = window(@nuttallwin,c); % TODO: nuttallwin or flattopwin?
wr = window(@nuttallwin,r);
[maskr,maskc] = meshgrid(wr,wc);
imageWindow = maskr.*maskc;
imageWindow(imageWindow<0)=0; % not really needed with nuttallwin

end

