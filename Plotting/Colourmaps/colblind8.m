function cmap = colblind8()
% 8 colour blind friendly colours: dark blue, dark green, light green, light
% blue, pale yellow, pink, purple/pink, plum/maroon
% From: https://davidmathlogic.com/colorblind/#%23332288-%23117733-%2344AA99-%2388CCEE-%23DDCC77-%23CC6677-%23AA4499-%23882255
% 07/03/2024

cmap = ones(5,3);
cmap(1,:) = [51,34,136]/256;    % dark blue
cmap(2,:) = [17,119,51]/256;    % dark green
cmap(3,:) = [68,170,153]/256;   % light green
cmap(4,:) = [136,204,238]/256;  % light blue
cmap(5,:) = [221,204,119]/256;  % pale yellow
cmap(6,:) = [204,102,119]/256;  % pink
cmap(7,:) = [170,68,153]/256;   % purple/pink
cmap(8,:) = [136,34,85]/256;    % plum/maroon
end