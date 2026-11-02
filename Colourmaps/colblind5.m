function cmap = colblind5()
% 5 colour blind friendly colours: blue, purple, pink, orange, yellow
% From: https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
% 06/03/2024

colblue = [100,143,255]/256;
colpurp = [120,94,240]/256;
colpink = [220,38,127]/256;
colorange = [254,97,0]/256;
colyell = [255,176,0]/256;

cmap = ones(5,3);
cmap(1,:) = colblue;
cmap(2,:) = colpurp;
cmap(3,:) = colpink;
cmap(4,:) = colorange;
cmap(5,:) = colyell;
end