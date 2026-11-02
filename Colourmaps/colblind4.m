function cmap = colblind4()
% 4 colour blind friendly colours: pink, blue, yellow, green
% From: https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
% 06/03/2024

colpink = [216,27,96]/256;
colblue = [30,136,229]/256;
colyell = [255,193,7]/256;
colgreen = [0,77,64]/256;

cmap = ones(4,3);
cmap(1,:) = colpink;
cmap(2,:) = colblue;
cmap(3,:) = colyell;
cmap(4,:) = colgreen;
end