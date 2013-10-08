file = '/Users/jearly/Desktop/QuasigeostrophyTracers.nc';
FramesFolder ='/Users/jearly/Desktop/QuasigeostrophyTracerFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder) == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x');
y = ncread(file, 'y');
% drifterID = ncread(file, 'drifterID');
t = ncread(file, 'time');

timeScale = 12;
distanceScale = 47;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%
stride = 4;
size = 2;

xposInitial = ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1], [length(y)/stride length(x)/stride 1], [stride stride 1]);
xposInitial = reshape(xposInitial, length(y)*length(x)/(stride*stride), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for iTime=1:length(t)
% for iTime=20:20
% 	tracer = double(ncread(file, 'x-tracer', [1 1 iTime], [length(y) length(x) 1], [1 1 1]));
	
	xpos = ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]);
	ypos = ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]);
	xpos = reshape(xpos, length(y)*length(x)/(stride*stride), 1);
	ypos = reshape(ypos, length(y)*length(x)/(stride*stride), 1);
% 	
% 	[X,Y]=meshgrid(x,y);
% 	tracerInterp = interp2( X, Y, tracer, xpos, ypos);

% 	xpos = ncread(file, 'x-position', [1 1 iTime], [length(y) length(x) 1]);
% 	ypos = ncread(file, 'y-position', [1 1 iTime], [length(y) length(x) 1]);
% 	xpos = reshape(xpos, length(y)*length(x), 1);
% 	ypos = reshape(ypos, length(y)*length(x), 1);

% 	xpos = ncread(file, 'x-position', [1 iTime], [length(drifterID) 1]);
% 	ypos = ncread(file, 'y-position', [1 iTime], [length(drifterID) 1]);
	
% 	[X,Y]=meshgrid(x,y);
% 	tracerInterp = interp2( X, Y, tracer, xpos, ypos);
	
	
	scatter(xpos, ypos, size*size, xposInitial, 'filled')	
	set( gca, 'TickDir', 'out');
	set(gca, 'Linewidth', 1.0);
	axis equal tight
	
	set( gca, 'xtick', [])
	
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
	
	title( sprintf('Floats advected by a Quasigeostrophic eddy, colored by initial position, day %d', round(t(iTime)*12)), 'fontsize', 16 );
	
	cb = colorbar( 'location', 'SouthOutside' );
	set(get(cb,'xlabel'),'String', 'Distance (Rossby radius)', 'FontSize', 12.0);
	set( gca, 'clim', [min(x) max(x)] );
	
% 	set(gca,'LooseInset',get(gca,'TightInset'))
	
	output = sprintf('%s/Day_%03d', FramesFolder,iTime-1);
	print('-depsc', output)
% 	print('-r144', '-djpeg', output );
end