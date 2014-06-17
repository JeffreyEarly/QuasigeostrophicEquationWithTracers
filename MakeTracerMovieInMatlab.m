file = '/Users/jearly/Desktop/QuasigeostrophyTracers.nc';
FramesFolder ='/Users/jearly/Desktop/QuasigeostrophyTracerFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x')/1000;
y = ncread(file, 'y')/1000;
t = ncread(file, 'time')/86400;

timeScale = ncreadatt(file, '/', 'time_scale');
distanceScale = ncreadatt(file, '/', 'length_scale');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%
stride = 2;
floatSize = 4;

% Read in the initial position of the floats.
% We will use this information to maintain a constant color on each float.
xposInitial = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1], [length(y)/stride length(x)/stride 1], [stride stride 1]))/1000;
xposInitial = reshape(xposInitial, length(y)*length(x)/(stride*stride), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

%for iTime=1:length(t)
for iTime=20:20

	% read in the position of the floats for the given time	
	xpos = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]))/1000;
	ypos = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]))/1000;
	
	% make everything a column vector
	xpos = reshape(xpos, length(y)*length(x)/(stride*stride), 1);
	ypos = reshape(ypos, length(y)*length(x)/(stride*stride), 1);
	
	% default color map is only 128 shades---we need more!
	colormap(jet(1024))	
	
	% now plot the floats, colored by initial position
	% Scatter works, but is substantially slower than using mesh.
	% scatter(xpos, ypos, floatSize*floatSize, xposInitial, 'filled')	
	mesh([xpos';xpos'],[ypos';ypos'],[xposInitial';xposInitial'],'mesh','column','marker','.','MarkerSize',floatSize*floatSize), view(2)
	grid off
	
	% make the axes look better
	set( gca, 'TickDir', 'out');
	set( gca, 'Linewidth', 1.0);
	axis equal tight
	
	% get rid of the xticks because we're going to put a colorbar below with the same info.
	set( gca, 'xtick', [])
	
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
	
	% label everything
	title( sprintf('Floats advected by a Quasigeostrophic eddy, colored by initial position, day %d', round(t(iTime))), 'fontsize', 28, 'FontName', 'Helvetica' );
	ylabel( 'distance (km)', 'FontSize', 24.0, 'FontName', 'Helvetica');
	
	% add a color bar
	cb = colorbar( 'location', 'SouthOutside' );
	set(get(cb,'xlabel'),'String', 'distance (km)', 'FontSize', 24.0, 'FontName', 'Helvetica');
	set( gca, 'clim', [min(x) max(x)] );
	
	% write everything out	
	output = sprintf('%s/Day_%03d', FramesFolder,iTime-1);
	print('-depsc2', output)
end