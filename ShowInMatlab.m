file = '/Users/jearly/Desktop/QuasigeostrophyTracers.nc';
trackfile = '/Users/jearly/Desktop/QuasigeostrophyTracers_tracks.mat';

shouldDisplayTracks = 1;
if exist(trackfile,'file')
	load(trackfile);
else
	shouldDisplayTracks = 0;
end

endTime = 201;


x = ncread(file, 'x')/1000;
y = ncread(file, 'y')/1000;
t = ncread(file, 'time')/86400;
ssh = double(ncread(file, 'SSH'));
tracer = double(ncread(file, 'x-tracer'));

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)/2])

subplot(2,2,1)
pcolor(x, y, ssh(:,:,1)), axis equal tight, shading interp
title('SSH of gaussian eddy, day 0')

subplot(2,2,3)
pcolor(x, y, ssh(:,:,endTime)), axis equal tight, shading interp
title(sprintf('SSH of gaussian eddy, day %d', round(t(endTime))))


subplot(2,2,2)
pcolor(x, y, tracer(:,:,1)), axis equal tight, shading interp
title('Tracer advected by a gaussian eddy, day 0')

subplot(2,2,4)
pcolor(x, y, tracer(:,:,endTime)), axis equal tight, shading interp
if (shouldDisplayTracks == 1)
	hold on
	plot([eddy.rv(1:endTime).x]/1000, [eddy.rv(1:endTime).y]/1000, 'LineWidth', 1, 'Color', 'black')
	plot(eddy.rv(endTime).xContour/1000, eddy.rv(endTime).yContour/1000, 'LineWidth', 2, 'Color', 'black')
	plot(eddy.trappedFluid(endTime).xContour/1000, eddy.trappedFluid(endTime).yContour/1000, 'LineWidth', 2, 'Color', [0.3 0.3 0.3])
end
title(sprintf('Tracer advected by a gaussian eddy, day %d', round(t(endTime))))