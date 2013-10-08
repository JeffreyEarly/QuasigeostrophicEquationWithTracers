file = '/Users/jearly/Desktop/QuasigeostrophyTracers.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time')*ncreadatt(file, '/', 'time_scale');
ssh = double(ncread(file, 'SSH'));
tracer = double(ncread(file, 'x-tracer'));

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(4) scrsz(4)/2])

subplot(2,2,1)
pcolor(x, y, ssh(:,:,1)), axis equal tight, shading interp
title('SSH of gaussian eddy, day 0')
subplot(2,2,3)
pcolor(x, y, ssh(:,:,end)), axis equal tight, shading interp
title(sprintf('SSH of gaussian eddy, day %d', round(t(end))))


subplot(2,2,2)
pcolor(x, y, tracer(:,:,1)), axis equal tight, shading interp
title('Tracer advected by a gaussian eddy, day 0')
subplot(2,2,4)
pcolor(x, y, tracer(:,:,end)), axis equal tight, shading interp
title(sprintf('Tracer advected by a gaussian eddy, day %d', round(t(end))))