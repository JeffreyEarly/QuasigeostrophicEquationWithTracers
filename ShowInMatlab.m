file = '/Users/jearly/Desktop/LinearQuasigeostrophy.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
ssh = double(ncread(file, 'SSH'));

figure
subplot(2,1,1)
pcolor(x, y, ssh(:,:,1)), axis equal tight, shading interp
title('Gaussian eddy, day 0')
subplot(2,1,2)
pcolor(x, y, ssh(:,:,end)), axis equal tight, shading interp
title(sprintf('Gaussian eddy, day %d', round(t(end)*12)))