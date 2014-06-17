%function eddy = TrappedFluid( NetCDFFile, TrackFile, TrappedFluidFile )

% First iteration divergence was calculated with shifted height -- RV stayed very constant
% Second iteration divergence was calculated with height
% Third iteration all values were multiplied by height (to get volume). Divergence versus time looks promising
% Fourth iteration, divergence with shifted height, but multiplied by height. Osillates, and/or stays near zero.
% Fifth iteration, shifted height, not multipled by height, by speed, VERY smoothed.
% Sixth iteration, I finally actually use the fluid depth, instead of the surface height
% Seventh iteration, whoops, really using the shifted divergence -- and computing the trapped volume.
% Eighth iteration I eliminated f. Which, doesn't require rerunning.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load the variables from the file
%

TrackFile = '/Users/jearly/Desktop/QuasigeostrophyTracers_tracks.mat';
NetCDFFile = '/Users/jearly/Desktop/QuasigeostrophyTracers.nc';

load(TrackFile);

[time, x, y] = FieldsFromTurbulenceFile(NetCDFFile, 1, 't', 'x', 'y');
time = time ./ 86400;

% Grab and create the global attributes
lengthScale = ncreadatt(NetCDFFile, '/', 'length_scale');
latitude = ncreadatt(NetCDFFile, '/', 'latitude');
EARTH_RADIUS = 6378100.0;
SIDEREAL_DAY = 86164.1;
g = 9.81;
coriolisParameter = 4.0 * pi * sin( latitude * pi / 180.0 ) / SIDEREAL_DAY;
equivalentDepth = (coriolisParameter * lengthScale).^2 / g;
betaParameter = 4.0 * pi * cos( latitude * pi / 180.0 ) / (SIDEREAL_DAY * EARTH_RADIUS);

if eddy.height(1).amplitude > 0.0
	isPositive = 1;
else
	isPositive = 0;
end

timeTracked = length(eddy.height);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the loop over time steps to find the trapped fluid region at each step.
%

a.contourAmplitude = 0.0;
a.contourArea = 0.0;
a.xContour = [];
a.yContour = [];

eddy.trappedFluid(timeTracked) = a;
eddy.AmbientVorticityTrappedSummed(timeTracked) = 0.0;

eddy.relativeVorticityTrappedSummed(timeTracked) = 0.0;
eddy.planetaryVorticityTrappedSummed(timeTracked) = 0.0;
eddy.surfaceHeightTrappedSummed(timeTracked)  = 0.0;
eddy.volumeTrapped(timeTracked) = 0.0;

eddy.relativeVorticityZeroContourSummed(timeTracked) = 0.0;
eddy.planetaryVorticityZeroContourSummed(timeTracked) = 0.0;
eddy.surfaceHeightZeroContourSummed(timeTracked)  = 0.0;

eddy.kineticEnergyTrappedSummed(timeTracked) = 0.0;
eddy.kineticEnergyZeroContourSummed(timeTracked) = 0.0;
eddy.potentialEnergyTrappedSummed(timeTracked) = 0.0;
eddy.potentialEnergyZeroContourSummed(timeTracked) = 0.0;

xSpeed = SmoothFiniteDifference(time*86400,[eddy.rv.x]);
ySpeed = SmoothFiniteDifference(time*86400,[eddy.rv.y]);


for timePoint= 1:25%length(timeTracked)

	if (mod(timePoint, 10) == 0)
		sprintf('Trapped Fluid Algorithm, time = %d', timePoint)
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Transform the surface height into moving coordinates
	%
    [x, y, surfaceHeight, relativeVorticityMag, u, v] = FieldsFromTurbulenceFile(NetCDFFile, timePoint, 'x', 'y', 'ssh', 'rv', 'u', 'v');
    
	fluidDepth = equivalentDepth + surfaceHeight;
	dx = x(2) - x(1);
    dy = y(2) - y(1);
    
    
	cy = xSpeed(timePoint) * repmat(y,[1,length(x)]);
	cx = ySpeed(timePoint) * repmat(x',[length(y),1]);
	
	heightShift = double(-cx * coriolisParameter/g + cy * coriolisParameter/g + surfaceHeight);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Subdomain Box
	%
	% Let's shrink our search domain for the largest contour to 3x the rv-zero-contour diameter.
	% Okay, the first time step doesn't converge if we don't set the domain to 5x.
	
	boxRadius = 3*sqrt(eddy.rv(timePoint).contourArea/pi);
	
	% Compute the range of indices to look for the height extremum.
	cellWidth = x(2) - x(1);
	indexDelta = ceil( boxRadius / cellWidth );
	xMin = eddy.rv(timePoint).xIndex - indexDelta; xMax = eddy.rv(timePoint).xIndex + indexDelta;
    
    % momentarily changing the lengh of the range to debug indices
	yMin = eddy.rv(timePoint).yIndex - ceil(indexDelta/1.2); yMax = eddy.rv(timePoint).yIndex + ceil(indexDelta/1.2);
	if (xMin < 1) xMin = 1; end; if (xMax > length(x)) xMax = length(x); end;
	if (yMin < 1) yMin = 1; end; if (yMax > length(y)) yMax = length(y); end;
	xRange = xMin:xMax; yRange = yMin:yMax;
	
	% Create a mesh grid of the subdomain, for use later
	[xGridSubdomain, yGridSubdomain] = meshgrid( x(xRange), y(yRange) );
    %[yGridSubdomain, xGridSubdomain] = meshgrid( y(yRange), x(xRange) );
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Transformed Height Extremum
	%
	% Estimate the transformed height extremum by looking at the unshifted extremum location
	mShiftedHeight = heightShift(eddy.height(timePoint).yIndex,eddy.height(timePoint).xIndex);
	if isPositive == 1
		minShiftedHeight = min(min(heightShift(yRange,xRange)));
	else
		minShiftedHeight = max(max(heightShift(yRange,xRange)));
	end
	
	searchHeight = mShiftedHeight;
	stillSearching = 1;

	while ( stillSearching == 1 )
		% Search by dropping the height 5% each iteration
		if isPositive == 1
			searchHeight = searchHeight - abs(mShiftedHeight-minShiftedHeight)*0.02;
		else
			searchHeight = searchHeight + abs(mShiftedHeight-minShiftedHeight)*0.02;
		end
		
		C = contourc(x(xRange), y(yRange), heightShift(yRange,xRange), [searchHeight searchHeight]);
		startPosition=1;
		while (startPosition < length( C(1,:)))
			theLength = C(2,startPosition);
			xContour = C(1, (startPosition+1):(startPosition+theLength));
			yContour = C(2, (startPosition+1):(startPosition+theLength));
			
			if ( xContour(1) == xContour(end) && yContour(1) == yContour(end) )
				in = inpolygon( eddy.height(timePoint).x, eddy.height(timePoint).y, xContour, yContour );
				if ( in(1) == 1)
					break;
				end
			else
				in(1) = 0;
			end
			startPosition = startPosition + theLength + 1;
		end
		
		% This is a candidate for the largest closed contour
		if ( (isPositive == 1 && (searchHeight < mShiftedHeight - abs(mShiftedHeight-minShiftedHeight))) || (isPositive == 0 && (searchHeight > mShiftedHeight + abs(mShiftedHeight-minShiftedHeight))) )
			stillSearching = 0;
		elseif ( in(1) == 1)
			lastContourMagnitude = searchHeight;
			lastXContour = xContour;
			lastYContour = yContour;
		else
			stillSearching = 0;
		end

	end
	
	% Record the trapped fluid contour
    eddy.trappedFluid(timePoint).contourAmplitude = lastContourMagnitude;
    eddy.trappedFluid(timePoint).contourArea = polyarea(lastXContour, lastYContour);
    eddy.trappedFluid(timePoint).xContour = lastXContour;
    eddy.trappedFluid(timePoint).yContour = lastYContour;
	
	% Create a matrix of 0s and 1s to identify the region inside the contour
	trappedFluidSubdomain = inpolygon( xGridSubdomain, yGridSubdomain, lastXContour, lastYContour );
	%figure, pcolor(xGridSubdomain,yGridSubdomain,trappedFluidSubdomain)
    
    
	% Convert this into an integration matrix, so multiply by the cell size.
	trappedFluidSubdomain = trappedFluidSubdomain * dx * dy;
	
	% Do the same thing for the relative vorticity zero contour
	rvZeroContourFluidSubdomain = inpolygon( xGridSubdomain, yGridSubdomain, eddy.rv(timePoint).xContour, eddy.rv(timePoint).yContour );
	rvZeroContourFluidSubdomain = rvZeroContourFluidSubdomain * dx * dy;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Relative vorticity
	%
	relativeVorticityMagTrappedSubdomain = relativeVorticityMag(yRange,xRange) .* trappedFluidSubdomain;
    
	eddy.relativeVorticityTrappedSummed(timePoint) = sum( sum( relativeVorticityMagTrappedSubdomain ));
	eddy.relativeVorticityZeroContourSummed(timePoint) = sum( sum( relativeVorticityMag(yRange, xRange) .* rvZeroContourFluidSubdomain ));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Planetary vorticity
	%
	[~, AmbientVorticity] = meshgrid(x, y);
	AmbientVorticity = betaParameter * AmbientVorticity;
	AmbientVorticityTrappedSubdomain = AmbientVorticity(yRange, xRange) .* trappedFluidSubdomain;
    
	eddy.planetaryVorticityTrappedSummed(timePoint) = sum( sum( AmbientVorticityTrappedSubdomain ));
	eddy.planetaryVorticityZeroContourSummed(timePoint) = sum( sum( AmbientVorticity(yRange,xRange) .* rvZeroContourFluidSubdomain ));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Vorticity scaled height
	%
	surfaceHeightTrappedSubdomain = - (coriolisParameter / equivalentDepth) * surfaceHeight(yRange,xRange) .* trappedFluidSubdomain; % .* fluidDepth(xRange, yRange)';
	eddy.surfaceHeightTrappedSummed(timePoint) = sum( sum( surfaceHeightTrappedSubdomain ));
	eddy.surfaceHeightZeroContourSummed(timePoint) = sum( sum( - (coriolisParameter / equivalentDepth) * surfaceHeight(yRange, xRange) .* rvZeroContourFluidSubdomain ));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Trapped fluid volume
	%
	volumeTrappedSubdomain = trappedFluidSubdomain .* fluidDepth(yRange, xRange);
	eddy.volumeTrapped(timePoint) = sum( sum( volumeTrappedSubdomain ));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Kinetic Energy
	%
	u2 = equivalentDepth * (u.^2 + v.^2);
	
	eddy.kineticEnergyTrappedSummed(timePoint) = sum(sum(u2(yRange, xRange) .* trappedFluidSubdomain));
	eddy.kineticEnergyZeroContourSummed(timePoint) = sum(sum(u2(yRange, xRange) .* rvZeroContourFluidSubdomain));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Potential Energy
	%
	PE = g * surfaceHeight(yRange, xRange) .* surfaceHeight(yRange, xRange);
	
	eddy.potentialEnergyTrappedSummed(timePoint) =  sum(sum(PE .* trappedFluidSubdomain));
	eddy.potentialEnergyZeroContourSummed(timePoint) = sum(sum(PE .* rvZeroContourFluidSubdomain));
end