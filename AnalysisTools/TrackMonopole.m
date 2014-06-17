function [eddy] = TrackMonopole( file, shouldSave )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load the variables from the file
%

[time, x, rv0] = FieldsFromTurbulenceFile(file, 1, 't', 'x', 'rv');
time = time ./ 86400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate the eddy propagation speed
%

lengthScale = ncreadatt(file, '/', 'length_scale');
timeScale = ncreadatt(file, '/', 'time_scale');
EddySpeed = lengthScale / timeScale;

% Convert the eddy speed from m/s into km/day
EddySpeed = EddySpeed * 86400 / 1000;

% Use this to put a cap on how far we'll allow the extremum to shift each time step.
maxDistanceAllowed = 2 * EddySpeed * ( time(2) - time(1) );
cellsPerTimeStep = ceil( maxDistanceAllowed / ( (x(2) - x(1) ) ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine whether the monopole's laplacian extremum is a min or max.
% Compute the laplacian to look for an extremum -- we don't care about exact value, just relative magnitude
%

maximumRV = max(max(rv0));
minimumRV = min(min(rv0));
if ( abs(maximumRV) > abs(minimumRV) )
	monopoleIsPositive = 1;
else
	monopoleIsPositive = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup the arrays -- necessary for some cases where we append.
%

a.amplitude=0.0;
a.x=0.0;
a.y=0.0;
a.xIndex=0.0;
a.yIndex=0.0;

a.contourAmplitude = 0.0;
a.contourArea = 0.0;
a.xContour = [];
a.yContour = [];

eddy.rv(length(time)) = a;
eddy.height(length(time)) = a;
eddy.uMax(length(time)) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Record the previous centroid location to provide additional criteria for eddy tracking.
%

previousXCentroid = 0;
previousYCentroid = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop through the file by time step, find the extremum, then determine if it shifted too far to be the same one.
%

% Matlab vectors start with index 1, NetCDF vectors with index 0.
for t=1:length(time)
	if mod(t, 10) == 0
		fprintf('Tracking Algorithm, time = %d\n', t)
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Compute the RV and fluid speed.
	%
	
    % Allow the possibility of a co-moving frame, so grab x and y again.
    [x, y, height, relativeVorticityMag, u, v] = FieldsFromTurbulenceFile(file, t, 'x', 'y', 'ssh', 'rv', 'u', 'v');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Create a mesh grid of the domain
	%
	
	[xGridSubdomain, yGridSubdomain] = meshgrid( x, y );
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Find the RV global extremum
	%
	
	% Compute the max relative vorticity
	if monopoleIsPositive == 1
		[C, I] = max(relativeVorticityMag);
		[~, J] = max(C);
	else
		[C, I] = min(relativeVorticityMag);
		[~, J] = min(C);
	end
	
	rvGlobalExtremum = relativeVorticityMag(I(J),J);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Use the RV global extremum to find candidate contours
	%
	contourHeight = rvGlobalExtremum/20; % RV *should* go through zero, so we're looking for a contour near zero.
	C = contourc(x, y, relativeVorticityMag, [contourHeight contourHeight]);
	startPosition = 1;
	
	candidateXCentroid = 0;
	candidateYCentroid = 0;	
	candidateXContour = 0;
	candidateYContour = 0;
	candidateLengthScale = 0;
	foundCandidate = 0;
	
	while (startPosition < length( C(1,:)))		
		%
		% Candidate contour
		%
		theLength = C(2,startPosition);
		xContour = C(1, (startPosition+1):(startPosition+theLength));
		yContour = C(2, (startPosition+1):(startPosition+theLength));
		startPosition = startPosition + theLength + 1;
		
		%
		% Skip if the radius is less than 5 km
		%
		if (  sqrt(polyarea(xContour, yContour)/3.14) < 5.0 )
			continue;
		end
		
		%
		% Determine the convex hull of the contour
		%
%		convexHullIndexes = convhull(xContour, yContour);
		convexHullIndexes = 1:length(xContour);
		
		%
		% Find the centroid of the convex hull of the contour
		%
		[ geom, ~, ~ ] = polygeom( xContour(convexHullIndexes), yContour(convexHullIndexes) );
		xCentroid = geom(2); yCentroid = geom(3);
		
		%
		% If the centroid moved more than the previous contour's diameter, bail out.
		% This is designed to be fairly generous...
		%
		if ( t > 1 )
			distanceTraveled = sqrt( (xCentroid-previousXCentroid)^2 + (yCentroid-previousYCentroid)^2 );
			previousRadius = sqrt(eddy.rv(t-1).contourArea/pi)/2;
			
			if (distanceTraveled > previousRadius)
				continue;
			end
		end
		
		%
		% We have a candidate, so lets compute its length scale.
		%
		ourLengthScale = sqrt(polyarea(xContour(convexHullIndexes), yContour(convexHullIndexes))/3.1415);
		
		%
		% If this candidate has a long lengther scale than the other candidate (or it's the first candidate), record it.
		% Maybe looking for the bigger height extremum would make more sense.
		if ( foundCandidate == 0 || ourLengthScale > candidateLengthScale)
			foundCandidate = 1;
			
			candidateXCentroid = xCentroid;
			candidateYCentroid = yCentroid;
			
			candidateXContour = xContour(convexHullIndexes);
			candidateYContour = yContour(convexHullIndexes);
			
			candidateLengthScale = ourLengthScale;		
		end

		
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	
	% This is our definition of the eddy, so let's record.
	%
	if ( foundCandidate == 1)
		eddy.rv(t).contourAmplitude = contourHeight;
        eddy.rv(t).xContour = candidateXContour;
        eddy.rv(t).yContour = candidateYContour;
		
		% Compute of the area of the enclosing contour in order to determine the length scale.
		eddy.rv(t).contourArea = polyarea(candidateXContour, candidateYContour);
		
		% Store the previous centroid's location
		previousXCentroid = candidateXCentroid;
		previousYCentroid = candidateYCentroid;
	else
		disp('TrackingAlgorithm breaking because foundCandidate != 1');
		break;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Save the location and magnitude of the relative vorticity extremum.
	%
	
	% Create a matrix of 0s and 1s to identify the region inside the contour
	rvZeroContourFluidSubdomain = inpolygon( xGridSubdomain, yGridSubdomain, candidateXContour, candidateYContour );
	
	% Compute the max relative vorticity
	if ( monopoleIsPositive == 1 )
        [C, I] = max(relativeVorticityMag .* rvZeroContourFluidSubdomain);
        [~, J] = max(C);
    else
        [C, I] = min(relativeVorticityMag .* rvZeroContourFluidSubdomain);
        [~, J] = min(C);
    end
	
    
	eddy.rv(t).amplitude = relativeVorticityMag(I(J),J);
	eddy.rv(t).x = x(J);
	eddy.rv(t).y = y(I(J));
	eddy.rv(t).xIndex = J;
	eddy.rv(t).yIndex = I(J);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Height
	%
	% The height extremum should be nearby the RV extremum. We we will look in
	% a box with diameter equal to the mRVZeroContourArea diameter.
	
	boxRadius = sqrt(eddy.rv(t).contourArea/pi);
	
	% Compute the range of indices to look for the height extremum.
	cellWidth = x(2) - x(1);
	indexDelta = ceil( boxRadius / cellWidth );
	xMin = eddy.rv(t).xIndex - indexDelta; xMax = eddy.rv(t).xIndex + indexDelta;
	yMin = eddy.rv(t).yIndex - indexDelta; yMax = eddy.rv(t).yIndex + indexDelta;
	if (xMin < 1), xMin = 1; end; if (xMax > length(x)), xMax = length(x); end;
	if (yMin < 1), yMin = 1; end; if (yMax > length(y)), yMax = length(y); end;
	xRange = xMin:xMax; yRange = yMin:yMax;
	
	% RV minimums correspond to height maximums, and vice-a-versa
	if monopoleIsPositive == 0
		[C, I] = max(height(yRange,xRange));
		[~, J] = max(C);
	else
		[C, I] = min(height(yRange,xRange));
		[~, J] = min(C);
	end
	
	yMin = yMin-1; xMin = xMin-1;
	% Save the location and magnitude of the height extremum.
	eddy.height(t).amplitude = height(yMin+I(J),xMin+J);
	eddy.height(t).x = x(xMin+J);
	eddy.height(t).y = y(yMin+I(J));
	eddy.height(t).xIndex = xMin+J;
	eddy.height(t).yIndex = yMin+I(J);
	
	% Determine the contour for the efold height
	efoldHeight = eddy.height(t).amplitude/exp(1);
	C = contourc(x, y, height, [efoldHeight efoldHeight]);
	startPosition=1;
	while (startPosition < length( C(1,:)))
		theLength = C(2,startPosition);
		xContour = C(1, (startPosition+1):(startPosition+theLength));
		yContour = C(2, (startPosition+1):(startPosition+theLength));
		in = inpolygon( eddy.height(t).x, eddy.height(t).y, xContour, yContour );
		if ( in(1) == 1)
			break;
		end
		startPosition = startPosition + theLength + 1;
	end
	
	% If this is true, we found a contour containing the extremum, so let's record it.
	if ( in(1) == 1)
		eddy.height(t).contourAmplitude = efoldHeight;
		eddy.height(t).xContour = xContour;
		eddy.height(t).yContour = yContour;
		
		% Compute of the area of the enclosing contour in order to determine the length scale.
		eddy.height(t).contourArea = polyarea(xContour, yContour);
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Find the maximum fluid speed within our box
	%
	eddy.uMax(t) = sqrt(max(max(u.^2 + v.^2)));
    lastRecordedTime = t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Did we bail early? If so, let's shorten the array.
%
if (lastRecordedTime < length(t))
    eddy.rv(lastRecordedTime:end) = [];
    eddy.height(lastRecordedTime:end) = [];
    eddy.uMax(lastRecordedTime:end) = [];
end

if shouldSave == 1
    [pathstr,name,~] = fileparts(file);
    outputFile = fullfile(pathstr,sprintf('%s_tracks.mat',name));
    save(outputFile, 'eddy');
end

return

% % Length Decay Rate in km/year
% LengthDecayRate = SmoothFiniteDifference(timeTracked, lengthRV) * 365;
% 
% % Height Decay Rate in cm/year
% HeightDecayRate = SmoothFiniteDifference(timeTracked, mHeight) * 100 * 365;
% 
% % Speeds, in m/s
% mRVXSpeed = SmoothFiniteDifference(timeTracked, mRVXPosition) * 1000 / 86400;
% mRVYSpeed = SmoothFiniteDifference(timeTracked, mRVYPosition) * 1000 / 86400;
% mRVSpeed = sqrt( mRVXSpeed.^2 + mRVYSpeed .^2 );
% 
% mHeightXSpeed = SmoothFiniteDifference(timeTracked, mHeightXPosition) * 1000 / 86400;
% mHeightYSpeed = SmoothFiniteDifference(timeTracked, mHeightYPosition) * 1000 / 86400;
% mHeightSpeed = sqrt( mHeightXSpeed.^2 + mHeightYSpeed .^2 );
% 
% save(outputFile, 'mRVZeroContourArea', 'mHeightEFoldContourArea', 'latitude', 'modeSpeed', 'EARTH_RADIUS', 'SIDEREAL_DAY', 'g', 'coriolisParameter', 'betaParameter', 'equivalentDepth', 'timeTracked', 'lengthRV', 'lengthEFold', 'mRVXSpeed', 'mRVYSpeed', 'mRVSpeed', 'mHeightXSpeed', 'mHeightYSpeed', 'mHeightSpeed', 'uMax', 'mRV', 'mRVXPosition', 'mRVYPosition', 'mRVXIndex', 'mRVYIndex', 'mRVContourMagnitude', 'mRVContourStartIndex', 'mRVContourLengthIndex', 'mRVXContour', 'mRVYContour', 'mHeight', 'mHeightXPosition', 'mHeightYPosition', 'mHeightXIndex', 'mHeightYIndex', 'mHeightContourMagnitude', 'mHeightContourStartIndex', 'mHeightContourLengthIndex', 'mHeightXContour', 'mHeightYContour', 'LengthDecayRate', 'HeightDecayRate')
