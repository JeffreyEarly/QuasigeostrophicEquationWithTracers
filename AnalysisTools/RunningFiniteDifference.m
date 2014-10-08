function [xnew smoothVec] = RunningFiniteDifference( x, y, smoothness )
%	Returns a running average of vec. Tries to handle end points well.

smoothnessHalfLength = (smoothness-1)/2;

smoothVec = [];

for index=2:(length(y)-1)
	% Determine the largest width we're allowed for our running average=
	startDistance = index-1;
	endDistance = length(y) - index;
	restrictionDistance = startDistance;
	if ( endDistance < restrictionDistance)
		restrictionDistance = endDistance;
	end
	
	if ( smoothnessHalfLength < restrictionDistance)
		restrictionDistance = smoothnessHalfLength;
	end
	
	final = index+restrictionDistance;
	initial = index-restrictionDistance;
	
	smoothVec = [smoothVec; (y(final) - y(initial))/(x(final)-x(initial))];
end

xnew = x(2:(length(y)-1));


