function cx = SmoothFiniteDifference( time, xPosition )
% Smooooooth speed, by simply taking the difference in position

sparseXPosition = xPosition(1);
sparseXTime = time(1);
for day = 1:length(xPosition)
	if (sparseXPosition(end) ~= xPosition(day))
		sparseXPosition = [sparseXPosition xPosition(day)];
		sparseXTime = [sparseXTime time(day)];
	end
end

interpXTime = time(find( time == sparseXTime(1), 1, 'first'):find( time == sparseXTime(end), 1, 'first') );

interpXPosition = interp1( sparseXTime, sparseXPosition, interpXTime, 'linear');

[x cx] = RunningFiniteDifference(interpXTime, interpXPosition, 51);
%cx = RunningAverage( cx, 7 );
cx = cx';

while length(cx) < length(time)
	cx = [cx cx(end)];
end
