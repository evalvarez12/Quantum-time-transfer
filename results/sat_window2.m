
close('all')

set(0, 'DefaultAxesFontName', 'Times');

startTime = datetime(2020,5,12,0,0,0);
stopTime = startTime + hours(24);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);


alt = 500e3 + 6371e3;

% sat = walkerDelta(sc, alt, 56, 240, 3, 1, ...
%        ArgumentOfLatitude=15, Name="Galileo");

sat = walkerDelta(sc, alt, 56, 680, 40, 3, ...
       ArgumentOfLatitude=15, Name="Galileo");

names = sat.Name + " Camera"; 
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",60);


name = "Ground station";
minElevationAngle = 40; % degrees
gs = groundStation(sc, ...
    "Name",name, ...
    "MinElevationAngle",minElevationAngle, ...
    "Latitude", -34.7201, ...
    "Longitude", 138.659);


ac = access(cam,gs);

% Properties of access analysis objects
ac(1)


% v = satelliteScenarioViewer(sc,"ShowDetails",true);
% sat(1).ShowLabel = true;
% gs.ShowLabel = true;
% show(sat(1));
% 
% fov = fieldOfView(cam([cam.Name] == "Galileo_1 Camera"));



pointAt(sat,gs);

ai = accessIntervals(ac);
sortrows(ai, 4)



for idx = 1:numel(ac)
    [s,time] = accessStatus(ac(idx));
    
    if idx == 1
        % Initialize system-wide access status vector in the first iteration
        systemWideAccessStatus = s;
    else
        % Update system-wide access status vector by performing a logical OR
        % with access status for the current camera-site access
        % analysis
        systemWideAccessStatus = or(systemWideAccessStatus,s);
    end
end


% plot(time,systemWideAccessStatus,"LineWidth",2);
% grid on;
% xlabel("Time");
% ylabel("System-Wide Access Status");

l = size(systemWideAccessStatus);
l = l(2);

s = sum(systemWideAccessStatus);

s/l

