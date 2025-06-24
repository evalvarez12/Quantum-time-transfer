
close('all')

set(0, 'DefaultAxesFontName', 'Times');

startTime = datetime(2020,5,12,0,0,0);
stopTime = startTime + hours(24);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

% tleFile = "leoSatelliteConstellation.tle";
tleFile = 'Satellite_test.tle';
sat = satellite(sc,tleFile);


names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",60);


name = "Ground station";
minElevationAngle = 30; % degrees
gs = groundStation(sc, ...
    "Name",name, ...
    "MinElevationAngle",minElevationAngle, ...
    "Latitude", -34.7201, ...
    "Longitude", 138.659);


ac = access(cam,gs);

% Properties of access analysis objects
ac(1)


v = satelliteScenarioViewer(sc,"ShowDetails",false);
sat(1).ShowLabel = true;
gs.ShowLabel = true;
show(sat(1));

fov = fieldOfView(cam([cam.Name] == "Satellite Camera"));



pointAt(sat,gs);

access_table = accessIntervals(ac);


% timerange(startTime)


disp('---------------------------');
% Altitudes of the satellite
pos = states(sat(1));
hs = sqrt(pos(1,:).^2 + pos(2,:).^2 + pos(3,:).^2);

[nPasses, ll] = size(access_table.StartTime);

passes = [];

timeStep = seconds(1);

az_arr = [0];
el_arr = [0];
t_arr = [NaT(1, 'TimeZone', 'UTC')];
pass_groups = [0];

for i = 1:nPasses 

    timeA = access_table.StartTime(i);
    timeB = access_table.EndTime(i);
    
    nTime = floor((timeB - timeA)/timeStep);
    
    azs = zeros(1, nTime);
    els = zeros(1, nTime);
    ts = NaT(1, nTime, 'TimeZone', timeA.TimeZone);
    
    for j = 1:nTime
        t = timeA + j*timeStep;
        [az, el, ~] = aer(gs,sat(1), t);
    
        azs(j) = az;
        els(j) = el;
        ts(j) = t;
    end   
    
    passes = [passes; {[azs; els]}];

    az_arr = cat(2, az_arr, azs);
    el_arr = cat(2, el_arr, els);
    t_arr = cat(2, t_arr, ts);
    pass_groups = cat(2, pass_groups, ones(1, nTime)*i);

end

az_arr = az_arr(2:end);
el_arr = el_arr(2:end);
pass_groups = pass_groups(2:end);
t_arr = t_arr(2:end);

skyplot(az_arr(1:10:end), el_arr(1:10:end), GroupData=categorical(string(pass_groups(1:10:end))))


figure();
plot(t_arr, el_arr, 'o');

figure();
plot(hs- 6371e3);
% 
% t1 = t_arr(pass_groups == 1);
% t2 = t_arr(pass_groups == 2);
% t3 = t_arr(pass_groups == 3);
% t4 = t_arr(pass_groups == 4);
% t5 = t_arr(pass_groups == 5);
% 
% 
% el1 = el_arr(pass_groups == 1);
% el2 = el_arr(pass_groups == 2);
% el3 = el_arr(pass_groups == 3);
% el4 = el_arr(pass_groups == 4);
% el5 = el_arr(pass_groups == 5);
% 
% tt = t1(el1 < 35);
% disp('------------Pass 1');
% disp('35 - 45');
% seconds(tt(end) - tt(1))
% 
% disp('------------Pass 2');
% [~, mi] = max(el2);
% 
% 
% % disp('30 - 45');
% % ii = el2 < 35;
% % ii(mi:end) = 0;
% % tt = t2(ii);
% % seconds(tt(end) - tt(1))
% 
% 
% disp('35 - 45');
% ii = el2 < 45;
% ii(mi:end) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 55');
% ii = el2 < 55;
% ii(mi:end) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 65');
% ii = el2 < 65;
% ii(mi:end) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% 
% disp('65 - 75');
% ii = el2 < 75;
% ii(mi:end) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('75 - 85');
% ii = el2 < 85;
% ii(mi:end) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('85 - 90 - 85');
% ii = el2 > 85;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('85 - 75');
% ii = el2 < 85;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('75 - 65');
% ii = el2 < 75;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('65 - 55');
% ii = el2 < 65;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 45');
% ii = el2 < 55;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 35');
% ii = el2 < 45;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% disp('35 - 30');
% ii = el2 < 35;
% ii(1:mi) = 0;
% tt = t2(ii);
% seconds(tt(end) - tt(1))
% 
% 
% disp('-----------------Pass 3');
% [~, mi] = max(el3);
% 
% disp('35 - 45');
% ii = el3 < 45;
% ii(mi:end) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 55');
% ii = el3 < 55;
% ii(mi:end) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 65');
% ii = el3 < 65;
% ii(mi:end) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% 
% disp('65 - 70 - 65');
% ii = el3 < 70;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('65 - 55');
% ii = el3 < 65;
% ii(1:mi) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 45');
% ii = el3 < 55;
% ii(1:mi) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 35');
% ii = el3 < 45;
% ii(1:mi) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% disp('35 - 30');
% ii = el3 < 35;
% ii(1:mi) = 0;
% tt = t3(ii);
% seconds(tt(end) - tt(1))
% 
% 
% 
% disp('-------------Pass 4');
% [~, mi] = max(el4);
% 
% 
% disp('30 - 45');
% ii = el4 < 35;
% ii(mi:end) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% 
% disp('35 - 45');
% ii = el4 < 45;
% ii(mi:end) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 55');
% ii = el4 < 55;
% ii(mi:end) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 65');
% ii = el4 < 65;
% ii(mi:end) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% 
% disp('65 - 75');
% ii = el4 < 75;
% ii(mi:end) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('75 - 85 - 75');
% ii = el4 < 85;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('75 - 65');
% ii = el4 < 75;
% ii(1:mi) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('65 - 55');
% ii = el4 < 65;
% ii(1:mi) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('55 - 45');
% ii = el4 < 55;
% ii(1:mi) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('45 - 35');
% ii = el4 < 45;
% ii(1:mi) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('35 - 30');
% ii = el4 < 35;
% ii(1:mi) = 0;
% tt = t4(ii);
% seconds(tt(end) - tt(1))
% 
% disp('--------Pass 5');
% disp('30 - 35 - 30');
% ii = el5 < 35;
% tt = t5(ii);
% seconds(tt(end) - tt(1))
