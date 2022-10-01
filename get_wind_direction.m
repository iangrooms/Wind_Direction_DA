function theta = get_wind_direction(u,v)
% Gets wind direction, then adds observation error and categorizes the result
% Angles run from 0 (east) counterclockwise through 8 (west) to 15 (east by
% southeast).
% theta = mod(round((8/pi)*(atan2(v,u) + (pi/8)*(rand(size(u)) - 0.5))),16);
theta = mod(round((16/pi)*(atan2(v,u) + (pi/8)*((rand(size(u))+rand(size(u)))/2 - 0.5))),32);