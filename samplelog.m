function [t_log,d_log] = samplelog(time, data)

num_points = 50;  % Number of points you want in the new logarithmic scale
t_log = logspace(log10(time(find(time>0,1))), log10(time(end)), num_points);

% Interpolate the data to match the new time vector
d_log = interp1(time, data, t_log, 'pchip');
end