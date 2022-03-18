function str = format_duration(duration)
%FORMAT_DURATION Summary of this function goes here
%   Detailed explanation goes here

days     = floor(duration / (60*60*24));
duration = mod(duration, 60*60*24);
hours    = floor(duration / (60*60));
duration = mod(duration, 60*60);
minutes  = floor(duration / 60);
duration = mod(duration, 60);
seconds  = floor(duration);
duration = duration - seconds;

if days > 0
    str = sprintf('%dd ', days);
else
    str = '';
end
if hours > 0
    str = [str, sprintf('%dh ', hours)];
end
if minutes > 0
    str = [str, sprintf('%dm ', minutes)];
end
if days == 0
    if seconds > 0
        str = [str, sprintf('%ds ', seconds)];
    end
    if minutes == 0
        str = [str, sprintf('%dms', floor(duration*1e3))];
    end
end
str = strtrim(str);
