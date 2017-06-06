function [amp, phase, offset, trace, stats] = fitsine(time, data, freq, ploton, useconstant)
% [amp, phase, offset, trace] = fitsine(time, data, freq, ploton)
% [amp, phase, offset, trace] = fitsine(dt, data, freq, ploton)
% Row-ize


if length(time)==1
    dt = time;
    time = (1:length(data))'*dt - dt;
end
time = time(:)';
data = data(:);

keep = ~isnan(data);

y1 = sin(2*pi*freq*time);
y2 = cos(2*pi*freq*time);
constant = ones(size(y1)); % 1/6/13 Only if you want to allow an offset

if exist('useconstant','var') && ~useconstant
    vars = [y1; y2];
else
    vars = [y1; y2; constant];
end
[fit, ~,residuals,~,stats]= regress(data(keep), vars(:,keep)');
% rsquare=stats(1);
trace = fit'*vars;
trace(~keep) = NaN;

if exist('useconstant','var') && ~useconstant
    offset  = 0;
else
    offset = fit(3);
end
% stdresid = std(residuals);

amp = sqrt(fit(1)^2+fit(2)^2);
phase = atan2d(fit(2), fit(1));
% trace2 = amp*sin(2*pi*freq*time + phase*pi/180);

% const = fit(3);
%% Debugging - plot
if exist('ploton','var') && ploton
    clf
    plot(time, data); hold on
    plot(time, trace,'r')
%     plot(time, trace2,'g--')
    ylim([min(data)-.1 max(data)+.1])
end
