% example - An example of how to use the RegionalDyssynchrony class to
% compute regional delay times
%
% Copyright (c) 2012, Jonathan Suever

% Create some regional delay curves that are shifted verions of a sinusoid
nCurves = 100;
shifts = round(randn(nCurves, 1) * 10);

% "True" reference curve
curve = sin(linspace(0, 2*pi, 101));
curve(end) = [];

% Create all of the shifted curves
curves = repmat(curve, [nCurves, 1]);
for k = 1:numel(shifts)
    curves(k,:) = circshift(curve, [0 -shifts(k)]);
end

% Now plot all of the regional displacement curves
figure('Position', [100 100 1423 977]);
subplot(2,2,1);
plot(curves');
title('Regional Displacement Curves');
xlabel('Time')
ylabel('Radial Displacement');

% Construct the RegionalDyssynchrony object
obj = RegionalDyssynchrony(curves);

% Set the threshhold for QTClustering (default value is 0.3265)
obj.ClusteringThreshold = 0.3265;

[delays, correlations] = obj.computeRegionalDelays();

% Look at the cluster sizes
subplot(2,2,2)
bar(obj.ClusterSizes);
xlabel('Cluster Index')
ylabel('Number in Cluster');
title('Cluster Sizes')

% Now we can look at the reference curve
subplot(2,2,3);
plot(obj.ReferenceCurve)
title('Reference Curve')
xlabel('Time')
ylabel('Radial Displacement');

% Now we can compare the delays that we got with the shifts that we
% prescribed
subplot(2,2,4)
plot(delays)
hold on
plot(shifts / nCurves, 'r')
title('Results')
xlabel('Samples')
ylabel('Delay Time (%)')
legend({'Computed Delays', 'Prescribed Delays'})

% The differences here are expected since the reference curve isn't exactly
% equal to our input curve with zero shift since some filtering is applied
% and QT Clustering is used.
