classdef RegionalDyssynchrony < handle
% RegionalDyssynchrony - Class for performing regional dyssynchrony analysis
%
%   This is a MATLAB implementation of the delay time methodology presented
%   in the following manuscript.
%
%       Suever, JD, et al. Method to create regional mechanical dyssynchrony
%       maps from short-axis cine steady-state free-precession images.
%       Journal of Magnetic Resonance Imaging, 2013;16(4):4.
%
% USAGE:
%   r = RegionalDyssynchrony(curves)
%
% INPUTS:
%   curves: [N x T] Matrix, Rows are the different sampled regions and the
%           columns are the displacement (or strain) values over time.
%
% OUTPUTS:
%   r:      RegionalDyssynchrony object
%

%   Copyright (c) 2013, Jonathan Suever
%   Last Modified: 06-26-2012
%   Modified By: Suever (suever@gmail.com)

    properties
        TimeSeriesData                  % Data used to compute regional delay times

        Correlations                    % Maximum correlation values for each curve
        Delays                          % Delay time corresponding to maximum correlation

        ClusteredCurves                 % A cell array grouping TimeSeriesData into clusters
        ClusteringThreshold = 0.3265    % QT Clustering threshold (RMSE)
        ClusterSizes                    % Number of data in each cluster

        MaximumDetection = 'global'     % Method to use for finding maximum correlation
    end

    properties (Hidden)
        FullCorrelations                % All correlations values for various shifts relative to the reference
        FullDelays                      % All delays corresponding to the correlation values in FullCorrelations
    end

    properties (Access = 'private')
        ReferenceCurve_                 % Private variable for gracefully handling reference curve
    end

    properties (Dependent)
        ReferenceCurve                  % Reference curve to which all data is compared via xcorr
    end

    properties (Dependent, Hidden)
        NormalizedTimeSeriesData        % Normalized version of TimeSeriesData for clustering

        nLocations                      % Number of Sample Locations
        nPhases                         % Number of Phases in the data
    end

    methods
        function self = RegionalDyssynchrony(regionalData)
            % RegionalDyssynchrony - Constructor for RegionalDyssychrony
            %
            % USAGE:
            %   r = RegionalDyssynchrony(curves)
            %
            % INPUTS:
            %   curves: [N x T] Matrix, Rows are the different sampled regions and the
            %           columns are the displacement (or strain) values over time.
            %
            % OUTPUTS:
            %   r:      RegionalDyssynchrony object

            if ~exist('regionalData', 'var') || ~ismatrix(regionalData)
                error('RegionalDyssynchrony:InvalidInput',...
                      'Data must be of the form [nLocations, nPhases]');
            end
            self.TimeSeriesData = regionalData;
        end

        function varargout = computeRegionalDelays(self)
            % computeRegionalDelays - Computes the delay times throughout the LV
            %
            %   This method coordinates all of the regional dyssynchrony
            %   assessment including prepping the data, generating the
            %   reference curve using QT clustering, and determining regional
            %   delay times.
            %
            % USAGE:
            %   D = computeRegionalDelays(self)
            %
            % INPUTS:
            %   self:       Object instance
            %
            % OUTPUTS:
            %   D:          [nSectors x 1] Array, delay times throughout
            %               the left ventricle
            %
            % Last Modified: 06-26-2012
            % Modified By: Suever (suever@gmail.com)

            curves = self.NormalizedTimeSeriesData;

            self.ReferenceCurve_ = self.computeReferenceCurve(curves);

            dim = [self.nLocations, 1];

            self.FullCorrelations   = cell(dim);
            self.FullDelays         = cell(dim);
            [delays, corrs]         = self.xcorr(self.ReferenceCurve, curves);
            self.Delays             = reshape(delays, dim);
            self.Correlations       = reshape(corrs, dim);

            if nargout
                varargout = {delays, corrs};
            end
        end
    end

    methods (Access = 'protected')
        function [delays,corrs] = xcorr(self,reference,curves)
            % xcorr - Perform cross correlation of curves relative to reference
            %
            %   This function actually preps the data and performs the
            %   cross-correlation analysis necessary to obtain delay times
            %
            % USAGE:
            %   [delays,corrs] = xcorr(self,ref,curves)
            %
            % INPUTS:
            %   ref:    [1 x M] Array, reference curve obtained from clustering
            %   curves: [N x M] Matrix, All other curves that have to be
            %           compared to the reference
            %
            % OUTPUTS:
            %   delays: [nSectors x 1] Thickmat representation of delay
            %           times throughout the ventricle
            %   corrs:  [nSectors x 1] Thickmat representation of the
            %           correlation values associated with the delay times
            %
            % Last Modified: 06-26-2012
            % Modified By: Suever (suever@gmail.com)

            nPoints = size(curves,1);

            delays  = zeros(nPoints,1);
            corrs   = zeros(nPoints,1);

            limits = [-floor(self.nPhases / 2) ceil(self.nPhases / 2)];

            for i = 1:nPoints
                [corr,lags] = xcorrPeriodic(reference,curves(i,:));

                inds = (lags >= limits(1) & lags <= limits(2));

                corr = corr(inds);
                lags = lags(inds);

                self.FullCorrelations{i} = corr;
                self.FullDelays{i} = lags;

                switch self.MaximumDetection
                    case 'global'
                        [corrs(i),ndx] = max(corr(:));
                        delays(i) = lags(ndx);
                    case 'local'
                        dc = conv(corr([end 1:end 1]),[-0.5 0 0.5],'valid');

                        % Find sign changes
                        tmp = diff(sign(dc));

                        % Find the maxima
                        ndx = find(tmp > 0);

                        if isempty(ndx)
                            % Then just use the global maximum (because it's
                            % probably at an end)
                            [corrs(i),ndx] = max(corr(:));
                            delays(i) = lags(ndx);
                            continue
                        end

                        % Now find the one closes to zero
                        maxlags = lags(ndx);
                        [~,mindx] = min(abs(maxlags));

                        bestndx = ndx(mindx);

                        delays(i) = lags(bestndx);
                        corrs(i) = corr(bestndx);
                    otherwise
                        error('RegionalDyssynchrony:InvalidParameter',...
                              'MaximumDetection must be global or local')
                end
            end

            % Now we need to convert these delays to percentages
            delays = delays / self.nPhases;
        end

        function curve = computeReferenceCurve(self,curves)
            % computeReferenceCurve - Compute the reference via QT clustering
            %
            %   using the input curves we can perform QT clustering between all
            %   curves and group them into clusters. The largest cluster is
            %   considered to be the most "synchronous" and is averaged to
            %   obtain the reference curve
            %
            % USAGE:
            %   C = computeReferenceCurve(self,curves)
            %
            % INPUTS:
            %   curves: [M x N] Array, Each row represents a different curve to
            %           be clustered.
            %
            % OUTPUTS:
            %   C:      [1 x N] Array, Reference curve from QT clustering
            %
            % Last Modified: 01-25-2012
            % Modified By: Suever (suever@gmail.com)

            [indices] = qtClustering(curves,self.ClusteringThreshold,1);

            self.ClusteredCurves = cell(1,max(indices));

            for i = 1:numel(self.ClusteredCurves)
                self.ClusteredCurves{i} = curves(indices == i,:);
            end

            % Find the most prominent index
            h = hist(indices,1:max(indices));

            self.ClusterSizes = h;
            [~,biggest] = max(h);

            % Average all curves in this group
            curve = mean(curves(indices == biggest,:),1);
        end
    end

    methods
        %--- Get Methods ---%
        function res = get.nPhases(self)
            res = size(self.TimeSeriesData,2);
        end

        function res = get.nLocations(self)
            res = size(self.TimeSeriesData,1);
        end

        function res = get.NormalizedTimeSeriesData(self)
            res = self.TimeSeriesData;
            res = bsxfun(@minus,res,min(res,[],2));
            res = bsxfun(@rdivide,res,max(res,[],2));
        end

        function res = get.ReferenceCurve(self)
            if isempty(self.ReferenceCurve_)
                error('RegionalDyssynchrony:NotComputed',...
                      'Reference Curve not computed yet');
            end
            res = self.ReferenceCurve_;
        end
    end
end
