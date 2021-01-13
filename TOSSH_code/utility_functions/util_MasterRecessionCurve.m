function [MRC, fig_handles] = util_MasterRecessionCurve(Q, flow_section, varargin)
%util_MasterRecessionCurve fits a master recession curve to recession segments.
%
%   INPUT
%   Q: streamflow [mm/timestep]
%   flow_section: n-by-2 array where n is the number of recession segments;
%       columns are the indices into the flow array of the start and end of
%       the recession segments
%   OPTIONAL
%   fit_method: 'exponential' (approximates each recession segment as an
%       exponential before stacking into MRC), 'nonparameteric' (fits 
%       horizontal time shift for minimum standard deviation at each lag
%       time, does not assume any form of the curve)
%   match_method: how to space points on the MRC used for alignment,
%       'linear' or 'log'
%   plot_results: whether to plot results, default = false
%
%   OUTPUT
%   MRC: two-column array of time and flow, specifying the MRC
%   fig_handles: figure handles to manipulate figures (empty if plotting is
%       not requested)
%
%   EXAMPLE
%   % load example data
%   data = load('example/example_data/33029_daily.mat');
%   Q = data.Q;
%   t = data.t;
%   flow_section = util_RecessionSegments(Q,t); % get recession segments
%   [mrc] = util_MasterRecessionCurve(Q, flow_section); % get MRC
%
%   References
%   Posavec, K., Bacani, A. and Nakic, Z., 2006. A visual basic spreadsheet
%   macro for recession curve analysis. Groundwater, 44(5), pp.764-767.
%
%   Copyright (C) 2020
%   This software is distributed under the GNU Public License Version 3.
%   See <https://www.gnu.org/licenses/gpl-3.0.en.html> for details.

if nargin < 2
    error('Not enough input arguments.')
end

ip = inputParser;
ip.CaseSensitive = true;

% required input arguments
% time series have to be numeric and either a (n,1) or a (1,n) vector
addRequired(ip, 'Q', @(Q) isnumeric(Q) && (size(Q,1)==1 || size(Q,2)==1))
addRequired(ip, 'flow_section', @(flow_section) isnumeric(flow_section) && (size(flow_section,2)==2))

addParameter(ip, 'fit_method', 'exponential', @ischar) % defines method for aligning flow segments
addParameter(ip, 'match_method', 'log', @ischar) % defines method for aligning flow segments
addParameter(ip, 'plot_results', false, @islogical) % whether to show plot of MRC

parse(ip, Q, flow_section, varargin{:})
fit_method = ip.Results.fit_method;
match_method = ip.Results.match_method;
plot_results = ip.Results.plot_results;

% create empty figure handle
fig_handles = [];

% calculate the master recession curve (MRC)
% MRC is an array [time after recession start, flow]

switch fit_method
    
    case 'exponential'
        % sort the flow sections according to highest starting value
        start_values = Q(flow_section(:,1).');
        start_values = sortrows([(1:length(start_values)).',start_values],-2);
        
        % start the MRC with the highest segment
        MRC = [(1:(flow_section(start_values(1,1),2)-flow_section(start_values(1,1),1)+1)).',...
            (Q(flow_section(start_values(1,1),1):flow_section(start_values(1,1),2)))];
        
        % loop adding segment to MRC each time
        for i = 2:size(flow_section,1)
            % fit an exponential to the mrc so far lny=ax+b
            % Intercept is mdl(1), slope is mdl(2)
            mdl = [(MRC(:,1).^0) (MRC(:,1))]\log(MRC(:,2));
            % calculate the time shift required to place the initial point
            % of the next recession segment on the first regression curve
            timeshift = ((log(start_values(i,2))-mdl(1))/mdl(2));
            % add the shifted segment to the master recession
            MRC = [MRC; ...
                [timeshift+(1:(flow_section(start_values(i,1),2)-flow_section(start_values(i,1),1)+1)).',...
                (Q(flow_section(start_values(i,1),1):flow_section(start_values(i,1),2)))]];
        end
        
    case 'nonparametric_analytic'
        % download all the flow segments, add jitter to avoid long constant
        % flow values that can't be interpolated, sort values to avoid
        % cases where flow is not decreasing, find min and max flow
        
        % constants
        jitter_size = 1e-8;
        % number of interpolated points in the MRC
        numflows = 500;
        % reset random number seed for reproducibility
        rng('default')
        
        % get number of flow segments
        numsegments = size(flow_section,1);
        
        % order flow sections starting with the largest initial flow value
        flow_init_value = Q(flow_section(:,1));
        [~,sortind]=sort(flow_init_value,'descend');
        % keep running tally of minimum
        running_min = max(flow_init_value);
        
        % create cell array of recession segments, starting with highest flow
        % add jitter to everything except the first value of each segment
        segments = [];
        for i = 1:size(flow_section,1)
            % retrieve the segment
            segment = Q(flow_section(sortind(i),1):flow_section(sortind(i),2));
            % add jitter
            segment(2:end) = segment(2:end) + ...
                normrnd(0,jitter_size,size(segment,1)-1,1);
            % avoid negative segment values
            segment = abs(segment)+1e-20;
            % sort the segment with jitter, in case eps parameter was used
            % and so thereare small increases during the recessions
            segment = sort(segment,'descend');
            % store in cell array
            segments{i} = segment.';
        end
        
        % get flow values where curves should be matched
        max_flow = max([segments{:}]);
        min_flow = min([segments{:}]);
        if min_flow <=0
            min_flow = jitter_size;
        end
        % get interpolated flow values where MRC will be evaluated
        switch match_method
            case 'linear'
                flow_vals = linspace(max_flow,min_flow,numflows);
            case 'log'
                frac_log = 0.2;
                gridspace = (max_flow - min_flow)/numflows;
                flow_vals = sort([linspace(max_flow-gridspace/2,min_flow+gridspace/2,numflows-floor(frac_log.*numflows)),...
                    logspace(log10(max_flow),log10(min_flow),floor(frac_log.*numflows))],'descend');
                flow_vals(end) = min_flow;
                flow_vals(1) = max_flow;
                flow_vals = sort(unique(flow_vals),'descend');
                numflows = numel(flow_vals);
            otherwise
                error('Match method for MRC not a recognised option.')
        end
        
        %Keep track of good segments
        short_segs = false(size(flow_section,1),1);
        
        % extract and interpolate each segment, and check validity; remove
        % invalid segments
        for i = 1:numsegments
                       
            % extract segment
            segment = segments{i};

            % find indices of max and min interpolated flow values for this segment
            fmax_index = find(segment(1) >= flow_vals,1,'first');
            if segment(end) <= flow_vals(end)
                fmin_index = numel(flow_vals);
            else
                fmin_index = find(segment(end) > flow_vals,1,'first')-1;
            end
                       
            % find number of interpolated values
            nf = fmin_index-fmax_index+1;
            
            % if no interpolated values (occurs when min and max of segment
            % are too close together, remove the segment
            if nf <= 1
                %Collect segment number
                short_segs(i) = true;

            end
                 
        end
        
        %If some segments were rejected, recalculate flow values for
        %interpolation and flow value initialisations and counts
        if sum(short_segs) > 0
            %Remove segments without interpolated values
            segments(short_segs) = [];
            numsegments = numel(segments);
            % keep running tally of minimum
            running_min = max(flow_init_value(sortind(~short_segs)));
            %Remove flow vals for interpolation if the reduced 'good' segment
            %set no longer cover those values
            max_flow = max([segments{:}]);
            min_flow = min([segments{:}]);
            flow_vals(flow_vals > max_flow)=[];
            flow_vals(flow_vals < min_flow)=[];
            numflows = numel(flow_vals);
        end
        
        % set up the optimisation matrix
        msp_matrix = zeros(numsegments*numflows*2,3);
        b_matrix = zeros(numsegments*numflows,1);
        % initialise count into that matrix
        mcount = 1;
        % initialise count into sparse matrix
        mspcount = 1;
        %Keep track of any segments with no interpolated values
        bad_segs = [];
        
        
        % extract and interpolate each segment
        for i = 1:numsegments
                       
            % extract segment
            segment = segments{i};
            % if there is a gap between previous segments and this one,
            % then interpolate with a vertical line
            if segment(1) < running_min
                segment = [running_min , segment];
            end
            % find indices of max and min interpolated flow values for this segment
            fmax_index = find(segment(1) >= flow_vals,1,'first');
            if segment(end) <= flow_vals(end)
                fmin_index = numel(flow_vals);
            else
                fmin_index = find(segment(end) > flow_vals,1,'first')-1;
            end
            % interpolate each segment onto the flow values
            interp_segment = interp1(segment,1:numel(segment),flow_vals(fmax_index:fmin_index));
            % keep running tally of minimum
            running_min = min(running_min,flow_vals(fmin_index));
            
            % find number of interpolated values
            nf = fmin_index-fmax_index+1;
            
            % if no interpolated values (occurs when min and max of segment
            % are too close together 
            if nf == 0
                %Collect segment number
                bad_segs = [bad_segs, i];
                %Don't add to minimisation matrix
                continue
            end
            
            % construct the minimisation matrix block for each segment            
            if i==1
                % lag of the first segment is set to zero
                msp_matrix(mspcount:mspcount+nf-1,:)=[[mspcount:mspcount+nf-1].',...
                    [numsegments+fmax_index-1:numsegments+fmin_index-1].',-ones(nf,1)];
                b_matrix(mcount:mcount+nf-1)=interp_segment(:);
            else
                % lags of other segments can be minimised, along with the
                % fitted MRC
                msp_matrix(mspcount:mspcount+2*nf-1,:) = [[mcount:mcount+nf-1,mcount:mcount+nf-1].',...
                    [(i-1)*ones(nf,1);...
                    [numsegments+fmax_index-1:numsegments+fmin_index-1].'],[ones(nf,1);-ones(nf,1)]];
                b_matrix(mcount:mcount+nf-1)=interp_segment(:);
            end
            % update count of rows in the optimisation matrix
            mcount = mcount + nf;
            if i==1
                mspcount = mspcount + nf;
            else
                mspcount = mspcount + 2*nf;
            end
        end
        
        % create sparse matrix
        msp_matrix = msp_matrix(1:mspcount-1,:);
        m_sparse = sparse(msp_matrix(:,1),msp_matrix(:,2),msp_matrix(:,3),mcount-1,numsegments-1+numflows);
        
        % cut off unused rows of optimisation matrix
        B_mat = -b_matrix(1:mcount-1);
        
        %Delete unused columns of minimimsation matrix
        seg_indices = [1:numsegments];
        m_sparse(:,bad_segs-1)=[];
        seg_indices(bad_segs)=[];
        numsegments = numsegments - length(bad_segs);
        
        % minimise the differences to a Master recession curve
        MRC_solve = m_sparse\B_mat;
        
        % extract the time lags and flow values for the MRC
        lags = [0; MRC_solve(1:numsegments-1)];
        mrc_time = MRC_solve(numsegments:end);
        % sort the MRC to avoid any places where not strictly decreasing
        mrc_time = sort(mrc_time,'ascend');
        %Have the MRC start at 0 time
        offset = min(mrc_time);
        mrc_time = mrc_time - offset;
        lags = lags - offset;
        
        % output
        MRC = [mrc_time(:),flow_vals(:)];
        
        % optional plotting
        if plot_results
            fig = figure('Position',[100 100 350 300]); hold on
            for i = 1:numsegments
                % extract segment
                segment = segments{seg_indices(i)};
                h1 = plot([1:length(segment)]+lags(i),segment,'b-');
            end
            h2 = plot(mrc_time,flow_vals,'g','linewidth',2);
            xlabel('Relative time')
            ylabel('Flow')
            legend([h1 h2],{'Recession Segments','Fitted MRC'})
            title('Nonparametric MRC fit')
            fig_handles.MRC_nonparametric = fig;
        end
        
    otherwise
        error('Fit method for MRC not a recognised option.')
end

end

function f = mrc_nonparameteric(offsets, segments, S)

[segments_aligned_locs,~] = offset_matrix(segments, S, offsets);

% sum the standard deviations of each column, ignoring zero values
segments_aligned_locs(segments_aligned_locs==0)=nan;
sd = var(segments_aligned_locs,1,'omitnan');
numvar = sum(segments_aligned_locs > 0);
sd = sd.*numvar;
sd(isnan(sd))=0;

f = sum(sd);

end

function [segments_aligned_locs,locations] = offset_matrix(segments, S, offsets)

offsets = [0, offsets];
num_segments = size(segments,2);
cellsz = cellfun(@length,segments,'uni',false);
max_length = ceil(max([(cellsz{:})]+offsets));

% line up all the segments, with offsets into an array
if ~(num_segments==floor(num_segments) && max_length==floor(max_length))
    error('non integer dimensions')
end

% offset each segment by the prescribed offset and add these locations to a
% master locations array
locations = [];
for i = 1:num_segments
    segment_i = segments{i};
    segment_loc = offsets(i)+1:offsets(i)+length(segment_i);
    locations = [locations, segment_loc];
end

% sort all the locations from all the segments into once increasing array
locations = sort(locations);
% set up an array where the segment values can be stored against their locations
segments_aligned_locs = zeros(num_segments,length(locations));

% for each segment
for i = 1:num_segments
    % retrieve segment values from master array
    segment_i = segments{i};
    % recalculate the segment locations
    segment_loc = 1:length(segment_i);
    locations_i = locations - offsets(i);
    
    % find the start and end point of the segment within the locations array
    locations_start = find(segment_loc(1)<=locations_i,1,'first');
    locations_end = find(segment_loc(end)>=locations_i,1,'last');
    
    % interpolate the segment onto all the intervening locations
    Si=S{i};
    segment_i_interp=Si(locations_i(locations_start:locations_end).');
    
    % write the interpolated values into the master array
    segments_aligned_locs(i,locations_start:locations_end) = segment_i_interp;
    
end

end

