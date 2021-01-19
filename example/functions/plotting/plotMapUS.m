function [] = plotMapUS(lat,lon,z,varargin)
%plotMapUS Plots US map with dots coloured according to an attribute.
%   Options:
%       - various plotting options, e.g. axes limits
%       - save plot as PDF
%
%	INPUT
%   lat: latitude
%   lon: longitude
%   z: attribute to be coloured in, e.g. BFI
%   attribute_name: name of attribute
%   ID: catchment ID
%   colour_scheme: name of colour scheme
%   flip_colour_scheme: flip colour scheme?
%   c_limits: limits of colour axis, e.g. [0 1]
%   c_lower_limit_open: is the lower limit open?
%   c_upper_limit_open: is the upper limit open?
%   nr_colours: nr of colours used for colourscale
%   figure_title: title of plot, e.g. '(a)'
%   figure_name: name for saving, e.g. CAMELS_BFI
%   save_figure: save plot true/false
%   figure_path: path to folder where figure should be saved
%   figure_type: figure type, e.g. -dpdf or -dmeta
%
%   OUTPUT
%   plot and saved figure
%
%   ---
%
%   Sebastian Gnann, sebastian.gnann@bristol.ac.uk (2020)

if nargin < 3
    error('Not enough input arguments.')
end

ip = inputParser;

addRequired(ip, 'latitude', ...
    @(lat) isnumeric(lat) && (size(lat,1)==1 || size(lat,2)==1))
addRequired(ip, 'longitude', ...
    @(lon) isnumeric(lon) && (size(lon,1)==1 || size(lon,2)==1))
addRequired(ip, 'attribute', ...
    @(z) isnumeric(z) || islogical(z))

addParameter(ip, 'attribute_name', char(), @ischar)
addParameter(ip, 'ID', NaN(size(z)), @isnumeric)
addParameter(ip, 'colour_scheme', 'parula', @ischar)
addParameter(ip, 'flip_colour_scheme', false, @islogical)
addParameter(ip, 'c_limits', [min(z) max(z)], @(x) isnumeric(x) && length(x)==2)
addParameter(ip, 'c_lower_limit_open', false, @islogical)
addParameter(ip, 'c_upper_limit_open', false, @islogical)
addParameter(ip, 'nr_colours', 10, @isnumeric)
addParameter(ip, 'figure_title', '', @ischar)
addParameter(ip, 'figure_name', 'no_name', @ischar)
addParameter(ip, 'save_figure', false, @islogical)
addParameter(ip, 'figure_path', '', @ischar)
addParameter(ip, 'figure_type', '-dpdf', @ischar)

parse(ip, lat, lon, z, varargin{:})

attribute_name = ip.Results.attribute_name;
ID = ip.Results.ID;
colour_scheme = ip.Results.colour_scheme;
flip_colour_scheme = ip.Results.flip_colour_scheme;
c_limits = ip.Results.c_limits;
c_lower_limit_open = ip.Results.c_lower_limit_open;
c_upper_limit_open = ip.Results.c_upper_limit_open;
nr_colours = ip.Results.nr_colours;
figure_title = ip.Results.figure_title;
figure_name = ip.Results.figure_name;
save_figure = ip.Results.save_figure;
figure_path = ip.Results.figure_path;
figure_type = ip.Results.figure_type;

%% plotting
index = [1:length(ID)]';

fig = figure('Name',figure_name,'NumberTitle','off','pos',[100 100 600 300]);
ax = axesm('MapProjection','mercator','MapLatLimit',[24 52],'MapLonLimit',[-130 -65]);
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax, states, ...
    'DisplayType','polygon','DefaultFaceColor','white','DefaultEdgeColor','black')
hold on
% grid on

% create colormap
% if flip_colour_scheme
%     colour_mat = flip(brewermap(nr_colours,colour_scheme));
% else
%     colour_mat = (brewermap(nr_colours,colour_scheme));
% end

% plot
scatterm(lat(isnan(z)),lon(isnan(z)),'x k','linewidth',1.2);
scatterm(lat,lon,25,z,'filled')
% xlabel('Latitude [km]'); ylabel('Longitude [km]')
set(gca,'Visible','off')
title(figure_title,'Visible','on')
axis equal
% colormap(colour_mat)
colormap(colour_scheme);
if flip_colour_scheme
    cmap = colormap;
    colormap(flipud(cmap));
end
c=colorbar;
title(c,attribute_name)
x1=get(gca,'position');
x=[0.8 0.2 0.01 0.25];
set(c,'Position',x)
set(gca,'position',x1)
caxis(c_limits)
if c_lower_limit_open
    c.TickLabels{1} = ['<' c.TickLabels{1}];
end
if c_upper_limit_open
    c.TickLabels{end} = ['>' c.TickLabels{end}];
end

% update cursor
dcm_obj = datacursormode(figure(fig));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,index})

%% save fig
if save_figure
    fig_name = strcat('map_US','_',figure_name);
    saveFig(fig,fig_name,figure_path,figure_type)
end

end

