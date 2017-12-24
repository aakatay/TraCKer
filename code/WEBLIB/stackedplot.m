function [h,hlabels,mode] = stackedplot(varargin)
%STACKEDPLOT   Non-overlapping plot
%   STACKEDPLOT(X,Y) plots vector Y versus vector X like built-in PLOT
%   function. If multiple lines are plotted, they are offset from each
%   other so that they do not overlap. If all line's YData are the same,
%   they are stacked horizontally; otherwise, they are stacked vertically.
%   All other standard PLOT calling conventions are supported with the same
%   stacking default behavior.
%
%   The baselines and rangelines (displayed along y-axis to indicate the
%   span of each plot) are automatically adjusted by event listeners when
%   the limits of x-axis is later changed.
%
%   STACKEDPLOT(AX,...) plots into the axes with handle AX.
%  
%   [HLINES,HLABELS] = STACKEDPLOT(...) returns a column vector of handles
%   HLINES to lineseries objects, one handle per plotted line as well as a
%   column vector of handles HLABELS of text objects for the plot labeling.
%  
%   STACKEDPLOT(..., Option1Name,Option1Value,...) to specify either
%   parameter/value pairs of the plotted lines or option/value pairs of the
%   plotting format.
%
%   Option Name      Description
%   -----------------------------------------------------------------------
%   PlotOrder        [{'topdown'}|'bottomup'|{'leftright'}|'rightleft']
%                    specifies the order of plot appearance. If not
%                    specified, 'topdown' for stacking in y-axis or
%                    'leftright' for stacking in x-axis.
%   PlotLabels       Cell array of strings.
%   PlotScaling      ['normalize'|vector] Pre-scale each plot by this
%                    factor. If set to 'normalize', plots are scaled to
%                    have unit standard deviation. To customize, provide a
%                    vector of scaling factors for plots. Displayed range
%                    (if shown) remains true to the original plot.
%   PlotSpacing      [finite vector {0.25}] space between neighboring plots
%                    specified in the fraction of equally allocated span
%                    for each plot. The value is recycled if not matching
%                    the number of plots. Value may be negative to overlap
%                    the plots.
%   PlotSpans         [{'equal'}|'dynamic'|vector] span of each plot. If
%                    'equal', all plot gets equal space. If 'dynamic',
%                    space fo each plot is scaled according to its range.
%                    If a vector is given, spans are specified to the given
%                    ratio.
%   BaseValue        [{'mean'}|finite vector] Locations of base lines of ,
%                    each plot. 'mean' places the base line through the
%                    middle of each plot. A numeric entry defines a custom
%                    base line location. If there are more plots than the
%                    value specified, STACKEDPLOT recycles the values.
%   PlotRangeMode    [{'all'}|'above'|'below'] specifies how to compute
%                    the range for each plot. If 'all', range is given by
%                    min([data,BaseValue]) to max([data,BaseValue]). If
%                    'above', range is given by BaseValue to
%                    max([data,BaseValue]). If 'below' range is given by
%                    min([data,BaseValue]) to BaseValue.                    
%   MinimumPlotSpans  [scalar {0.36}] used for dynamic PlotRangeMode to keep
%                    the range of minimally varying plot to a reasonable
%                    size. The scalar value defines the fraction of equally
%                    allocated span for each plot.
%   ShowBaseLine     [logical scalar {false}] true to plot the base line 
%                    along the common axis
%   ShowDataRange    [logical scalar {true}] true to plot a vertical line
%                    indicating the data range specified by PlotRangeMode.
%                    RangeLines are only shown if PlotSpacing>0.
%   DataTickMode     [{'none'}|'baseline'|'range'|'all'] Data axis tick
%                    spacing mode. 'baseline' to place a tick at the
%                    baseline, 'range' to place ticks at the extremities of
%                    the plot range, and 'all' for both.
%   LabelLocation    [{'plot'}|'baseline'] Data label placement. 'plot' to
%                    place it next to the first/last sample or 'baseline'
%                    to place it at the baseline location
%   XAxisScale       ['linear'|'log'] sets x-axis scaling of parent axes
%   YAxisScale       ['linear'|'log'] sets y-axis scaling of parent axes
%   XAxisLocation    ['bottom'|'top'] sets x-axis location of parent axes
%   YAxisLocation    ['left'|'right'] sets y-axis locatino of parent axes
%  
%   Example
%      x = (-pi:pi/10:pi)';
%      y = [tan(sin(x)) sin(tan(x))];
%      stackedplot(x,y,'--rs','LineWidth',2,...
%                      'MarkerEdgeColor','k',...
%                      'MarkerFaceColor','g',...
%                      'MarkerSize',10)
%  

% Copyright 2015 Takeshi Ikuma
% History:
% rev. 4: (11-04-2015)
%   * Fixed the error when calling stackedplot(Y,S)
%   * Adjusted behavior when plotting a single line (default to PlotOrder='topdown')
% rev. 3: (11-02-2015)
%   * Removed nagging warning in R2015b by the use of feature('UseHG2')
% rev. 2: (10-25-2015) 
%   * Added XLim & YLim monitoring for automatic placement of baselines and
%     rangelines
%   * Bug fix on PlotOrder option handling
%   * Bug fix on RangeTick display
% rev. 1: (10-23-2015) 
%   * Added 'normalize' option to PlotScaling option
%   * Fixed DataTick display bug when data is scaled
% rev. -: (10-21-2015) original release

narginchk(1,inf);
[plotargs,varargin] = parse_plotargs(varargin);

% do the plot
try
   h = handle(plot(plotargs{:}));
catch ME
   ME.throwAsCaller();
end
if isempty(h)
   return;
end

% if plotted on axes previously used, delete all
axescleanup(handle(h(1).Parent));

% parse the rest of the data
[opts,props] = parse_input(varargin,numel(h),fieldnames(set(h(1))));

% apply the line properties first
try
   set(h,props);
catch ME
   ME.throwAsCaller();
end

% then get started on moving stuff around
mode = stackplots(h,opts);

% format axes
hlabels = formataxes(h,mode,opts);

% prep output
if nargout>0
   if verLessThan('matlab','8.4')
      h = double(h);
   end
else
   clear h
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hlabels = formataxes(h,mode,opts)

ax = handle(h(1).Parent);
N  = numel(h);
color = get(h,{'Color'});

if ~isempty(opts.XAxisScale)
   ax.XScale = opts.XAxisScale;
end
if ~isempty(opts.YAxisScale)
	ax.YScale = opts.YAxisScale;
end
if ~isempty(opts.XAxisLocation)
   ax.XAxisLocation = opts.XAxisLocation;
end
if ~isempty(opts.YAxisLocation)
   ax.YAxisLocation = opts.YAxisLocation;
end

% set axis limit
if mode.data(1)=='X'
   ax.XLim = mode.axlim;
   if ax.XAxisLocation(1)=='b'
      axloc = ax.YLim(1);
   else
      axloc = ax.YLim(2);
   end
else
   ax.YLim = mode.axlim;
   if ax.YAxisLocation(1)=='l'
      axloc = ax.XLim(1);
   else
      axloc = ax.XLim(2);
   end
end

% set options for baseline object
bline_props = {[] [] 'Visible','','Color','k','Tag','StackedPlotBaselines','HandleVisibility','off'};
if mode.data(1)=='X'
   lim = ax.YLim;
   bline_props{1} = repmat(mode.pos(:,2)',2,1);
   bline_props{2} = repmat(lim',1,N);
else
   lim = ax.XLim;
   bline_props{1} = repmat(lim',1,N);
   bline_props{2} = repmat(mode.pos(:,2)',2,1);
end
if opts.ShowBaseLine
   bline_props{4} = 'on';
else
   bline_props{4} = 'off';
end

% set properties for rangeline object
rline_props = {[],[],'Visible','','Tag','StackedPlotRangelines','HandleVisibility','off'};
rline_color = cell(N,1);
if mode.data(1)=='X'
   if strcmp(opts.DataTickMode,'none')
      if verLessThan('matlab','8.4')
         ax.XColor = ax.Color;
      else
         ax.XColor = 'none';
      end
      [rline_color{:}] = deal('k');
   else
      rline_color = color;
   end
   rline_props{1} = mode.pos(:,[1 3])';
   rline_props{2} = repmat(axloc,2,N);
else
   if strcmp(opts.DataTickMode,'none')
      if verLessThan('matlab','8.4')
         ax.YColor = ax.Color;
      else
         ax.YColor = 'none';
      end
      [rline_color{:}] = deal('k');
   else
      rline_color = color;
   end
   rline_props{1} = repmat(axloc,2,N);
   rline_props{2} = mode.pos(:,[1 3])';
end
if opts.ShowDataRange
   rline_props{4} = 'on';
else
   rline_props{4} = 'off';
end

if strcmp(opts.DataTickMode,'none')
   if mode.data(1)=='X'
      ax.XTick = [];
   else
      ax.YTick = [];
   end
else
   tick = mode.pos';
   ticklabel = nan(size(tick));
   if ~strcmp(opts.DataTickMode,'range')
      ticklabel(2,:) = mode.bases;
   end
   if ~strcmp(opts.DataTickMode,'baseline') && all(opts.PlotSpacing>0)
      ticklabel([1 3],:) = mode.ranges';
   end
   
   ticklabel = num2cell(ticklabel);
   for n = 1:N
      vals = [ticklabel{:,n}]';
      idx = isnan(vals);
      [ticklabel{idx,n}] = deal('');
      idx = ~idx;
      idx0 = idx & (vals~=0);
      if any(idx0)
         mag = 10.^ceil(log10(max(abs(vals(idx0)))));
         vals(idx0) = round(vals(idx0)*mag)/mag;
      end
      if mod(vals(idx),1)==0
         fmt = '%d';
      else
         fmt = '%2.1f';
      end
      ticklabel(idx,n) = strtrim(cellstr(num2str(vals(idx),fmt)));
   end
   
   I = [(find(tick(2,:)==tick(1,:))-1)*3+1 (find(tick(2,:)==tick(3,:))-1)*3+3];
   tick(I) = [];
   ticklabel(I) = [];
   
   [tick,I] = sort(tick(:));
   ticklabel = ticklabel(I);
   
   if mode.data(1)=='X'
      ax.XTick = tick;
      ax.XTickLabel = ticklabel;
   else
      ax.YTick = tick;
      ax.YTickLabel = ticklabel;
   end
end

% if PlotLabels not given, default to plot's DisplayName if given
if isnumeric(opts.PlotLabels) && isnan(opts.PlotLabels) % default
   opts.PlotLabels = get(h,{'DisplayName'});
   if all(cellfun(@isempty,opts.PlotLabels))
      opts.PlotLabels = strtrim(cellstr(num2str((1:numel(h))')));
   end
end

if mode.data(1)=='X'
   hlab = handle(ax.XLabel);
   hlab.Units = 'data';
   pos = hlab.Position;
   y = repmat(pos(2),N,1);
   if opts.LabelLocation(1)=='b'%aseline
      x = mode.pos(:,2);
   elseif ax.XAxisLocation(1)=='b'%ottom | 'plot'
      x = mode.endpoints(:,1);
   else% 'top'|'plot'
      x = mode.endpoints(:,2);
   end
else
   hlab = handle(ax.YLabel);
   hlab.Units = 'data';
   pos = hlab.Position;
   x = repmat(pos(1),N,1);
   if opts.LabelLocation(1)=='b'%aseline
      y = mode.pos(:,2);
   elseif ax.XAxisLocation(1)=='b'%ottom | 'plot'
      y = mode.endpoints(:,1);
   else% 'top'|'plot'
      y = mode.endpoints(:,2);
   end
end
halign = hlab.HorizontalAlignment;
valign = hlab.VerticalAlignment;
rot = hlab.Rotation;
delete(hlab)

hlabels = text(x,y,opts.PlotLabels,...
   'HorizontalAlignment',halign,'VerticalAlignment',valign,...
   'Rotation',rot);
set(hlabels,{'Color'},color,'Units','normalized');

hb = handle(line(bline_props{:}));
hr = handle(line(rline_props{:}));
set(hr,{'Color'},rline_color);

% configure axes for event listeners
axessetup(ax,h,hlabels,hb,hr,opts,mode);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mode = stackplots(h,opts)

% get plot data
N = numel(h);
ax = handle(ancestor(h(1),'axes'));
data = cell(N,1);

% determine the plot mode
if isempty(opts.PlotOrder) % auto detect
   % get YData 
   for n = 1:N
      data{n} = h(n).YData;
   end
   
   % if all YData are equal, offset X axis
   if ~isscalar(data) && all(cellfun(@(vals)isequal(vals,data{1}),data(2:end))) % 'leftright'
      reverse = false;
      isX = true; % offset x, eliminate y
   else %'topdown'
      reverse = true;
      isX = false; % offset y, eliminate x
   end
else % manual
   reverse = any(opts.PlotOrder(1)=='tr'); % topdown/rightleft stacks backwards wrt the axes origin
   isX = any(opts.PlotOrder(1)=='lr');
   
   % grab the YData
   if ~isX
      for n = 1:N
         data{n} = h(n).YData;
      end
   end
end

% eliminate as-is data
if isX
   mode.data = 'XData';
   for n = 1:N
      data{n} = h(n).XData;
   end
   axisscale = 'XScale';
else
   mode.data = 'YData';
   axisscale = 'YScale';
end

% set axis scale mode
if ~isempty(opts.XAxisScale)
   ax.XScale = opts.XAxisScale;
end
if ~isempty(opts.YAxisScale)
   ax.YScale = opts.YAxisScale;
end
if ax.(axisscale)(2)=='o' % 'log'
   for n = 1:N
      data{n}(:) = log10(data{n});
   end
end

% set base line values
mode.bases = zeros(N,1);
if ischar(opts.BaseValue) % 'mean'
   for n = 1:N
      mode.bases(n) = mean(data{n});
   end
else
   basevals = repmat(opts.BaseValue(:),ceil(N/numel(opts.BaseValue)),1);
   mode.bases(:) = basevals(1:N);
end

% apply data scaling
datascaled = ischar(opts.PlotScaling) || ~isscalar(unique(opts.PlotScaling));
if datascaled
   if ischar(opts.PlotScaling) % 'normalize'
      scale = cellfun(@(x)std(x),data);
      scale(:) = max(scale)./scale;
   else
      scale = opts.PlotScaling(:);
   end
   % scale both data & baseline
   for n = 1:N
      data{n}(:) = data{n}.*scale(n);
   end
   mode.bases(:) = mode.bases.*scale;
end

% set data range
mode.ranges = zeros(N,2);
if strcmp(opts.PlotRangeMode,'all')
   fcn = @(data,base)[min(min(data),base) max(max(data),base)];
elseif strcmp(opts.PlotRangeMode,'above')
   fcn = @(data,base)[base max(max(data),base)];
else % 'below'
   fcn = @(data,base)[min(min(data),base),base];
end
for n = 1:N
   vals = data{n};
   vals(isnan(vals)|isinf(vals)) = [];
   if isempty(vals)
      mode.ranges(n,:) = [0 0];
   else
      mode.ranges(n,:) = fcn(vals,mode.bases(n));
   end
   data{n}(:) = data{n} - mode.ranges(n,1); % data=0 == low lim
end

% ranges & baselines are translated to a new coordinate system
mode.pos = [zeros(N,1) (mode.bases-mode.ranges(:,1)) diff(mode.ranges,[],2)];

% set plot position
spaces = repmat(opts.PlotSpacing(:),ceil((N-1)/numel(opts.PlotSpacing)),1);
spaces(N:end) = [];
spaces(:) = spaces*mean(mode.pos(:,3));
if reverse
   pos0 = flipud(cumsum([0;mode.pos(end:-1:2,3)+spaces]));
else
   pos0 = cumsum([0;mode.pos(1:end-1,3)+spaces]);
end
mode.pos(:) = bsxfun(@plus,mode.pos,pos0);
mode.axlim = [min(mode.pos(:,1)) max(mode.pos(:,3))];

% adjust the spans if PlotSpans = 'equal' or 'dynamic'
span = diff(mode.pos(:,[1 3]),[],2);
zerospan = span==0;
if ~all(zerospan) % if all zero, nothing to adjust
   if strncmp(opts.PlotSpans,'d',1) %ynamic
      spanavg = mean(span);
      spanmin = spanavg*opts.MinimumPlotSpans;
      toothin = span<spanmin;
      while any(toothin)
         % fatten the too thin ones
         for n = find(toothin)'
            dpos = spanavg - span(n);
            if span(n)==0
               delta = repmat(dpos/2,1,2);
            else
               delta = diff(mode.pos(n,:))/span(n)*dpos;
            end
            if reverse
               % raise the plots above
               mode.pos(1:n-1,:) = mode.pos(1:n-1,:) + delta(2);
               % lower the plots below
               mode.pos(n+1:end,:) = mode.pos(n+1:end,:) - delta(1);
            else
               % raise the plots above
               mode.pos(n+1:end,:) = mode.pos(n+1:end,:) + delta(2);
               % lower the plots below
               mode.pos(1:n-1,:) = mode.pos(1:n-1,:) - delta(1);
            end
            
            mode.axlim(1) = mode.axlim(1) - delta(1);
            mode.axlim(2) = mode.axlim(2) + delta(2);
            span(n) = span(n) + dpos;
         end
         
         spanavg = mean(span);
         spanmin = spanavg*opts.MinimumPlotSpans;
         toothin = span<spanmin;
      end
   else
      if strncmp(opts.PlotSpans,'e',1)%qual -> need to further scale data
         dpos = max(span)-span;
         %delta = bsxfun(@times,diff(mode.pos,[],2),dpos./span);
         %delta(zerospan,:) = dpos(zerospan)/2;
         delta = repmat(dpos/2,1,2);
         
      else %if isnumeric(opts.PlotSpans) % manual plot spans given
         
         spanscale = repmat(opts.PlotSpans(:),ceil(N/numel(opts.PlotSpans)),1);
         spanscale(N+1:end) = [];
         spanscale(:) = spanscale/min(spanscale);
         
         % if zero-span, set it to the mean of non-zero spans
         meanspan = mean(span(~zerospan));
         
         % if plot range is zero (i.e., constant value) set span according to the
         % non-zero-range plots
         dpos = span.*spanscale - span;
         delta = bsxfun(@times,diff(mode.pos,[],2),dpos./span);
         delta(zerospan,:) = meanspan*spanscale(zerospan)/2;
      end
         
      posdelta = zeros(N,1);
      if reverse
         % raise the plots above
         posdelta(1:end-1) = flipud(cumsum(delta(end:-1:2,2)));
         % lower the plots below
         posdelta(2:end) = posdelta(2:end) - cumsum(delta(1:end-1,1));
      else
         % raise the plots above
         posdelta(2:end) = cumsum(delta(1:end-1,2));
         % lower the plots below
         posdelta(1:end-1) = posdelta(1:end-1) - flipud(cumsum(delta(end:-1:2,1)));
      end
      mode.pos(:) = bsxfun(@plus,mode.pos,posdelta);
      mode.axlim(1) = mode.axlim(1)-sum(delta(:,1));
      mode.axlim(2) = mode.axlim(2)+sum(delta(:,2));
   end
end

% offset data to the final position
for n = 1:N
   data{n} = data{n} + mode.pos(n,1);
end

if datascaled
   mode.ranges(:) = bsxfun(@rdivide,mode.ranges,scale);
   mode.bases(:) = bsxfun(@rdivide,mode.bases,scale);
end

% if in log, revert to linear scale
if ax.(axisscale)(2)=='o' % 'log'
   val = min(mode.pos(:));
   if val<1
      val = 1-val;
   else
      val = 0;
   end
   mode.pos(:) = 10.^(mode.pos+val);
   for n = 1:N
      data{n}(:) = 10.^(data{n}+val);
   end
   mode.ranges(:) = 10.^(mode.ranges);
   mode.bases(:) = 10.^(mode.bases);
   mode.axlim = 10.^(mode.axlim+val);
end

% relocate the plots
mode.endpoints = zeros(N,2);
for n = 1:N
   h(n).(mode.data) = data{n};
   mode.endpoints(n,:) = data{n}([1 end]);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [opts,props] = parse_input(argin,N,pnames)

opts_list = {
   'PlotOrder'       ''        @ischar {'','topdown','bottomup','leftright','rightleft'}
   'PlotLabels'      nan       @(v)(iscellstr(v) && numel(v)==N) ''
	'PlotScaling'     ones(1,N) @(v)validatePlotScaling(v,N) {'normalize'}
   'PlotSpacing'     0.25      @(v)validateattributes(v,{'numeric'},{'vector','finite'}) ''
   'PlotSpans'        'dynamic' @validatePlotSpans {'equal','dynamic'}
   'BaseValue'       'mean'    @validatebasevalue {'mean'}
   'PlotRangeMode'   'all'     @ischar {'all','above','below'}
   'MinimumPlotSpans' 0.36      @(v)validateattributes(v,{'numeric'},{'scalar','positive','finite'}) ''
   'ShowBaseLine'    true      @(v)validateattributes(v,{'logical'},{'scalar'}) ''
   'ShowDataRange'   true      @(v)validateattributes(v,{'logical'},{'scalar'}) ''
	'DataTickMode'    'none'    @ischar {'none','baseline','range','all'}
	'LabelLocation'   'plot'    @ischar {'plot','baseline'}
   'XAxisScale'      ''        @ischar {'','linear','log'}
	'YAxisScale'      ''        @ischar {'','linear','log'}
   'XAxisLocation'   ''        @ischar {'','bottom','top'}
   'YAxisLocation'   ''        @ischar {'','left','right'}
};

p = inputParser;
for n = 1:size(opts_list,1)
   p.addParameter(opts_list{n,1:3});
end
for n = 1:numel(pnames)
   try % if same as an option name, skip
      p.addParameter(pnames{n},[]);
   catch
   end
end
p.parse(argin{:});

rnames = fieldnames(p.Results);
props = struct([]);
opts = struct([]);
for n = 1:numel(rnames)
   tf = strcmp(rnames{n},opts_list(:,1));
   if any(tf)
      opts(1).(rnames{n}) = p.Results.(rnames{n});
      if ischar(opts.(rnames{n})) && ~isempty(opts_list{tf,4})
         opts.(rnames{n}) = validatestring(opts.(rnames{n}),opts_list{tf,4});
      end
   elseif ~isempty(p.Results.(rnames{n}))
      props(1).(rnames{n}) = p.Results.(rnames{n});
   end
end

end

function validatePlotScaling(v,N)
if ~ischar(v)
   validateattributes(v,{'numeric'},{'vector','numel',N,'positive','finite'})
end
end

function validatePlotSpans(v)
if ~ischar(v)
   validateattributes(v,{'numeric'},{'vector','positive','finite'});
end
end

function validatebasevalue(v)
if ~ischar(v)
   validateattributes(v,{'numeric'},{'vector','finite'});
end
end

function [plotargs,argin] = parse_plotargs(argin)
%   STACKEDPLOT(X,Y) plots vector Y versus vector X like built-in PLOT
%   STACKEDPLOT(AX,...) plots into the axes with handle AX.
%   STACKEDPLOT(Y)
%   STACKEDPLOT(X1,Y1,S1,X2,Y2,S2,X3,Y3,S3,...)

if isscalar(argin{1}) && ishghandle(argin{1},'axes')
   idx = 1;
else
   idx = 0;
end

% first argument must be numeric
tf = cellfun(@isnumeric,argin(idx+1:end));
if ~tf(1)
   error('The first argument of STACKEDPLOT must be numeric data.');
end

% find the first non-numeric input
Ich = find(~tf,1);
if isempty(Ich) % that's all we got
   plotargs = argin;
   argin(:) = [];
else
   % account for all numeric
   idx = idx + Ich-1;
   tf(1:Ich-1) = [];
   
   % check the first non-numeric argument
   while isplotspec(idx+1,argin) % if plotspec
      % find the next non-numeric argument
      idx = idx + 1;
      tf(1) = [];

      % get next non-numeric input
      Ich = find(~tf,1);
      if isempty(Ich) % only plot arguments
         break;
      else
         idx = idx + Ich-1;
         tf(1:Ich-1) = [];
      end
   end
   plotargs = argin(1:idx);
   argin(1:idx) = [];
end
end

function tf = isplotspec(idx,argin)

Schoice = '*+-.:<>^bcdghkmoprsvwxy'; % possible characters for S argument

if isempty(idx) || idx>numel(argin)
   tf = false;
else
   S = argin{idx};
   tf = ischar(S) && isrow(S) && all(ismember(S,Schoice));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axes property listeners to keep the plots right

function axessetup(ax,h,hl,hb,hr,opts,mode)
% call this to set up all the listeners and necessary appdata

% store data
setappdata(ax,'StackedPlotMode',mode);
setappdata(ax,'StackedPlotOpts',opts);
setappdata(ax,'StackedPlotHandles',h);
setappdata(ax,'StackedPlotBaseLines',hb);
setappdata(ax,'StackedPlotRangeLines',hr);
setappdata(ax,'StackedPlotLabels',hl);

% store current state
pnames = {...
   'XLim','XDir','XScale',...
   'YLim','YDir','YScale'};
states = cell2struct(get(ax,pnames),pnames,2);
setappdata(ax,'StackedPlotAxesStates',states);

try
   setappdata(ax,'StackedPlotMarkedClean',...
      addlistener(ax,'MarkedClean',@axesmarkedclean));
catch
   pfcns = {
      @xlimchanged
      @ylimchanged
      @xdirchanged
      @ydirchanged
      @xscalechanged
      @yscalechanged};
   
   % set up new listeners
   for n = 1:numel(pnames)
      setappdata(ax,['StackedPlotPostSet' pnames{n}],...
         addlistener(ax,'XLim','PostSet',pfcns{n}));
   end
end
end

function axescleanup(ax)
% called to explicitly remove all stackedplot objects

if isappdata(ax,'StackedPlotOpts')
   
   if verLessThan('matlab','8.4')
      pnames = {...
         'XLim','XDir','XScale',...
         'YLim','YDir','YScale'};
      
      for n = 1:numel(pnames)
         adname = ['StackedPlotPostSet' pnames{n}];
         delete(getappdata(ax,adname));
         rmappdata(ax,adname);
      end
   else
      delete(getappdata(ax,'StackedPlotMarkedClean'));
      rmappdata(ax,'StackedPlotMarkedClean');
   end
   
   try
      delete(getappdata(ax,'StackedPlotHandles'));
   catch
   end
   try
      delete(getappdata(ax,'StackedPlotBaseLines'));
   catch
   end
   try
      delete(getappdata(ax,'StackedPlotRangeLines'));
   catch
   end
   try
      delete(getappdata(ax,'StackedPlotLabels'));
   catch
   end
   
   rmappdata(ax,'StackedPlotMode');
   rmappdata(ax,'StackedPlotOpts');
   rmappdata(ax,'StackedPlotAxesStates');
   rmappdata(ax,'StackedPlotHandles');
   rmappdata(ax,'StackedPlotBaseLines');
   rmappdata(ax,'StackedPlotRangeLines');
   rmappdata(ax,'StackedPlotLabels');
end

end

function setlisena(lis,ena)
try
   lis.Enabled = ena;
catch
   if ena
      set(lis,'Enabled','on');
   else
      set(lis,'Enabled','off');
   end
end
end

function axesmarkedclean(ax,evt)
try
   ax = evt.AffectedObject;
catch
end

lis = getappdata(ax,'StackedPlotMarkedClean');
setlisena(lis,false);
C = onCleanup(@()setlisena(lis,true));

states = getappdata(ax,'StackedPlotAxesStates');
changed = false;
if ~isequal(states.XLim,ax.XLim)
   xlimchanged(ax,evt);
   changed = true;
   states.XLim = ax.XLim;
end
if ~isequal(states.YLim,ax.YLim)
   ylimchanged(ax,evt);
   changed = true;
   states.YLim = ax.YLim;
end

if changed
   setappdata(ax,'StackedPlotAxesStates',states);
end

end

function xlimchanged(ax,evt)
% axes' XLim property changed

try
   ax = evt.AffectedObject;
catch
end

lis = getappdata(ax,'StackedPlotPostSetXLim');
setlisena(lis,false);
C = onCleanup(@()setlisena(lis,true));
mode = getappdata(ax,'StackedPlotMode');
if any(mode.data(1)=='Y') % topdown or bottomup
   xdata = ax.XLim;
   hb = getappdata(ax,'StackedPlotBaseLines');
   hr = getappdata(ax,'StackedPlotRangeLines');
   if all(ishghandle([hb;hr]))
      for n = 1:numel(hb)
         hb(n).XData = xdata;
         hr(n).XData(:) = xdata(1);
      end
   end
end

end

function ylimchanged(ax,evt)
% axes' YLim property changed

try
   ax = evt.AffectedObject;
catch
end

lis = getappdata(ax,'StackedPlotPostSetYLim');
setlisena(lis,false);
C = onCleanup(@()setlisena(lis,true));
mode = getappdata(ax,'StackedPlotMode');
if any(mode.data(1)=='X') % leftright or rightleft
   ydata = ax.YLim;
   hb = getappdata(ax,'StackedPlotBaseLines');
   hr = getappdata(ax,'StackedPlotRangeLines');
   if all(ishghandle([hb;hr]))
      for n = 1:numel(hb)
         hb(n).YData = ydata;
         hr(n).YData(:) = ydata(1);
      end
   end
end

end

% To be supported
function xscalechanged(h,evt)
% axes' XScale property changed

% to-be-supported

end

function yscalechanged(h,evt)
% axes' YScale property changed

% to-be-supported

end

% TO be supported
function xdirchanged(h,evt)
% axes' XScale property changed

% to-be-supported

end

function ydirchanged(h,evt)
% axes' YScale property changed

% to-be-supported

end
