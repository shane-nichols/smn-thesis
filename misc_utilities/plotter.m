function handle = plotter(X,Y,varargin)
% this program just makes line-plots easier. Documentation is similar to
% the MMplot program, except that this only makes 1 plot not a 4x4 plot array.
% EXAMPLE:
% 
% plotStuff = {...
%     'size',[700,500],...
%     'fontsize',16,...
%     'title','Title of Graph',...
%     'xLabel','X Axis',...
%     'yLabel','Y Axis',...
%     'limy',0.1,...
%     'lineNV',{'lineWidth',2},...
%     'axNV',{'XGrid','on','YGrid','on'}...
%     };
% 
% h = plotter(Lam,MMgetp(MM1,'ld'),'b',plotStuff{:});
% plotter(Lam,MMgetp(MM1,'ldp'),'r',plotStuff{:},'handle',h);
%
% or
%
% h = plotter(Lam,[MMgetp(MM1,'ld') ; MMgetp(MM1,'ldp')],plotStuff{:});

p = inputParser;
% input validation functions
valFun1 = @(x) ischar(x) && ...
    all(~strcmpi(x,...
    {'handle','lineNV','limY','fontsize','axNV','size','title','xLabel','yLabel','legend'}));
valFun2 = @(x) isscalar(x)&&isnumeric(x);
% setup input scheme
addRequired(p,'X',@isnumeric);
addRequired(p,'Y',@isnumeric);
addOptional(p,'LineSpec','-',valFun1)
addParameter(p,'handle',gobjects(1), @ishandle);
addParameter(p,'limY',0,valFun2)
addParameter(p,'fontsize',12,valFun2)
addParameter(p,'axNV',{},@iscell)
addParameter(p,'lineNV',{},@iscell)
addParameter(p,'size',[700 500],@(x) length(x) == 2 && isnumeric(x))
addParameter(p,'title','',@ischar)
addParameter(p,'xLabel','',@ischar)
addParameter(p,'yLabel','',@ischar)
addParameter(p,'legend',{},@iscell)
parse(p,X,Y,varargin{:}) %parse inputs

% create new figure if no valid handles were given
if any(strcmpi('handle',p.UsingDefaults))
    % Determine how large to make the figure window, according to the screensize.
    scrsz = get(0,'screensize');
    figPos = [1 5 p.Results.size];
    if figPos(3) > scrsz(3)
        figPos(3) = scrsz(3);
        warning(['Figure horizontal dimension set to the maximum value of ',...
            num2str(figPos(3)),' pixels.'])
    end
    if figPos(4) > (scrsz(4) - 99)   % 99 pixels is the height of the OSX status bar on my machine
        figPos(4) = (scrsz(4) - 99);
        warning(['Figure vertical dimension set to the maximum value of ',...
            num2str(figPos(4)),' pixels.'])
    end
    h_fig = figure('position',figPos,'units','pixels'); %create figure
    handle = axes;
    hold(handle,'on')
    box(handle,'on')
else
    handle = p.Results.handle;
    h_fig = get(handle,'parent');
    figPos = get(h_fig,'Position');
end
% plot line and set Line Properties
plot(handle,X,Y(:,:),p.Results.LineSpec,p.Results.lineNV{:})
% set Axes properties
axis(handle,'tight'); % first, axes are set to tight
if ~isempty(p.Results.axNV)
    set(handle,p.Results.axNV{:});
end
if p.Results.limY ~= 0 % modify axes bounds if limY is set
    lim = p.Results.limY;
    Ylim = get(handle,'YLim');
    if (Ylim(2) - Ylim(1)) < lim
        avg = (Ylim(2) + Ylim(1))./2;
        Ylim(2) = avg + lim/2;
        Ylim(1) = avg - lim/2;
        set(handle,'Ylim',Ylim);
    end
end
% Adjust plot limits so that lines do not overlap axis borders.
lineHandle = get(handle,'children');
lineWidth = zeros(size(lineHandle));
for j = 1:length(lineHandle)
    if strcmp(get(lineHandle(j),'Marker'),'none')
        lineWidth(j) = get(lineHandle(j),'LineWidth');
    else
        lineWidth(j) = get(lineHandle(j),'MarkerSize');
    end
end
lineWidth = max(lineWidth);
plotPos = get(handle,'Position');
xlim = get(handle,'xLim');
ylim = get(handle,'yLim');
xStep = (xlim(2) - xlim(1))/plotPos(3)/figPos(3)*lineWidth/2;
yStep = (ylim(2) - ylim(1))/plotPos(4)/figPos(3)*lineWidth;
set(handle,'XLim',[xlim(1)-xStep,xlim(2)+xStep]);
set(handle,'YLim',[ylim(1)-yStep,ylim(2)+yStep]);
% add the labels if passed
if ~any(strcmpi('title',p.UsingDefaults))
    title(p.Results.title,'FontSize',p.Results.fontsize,'FontWeight','normal');
end
if ~any(strcmpi('xLabel',p.UsingDefaults))
    xlabel(p.Results.xLabel,'FontSize',p.Results.fontsize);
end
if ~any(strcmpi('yLabel',p.UsingDefaults))
    ylabel(p.Results.yLabel,'FontSize',p.Results.fontsize);
end
% set font size of all graphics objects if fontsize was passed
if ~any(strcmpi('fontsize',p.UsingDefaults))
    set(get(gcf,'children'),'FontSize',p.Results.fontsize);
end
% optionally create legend (this will increase the width of the figure!)
if ~any(strcmpi('legend',p.UsingDefaults))
    Labels = p.Results.legend{2};
    if isnumeric(Labels); Labels = strread(num2str(Labels),'%s'); end
    set(handle,'units','pixels');
    pos = get(handle,'Position');
    lgd = legend(handle,Labels,'location','northeastoutside');
    set(lgd,'units','pixels','fontsize',p.Results.fontsize);
    title(lgd,p.Results.legend{1},'FontSize',p.Results.fontsize);
    lgd_pos = get(lgd,'Position');
    h_fig.Position = h_fig.Position + [0 0 lgd_pos(3) 0];
    set(handle,'Position',pos);
end
end