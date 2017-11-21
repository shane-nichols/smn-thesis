classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MPlot < handle
    
    properties
        
        size = [1000 700]
        fontsize = 12
        limy = 1e-4
        title = ''
        legend = {}
        axNV = {}
        lineNV = {}
        ev = false
        vSpace = 0
        borderFactor = 0
        
    end
    
    properties (SetAccess = protected)
        figHandle
        axesHandles = gobjects(4);
    end
    
    methods
        
        function obj = MPlot(varargin)
            obj.figHandle = figure;
            makePlots(obj,varargin{:});
        end
        
        function add(obj, wavelengths, mmData, varargin)
            [mmError, lineSpec] = parseVariableInputs(obj, varargin{:});
            plot(obj, wavelengths, mmData, mmError, lineSpec);
        end
        
        function print(obj,filepath)
            print(obj.figHandle, filepath, '-depsc');
        end
        
    end
    
    methods (Hidden)
        
        function [mmError, lineSpec] = parseVariableInputs(obj, varargin)

            % Optional positional inputs:
            %   mmError:  array the same size as MMdata representing error bar values;
            %   lineSpec: string containing a valid lineSpec. Type "doc LineSpec" in
            %       command window for more info. Default is "-", a solid line.
            % Optional Name-Value pairs inputs:
            %   ev: bool. converts X axis from nanometers to eV.
            %   limY: scalar numeric. limits how small the range of the y-axes can be.
            %   fontsize: sets font-size. Default is 12 pts. Changing the fontsize
            %       of existing plots is not recommended. (Set on first call).
            %   lineNV: a 1D cell array containing Name-Value pair arguments valid for
            %       Chart Line Properties.
            %   axNV: a 1D cell array containing Name-Value pairs arguments valid for
            %       Axes Properties.
            %   size: Size of the figure in pixels given as a two element vector [X Y].
            %       A warning is issued if the requested size is larger than the screen
            %       size minus the height of the OSX status bar (on my machine).
            %       Default size is [1000 700].
            %   title: string containing a title to place at the top of the figure.
            %   legend: two-element cell array. First element is a string to use for
            %       title of the legend. Second element is either a numeric array
            %       containing values to use for labels of each plot, or a cell array
            %       of strings to use as labels. Only set legend on last call, or just
            %       write all plots at once (better).
            %   vSpace: Adds extra space vertical between plots, in pixels
            %   borderFactor: Increases white space around plots. This value is a
            %       multiple of the largest line width on the plots.
    
            p = inputParser;
            p.CaseSensitive = false;
            % input validation functions
            valFun1 = @(x) ischar(x) && ...
                all(~strcmpi(x,{'ev','lineNV','limy','fontsize','axNV','size',...
                'title','legend','vSpace','borderFactor'}));
            valFun2 = @(x) isscalar(x)&&isnumeric(x);
            valFun3 = @(x) isnumeric(x);
            % setup input scheme
            addOptional(p,'mmError',[],valFun3)
            addOptional(p,'lineSpec','-',valFun1)
            addParameter(p,'ev',obj.ev,@islogical)
            addParameter(p,'limy',obj.limy,valFun2)
            addParameter(p,'fontsize',obj.fontsize,valFun2)
            addParameter(p,'axNV',obj.axNV,@iscell)
            addParameter(p,'lineNV',obj.lineNV,@iscell)
            addParameter(p,'size',obj.size,@(x) length(x) == 2 && isnumeric(x))
            addParameter(p,'title',obj.title,@ischar)
            addParameter(p,'legend',obj.legend,@(x) iscell(x) || strcmp(x,'none'))
            addParameter(p,'vSpace',obj.vSpace,@isscalar)
            addParameter(p,'borderFactor',obj.borderFactor,@isscalar)
            parse(p, varargin{:}) %parse inputs
            obj.ev = p.Results.ev;
            obj.limy = p.Results.limy;
            obj.fontsize = p.Results.fontsize;
            obj.axNV = p.Results.axNV;
            obj.lineNV = p.Results.lineNV;
            obj.size = p.Results.size;
            obj.title = p.Results.title;
            obj.legend = p.Results.legend;
            obj.vSpace = p.Results.vSpace;
            obj.borderFactor = p.Results.borderFactor;
            mmError = p.Results.mmError;
            lineSpec = p.Results.lineSpec; 
        end
        
        function makePlots(obj, wavelengths ,mmData, varargin)
            [mmError, lineSpec] = parseVariableInputs(obj, varargin{:});
            % Determine if figure fits on screen and resize
            scrsz = get(0,'screensize');
            figPos = [1 5 obj.size];
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
            set(obj.figHandle,'position',figPos,'units','pixels'); %size figure
            obj.size = figPos(3:4);
            xLabel = uicontrol('style','text','BackgroundColor','w',...
                'units','pixels','FontSize',obj.fontsize,...
                'tag','xLabelObject'); % create x-label
            if obj.ev == true
                set(xLabel,'String','Energy (eV)');
            else
                set(xLabel,'String','Wavelength (nm)');
            end
            xLabel_sz = get(xLabel,'extent');
            set(xLabel,'Position',[(figPos(3) - xLabel_sz(3) )./2, 0, xLabel_sz(3), xLabel_sz(4)]);
            
            if ~isempty(obj.title) % create title if given
                figTitle = uicontrol('style','text','BackgroundColor','w',...
                    'units','pixels','FontSize',obj.fontsize,...
                    'tag','titleObject');
                set(figTitle,'String',obj.title)
                figTitle_sz = get(figTitle,'extent');
                set(figTitle,'Position',[( figPos(3) - figTitle_sz(3) )./2,...
                    ( figPos(4) - figTitle_sz(4) ), figTitle_sz(3), figTitle_sz(4)]);
            end
            % determine the horizontal extent of y-axis marker labels
            dummy = uicontrol('style','text','fontsize',obj.fontsize,'units','pixels');
            set(dummy,'String','-0.000');
            yAxSz = get(dummy,'extent');
            delete(dummy)
            
            plotSzX = figPos(3)/4 - yAxSz(3) - yAxSz(3)./5; % X size of plot area in pixels
            plotSzY = ( figPos(4) - 4*yAxSz(4) )/4 - 6 - obj.vSpace; % Y size of plot area in pixels
            for i=1:4
                for j=1:4
                    plotPos = [ ( (plotSzX + yAxSz(3) + 3)*(j-1) + yAxSz(3) +5)./figPos(3) ,...
                        ((plotSzY + yAxSz(4)./2 + obj.vSpace)*(4-i)+yAxSz(4)*2 + 3)./figPos(4),...
                        plotSzX./figPos(3), plotSzY./figPos(4)];
                    hand = subplot('Position',plotPos);
                    hold(hand,'on')
                    box(hand,'on')
                    if i ~= 4
                        set(hand,'XTickLabel',[]) % keep X lables only for bottom row
                    end
                    obj.axesHandles(i, j) = hand;
                end
            end
        plot(obj, wavelengths, mmData, mmError, lineSpec);
        end
        
        function plot(obj, wavelengths, mmData, mmError, lineSpec)
            % set Axes properties
            axis(obj.axesHandles(:),'tight'); % first, axes are set to tight
            if ~isempty(obj.axNV)
                for j=1:4
                    for k=1:4
                        set(obj.axesHandles(j,k),obj.axNV{:});
                    end
                end
            end
            %plot data and set Line properties.
            if obj.ev; wavelengths = 1239.8./wavelengths; end
            if isempty(mmError)
                for j = 1:4
                    for k = 1:4
                        plot(obj.axesHandles(j, k),wavelengths,squeeze(mmData(j,k,:,:)),...
                            lineSpec, obj.lineNV{:})
                    end
                end
            else
                for j = 1:4
                    for k = 1:4
                        errorbar(obj.axesHandles(j,k),wavelengths,squeeze(mmData(j,k,:,:)),...
                            squeeze(mmError(j,k,:,:)),...
                            LineSpec,'CapSize',0,obj.lineNV{:})
                    end
                end
            end
            if obj.limy ~= 0 % modify axes bounds if limY is set
                lim = obj.limy;
                for j=1:4
                    for k=1:4
                        Ylim = get(obj.axesHandles(j,k),'YLim');
                        if (Ylim(2) - Ylim(1)) < lim
                            avg = (Ylim(2) + Ylim(1))./2;
                            Ylim(2) = avg + lim/2;
                            Ylim(1) = avg - lim/2;
                            set(obj.axesHandles(j,k),'Ylim',Ylim);
                        end
                    end
                end
            end
            % Adjust plot limits so that lines do not overlap axis borders.
            % *** If you like to use Markers, then perhaps change 'lineWidth' to 'MarkerSize'
            lineHandle = get(obj.axesHandles(1,1),'children');
            lineWidth = zeros(size(lineHandle));
            for j = 1:length(lineHandle)
                lineWidth(j) = get(lineHandle(j),'lineWidth');
            end
            lineWidth = max(lineWidth)*obj.borderFactor;
            plotPos = get(obj.axesHandles(1,1),'Position');
            for j=1:4
                for k=1:4
                    xlim = get(obj.axesHandles(j,k),'xLim');
                    ylim = get(obj.axesHandles(j,k),'yLim');
                    xStep = (xlim(2) - xlim(1))/plotPos(3)/obj.size(1)*lineWidth/2;
                    yStep = (ylim(2) - ylim(1))/plotPos(4)/obj.size(1)*lineWidth;
                    set(obj.axesHandles(j,k),'XLim',[xlim(1)-xStep,xlim(2)+xStep]);
                    set(obj.axesHandles(j,k),'YLim',[ylim(1)-yStep,ylim(2)+yStep]);
                end
            end
            % set font size of all graphics objects
            set(get(obj.figHandle,'children'),'FontSize',obj.fontsize);
            % optionally create legend (this will increase the width of the figure!)
            if ~isempty(obj.legend)
                if iscell(obj.legend)
                    Labels = obj.legend{2};
                    if isnumeric(Labels)
                        Labels = cellfun(@(x) num2str(x),num2cell(Labels),'uniformoutput',0);
                    end
                    pos = zeros(4,4,4);
                    for j=1:4
                        for k=1:4
                            set(obj.axesHandles(j,k),'units','pixels');
                            pos(:,j,k) = get(obj.axesHandles(j,k),'Position');
                        end
                    end
                    lgd = legend(obj.axesHandles(1,4),Labels,'location','northeastoutside'); %#ok<*CPROPLC>
                    set(lgd,'units','pixels','fontsize',obj.fontsize);
                    title(lgd,obj.legend{1},'FontSize',obj.fontsize);
                    lgd_pos = get(lgd,'Position');
                    obj.figHandle.Position = obj.figHandle.Position + [0 0 lgd_pos(3) 0];
                    for j=1:4
                        for k=1:4
                            set(obj.axesHandles(j,k),'Position',pos(:,j,k));
                        end
                    end
                end
            end
        end
        
    end
    
end
