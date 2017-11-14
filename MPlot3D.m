classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MPlot3D < handle
    
    properties
       
        uniquezero = true
        palette =  'HotCold Bright'
        gs = 0
        width
        fontsize = 14
        limz = 1e-3
        norm = true
        hSpacing = 3; 
        vSpacing = 3; 
        cbw = 10;
    end
    
    properties (SetAccess = protected)
        figHandle
        axesHandles = gobjects(4);
        colorbarHandles = gobjects(4);
    end
    
    properties (Hidden)
        maskHandles
    end
    
    methods
        
        function obj = MPlot3D(varargin)
            obj.figHandle = figure;
            obj.width = getFigWidth;
            plot(obj,varargin{:});
        end
        
        function plot(obj,data,varargin)
            % \\ Required positional input:
            
            % data: [4,4,X,Y] array. X and Y are horizontal and vertical plot
            %     dimensions.
            
            % \\ Optional Name-Value pairs that set object properties
            
            %   'uniquezero', logical:   make zero white or black in colormaps
            %       Default is true.
            %   'palette', string: name of a colormap, including custom
            %       ones in local function colPalette
            %       Default is 'Fireice'
            %   'gs', [min max]: GlobalScale. plot limits of all Z-scales between min, max.
            %       If not given, each MM element maps to its own min and max value.
            %       Only 1 colorbar is drawn with GlobalScale is set
            %   'fontsize', scalar: Size of font in colorbars
            %   'width', scalar: Width of figure in inches. Height is
            %       computed automatically to ensure no streching of plots (figure can go
            %       off page, in which case, reduce value of 'width'. Default is %60 of
            %       the monitor or with dual displays of different size, who knowns...
            %   'limz', scalar: limits how small the range of the z-axes can be.
            %   'hSpacing', scalar: sets the horizontal space between plots in pixels
            %   'vSpacing', scalar: sets the vertical space between plots in pixels
            %   'cbw', scalar: Colorbar width in pixels.
            
            p = inputParser;
            % setup input scheme
            addRequired(p,'obj',@(x) isa(x,'MPlot3D'))
            addRequired(p,'data',@(x) isnumeric(x) && ndims(x) == 4)
            addParameter(p,'norm',obj.norm,@(x) x == 1 || x == 0)
            addParameter(p,'uniquezero',obj.uniquezero,@(x) x == 1 || x == 0)
            addParameter(p,'palette',obj.palette,@(x) ischar(x) )
            addParameter(p,'limz',obj.limz,@(x) isscalar(x)&&isnumeric(x))
            addParameter(p,'fontsize',obj.fontsize,@(x) isscalar(x)&&isnumeric(x))
            addParameter(p,'width',obj.width,@(x) isscalar(x) && isnumeric(x)) % inches
            addParameter(p,'gs',obj.gs,@(x) length(x) == 2 && isnumeric(x))
            addParameter(p,'hSpacing',obj.hSpacing,@isscalar)
            addParameter(p,'vSpacing',obj.vSpacing,@isscalar)
            addParameter(p,'cbw',obj.cbw,@isscalar)
            parse(p,obj,data,varargin{:}) %parse inputs
            sz = size(data);
            obj.norm = p.Results.norm;
            obj.uniquezero = p.Results.uniquezero;
            obj.palette = p.Results.palette;
            obj.gs = p.Results.gs;
            obj.limz = p.Results.limz;
            obj.fontsize = p.Results.fontsize;
            obj.width = p.Results.width;
            obj.hSpacing = p.Results.hSpacing;
            obj.vSpacing = p.Results.vSpacing;
            obj.cbw = p.Results.cbw;
            
            % normalize and replace NaN with 0 if obj.norm is set
            if obj.norm
                data = data ./ data(1,1,:,:);
                data(isnan(data)) = 0;
            end
            
            dummy = uicontrol('style', 'text', 'fontsize', obj.fontsize, 'units', 'pixels');
            set(dummy,'String', '-0.000');
            cblbextents = get(dummy, 'extent');
            cblbsz = cblbextents(3); % colorbar label size
            delete(dummy)
            figWidth = (obj.width) * obj.figHandle.Parent.ScreenPixelsPerInch;

            if obj.gs==0
                plotW = (figWidth - 9*obj.vSpacing-4*(obj.cbw + cblbsz))/4;
                plotH = sz(3)/sz(4)*plotW;
                figHeight = plotH*4+5*obj.hSpacing;
                
                totalPlotWidth = obj.vSpacing*2+obj.cbw+cblbsz+plotW;
                plotPosFun = @(j,k) [ (obj.vSpacing+(k-1)*totalPlotWidth)/figWidth...
                    ,(obj.hSpacing+(4-j)*(plotH+obj.hSpacing))/figHeight,...
                    plotW/figWidth,...
                    plotH/figHeight];
                set(obj.figHandle,'Position',[0,0,figWidth,figHeight],'units','pixels');
                for j=1:4
                    for k=1:4
                        if isgraphics(obj.axesHandles(j,k))
                            subplot(obj.axesHandles(j,k), ...
                                'position',plotPosFun(j,k),'units','pixels');
                        else
                            obj.axesHandles(j,k) = ...
                                subplot('position',plotPosFun(j,k),'units','pixels');
                        end
                        clim = [min(min(data(j,k,:,:))),max(max(data(j,k,:,:)))];
                        if obj.limz ~= 0 % modify axes bounds if limz is set
                            if (clim(2) - clim(1)) < obj.limz
                                avg = (clim(2) + clim(1))./2;
                                clim(2) = avg + obj.limz/2;
                                clim(1) = avg - obj.limz/2;
                            end
                        end
                        pos = get(obj.axesHandles(j,k),'Position');
                        imagesc(squeeze(data(j,k,:,:)),'Parent',obj.axesHandles(j,k),clim)
                        axis(obj.axesHandles(j,k),'off')
                        colormap(obj.axesHandles(j,k),makeColormap(clim,obj.uniquezero,obj.palette))
                        obj.colorbarHandles(j,k) = colorbar(obj.axesHandles(j,k),'units','pixels',...
                            'Position',[pos(1)+pos(3)+obj.vSpacing,pos(2)+cblbextents(4)/4,...
                            obj.cbw,pos(4)-cblbextents(4)/2],...
                            'fontsize',obj.fontsize);
                    end
                end
                if any(strcmp('nonorm', p.UsingDefaults))
                    obj.axesHandles(1,1).CLim = [0 1];
                end
            else
                plotW = (figWidth - 6*obj.vSpacing - 2*obj.cbw - cblbsz)/4;
                plotH = sz(3)/sz(4)*plotW;
                figHeight = plotH*4+5*obj.hSpacing;
                plotPosFun = @(j,k) [ (obj.vSpacing+(k-1)*(plotW+obj.vSpacing))/figWidth,...
                    (obj.hSpacing+(4-j)*(plotH+obj.hSpacing))/figHeight,...
                    plotW/figWidth,...
                    plotH/figHeight];
                set(obj.figHandle,'Position',[0,0,figWidth,figHeight],'units','pixels');
                for j=1:4
                    for k=1:4
                        if isgraphics(obj.axesHandles(j,k))
                            subplot(obj.axesHandles(j,k),...
                                'position',plotPosFun(j,k),'units','pixels');
                        else
                            obj.axesHandles(j,k) = ...
                                subplot('position',plotPosFun(j,k),'units','pixels');
                        end
                        pos = get(obj.axesHandles(j,k),'Position');
                        imagesc(squeeze(data(j,k,:,:)),'Parent',obj.axesHandles(j,k),obj.gs)
                        colormap(obj.axesHandles(j,k),makeColormap(obj.gs,obj.uniquezero,obj.palette))
                        axis(obj.axesHandles(j,k),'off')
                    end
                end
                obj.colorbarHandles(1,4) = colorbar(obj.axesHandles(1,4),'units','pixels',...
                    'Position',[pos(1)+pos(3)+obj.vSpacing,cblbextents(4)/4+6,...
                    obj.cbw,figHeight-cblbextents(4)/2-12],...
                    'fontsize',obj.fontsize);
            end
        end
        
        function mmdata = getPlotData(obj)
            % h: [4,4] array of axis handles
            mmdata = zeros([4, 4, size(obj.axesHandles(1,1).Children.CData)], ...
                class(obj.axesHandles(1,1).Children.CData));
            for j=1:4
                for k=1:4
                    mmdata(j,k,:,:) = obj.axesHandles(j,k).Children.CData;
                end
            end
        end

        function replacePlotData(obj,mmdata)
            % MMreplace3DplotData  replaces the data in 4x4 intensity plots.
            % h is a [4,4] array of axis handles
            % Data is a 4x4xNxM array. Data size should not be different than data in
            % plots.
            for j=1:4
                for k=1:4
                    obj.axesHandles(j,k).Children.CData = squeeze(mmdata(j,k,:,:));
                end
            end
        end
        
        function update(obj,varargin)
            obj.figHandle.Visible = 'off';
            data = getPlotData(obj);
            delete(obj.colorbarHandles)
            obj.colorbarHandles = gobjects(4);
            plot(obj,data,varargin{:});
            obj.figHandle.Visible = 'on';
        end
        
        function drawMask(obj, i, j)
            sz = size(obj.axesHandles(1,1).Children.CData);
            h_im = obj.axesHandles(i,j).Children;
            e = imellipse(obj.axesHandles(i,j),...
                [sz(1)*0.1,sz(2)*0.1,0.8*sz(1),0.8*sz(1)]);
            obj.maskHandles = {e,h_im};
        end
        
        function setElipse(obj, position)
            setPosition(obj.maskHandles{1}, position);
        end
        
        function applyMask(obj)
            mask = createMask(obj.maskHandles{1}, obj.maskHandles{2});
            delete(obj.maskHandles{1})
            for j=1:4
                for k=1:4
                    obj.axesHandles(j,k).Children.CData = ...
                        obj.axesHandles(j,k).Children.CData.*mask;
                end
            end
        end
        
        function applyMaskWithTrim(obj)
            mask = createMask(obj.maskHandles{1},obj.maskHandles{2});
            pos = obj.maskHandles{1}.getPosition;
            pos = [floor(pos(1:2)),ceil(pos(3:4))];
            delete(obj.maskHandles{1})
            data = zeros(4,4,pos(4),pos(3),'single');
            idx = {pos(2):(pos(2)+pos(4)-1),pos(1):(pos(1)+pos(3)-1)};
            for j=1:4
                for k=1:4
                    data(j,k,:,:) = ...
                        obj.axesHandles(j,k).Children.CData(idx{:}).*mask(idx{:});
                end
            end
            plot(obj,data)
        end
        
        function print(obj,filepath)
            print(obj.figHandle,filepath,'-depsc');
        end
        
        function flipX(obj)
            replacePlotData(obj, flip(getPlotData(obj), 4));
        end
    end
end

function width = getFigWidth  % sets default width to 60% of display width
            scrsz = get(0,'screensize');
            width = 0.6*scrsz(3)/get(0,'ScreenPixelsPerInch');
end

function colAr = colPalette(palette)
% these are custom colormaps. A colormap is just a Nx4 matrix. The
% first column are values between 0 and 256 that position a color marker.
% The 2nd, 3rd, and 4th columns are RGB color values. Names of matlab
% colormaps can also be handed to this function.
switch palette
        
    case 'HotCold Bright'
        colAr = ...
            [0	0	65	220;...
            36	0	90	240;...
            76	0	253	253;...
            128	250	250	250;...
            182	255	242	0;...
            224	255	127	0;...
            256	255	0	0];
        
    case 'HotCold Dark'
        colAr = ...
            [0	0	253	253;...
            36	1	114	239;...
            76	0	90	240;...
            128	0	0	0;...
            182	255	0	0;...
            224	255	127	0;...
            256	255	242	0];
        
    case 'TwoTone Bright'
        colAr = ...
            [0	0	0	255;...
            128	255	255	255;...
            256	255	0	0];
        
    case 'TwoTone Dark'
        colAr = ...
            [0	0	0	255;...
            128	0	0	0;...
            256	255	0	0];
        
    case 'Fireice'        
        %Copyright (c) 2009, Joseph Kirk 
        %All rights reserved.
        clrs = [0.75 1 1; 0 1 1; 0 0 1;...
            0 0 0; 1 0 0; 1 1 0; 1 1 0.75];
        
        y = -3:3;
        m = 64;
        if mod(m,2)
            delta = min(1,6/(m-1));
            half = (m-1)/2;
            yi = delta*(-half:half)';
        else
            delta = min(1,6/m);
            half = m/2;
            yi = delta*nonzeros(-half:half);
        end
        colAr = cat(2,(0:4:255).',255*interp2(1:3,y,clrs,1:3,yi));
        
    case 'Spectral'
        colAr = cbrewer('div', 'Spectral', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'RdYlGn'
        colAr = cbrewer('div', 'RdYlGn', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'RdYlBu'
        colAr = cbrewer('div', 'RdYlBu', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'RdBu'
        colAr = cbrewer('div', 'RdBu', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'RdGy'
        colAr = cbrewer('div', 'RdGy', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'PuOr'
        colAr = cbrewer('div', 'PuOr', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'PRGn'
        colAr = cbrewer('div', 'PRGn', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'PiYG'
        colAr = cbrewer('div', 'PiYG', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    case 'BrBG'
        colAr = cbrewer('div', 'BrBG', 11) .* 255;
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
        
    otherwise
        colAr = colormap(palette) .* 255; % to use other colormaps
        t1 = linspace(0,256,size(colAr, 1)).';
        colAr = [t1, colAr];
end
end

function fracIndx = fracIndex(array,x)

fracIndx = zeros(1,length(x));
for idx = 1:length(x)
    if x >= array(end)
        fracIndx(idx) = length(array);
    elseif x(idx) <= array(1)
        fracIndx(idx) = 1;
    else
        a = find(array <= x(idx));
        a = a(length(a));
        b = find(array > x(idx));
        b = b(1);
        fracIndx(idx) = a+(x(idx)-array(a))/(array(b)-array(a));
    end
    
end
end

function cm = makeColormap(clim,b_uniqueZero,palette)
dmin=clim(1);
dmax=clim(2);
if dmax == dmin
    dmax=1;
    dmin=0;
end
if b_uniqueZero == true
    Zscale = zeros(1,256);
    if abs(dmin) < abs(dmax)
        didx = (dmax - dmin)/(2*dmax);
        for idx = 0:255
            Zscale(idx+1) = 256 - didx*idx;
        end
    else
        didx = (dmin-dmax)/(2*dmin);
        for idx = 0:255
            Zscale(idx+1) = idx*didx;
        end
        Zscale = flip(Zscale);
    end
else
    Zscale = flip(1:256);
end
colAr = colPalette(palette);
cm = zeros(256,3);
for n = 1:256
    x = fracIndex(colAr(:,1),Zscale(n));
    cm(n,1) = interp1(colAr(:,2),x);
    cm(n,2) = interp1(colAr(:,3),x);
    cm(n,3) = interp1(colAr(:,4),x);
end
cm = cm./255;
cm = flip(cm,1);
end
