classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MPlot4D < handle
    % this is a more powerful but less polished version of MPlot3D. It can
    % accept arrays of dimension 5 and make videos that run over the 5th
    % dimension. The constructor requires an array xData, which is the physical
    % values ascribed to the 5th dimension of Data (usually wavelength).
    % These values are put into a text box above the MM plots. I never got
    % around to documenting this class...
    
    properties
        
        Data
        xData
        uniquezero = true
        palette =  'Fireice'
        gs = 0
        width
        fontsize = 14
        limz = 1e-3
        norm = true
        hSpacing = 3 %vertical spacing of subplots
        vSpacing = 3 % horizonal spacing of subplots
        cbw = 10
        
    end
    
    properties (SetAccess = protected)
        
        figHandle
        axesHandles = gobjects(4)
        colorbarHandles = gobjects(4)
        xDataTextBox
        
    end
    
    properties (Hidden)
        
        maskHandles
        
    end
    
    methods
        
        function obj = MPlot4D(data, xData, varargin)
            obj.figHandle = figure;
            obj.width = getFigWidth;
            obj.Data = data;
            obj.xData = xData;
            plot(obj, xData(1), varargin{:});
        end
        
        function parseProperties(obj, varargin)
            %   Optional Name-Value pairs
            
            %   'uniquezero',logical:   make 0 white or black in color palettes
            %       Default is true.
            %   'palette','colPalette': string giving name of a case in colPalette.m
            %       Default is 'Fireice'
            %   'gs',[min max]: GlobalScale, plot limits of all Z-scales between min, max.
            %       If not given, each MM element maps to its own min and max value.
            %       Only 1 colorbar is drawn with GlobalScale is set
            %   'fontsize',scalar: Size of font in colorbars
            %   'width',scalar: Width of figure in inches (but probably not inches). Height is
            %       computed automatically to ensure no streching of plots (figure can go
            %       off page, in which case, reduce value of 'width'. Default is %60 of
            %       monitor on a Mac.
            %   'limz',scalar: limits how small the range of the z-axes can be.
            %   'hSpacing',scalar: horizontal space between axes, in pixels.
            %   'vSpacing',scalar: vertical space between axes, in pixels.
            %   'cbw',scalar: Width of the colorbars, in pixels. 
            p = inputParser;
            % setup input scheme
            addRequired(p,'obj',@(x) isa(x,'MPlot4D'))
            addParameter(p,'norm',obj.norm,@(x) strcmp(x,'nonorm'))
            addParameter(p,'uniquezero',obj.uniquezero,@islogical)
            addParameter(p,'palette',obj.palette,@ischar)
            addParameter(p,'limz',obj.limz,@(x) isscalar(x)&&isnumeric(x))
            addParameter(p,'fontsize',obj.fontsize,@(x) isscalar(x)&&isnumeric(x))
            addParameter(p,'width',obj.width,@(x) isscalar(x) && isnumeric(x)) % inches
            addParameter(p,'gs',obj.gs,@(x) length(x) == 2 && isnumeric(x))
            addParameter(p,'hSpacing',obj.hSpacing,@isscalar)
            addParameter(p,'vSpacing',obj.vSpacing,@isscalar)
            addParameter(p,'cbw',obj.cbw,@isscalar)
            parse(p, obj, varargin{:}) %parse inputs
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
            
        end
        
        function normalize(obj)
            obj.Data = obj.Data ./ obj.Data(1,1,:,:,:);
            obj.Data(isnan(obj.Data)) = 0;
        end
        
        function plot(obj, xVal, varargin)
            if ~isempty(varargin)
                parseProperties(obj, varargin{:})
            end
            dataIndex = round(fracIndex(obj.xData,xVal));
            
            sz = size(obj.Data);
            dummy = uicontrol('style', 'text', 'fontsize', obj.fontsize, 'units', 'pixels');
            set(dummy,'String', '-0.000');
            cblbextents = get(dummy, 'extent');
            cblbsz = cblbextents(3); % colorbar label size
            delete(dummy)
            figWidth = (obj.width) * obj.figHandle.Parent.ScreenPixelsPerInch;

            if obj.gs==0
                plotW = (figWidth - 9*obj.hSpacing-4*(obj.cbw + cblbsz))/4;
                plotH = sz(3)/sz(4)*plotW;
                figHeight = plotH*4+5*obj.vSpacing + cblbextents(4);
                
                totalPlotWidth = obj.hSpacing*2+obj.cbw+cblbsz+plotW;
                plotPosFun = @(j,k) [ (obj.hSpacing+(k-1)*totalPlotWidth)/figWidth...
                    ,(obj.vSpacing+(4-j)*(plotH+obj.vSpacing))/figHeight,...
                    plotW/figWidth,...
                    plotH/figHeight];
                set(obj.figHandle,'Position',[0,0,figWidth,figHeight],'units','pixels');
                for j=1:4
                    for k=1:4
                        obj.axesHandles(j,k) = ...
                            subplot('position',plotPosFun(j,k),'units','pixels');
                        clim = [min(min(obj.Data(j,k,:,:,dataIndex))),max(max(obj.Data(j,k,:,:,dataIndex)))];
                        if obj.limz ~= 0 % modify axes bounds if limz is set
                            if (clim(2) - clim(1)) < obj.limz
                                avg = (clim(2) + clim(1))./2;
                                clim(2) = avg + obj.limz/2;
                                clim(1) = avg - obj.limz/2;
                            end
                        end
                        pos = get(obj.axesHandles(j,k),'Position');
                        imagesc(squeeze(obj.Data(j,k,:,:,dataIndex)),'Parent',obj.axesHandles(j,k),clim)
                        axis(obj.axesHandles(j,k),'off')
                        colormap(obj.axesHandles(j,k),makeColormap(clim,obj.uniquezero,obj.palette))
                        obj.colorbarHandles(j,k) = colorbar(obj.axesHandles(j,k),'units','pixels',...
                            'Position',[pos(1)+pos(3)+obj.hSpacing,pos(2)+cblbextents(4)/4,...
                            obj.cbw,pos(4)-cblbextents(4)/2],...
                            'fontsize',obj.fontsize);
                    end
                end
%                 if any(strcmp('nonorm', p.UsingDefaults))
%                     obj.axesHandles(1,1).CLim = [0 1];
%                 end
            else
                plotW = (figWidth - 6*obj.hSpacing - 2*obj.cbw - cblbsz)/4;
                plotH = sz(3)/sz(4)*plotW;
                figHeight = plotH*4+5*obj.vSpacing + cblbextents(4);
                plotPosFun = @(j,k) [ (obj.hSpacing+(k-1)*(plotW+obj.hSpacing))/figWidth,...
                    (obj.vSpacing+(4-j)*(plotH+obj.vSpacing))/figHeight,...
                    plotW/figWidth,...
                    plotH/figHeight];
                set(obj.figHandle,'Position',[0,0,figWidth,figHeight],'units','pixels');
                for j=1:4
                    for k=1:4
                        obj.axesHandles(j,k) = ...
                            subplot('position',plotPosFun(j,k),'units','pixels');
                        pos = get(obj.axesHandles(j,k),'Position');
                        imagesc(squeeze(obj.Data(j,k,:,:,dataIndex)),'Parent',obj.axesHandles(j,k),obj.gs)
                        colormap(obj.axesHandles(j,k),makeColormap(obj.gs,obj.uniquezero,obj.palette))
                        axis(obj.axesHandles(j,k),'off')
                    end
                end
                obj.colorbarHandles(1,4) = colorbar(obj.axesHandles(1,4),'units','pixels',...
                    'Position',[pos(1)+pos(3)+obj.hSpacing, cblbextents(4)/4+6,...
                    obj.cbw,figHeight-3*cblbextents(4)/2-12],...
                    'fontsize',obj.fontsize);
            end
            obj.xDataTextBox = ...
                uicontrol('style', 'text', 'fontsize', obj.fontsize, 'units', 'pixels', ...
                'position', [figWidth/2-cblbextents(3), figHeight-cblbextents(4), cblbextents(3:4)]);
            set(obj.xDataTextBox, 'String', num2str(obj.xData(dataIndex)));
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

        function replacePlotData(obj, idx)
            % MMreplace3DplotData  replaces the data in 4x4 intensity plots.
            % h is a [4,4] array of axis handles
            % Data is a 4x4xNxM array. Data size should not be different than data in
            % plots.
            if obj.gs == 0
                for j=1:4
                    for k=1:4
                        obj.axesHandles(j,k).Children.CData = squeeze(obj.Data(j,k,:,:,idx));
                        clim = [min(min(obj.Data(j,k,:,:,idx))),max(max(obj.Data(j,k,:,:,idx)))];
                        if obj.limz ~= 0 % modify axes bounds if limz is set
                            if (clim(2) - clim(1)) < obj.limz
                                avg = (clim(2) + clim(1))./2;
                                clim(2) = avg + obj.limz/2;
                                clim(1) = avg - obj.limz/2;
                            end
                        end
                        obj.axesHandles(j,k).CLim = clim;
                        colormap(obj.axesHandles(j,k),makeColormap(clim,obj.uniquezero,obj.palette))
                    end
                end
            else
                for j=1:4
                    for k=1:4
                        obj.axesHandles(j,k).Children.CData = squeeze(obj.Data(j,k,:,:,idx));
                    end
                end
            end
%                     
        end
        
        function makeAVI(obj, xRange, AVIfilename)
            if isempty(obj.xData)
                xVals = xRange;
            else
                [X,I] = sort(obj.xData); % added this to allow unsorted xData
                indices = unique(round(fracIndex(X,xRange)),'first');
                xVals = I(indices);
            end
            
            v = VideoWriter(AVIfilename);
            v.FrameRate = 10;
            open(v);
            for i=xVals
                replacePlotData(obj, i)
                set(obj.xDataTextBox, 'String', num2str(obj.xData(i)));
                writeVideo(v, getframe(obj.figHandle));
            end
            close(v);

        end
        
        function update(obj, varargin)
            obj.figHandle.Visible = 'off';
            data = getPlotData(obj);
            delete(obj.axesHandles);
            delete(obj.colorbarHandles)
            obj.axesHandles = gobjects(4);
            obj.colorbarHandles = gobjects(4);
            plot(obj,data,varargin{:});
            obj.figHandle.Visible = 'on';
        end
        
    end
end


function width = getFigWidth  % sets default width to 60% of display width
            scrsz = get(0,'screensize');
            width = 0.6*scrsz(3)/get(0,'ScreenPixelsPerInch');
end

function colAr = colPalette(palette)
% these are custom color palettes. A palette is just a Nx4 matrix. The
% first column are values between 0 and 256 that position a color marker.
% The 2nd, 3rd, and 4th columns are RGB color values.
switch palette
    case 'Rainbow'
        colAr = ...
            [0	255	0	241;...
            36	0	65	220;...
            86	0	253	253;...
            128	0	255	15;...
            171	255	242	0;...
            234	255	127	0;...
            256	255	0	0];
        
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
end
end

function fracIndx = fracIndex(array,x)

fracIndx = zeros(1,length(x));
for idx = 1:length(x)
    if x(idx) >= array(end)
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
