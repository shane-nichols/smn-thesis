classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) CRRmeasurement
    
    properties
        
        N
        periods
        dropped_frames
        frame_rate
        motor_speeds
        binning
        roi
        exposure
        dark_value
        time
        wavelength
        jpeg_multiplier
        mmtype
        p1
        p2
        d1
        d2
        pol1
        pol2
        illumination
        specimen
        mmdata
        frames
        path
        
    end
    
    methods
        
        function obj = CRRmeasurement(varargin)
            obj.illumination = obj.imageplane;  % populate structure fields
            obj.specimen = obj.imageplane;
            if nargin ~= 0 && isdir(varargin{1})
                obj.path = varargin{1};
                [~,basepath] = fileparts(obj.path);
                t = readtable(fullfile(obj.path,[basepath,'_parameters.dat']),'HeaderLines',0,'Delimiter',',');
                t = table2cell(t);
                t = cell2struct({t{:,2}},{t{:,1}},2);
                obj.N = str2double(t.Frames);
                obj.dropped_frames = double(t.DroppedFrames);
                obj.frame_rate = str2double(t.FrameRate);
                temp = str2double(t.FastSpeedMotor);
                obj.motor_speeds = [temp, temp./str2double(t.MotorRatio)];
                obj.binning = int32(str2double(t.BinningFactor));
                obj.roi = str2num(t.ROI); %#ok<ST2NM>
                obj.exposure = str2double(t.ExposureTime);
                obj.dark_value = str2double(t.DarkValue);
                obj.time = t.Time;
                obj.wavelength = str2double(t.Wavelength);
                obj.jpeg_multiplier = str2double(t.jpegMultiplier);
                if isfield(t,'MatrixType'); obj.mmtype = t.MatrixType; end
                obj.periods = abs(obj.N ./ obj.frame_rate .* obj.motor_speeds(2) ./ 2);
                opticsPhases = num2cell(str2num(t.p1p2d1d2pol1pol2)); %#ok<ST2NM>
                [obj.p1, obj.p2, obj.d1, obj.d2, obj.pol1, obj.pol2] = opticsPhases{:};
                obj.mmdata = readMMI(obj.path);
                if ~strcmp('N/A', t.PrecessionAmplitude)
                    obj.specimen.amplitude = str2double(t.PrecessionAmplitude);
                    obj.specimen.phase = str2double(t.PrecessionPhase);
                end
                if isfield(t,'IlluminationImagePath') && ~strcmp('N/A',t.IlluminationImagePath)
                    obj.illumination.image = imread(t.IlluminationImagePath);
                end
            end
        end
          
        function obj = setInstrument(obj, N, periods, p1, p2, d1, d2, pol1, pol2)
            inputs = {N, periods, p1, p2, d1, d2, pol1, pol2};
            [obj.N, obj.periods, obj.p1, obj.p2, obj.d1, obj.d2, obj.pol1, obj.pol2]...
                = inputs{:};
        end
           
        function obj = loadFrames(obj)
                filePattern = fullfile(obj.path,'frames', '*.jpg');
                jpgFiles   = dir(filePattern);
                for j=1:length(jpgFiles)
                    baseFileName = jpgFiles(j).name;
                    fullFileName = fullfile(obj.path,'frames', baseFileName);
                    if  j == 1
                        temp = imfinfo(fullFileName);
                        obj.frames = zeros(temp.Height, temp.Width,length(jpgFiles), 'uint16');
                    end
                    obj.frames(:,:,j) = imread(fullFileName);
                end
        end
        
        function obj = getSpecimenPrecession(obj,varargin)
            switch nargin
                case 1
                    framestep = 1;
                    useCV = false;
                case 2
                    if strcmp(varargin{1},'CV')
                        useCV = true;
                        framestep = 1;
                    else
                        useCV = false;
                        framestep = varargin{1};
                    end
                case 3
                    useCV = any(strcmp(varargin,'CV'));
                    framestep = varargin{cellfun(@isnumeric,varargin)};
            end
            X = 1:framestep:obj.N;  % subset of frames to process
            
            j = 1;
            if useCV
                frame1 = imadjust(obj.frames(:,:,1));
                T = zeros(2,length(X));
                for i=X(2:end)
                    j = j + 1;
                    % Estimate transform from frame A to frame B, and fit as an s-R-t
                    H = cvexEstStabilizationTform(frame1,imadjust(obj.frames(:,:,i)),0.01);
                    HsRt = cvexTformToSRT(H);
                    T(:,j) = HsRt(3,1:2);
                end
            else
                frame1 = obj.frames(:,:,1);
                tform(1,length(X)) = affine2d;
                transformation = 'translation';
                [optimizer,metric] = imregconfig('multimodal');
                optimizer.InitialRadius = 2e-4;
                optimizer.GrowthFactor = 1.02;
                for i=X(2:end)
                    j = j + 1;
                    tform(j) = imregtform(obj.frames(:,:,i),frame1,transformation,optimizer,metric,...
                        'InitialTransformation',tform(j-1));
                end
                % exract the translation components from the affine transformation
                T = {tform.T};
                T = cell2mat(cellfun(@(x) x(3,1:2).',T,'uniformoutput',0));
                X = X - 1; % move origin of time to zero
            end
            % plotting and fitting
            p_1 = lsqcurvefit(@(p,X) p(1)*cos(-4*pi*obj.periods*X/obj.N - p(2)) + p(3),[6,0,0],X,T(1,:));
            p_2 = lsqcurvefit(@(p,X) p(1)*sin(-4*pi*obj.periods*X/obj.N - p(2)) + p(3),[6,0,0],X,T(2,:));
            figure
            h=axes;
            h.NextPlot = 'add';
            plot(h,X,T.','o');
            plot(h,X,p_1(1)*cos(-4*pi*obj.periods*X/obj.N - p_1(2)) + p_1(3),'LineWidth',3)
            plot(h,X,p_2(1)*sin(-4*pi*obj.periods*X/obj.N - p_2(2)) + p_2(3),'LineWidth',3)
            legend({'\delta_X','\delta_Y','\delta_X fitted','\delta_Y fitted'});
            h.XLabel.String = 'Frame #';
            h.YLabel.String = 'Pixels';
            if p_1(1) > 0
                p_1(2) = wrapTo2Pi(p_1(2) + pi);
            else
                p_1(2) = wrapTo2Pi(p_1(2));
            end
            if p_2(1) > 0
                p_2(2) = wrapTo2Pi(p_2(2) + pi);
            else
                p_2(2) = wrapTo2Pi(p_2(2));
            end
            obj.specimen.amplitude = (abs(p_1(1)) + abs(p_2(1)))/2;
            obj.specimen.phase = (p_1(2) + p_2(2))/2;
        end
        
        function obj = simulateRawData(obj,varargin)
            switch nargin
                case 1
                    wf = obj.waveform(eye(4));
                case 2
                    wf = obj.waveform(varargin{:});
            end
            %sz1 = size(obj.illumination.image);
            sz2 = size(obj.specimen.image);
            
            X = (0:(obj.N - 1)) / obj.N;
            illuminDisp = obj.illumination.amplitude * ...
                exp(-5i * obj.periods * X - 1i*obj.illumination.phase );
            illuminDisp = [real(illuminDisp); imag(illuminDisp)].';
            specimenDisp = obj.specimen.amplitude * ...
                exp(-4i * obj.periods * X - 1i*obj.specimen.phase );
            specimenDisp = [real(specimenDisp); imag(specimenDisp)].';
            
            obj.frames = zeros(sz2(1), sz2(2), obj.N, class(obj.specimen.image));
            for i=1:obj.N
                obj.frames(:,:,i) = wf(i).* ...
                    imtranslate( obj.specimen.image,     specimenDisp(i,:)) .* ...
                    imtranslate( obj.illumination.image, illuminDisp(i,:)) ;
            end
            A = ceil( max( obj.specimen.amplitude, obj.illumination.amplitude ));
            obj.frames = obj.frames(A:(end-A),A:(end-A),:);
        end
        
        function obj = correctFrames(obj)
            % Correct frames if necessary. Frames of an integer class are
            % converted to single precission float.
          
            ck1 = ~isempty(obj.illumination.amplitude) && ...
                  ~isempty(obj.illumination.phase) && ...
                  ~isempty(obj.illumination.image);
            ck2 = ~isempty(obj.specimen.amplitude) && ...
                  ~isempty(obj.specimen.phase);
            if bitand(~ck1,~ck2) && isempty(obj.illumination.image)
                error('No precession or illumination information is set.')
            end
              
            if isempty(obj.frames) || ndims(obj.frames) ~= 3
                error('The frames property must be a 3D array');
            end
            if isinteger(obj.frames)
                obj.frames = single(obj.frames);
            end
              
            if ck1 && ck2
                X = (0:(obj.N - 1)) / obj.N;
                illumin_disp = obj.illumination.amplitude * ...
                    exp(-5i*pi * obj.periods * X - 1i*obj.illumination.phase );
                illumin_disp = [real(illumin_disp); imag(illumin_disp)].';
                specimen_disp = obj.specimen.amplitude * ...
                    exp(-4i*pi * obj.periods * X - 1i*obj.specimen.phase );
                specimen_disp = -[real(specimen_disp); imag(specimen_disp)].';
                ill_frames = zeros([size(obj.illumination.image), obj.N], class(obj.frames));
                for i=1:obj.N
                    ill_frames(:,:,i) = ...
                        imtranslate( obj.illumination.image, illumin_disp(i,:) + specimen_disp(i,:));
                    obj.frames(:,:,i) = imtranslate( obj.frames(:,:,i), specimen_disp(i,:)) ;
                end
                A = ceil( max( obj.specimen.amplitude, obj.illumination.amplitude ));
                ill_frames = ill_frames(A:(end-A), A:(end-A),:);
                obj.frames = obj.frames ./ ill_frames;
                obj.frames = obj.frames(A:(end-A), A:(end-A),:);
                
            elseif ck2 && ~isempty(obj.illumination.image)
                X = (0:(obj.N - 1)) / obj.N;
                A = ceil(obj.specimen.amplitude);
                illumin_image = obj.illumination.image(A:(end-A), A:(end-A));
                obj.frames = bsxfun(@rdivide, obj.frames, illumin_image);
                specimen_disp = obj.specimen.amplitude * ...
                    exp(-4i*pi * obj.periods * X - 1i*obj.specimen.phase );
                specimen_disp = -[real(specimen_disp); imag(specimen_disp)].';
                for i=1:obj.N
                    obj.frames(:,:,i) = imtranslate( obj.frames(:,:,i), specimen_disp(i,:)) ;
                end
                obj.frames = obj.frames(A:(end-A), A:(end-A),:);
                
            elseif ck2
                X = (0:(obj.N - 1)) / obj.N;
                A = ceil(obj.specimen.amplitude);
                specimen_disp = obj.specimen.amplitude * ...
                    exp(-4i*pi * obj.periods * X - 1i*obj.specimen.phase );
                specimen_disp = -[real(specimen_disp); imag(specimen_disp)].';
                for i=1:obj.N
                    obj.frames(:,:,i) = imtranslate( obj.frames(:,:,i), specimen_disp(i,:)) ;
                end
                obj.frames = obj.frames(A:(end-A), A:(end-A),:);
            end
        end
        
        function obj = calculateMM(obj)
            if isinteger(obj.frames)
                obj.frames = single(obj.frames);
            end
            a = -sin(obj.d1);
            b = sin(obj.d2);
            c = (1+cos(obj.d1))/2;
            d = (1+cos(obj.d2))/2;
            e = (1-cos(obj.d1))/2;
            f = (1-cos(obj.d2))/2;
            % demodulation: determine 10 complex harmonic coefficients with retarder
            % phases. This is basically a discrete Fourier Transform
            t = 0 : 1 / obj.frame_rate : (obj.N - 1) / obj.frame_rate;
            t = t(:); % make sure that t is a column vector
            alpha = 2*(2*pi * obj.motor_speeds(1) .* t + obj.p1);
            beta = 2*(2*pi * obj.motor_speeds(2) .* t + obj.p2);
            alpha2 = 2*alpha;
            beta2 = 2*beta;
            B = exp(1i*...  % construct a basis matrix
                [alpha,...
                beta,...
                alpha2,...
                beta2,...
                alpha-beta,...
                alpha-beta2,...
                alpha2-beta,...
                alpha2+beta2,...
                alpha2-beta2]).';
            % demodulation, multiprod(B,I,[1 2],3). Dark value only effects
            % DC component. // should test DV correction on weak data. //
            B = cat(3, (sum(obj.frames, 3) - (obj.N * obj.dark_value))./ 2, crrMultiprod(B, obj.frames));
            
            % precalculate some exponentials of polarizer angles
            ep1 = exp(-2i * obj.pol1);
            ep2 = exp(-2i * obj.pol2);
            
            % take linear combinations of the Fourier coefficients
            M(1,1,:,:) = real(B(:,:,1) - c/e*ep1^2*B(:,:,4) - d/f*ep2^2*B(:,:,5) + ...
                c*d/(e*f)*(ep1^2*ep2^2*B(:,:,9) + ep1^2*ep2^2'*B(:,:,10)));
            temp = ep1*B(:,:,4)/e - d/(e*f)*(ep1*ep2^2*B(:,:,9) + ep1*ep2^2'*B(:,:,10));
            M(1,2,:,:) = real(temp);
            M(1,3,:,:) = imag(temp);
            M(1,4,:,:) = imag(ep1*B(:,:,2)/a - 2*d*ep1*ep2^2'*B(:,:,7)/(a*f));
            temp = ep2*B(:,:,5)/f - c/(e*f)*(ep1^2*ep2*B(:,:,9) + conj(ep1^2*ep2'*B(:,:,10)));
            M(2,1,:,:) = real(temp);
            M(3,1,:,:) = imag(temp);
            temp = (ep1*ep2*B(:,:,9) + ep1*ep2'*B(:,:,10))/(e*f);
            M(2,2,:,:) = real(temp);
            M(2,3,:,:) = imag(temp);
            temp = 2*ep1*ep2'*B(:,:,7)/(a*f);
            M(2,4,:,:) = imag(temp);
            M(3,4,:,:) = real(temp);
            temp = (-ep1*ep2*B(:,:,9) + ep1*ep2'*B(:,:,10))/(e*f);
            M(3,2,:,:) = -imag(temp);
            M(3,3,:,:) = real(temp);
            M(4,1,:,:) = imag(ep2*B(:,:,3)/b + 2*c*ep1^2*ep2'*B(:,:,8)/(b*e));
            temp =  2*ep1*ep2'*B(:,:,8)/(b*e);
            M(4,2,:,:) = -imag(temp);
            M(4,3,:,:) = real(temp);
            M(4,4,:,:) = real(2*ep1*ep2'*B(:,:,6)/(a*b));
            obj.mmdata = flip(permute(M./length(t)*8, [1,2,4,3]),4);
        end
        
        function obj = flipX(obj)
            obj.mmdata = flip(obj.mmdata,4);
        end
        
        function MPlot3Dobj = plot(obj, varargin)
            % plot can be overloaded because Axes class is inferior 
            MPlot3Dobj = MPlot3D(obj.mmdata, varargin{:});
        end
            
        function intensities = waveform(obj, mueller_matrix)
            t = linspace(0,2 * obj.periods, obj.N + 1);
            t = t(1:(end-1));
            intensities = 2*CRRmakeI2(...
                1.25, 1, obj.p1, obj.p2, obj.d1, obj.d2,...
                mueller_matrix, t, obj.pol1, obj.pol2);
        end
        
        function length = getScale(obj, magnification)
            % size of bottom edge in microns
            length = obj.roi(4)*obj.binning*6.5/magnification;
        end
        
    end
    
    methods(Static)
        
        function out = imageplane(varargin)
            out = struct('image',[],'amplitude',[],'phase',[]);
            switch nargin
                case 1
                    out.image = varargin{:};
                    
                case 2
                    [out.amplitude,out.phase] = varargin{:};
                    
                case 3
                    [out.image, out.amplitude, out.phase] = varargin{:};
            end
        end
        
    end
    
end

function I = CRRmakeI2(f1,f2,p1,p2,d1,d2,M,t,polAng1,polAng2)
% simulation of the continuous rotating retarder light intensity. Here, we
% are able to set the polarizer angles to any value.

% I = array of light intensity values at detector
% f0 = frequency of first retarder in Hz
% f1 = frequency of second retarder in Hz
% p0 = phase of first retarder in radians
% p1 = phase of second retarder in radians
% d0 = retardance of first retarder in radians
% d1 = retardance of second retarder in radians
% M = 4x4 test Mueller matrix
% t = time values to compute I
% polAng0 = angle of the input polarizer, in radians
% polAng1 = angle of the output polarizer, in radians

polAng1 = polAng1*2;
polAng2 = polAng2*2;
Cd1 = cos(d1);
Sd1 = sin(d1);
Cd2 = cos(d2);
Sd2 = sin(d2);
P1 = [1,cos(polAng1),sin(polAng1),0;...
    cos(polAng1),cos(polAng1)^2,cos(polAng1)*sin(polAng1),0;...
    sin(polAng1),cos(polAng1)*sin(polAng1),sin(polAng1)^2,0;...
    0,0,0,0]/2;
P2 = [1,cos(polAng2),sin(polAng2),0;...
    cos(polAng2),cos(polAng2)^2,cos(polAng2)*sin(polAng2),0;...
    sin(polAng2),cos(polAng2)*sin(polAng2),sin(polAng2)^2,0;...
    0,0,0,0]/2;
I = zeros(size(t));
for index = 1:length(t)
    arg1 = 2*(2*pi*f1.*t(index) + p1);
    arg2 = 2*(2*pi*f2.*t(index) + p2);
    Cr1 = cos(arg1);
    Sr1 = sin(arg1);
    Cr2 = cos(arg2);
    Sr2 = sin(arg2);
    Ret1 = [1,0,0,0;...
        0,Cr1.^2+Sr1.^2.*Cd1,Cr1.*Sr1.*(1-Cd1),-Sr1.*Sd1;...
        0,Cr1.*Sr1.*(1-Cd1),Sr1.^2+Cr1^2.*Cd1,Cr1.*Sd1;...
        0,Sr1.*Sd1,-Cr1.*Sd1,Sd1].';
    Ret2 = [1,0,0,0;...
        0,Cr2.^2+Sr2.^2.*Cd2,Cr2.*Sr2.*(1-Cd2),-Sr2.*Sd2;...
        0,Cr2.*Sr2.*(1-Cd2),Sr2.^2+Cr2^2.*Cd2,Cr2.*Sd2;...
        0,Sr2.*Sd2,-Cr2.*Sd2,Sd2].';
    I(index) = [1,0,0,0]*P2*Ret2*M*Ret1*P1*[1;0;0;0];
end
end

function mm = readMMI(Path,varargin)
% Supply the path to the directory containing the 16-jpg2000 images
% To normalize data, put input 'norm', as in: readMMI(Path,'norm')
filePattern = fullfile(Path, '*.jpg');
%where 'fullfile' builds the base file name
jpgFiles   = dir(filePattern);
fileIndex = 1;
b_norm =~ isempty(find(strcmpi(varargin,'norm'), 1));
for j = 1:4
    for k = 1:4
        baseFileName = jpgFiles(fileIndex).name;
        fullFileName = fullfile(Path, baseFileName);
        if  fileIndex == 1
            temp = imfinfo(fullFileName);
            mm = zeros(4, 4, temp.Height, temp.Width);
        end
        mm(j,k,:,:) = imread(fullFileName);
        if b_norm && fileIndex ~= 1
            mm(j,k,:,:) = mm(j,k,:,:)./mm(1,1,:,:);
        end
        fileIndex = fileIndex +1;
    end
end
end

function c = crrMultiprod(a,b)

sizeA = size(a);
sizeB = size(b);

% STEP 1 - Moving IDB(1) to first dimension
nd = length(sizeB);
order = [3 1:2 4:nd]; % Partial shifting
b = permute(b, order); % Q×...

% STEP 2 - Squashing B from N-D to 2-D
p = sizeA(1);
q = sizeA(2);
lengthorder = length(order);
collapsedsize = sizeB(order(2:lengthorder));
n = prod(collapsedsize);
b = reshape(b, [q, n]); % Q×N
fullsize = [p collapsedsize]; % Size to reshape C back to N-D

% FINAL STEPS - Multiplication, reshape to N-D, inverse permutation
invorder(order) = 1 : lengthorder;
c = permute (reshape(a*b, fullsize), invorder);
end
