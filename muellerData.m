classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) muellerData
    
    properties
        Label    % string
        Value    % 4,4,M,N,... array of Mueller matrix values
        ErValue  % 4,4,M,N,... array of Mueller matrix error values
        Size     % size of Value
        Dims     % cell array of length ndims(Value)-2 containing arrays of length M,N,...
        DimNames % cell array of strings with names of a dimensions M,N,...
        HV       % M,N,... array of detector high voltage values (4PEM specific)
        DC       % M,N,... array of waveform DC values (4PEM specific)
        reflection
    end
    
    methods
        function obj = muellerData(value) % Class Constructor
            obj.Size = size(value);
            obj.Value = value;
            obj.Label = '';
        end
        function varargout = subsref(obj,s)  % overload subsref for custom indexing
            switch s(1).type
                case '()'
                    if length(obj) == 1  % positional indexing of object properties
                        if length(s(1).subs) ~= length(obj.Size)
                            error('Error. Size of object and requested index are not equal');
                        end
                        if length(s) == 1
                            varargout = {objSubset(obj,s)};
                        else
                            varargout = {builtin('subsref',objSubset(obj,s(1)),s(2:end))};
                        end
                    else
                        if length(s) == 1
                            varargout = {builtin('subsref',obj,s)}; % index object array
                        else
                            obj = builtin('subsref',obj,s(1));
                            if numel(obj) == 1
                                varargout = {builtin('subsref',obj,s(2:end))};
                            else
                                temp = builtin('subsref',obj(1),s(2:end));
                                if isa(temp,'muellerData')
                                    for k=2:numel(obj)
                                        temp(k) = builtin('subsref',obj(k),s(2:end));
                                    end
                                else
                                    temp = {temp};
                                    for k=2:numel(obj)
                                        temp{k} = builtin('subsref',obj(k),s(2:end));
                                    end
                                end
                                varargout = {temp};
                            end
                        end
                    end
                    
                case '{}'
                    if length(obj) == 1
                        if length(s(1).subs) ~= length(obj.Size)
                            error('Error. Size of object and requested index are not equal');
                        end
                        if length(s) == 1
                            s = dims2index(obj,s);
                            varargout = {objSubset(obj,s)};
                        else
                            s(1) = dims2index(obj,s(1));
                            varargout = {builtin('subsref',objSubset(obj,s(1)),s(2:end))};
                        end
                    else
                        if any(arrayfun(@(x) length(s(1).subs) ~= length(x.Size),obj))
                            error('Error. Size of object and requested index are not equal');
                        end
                        if length(s) == 1
                            temp = obj;
                            for k=1:numel(obj)
                                subs = dims2index(obj(k),s);
                                temp(k) = objSubset(obj(k),subs);
                                varargout = {temp};
                            end
                        else
                            subs = dims2index(obj(1),s(1));
                            temp = builtin('subsref',objSubset(obj(1),subs),s(2:end));
                            if isa(temp,'muellerData')
                                for k=2:numel(obj)
                                    subs = dims2index(obj(k),s(1));
                                    temp(k) = builtin('subsref',objSubset(obj(k),subs),s(2:end));
                                end
                            else
                                temp = {temp};
                                for k=2:numel(obj)
                                    subs = dims2index(obj(k),s(1));
                                    temp{k} = builtin('subsref',objSubset(obj(k),subs),s(2:end));
                                end
                            end
                            varargout = {temp};
                        end
                    end
                    
                case '.'
                    if length(obj) > 1
                        temp = builtin('subsref',obj(1),s);
                        if isa(temp,'muellerData')
                            for k=2:numel(obj)
                                temp(k) = builtin('subsref',obj(k),s);
                            end
                        else
                            temp = {temp};
                            for k=2:numel(obj)
                                temp{k} = builtin('subsref',obj(k),s);
                            end
                        end
                        varargout = {temp};
                    else
                        varargout = {builtin('subsref',obj,s)};
                    end
            end
        end
        function n = numArgumentsFromSubscript(~,~,~)
            n = 1; % I don't like multiple outputs =P
        end
        function obj = merge(obj1,obj2) % merge two objects
            if ~(length(obj1.Size) == length(obj2.Size))
                error(['Objects not compatible with merge.'....
                    ' Length of obj.Size must be equal for objects.'])
            end
            if isempty(obj1.Dims) || isempty(obj2.Dims)
                error('Objects not compatible with merge. Dims must be defined.')
            end
            idx = find(cell2mat(cellfun(@isequal,obj1.Dims,obj2.Dims,'uniformoutput',0))==0);
            if length(idx) > 1 || ~isempty(intersect(obj1.Dims{idx},obj2.Dims{idx}))
                error('Objects not compatible with merge. Dims must differ in 1 element only.')
            end
            idx2 = length(obj1.Size) - length(obj1.Dims) + idx;
            obj = muellerData(cat(idx2,obj1.Value,obj2.Value));
            if ~isempty(obj1.ErValue) && ~isempty(obj2.ErValue)
                obj.ErValue = cat(idx2,obj1.ErValue,obj2.ErValue);
            end
            if ~isempty(obj1.HV) && ~isempty(obj2.HV)
                obj.HV = cat(idx,obj1.HV,obj2.HV);
            end
            if ~isempty(obj1.DC) && ~isempty(obj2.DC)
                obj.DC = cat(idx,obj1.DC,obj2.DC);
            end
            obj.Dims = obj1.Dims;
            obj.Dims{idx} = [obj1.Dims{idx} , obj2.Dims{idx}];
            obj.DimNames = obj1.DimNames;
            obj.reflection = obj1.reflection;
        end
        function obj = squeeze(obj)
            obj.Value = squeeze(obj.Value);
            obj.ErValue = squeeze(obj.ErValue);
            obj.Size = size(obj.Value);
            if ~isempty(obj.Dims)
                logicalIdx = cellfun(@(x) ~isscalar(x),obj.Dims);
                obj.Dims = obj.Dims(logicalIdx);
                if ~isempty(obj.DimNames)
                    obj.DimNames = obj.DimNames(logicalIdx);
                end
            end
            obj.HV = squeeze(obj.HV);
            obj.DC = squeeze(obj.DC);
        end
        function obj = plus(obj1,obj2) % overloading of + for muellerData.
            % to call, use: obj1 + obj2
            % Dims and DimNames and reflection are copied from obj1
            % It doesn't make sense to define HV and DC
            if isa(obj1,'muellerData') && isa(obj2,'muellerData')
                if isequal(obj1.Size,obj2.Size)
                    obj = muellerData(obj1.Value + obj2.Value);
                    obj.Dims = obj1.Dims;
                    obj.DimNames = obj1.DimNames;
                    obj.reflection = obj1.reflection;
                else
                    error('Error in obj1 + obj2 for muellerData. obj.Size must be equal for objects.')
                end
            elseif isa(obj1,'muellerData') && isscalar(obj2)
                obj = obj1;
                obj.Value = obj.Value + obj2;
            elseif isa(obj2,'muellerData') && isscalar(obj1)
                obj = obj2;
                obj.Value = obj.Value + obj1;
            end
        end
        function obj = minus(obj1,obj2) % overloading of - for muellerData.
            if isequal(obj1.Size,obj2.Size)
                obj = muellerData(obj1.Value - obj2.Value);
                obj.Dims = obj1.Dims;
                obj.DimNames = obj1.DimNames;
                obj.reflection = obj1.reflection;
            else
                error('Error in obj1 - obj2 for muellerData. obj.Size must be equal for objects.')
            end
        end
        function obj = times(obj1,obj2) % overloading of .* for muellerData.
            if isequal(obj1.Size,obj2.Size)
                obj = muellerData(obj1.Value .* obj2.Value);
                obj.Dims = obj1.Dims;
                obj.DimNames = obj1.DimNames;
                obj.reflection = obj1.reflection;
            else
                error('Error in obj1 .* obj2 for muellerData. obj.Size must be equal for objects.')
            end
        end
        function obj = rdivide(obj1,obj2) % overloading of ./ for muellerData.
            if isequal(obj1.Size,obj2.Size)
                obj = muellerData(obj1.Value ./ obj2.Value);
                obj.Dims = obj1.Dims;
                obj.DimNames = obj1.DimNames;
                obj.reflection = obj1.reflection;
            else
                error('Error in obj1 ./ obj2 for muellerData. obj.Size must be equal for objects.')
            end
        end
        function obj = mtimes(obj1,obj2) % overloading of * for muellerData.
            ck1 = isa(obj1, 'muellerData');
            ck2 = isa(obj2, 'muellerData');
            if ck1 && ck2
                if ndims(obj2.Value) > ndims(obj1.Value)
                    obj = obj2;
                else
                    obj = obj1;
                end
                    obj.Value = multiprod(obj1.Value, obj2.Value, [1 2], [1 2]);
            elseif ck1
                obj = obj1;
                obj.Value = multiprod(obj1.Value, obj2, [1 2], [1 2]);
            else
                obj = obj2;
                obj.Value = multiprod(obj1, obj2.Value, [1 2], [1 2]);
            end
%             if isequal(obj1.Size,obj2.Size)
%                 val1 = shapeDown(obj1.Value);
%                 val2 = shapeDown(obj2.Value);
%                 for i=1:size(val1,3); val1(:,:,i) = val1(:,:,i)*val2(:,:,i); end
%                 obj = muellerData(shapeUp(val1,obj1.Size));
%                 obj.Dims = obj1.Dims;
%                 obj.DimNames = obj1.DimNames;
%                 obj.reflection = obj1.reflection;
%             else
%                 error('Error in obj1 ./ obj2 for muellerData. obj.Size must be equal for objects.')
%             end
        end
        function obj = mrdivide(obj1,obj2) % overloading of / for muellerData.
            if isequal(obj1.Size,obj2.Size)
                val1 = shapeDown(obj1.Value);
                val2 = shapeDown(obj2.Value);
                for i=1:size(val1,3); val1(:,:,i) = val1(:,:,i)/val2(:,:,i); end
                obj = muellerData(shapeUp(val1,obj1.Size));
                obj.Dims = obj1.Dims;
                obj.DimNames = obj1.DimNames;
                obj.reflection = obj1.reflection;
            else
                error('Error in obj1 ./ obj2 for muellerData. obj.Size must be equal for objects.')
            end
        end
        function obj = mldivide(obj1,obj2) % overloading of \ for muellerData.
            if isequal(obj1.Size,obj2.Size)
                val1 = shapeDown(obj1.Value);
                val2 = shapeDown(obj2.Value);
                for i=1:size(val1,3); val1(:,:,i) = val1(:,:,i) \ val2(:,:,i); end
                obj = muellerData(shapeUp(val1,obj1.Size));
                obj.Dims = obj1.Dims;
                obj.DimNames = obj1.DimNames;
                obj.reflection = obj1.reflection;
            else
                error('Error in obj1 ./ obj2 for muellerData. obj.Size must be equal for objects.')
            end
        end
        function handles = plot(varargin)
            handles = prePlot(varargin{:});
        end
        function handles = subplot(varargin)
            % Example: % obj.subplot( {'lb','lbp','cb';'ld','ldp','cd'} , 'legend','none' )
            [obj,funcs] = varargin{:};
            figure
            M = size(funcs,1);
            N = size(funcs,2);
            funcs = funcs(:);
            handles = gobjects(1,M*N);
            for idx=1:M*N
                ax = subplot(M,N,idx);
                fn = str2func(funcs{idx});
                handles(idx) = plot(fn(obj),'handle',ax,varargin{3:end},...
                    'title',[', ',upper(funcs{idx})]);
            end
        end
        function handles = print(varargin)
            filePath = varargin{2}; % extract the filepath
            [pathStr,name] = fileparts(filePath);
            filePath = [pathStr,'/',varargin{1}.Label,name];
            handles = prePlot(varargin{[1,3:end]}); % make the figure
            print(gcf,filePath,'-depsc'); % print figure as .eps file
        end
        % Calls to static methods on obj.Value, returns new class instance %
        function obj = optProp(obj)
            obj.Value = obj.s_optProp(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lm(varargin)
            obj = varargin{1};
            if nargin == 1
                obj.Value = obj.s_lm(obj.Value);
            else
                obj.Value = obj.s_lm(obj.Value,varargin{2});
            end
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = logm(obj)
            obj.Value = obj.s_logm(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lu(obj)
            obj = obj.logm;
            g = diag([-1 1 1 1]);
            for n=1:size(obj.Value,3)
                obj.Value(:,:,n) = (obj.Value(:,:,n) + g*obj.Value(:,:,n).'*g)/2;
            end
        end
        function obj = lm2(obj)
            obj = obj.logm;
            g = diag([-1 1 1 1]);
            for n=1:size(obj.Value,3)
                obj.Value(:,:,n) = (obj.Value(:,:,n) - g*obj.Value(:,:,n).'*g)/2;
            end
        end
        function obj = expm(obj)
            obj.Value = obj.s_expm(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lb(obj)
            obj.Value = obj.s_lb(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = ld(obj)
            obj.Value = obj.s_ld(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lbp(obj)
            obj.Value = obj.s_lbp(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = ldp(obj)
            obj.Value = obj.s_ldp(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = cb(obj)
            obj.Value = obj.s_cb(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = cd(obj)
            obj.Value = obj.s_cd(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = a(obj)
            obj.Value = obj.s_a(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = a_aniso(obj)
            obj.Value = obj.s_a_aniso(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = a_iso(obj)
            obj.Value = obj.s_a_iso(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = ldmag(obj)
            obj.Value = obj.s_ldmag(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = ldang(obj)
            obj.Value = obj.s_ldang(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lbang(obj)
            obj.Value = obj.s_lbang(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lbmag(obj)
            obj.Value = obj.s_lbmag(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = di(obj)
            obj.Value = obj.s_di(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = jones(obj)
            obj.Value = obj.s_jones(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = nearestjones(obj)
            obj.Value = obj.s_nearestjones(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = mfilter(obj)
            obj.Value = obj.s_mfilter(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = covar(obj)
            obj.Value = obj.s_covar(obj.Value);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = mrotate(obj,angle_rad)
            obj.Value = obj.s_mrotate(obj.Value,angle_rad);
            obj.ErValue = [];
            obj.Size = size(obj.Value);
        end
        function obj = lm2optProp(obj)
            % [LB;LD;LBp;LDp;CB;CD;A]
            lm = obj.Value;
            sz = size(lm);
            lm = shapeDown(lm);
            val(1,:) = lm(4,3,:);
            val(2,:) = -lm(1,2,:);
            val(3,:) = lm(2,4,:);
            val(4,:) = -lm(1,3,:);
            val(5,:) = lm(2,3,:);
            val(6,:) = lm(1,4,:);
            val(7,:) = -lm(1,1,:);
            obj.Value = shapeUp(val, sz);
        end
    end
    
    methods(Static)
        % value = obj.Value
        function r = s_optProp(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
            T = acos( K.*( J(1,1,:) + J(2,2,:) )./2); % 2*T = sqrt(L.^2 + Lp.^2 + C.^2)
            O = (T.*K)./(sin(T));
            L=1i.*O.*( J(1,1,:) - J(2,2,:) );
            Lp=1i.*O.*( J(1,2,:) + J(2,1,:) );
            C=O.*( J(1,2,:) - J(2,1,:) );
            LB=real(L);
            LD=-imag(L);
            LBp=real(Lp);
            LDp=-imag(Lp);
            CB=real(C);
            CD=-imag(C);
            A = -2*real(log(1./K)); % mean absorption
            r =  shapeUp(squeeze([LB;LD;LBp;LDp;CB;CD;A]),sz);
        end
        function value = s_lm(varargin)
            value = varargin{1};
            sz = size(value);
            if nargin == 1
                value = shapeDown(value);
                %J = nearestJones(value);
                J = MJ2J(value);
                K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
                T = acos( K.*( J(1,1,:) + J(2,2,:) )./2);
                O = (T.*K)./(sin(T));
                L=1i.*O.*( J(1,1,:) - J(2,2,:) );
                Lp=1i.*O.*( J(1,2,:) + J(2,1,:) );
                C=O.*( J(1,2,:) - J(2,1,:) );
                LB=real(L);
                LD=-imag(L);
                LBp=real(Lp);
                LDp=-imag(Lp);
                CB=real(C);
                CD=-imag(C);
                A = 2*real(log(1./K)); % mean absorption
                value =  shapeUp([A,-LD,-LDp,CD ; -LD,A,CB,LBp ; -LDp,-CB,A,-LB ; CD,-LBp,LB,A],sz);
            else
                n_int = varargin{2};
                value = reshape(value,4,4,size(value,3),[]);
                for j = 1:size(value,4)
                    M = value(:,:,:,j);
                    M = flip(M,3);
                    J = nearestJones(M);
                    K=(J(1,1,1).*J(2,2,1) - J(1,2,1)*J(2,1,1)).^(-1/2);
                    T=2*acos((K.*(J(1,1,1) + J(2,2,1)))./2);
                    O=(T+2*pi*n_int).*K./(sin(T/2)*2);
                    
                    N = size(J,3);
                    L = zeros(1,N);
                    Lp = zeros(1,N);
                    C = zeros(1,N);
                    A = zeros(1,N);
                    
                    L(1) = 1i.*O.*(J(1,1,1) - J(2,2,1));
                    Lp(1) = 1i.*O.*(J(1,2,1) + J(2,1,1));
                    C(1) = O.*(J(1,2,1) - J(2,1,1));
                    A(1) = 2*real(log(1./K));
                    
                    n = n_int;
                    
                    for i = 2:N
                        if n==0 || n==-1
                            n_ar = [0,-1,1,-2,2];
                        else
                            n_ar = [n-1,-n,n,-(n+1),n+1];
                        end
                        K=(J(1,1,i).*J(2,2,i) - J(1,2,i)*J(2,1,i)).^(-1/2);
                        T=2*acos((K.*(J(1,1,i) + J(2,2,i)))./2);
                        O=(T+2*pi*n_ar).*K./(sin(T/2)*2);
                        l = 1i.*O.*(J(1,1,i) - J(2,2,i));
                        lp = 1i.*O.*(J(1,2,i) + J(2,1,i));
                        c = O.*(J(1,2,i) - J(2,1,i));
                        diffs = sum([L(i-1)-l;Lp(i-1)-lp;C(i-1)-c],1);
                        [~,I] = min(diffs);
                        L(i) = l(I);
                        Lp(i) = lp(I);
                        C(i) = c(I);
                        n = n_ar(I);
                        A(i) = 2*real(log(1./K));
                    end
                    
                    LB=reshape(real(L),1,1,[]);
                    LD=reshape(-imag(L),1,1,[]);
                    LBp=reshape(real(Lp),1,1,[]);
                    LDp=reshape(-imag(Lp),1,1,[]);
                    CB=reshape(real(C),1,1,[]);
                    CD=reshape(-imag(C),1,1,[]);
                    A = reshape(A,1,1,[]);
                    value(:,:,:,j) =  ...
                        flip([A,-LD,-LDp,CD ; -LD,A,CB,LBp ; -LDp,-CB,A,-LB ; CD,-LBp,LB,A],3);
                end
                value = reshape(value,sz);
            end
        end
        function r = s_logm(value)  % log of Mueller matrix with filtering
            sz = size(value);
            value = shapeDown(value);
            Mfiltered = filterM(value);
            r =  shapeUp(zeros(size(Mfiltered)),sz);
            for n=1:size(value,3); r(:,:,n) = logm(Mfiltered(:,:,n)); end
        end
        function r = s_expm(r)  % log of Mueller matrix with filtering
            sz = size(r);
            r = shapeDown(r);
            for n=1:size(r,3); r(:,:,n) = expm(r(:,:,n)); end
            r = shapeUp(r,sz);
        end
        function r = s_lb(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = real(1i.*r.*( J(1,1,:) - J(2,2,:) ));
            r =  shapeUp(r,sz);
        end % 0,90 linear retardance
        function r = s_ld(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = -imag(1i.*r.*( J(1,1,:) - J(2,2,:) ));
            r =  shapeUp(r,sz);
        end  % 0,90 linear extinction
        function r = s_lbp(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = real(1i.*r.*( J(1,2,:) + J(2,1,:) ));
            r =  shapeUp(r,sz);
        end  % 45,-45 linear retardance
        function r = s_ldp(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = -imag(1i.*r.*( J(1,2,:) + J(2,1,:) ));
            r =  shapeUp(r,sz);
        end  % 45,-45 linear extinction
        function r = s_cb(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = real(r.*( J(1,2,:) - J(2,1,:) ));
            r =  shapeUp(r,sz);
        end  % circular retardance
        function r = s_cd(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = jonesAnisotropy(J);
            r = -imag(r.*( J(1,2,:) - J(2,1,:) ));
            r =  shapeUp(r,sz);
        end  % circular extinction
        function r = s_a(value) % total mean extinction
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            r = -2*real(log( ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(1/2) ));
            r =  shapeUp(r,sz);
        end
        function r = s_a_aniso(value)  % anisotropic part of the mean extinction
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
            T = acos( K.*( J(1,1,:) + J(2,2,:) )./2); % 2*T = sqrt(L.^2 + Lp.^2 + C.^2)
            O = (T.*K)./(sin(T));
            LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            CD = -imag(O.*( J(1,2,:) - J(2,1,:) ));
            r =  shapeUp(sqrt(LD.^2 + LDp.^2 + CD.^2),sz);  % not same as imag(2*T) !
        end
        function r = s_a_iso(value) % isotropic part of the mean extinction
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
            T = acos( K.*( J(1,1,:) + J(2,2,:) )./2); % 2*T = sqrt(L.^2 + Lp.^2 + C.^2)
            O = (T.*K)./(sin(T));
            LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            CD = -imag(O.*( J(1,2,:) - J(2,1,:) ));
            r =  shapeUp(-2*real(log(1./K)) - sqrt(LD.^2 + LDp.^2 + CD.^2),sz);
        end
        function r = s_ldmag(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            O = jonesAnisotropy(J);
            LD = imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LDp = imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            r =  shapeUp(sqrt(LD.^2 + LDp.^2),sz);
        end
        function r = s_ldang(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            O = jonesAnisotropy(J);
            LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            r =  shapeUp(atan2(LDp , LD)./2,sz);
            %out = out + pi*(out < 0);
        end
        function r = s_lbang(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            O = jonesAnisotropy(J);
            LB = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LBp = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            r =  atan2(LBp , LB)./2;
            r =  shapeUp(r + pi*(r < 0),sz);
        end
        function r = s_lbmag(value)
            sz = size(value);
            value = shapeDown(value);
            J = nearestJones(value);
            O = jonesAnisotropy(J);
            LB = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
            LBp = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
            r =  shapeUp(sqrt(LB.^2 + LBp.^2),sz);
        end
        function r = s_di(value) % Depolarization Index
            sz = size(value);
            value = shapeDown(value);
            r =  shapeUp((sqrt(squeeze(sum(sum(value.^2,1),2))./squeeze(value(1,1,:)).^2-1)./sqrt(3)).',sz);
        end
        function r = s_jones(value) % Jones matrix of a Mueller-Jones matrix
            sz = size(value);
            value = shapeDown(value);
            r =  shapeUp(MJ2J(value),sz);
        end
        function r = s_nearestjones(value)
            sz = size(value);
            value = shapeDown(value);
            r =  nearestJones(value);  % Jones matrix
            % next line just phases the Jones matrix so that the
            % imaginary part of J(1,1) = 0. i.e., it matches case 'jones'
            for n=1:size(r,3); r(:,:,n) = exp( -1i*angle(r(1,1,n)) ) * r(:,:,n); end
            r = shapeUp(r,sz);
        end
        function r = s_mfilter(value) % closest physical Mueller matrix
            sz = size(value);
            value = shapeDown(value);
            r =  shapeUp(filterM(value),sz);
        end
        function r = s_covar(value) % Mueller to Cloude covariance
            sz = size(value);
            value = shapeDown(value);
            r =  shapeUp(M2Cov(value),sz);
        end
        function r = plotter(varargin)
            r = linePlot(varargin{:});
        end
        function r = s_mrotate(M,theta)
            % M is a Mueller matrix array of any dimension. The first two dimension
            % must be the Mueller matrix elements. MMout is a Mueller array with the
            % same dimension as the input array.
            
            % October 17, 2016: sign of theta changed so +LB transforms to +LB' with
            % theta = pi/4.
            sz = size(M);
            M = shapeDown(M);
            r = M;
            theta=-2*theta;
            C2=cos(theta);
            S2=sin(theta);
            r(1,2,:) = M(1,2,:)*C2 + M(1,3,:)*S2;
            r(1,3,:) = M(1,3,:)*C2 - M(1,2,:)*S2;
            r(2,1,:) = M(2,1,:)*C2 + M(3,1,:)*S2;
            r(3,1,:) = M(3,1,:)*C2 - M(2,1,:)*S2;
            r(2,4,:) = M(2,4,:)*C2 + M(3,4,:)*S2;
            r(3,4,:) = M(3,4,:)*C2 - M(2,4,:)*S2;
            r(4,2,:) = M(4,2,:)*C2 + M(4,3,:)*S2;
            r(4,3,:) = M(4,3,:)*C2 - M(4,2,:)*S2;
            r(2,2,:) = C2*(M(3,2,:)*S2 + M(2,2,:)*C2) + S2*(M(3,3,:)*S2 + M(2,3,:)*C2);
            r(2,3,:) = C2*(M(3,3,:)*S2 + M(2,3,:)*C2) - S2*(M(3,2,:)*S2 + M(2,2,:)*C2);
            r(3,2,:) = -C2*(M(2,2,:)*S2 - M(3,2,:)*C2) - S2*(M(2,3,:)*S2 - M(3,3,:)*C2);
            r(3,3,:) = S2*(M(2,2,:)*S2 - M(3,2,:)*C2) - C2*(M(2,3,:)*S2 - M(3,3,:)*C2);
            r = shapeUp(r,sz);
        end
        function fig = mergeAxes(h,sz)
            h = h(:);
            set(h,'Units','Pixels');
            p = get(h,'Position');
            ti = get(h,'TightInset');
            extents = ...
                cellfun(@(p,ti) [ti(1) + ti(3) + p(3) , ti(2) + ti(4) + p(4)],p,ti,'uniformoutput',0);
            extents = max(cell2mat(extents));
            [I,J] = ind2sub(sz,1:length(h));
            hspace = 10;
            vspace = 10;
            figSz = (flip(sz)).*[hspace,vspace] + flip(sz).*extents ;
            
            fig =  figure('Units','Pixels','Position',[0, 0, figSz(1), figSz(2)] );
            for i=1:length(h)
                os1 = p{i}(1) - ti{i}(1);
                os2 = p{i}(2) - ti{i}(2);
                obj = h(i).Parent.Children;
                set(obj,'Units','Pixels');
                pos = get(obj,'Position');
                obj = copyobj(obj,fig);
                if length(obj) == 1
                    pos = pos + [J(i) * hspace + (J(i) - 1) * extents(1) - os1 ,...
                        (sz(1)-I(i)) * vspace + (sz(1)-I(i)) * extents(2) - os2 ,...
                        0,0];
                    obj.Position = pos;
                else
                    for j=1:length(obj)
                        temp = pos{j} + ...
                            [(J(i)-1) * hspace + (J(i) - 1) * extents(1) - os1 ,...
                            (sz(1)-I(i)) * vspace + (sz(1)-I(i)) * extents(2) - os2 ,...
                            0,0];
                        obj(j).Position = temp;
                    end
                end
            end
        end
    end
end

% LOCAL FUNCTIONS
% =========================================================================
function s = dims2index(obj,s) % for indexing with Dims
if isempty(obj.Dims)
    error('Error. obj.Dims not defined.');
end
sz = length(s.subs) - length(obj.Dims);
for i=1:length(obj.Dims)
    if s.subs{i+sz} ~= ':'
        [X,I] = sort(obj.Dims{i}); % added this to allow unsorted Dims
        indices = unique(round(fracIndex(X,s.subs{i+sz})),'first');
        s.subs{i+sz} = I(indices);
    end
end
end
function obj = objSubset(obj,s) % obj parsing
obj.Value = obj.Value(s.subs{:});
obj.Size = size(obj.Value);
if ~isempty(obj.ErValue)
    obj.ErValue = obj.ErValue(s.subs{:});
end
obj.DimNames = obj.DimNames;
lsubs = length(s.subs) + 1;
if ~isempty(obj.HV)
    obj.HV = obj.HV(s.subs{(lsubs-sum(size(obj.HV) ~= 1)):end});
end
if ~isempty(obj.DC)
    obj.DC = obj.DC(s.subs{(lsubs-sum(size(obj.DC) ~= 1)):end});
end
if ~isempty(obj.Dims)
    sz = lsubs - length(obj.Dims) - 1;
    for i=1:length(obj.Dims)
        obj.Dims{i} = obj.Dims{i}(s.subs{i+sz});
    end
end
end
function out = shapeDown(out)
if ndims(out) > 3   % reshape array into 4,4,N
    out = reshape(out,4,4,[]);
end
end % reshape
function out = shapeUp(out,sz) % overly complicated reshaping
sz2 = size(out);
if  length(sz)>=3  %  reshape to match input dimensions
    out = reshape(out,[sz2(1:(length(sz2)-1)),sz(3:length(sz))]);
end
sz2 = size(out);
if sz2(1) == 1 % remove leading singletons if necessary
    if sz2(2) == 1
        out = shiftdim(out,2); % out = reshape(out,sz2(3:end));
    else
        out = shiftdim(out,1); %out = reshape(out,sz2(2:end));
    end
end
end
function J = MJ2J(M)  % Mueller-Jones to Jones
J(1,1,:) = ((M(1,1,:)+M(1,2,:)+M(2,1,:)+M(2,2,:))/2).^(1/2);
k = 1./(2.*J(1,1,:));
J(1,2,:) = k.*(M(1,3,:)+M(2,3,:)-1i.*(M(1,4,:)+M(2,4,:)));
J(2,1,:) = k.*(M(3,1,:)+M(3,2,:)+1i.*(M(4,1,:)+M(4,2,:)));
J(2,2,:) = k.*(M(3,3,:)+M(4,4,:)+1i.*(M(4,3,:)-M(3,4,:)));
end
function C = M2Cov(M) % Mueller to Cloude covariance
C(1,1,:) = M(1,1,:) + M(1,2,:) + M(2,1,:) + M(2,2,:);
C(1,2,:) = M(1,3,:) + M(1,4,:)*1i + M(2,3,:) + M(2,4,:)*1i;
C(1,3,:) = M(3,1,:) + M(3,2,:) - M(4,1,:)*1i - M(4,2,:)*1i;
C(1,4,:) = M(3,3,:) + M(3,4,:)*1i - M(4,3,:)*1i + M(4,4,:);
C(2,1,:) = M(1,3,:) - M(1,4,:)*1i + M(2,3,:) - M(2,4,:)*1i;
C(2,2,:) = M(1,1,:) - M(1,2,:) + M(2,1,:) - M(2,2,:);
C(2,3,:) = M(3,3,:) - M(3,4,:)*1i - M(4,3,:)*1i - M(4,4,:);
C(2,4,:) = M(3,1,:) - M(3,2,:) - M(4,1,:)*1i + M(4,2,:)*1i;
C(3,1,:) = M(3,1,:) + M(3,2,:) + M(4,1,:)*1i + M(4,2,:)*1i;
C(3,2,:) = M(3,3,:) + M(3,4,:)*1i + M(4,3,:)*1i - M(4,4,:);
C(3,3,:) = M(1,1,:) + M(1,2,:) - M(2,1,:) - M(2,2,:);
C(3,4,:) = M(1,3,:) + M(1,4,:)*1i - M(2,3,:) - M(2,4,:)*1i;
C(4,1,:) = M(3,3,:) - M(3,4,:)*1i + M(4,3,:)*1i + M(4,4,:);
C(4,2,:) = M(3,1,:) - M(3,2,:) + M(4,1,:)*1i - M(4,2,:)*1i;
C(4,3,:) = M(1,3,:) - M(1,4,:)*1i - M(2,3,:) + M(2,4,:)*1i;
C(4,4,:) = M(1,1,:) - M(1,2,:) - M(2,1,:) + M(2,2,:);
C = C./2;
end
function M = Cov2M(C) % Cloude covariance to Mueller
M(1,1,:) = C(1,1,:) + C(2,2,:) + C(3,3,:) + C(4,4,:);
M(1,2,:) = C(1,1,:) - C(2,2,:) + C(3,3,:) - C(4,4,:);
M(1,3,:) = C(1,2,:) + C(2,1,:) + C(3,4,:) + C(4,3,:);
M(1,4,:) = ( -C(1,2,:) + C(2,1,:) - C(3,4,:) + C(4,3,:) )*1i;
M(2,1,:) = C(1,1,:) + C(2,2,:) - C(3,3,:) - C(4,4,:);
M(2,2,:) = C(1,1,:) - C(2,2,:) - C(3,3,:) + C(4,4,:);
M(2,3,:) = C(1,2,:) + C(2,1,:) - C(3,4,:) - C(4,3,:);
M(2,4,:) = ( -C(1,2,:) + C(2,1,:) + C(3,4,:) - C(4,3,:) )*1i;
M(3,1,:) = C(1,3,:) + C(2,4,:) + C(3,1,:) + C(4,2,:);
M(3,2,:) = C(1,3,:) - C(2,4,:) + C(3,1,:) - C(4,2,:);
M(3,3,:) = C(1,4,:) + C(2,3,:) + C(3,2,:) + C(4,1,:);
M(3,4,:) = ( -C(1,4,:) + C(2,3,:) - C(3,2,:) + C(4,1,:) )*1i;
M(4,1,:) = ( C(1,3,:) + C(2,4,:) - C(3,1,:) - C(4,2,:) )*1i;
M(4,2,:) = ( C(1,3,:) - C(2,4,:) - C(3,1,:) + C(4,2,:) )*1i;
M(4,3,:) = ( C(1,4,:) + C(2,3,:) - C(3,2,:) - C(4,1,:) )*1i;
M(4,4,:) = C(1,4,:) - C(2,3,:) - C(3,2,:) + C(4,1,:);
M = real(M)./2;
end
function J = nearestJones(M)
C = M2Cov(M);
J = zeros(2,2,size(C,3));
for n=1:size(C,3)
    [V,D] = eig(C(:,:,n),'vector');
    [~,mx] = max(D);
    J(:,:,n) = sqrt(D(mx))*reshape(V(:,mx),2,2).';
end
end
function M = filterM(M)  % M to nearest physical M
C_raw = M2Cov(M);
C = zeros(size(C_raw));
for n=1:size(C_raw,3)
    [V,D] = eig(C_raw(:,:,n),'vector');
    list = find(D > 0.00001).';
    idx = 0;
    temp = zeros(4,4,length(list));
    for j = list
        idx = idx + 1;
        temp(:,:,idx) = D(j)*V(:,j)*V(:,j)';
    end
    C(:,:,n) = sum(temp,3);
end
M = Cov2M(C);
end
function O = jonesAnisotropy(J)
K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
T = acos( K.*( J(1,1,:) + J(2,2,:) )./2);
O = (T.*K)./(sin(T));
end
function fracIndx = fracIndex(X,y) %fractional index
% X: 1xN array of increasing values
% y: array of values in the range of X
% fracIndx is an array the length of y that contains the fractional
% index of the y values in array X.
% e.g., X = [2,4,6]; y = [4,5]; gives, fracIndx = [2,2.5];
fracIndx = zeros(1,length(y));
for idx = 1:length(y)
    if y(idx) >= X(length(X))
        fracIndx(idx) = length(X);
    elseif y(idx) <= X(1)
        fracIndx(idx) = 1;
    else
        a = find(X <= y(idx));
        a = a(length(a));
        b = find(X > y(idx));
        b = b(1);
        fracIndx(idx) = a+(y(idx)-X(a))/(X(b)-X(a));
    end
end
end
function handles = prePlot(varargin)
obj = varargin{1};
if all(obj.Size(1:2) == 4)
    plotTool = @MMplot;
else
    plotTool = @linePlot;
end
if ~isempty(obj.Label)
    if any(strcmpi('title',varargin))
        idx = find(strcmpi('title',varargin)) + 1;
        varargin{idx} = [obj.Label, ' ',varargin{idx}];
    else
        sz = length(varargin);
        varargin{sz+1} = 'title';
        varargin{sz+2} = obj.Label;
    end
end
if ~any(strcmpi('legend',varargin))
    if length(obj.Dims) >= 2 && ~isempty(obj.Dims{2})
        if length(obj.Dims) >= 3 && ~isempty(obj.Dims{3})
            idx = 1;
            Labels = cell(1,length(obj.Dims{2})*length(obj.Dims{3}));
            for i=1:length(obj.Dims{2})
                for j=1:length(obj.Dims{3})
                    Labels{idx} = [num2str(obj.Dims{2}(i)),' ; ',num2str(obj.Dims{3}(j))];
                    idx = idx + 1;
                end
            end
            LabelNames = [obj.DimNames{2},' ; ',obj.DimNames{3}];
        else
            Labels = obj.Dims{2};
            LabelNames = obj.DimNames{2};
        end
        sz = length(varargin);
        varargin{sz+1} = 'legend';
        varargin{sz+2} = {LabelNames,Labels};
    end
end
handles = plotTool(obj.Dims{1},obj.Value,obj.ErValue,varargin{2:end});
end
function handles = MMplot(Lam,MMdata,MMerror,varargin)
% Mueller matrix 2D plotting utility
% Makes a 4 x 4 array of 2-D line plots with full control over line and
% axes properties.
% Outputs: [1 x 16] array of axis handles
%
% Required positional inputs:
%   Lam: [1 x n] array of wavelengths (X-axis)
%   MMdata: [4 x 4 x n x ...] Mueller matrix array
% Optional positional inputs:
%   LineSpec: string containing a valid lineSpec. Type "doc LineSpec" in
%       command window for more info. Default is "-", a solid line.
% Optional Name-Value pairs inputs:
%   ev: bool. converts X axis to eV. e.g., 'ev',true
%   handles: [1 x 16] array of plot handles. New handles are created if not given.
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
%   multiple of the largest line width on the plots.

p = inputParser;
% input validation functions
valFun1 = @(x) ischar(x) && ...
    all(~strcmpi(x,{'ev','handles','lineNV','limY','fontsize','axNV','size',...
    'title','legend','vSpace','borderFactor'}));
valFun2 = @(x) isscalar(x)&&isnumeric(x);
% setup input scheme
addRequired(p,'Lam',@isnumeric);
addRequired(p,'MMdata',@isnumeric);
addRequired(p,'MMerror',@isnumeric);
addOptional(p,'LineSpec','-',valFun1)
addParameter(p,'ev',false,@islogical)
addParameter(p,'handles',gobjects(1,16), @(x) all(ishandle(x)))
addParameter(p,'limY',0,valFun2)
addParameter(p,'fontsize',12,valFun2)
addParameter(p,'axNV',{},@iscell)
addParameter(p,'lineNV',{},@iscell)
addParameter(p,'size',[1000 700],@(x) length(x) == 2 && isnumeric(x))
addParameter(p,'title','',@ischar)
addParameter(p,'legend',{},@(x) iscell(x) || strcmp(x,'none'))
addParameter(p,'vSpace',0,@isscalar)
addParameter(p,'borderFactor',0,@isscalar)
parse(p,Lam,MMdata,MMerror,varargin{:}) %parse inputs

% create new figure if no valid handles were given
handles = p.Results.handles;
if any(strcmpi('handles',p.UsingDefaults))
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
    xLabel = uicontrol('style','text','BackgroundColor','w',...
        'units','pixels','FontSize',p.Results.fontsize,...
        'tag','xLabelObject'); % create x-label
    if p.Results.ev == true
        set(xLabel,'String','Energy (eV)');
    else
        set(xLabel,'String','Wavelength (nm)');
    end
    xLabel_sz = get(xLabel,'extent');
    set(xLabel,'Position',[(figPos(3) - xLabel_sz(3) )./2, 0, xLabel_sz(3), xLabel_sz(4)]);
    
    if ~isempty(p.Results.title) % create title if given
        figTitle = uicontrol('style','text','BackgroundColor','w',...
            'units','pixels','FontSize',p.Results.fontsize,...
            'tag','titleObject');
        set(figTitle,'String',p.Results.title)
        figTitle_sz = get(figTitle,'extent');
        set(figTitle,'Position',[( figPos(3) - figTitle_sz(3) )./2,...
            ( figPos(4) - figTitle_sz(4) ), figTitle_sz(3), figTitle_sz(4)]);
    end
    % determine the horizontal extent of y-axis marker labels
    dummy = uicontrol('style','text','fontsize',p.Results.fontsize,'units','pixels');
    set(dummy,'String','-0.000');
    yAxSz = get(dummy,'extent');
    delete(dummy)
    
    plotSzX = figPos(3)/4 - yAxSz(3) - yAxSz(3)./5; % X size of plot area in pixels
    plotSzY = ( figPos(4) - 4*yAxSz(4) )/4 - 6 - p.Results.vSpace; % Y size of plot area in pixels
    for i=1:4
        for j=1:4
            plotPos = [ ( (plotSzX + yAxSz(3) + 3)*(j-1) + yAxSz(3) +5)./figPos(3) ,...
                ((plotSzY + yAxSz(4)./2 + p.Results.vSpace)*(4-i)+yAxSz(4)*2 + 3)./figPos(4),...
                plotSzX./figPos(3), plotSzY./figPos(4)];
            hand = subplot('Position',plotPos);
            hold(hand,'on')
            box(hand,'on')
            if i ~= 4
                set(hand,'XTickLabel',[]) % keep X lables only for bottom row
            end
            handles(j+4*(i-1)) = hand;
        end
    end
else
    h_fig = get(handles(1),'parent');
    figPos = get(h_fig,'Position');
end

%plot data and set Line properties.
if p.Results.ev == true; Lam = 1239.8./Lam; end
if isempty(MMerror)
    for j = 1:4
        for k = 1:4
            plot(handles(k+4*(j-1)),Lam,squeeze(MMdata(j,k,:,:)),...
                p.Results.LineSpec,p.Results.lineNV{:})
        end
    end
else
    for j = 1:4
        for k = 1:4
            errorbar(handles(k+4*(j-1)),Lam,squeeze(MMdata(j,k,:,:)),...
                squeeze(MMerror(j,k,:,:)),...
                p.Results.LineSpec,'CapSize',0,p.Results.lineNV{:})
        end
    end
end
% set Axes properties
axis(handles,'tight'); % first, axes are set to tight
if ~isempty(p.Results.axNV)
    for j=1:16; set(handles(j),p.Results.axNV{:}); end
end
if p.Results.limY ~= 0 % modify axes bounds if limY is set
    lim = p.Results.limY;
    for j=1:16
        Ylim = get(handles(j),'YLim');
        if (Ylim(2) - Ylim(1)) < lim
            avg = (Ylim(2) + Ylim(1))./2;
            Ylim(2) = avg + lim/2;
            Ylim(1) = avg - lim/2;
            set(handles(j),'Ylim',Ylim);
        end
    end
end
% Adjust plot limits so that lines do not overlap axis borders.
% *** If you like to use Markers, then perhaps change 'lineWidth' to 'MarkerSize'
lineHandle = get(handles(1),'children');
lineWidth = zeros(size(lineHandle));
for j = 1:length(lineHandle)
    lineWidth(j) = get(lineHandle(j),'lineWidth');
end
lineWidth = max(lineWidth)*p.Results.borderFactor;
plotPos = get(handles(1),'Position');
for j=1:16
    xlim = get(handles(j),'xLim');
    ylim = get(handles(j),'yLim');
    xStep = (xlim(2) - xlim(1))/plotPos(3)/figPos(3)*lineWidth/2;
    yStep = (ylim(2) - ylim(1))/plotPos(4)/figPos(3)*lineWidth;
    set(handles(j),'XLim',[xlim(1)-xStep,xlim(2)+xStep]);
    set(handles(j),'YLim',[ylim(1)-yStep,ylim(2)+yStep]);
end
% set font size of all graphics objects if fontsize was passed
if ~any(strcmpi('fontsize',p.UsingDefaults))
    set(get(gcf,'children'),'FontSize',p.Results.fontsize);
end
% optionally create legend (this will increase the width of the figure!)
if ~any(strcmpi('legend',p.UsingDefaults))
    if iscell(p.Results.legend)
        Labels = p.Results.legend{2};
        if isnumeric(Labels)
            Labels = cellfun(@(x) num2str(x),num2cell(Labels),'uniformoutput',0);
        end
        pos = zeros(4,16);
        for i=1:16
            set(handles(i),'units','pixels');
            pos(:,i) = get(handles(i),'Position');
        end
        lgd = legend(handles(4),Labels,'location','northeastoutside');
        set(lgd,'units','pixels','fontsize',p.Results.fontsize);
        title(lgd,p.Results.legend{1},'FontSize',p.Results.fontsize);
        lgd_pos = get(lgd,'Position');
        h_fig.Position = h_fig.Position + [0 0 lgd_pos(3) 0];
        for i=1:16
            set(handles(i),'Position',pos(:,i));
        end
    end
end

end
function handle = linePlot(X,Y,YEr,varargin)
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
    {'handle','lineNV','limY','fontsize','axNV','size','title','xLabel',...
    'yLabel','legend','legendLocation'}));
valFun2 = @(x) isscalar(x)&&isnumeric(x);
% setup input scheme
addRequired(p,'X',@isnumeric);
addRequired(p,'Y',@isnumeric);
addRequired(p,'YEr',@isnumeric);
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
addParameter(p,'legend',{},@(x) iscell(x) || strcmp(x,'none'))
addParameter(p,'legendLocation','northeastoutside',@ischar)
parse(p,X,Y,YEr,varargin{:}) %parse inputs

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
    if iscell(p.Results.legend)
        Labels = p.Results.legend{2};
        if isnumeric(Labels)
            Labels = cellfun(@(x) num2str(x),num2cell(Labels),'uniformoutput',0);
        end
        set(handle,'units','pixels');
        pos = get(handle,'Position');
        lgd = legend(handle,Labels,'location',p.Results.legendLocation);
        set(lgd,'units','pixels','fontsize',p.Results.fontsize);
        title(lgd,p.Results.legend{1},'FontSize',p.Results.fontsize);
        if ~isempty(regexp(p.Results.legendLocation,'.outside','ONCE'))
            lgd_pos = get(lgd,'Position');
            h_fig.Position = h_fig.Position + [0 0 lgd_pos(3) 0];
            set(handle,'Position',pos);
        end
    end
end
end
% =========================================================================