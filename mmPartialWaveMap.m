function [MM, MPlot_object] = mmPartialWaveMap(layerArray, wavelengths, Npts, maxAOI, bReflect, bNorm, bConoscopic, varargin) 

% // begin calclation //
Nlayers = size(layerArray,2);
thick = zeros(Nlayers,1);

for k=1:Nlayers
    thick(k) = layerArray{k}{4};
end
g = find(thick == true);

kmax = sin(maxAOI*pi/180);
kx_t = [0;0];
idx =  [0;0];
kmax = 2*kmax/Npts; % convert to pixels
for X = 1:Npts
    for Y = 1:Npts
        r = sqrt((X-Npts/2).^2 + (Y-Npts/2).^2);
        if r <= Npts/2
           kx_t = [kx_t , [r*kmax; atan2(X-Npts/2, Y-Npts/2)]]; %#ok<AGROW>
           idx = [idx, [X;Y]]; %#ok<AGROW>
        end
    end
end
kx_t = kx_t(:, 2:end);
idx = idx(:, 2:end);

[Psi0, n0] = psiAmbientMap(layerArray{1}, kx_t, wavelengths);

if layerArray{Nlayers}{5} == 0
    Psi2 = psiAnisoMap(layerArray{Nlayers}, wavelengths, kx_t, n0);
elseif layerArray{Nlayers}{2} == Inf
    Psi2 = psiIsoMap(layerArray{Nlayers}, wavelengths, kx_t, n0);
else
    Psi2 = psiAmbientMap(layerArray{Nlayers}, kx_t, wavelengths);
end

if g > 2 %if there are thin layers before thick one, compute layer matrices
    for m = 2:(g-1)
        layerArray{m}{2} = -layerArray{m}{2}; %change sign of d to invert layer matrix
        if layerArray{m}{5} == 0  %check if layer is anisotropic
            Psi0 = multiprod(layerBerremanMap(layerArray{m}, wavelengths, kx_t, n0), Psi0);
        else
            Psi0 = multiprod(layerBerremanIsoMap(layerArray{m}, wavelengths, kx_t, n0), Psi0);
        end
    end
end

if Nlayers - g > 1   %if there are thin layers after thick one, compute layer matrices
    for m = (Nlayers-1):-1:(g+1)
        if layerArray{m}{5} == 0  %check if layer is anisotropic
            Psi2 = multiprod(layerBerremanMap(layerArray{m}, wavelengths, kx_t, n0), Psi2);
        else
            Psi2 = multiprod(layerBerremanIsoMap(layerArray{m}, wavelengths, kx_t, n0), Psi2);
        end
    end
end

temp = partialWaveMap(Psi0, Psi2, layerArray{g}, wavelengths, kx_t, n0, bReflect);

if bNorm
    temp = temp ./ temp(1,1,:,:);
end

if bConoscopic
    if bReflect
        for i=1:length(wavelengths)
            temp(:,:,:,i) = mmRotateRefl(temp(:,:,:,i), -kx_t(2,:));
        end
    else
        for i=1:length(wavelengths)
            temp(:,:,:,i) = mmRotate(temp(:,:,:,i), -kx_t(2,:));
        end
    end
end

MM = zeros(4, 4, Npts, Npts, length(wavelengths));
for i=1:size(kx_t, 2)
    MM(:, :, idx(1,i), idx(2,i), :) = temp(:,:,i,:);
end

if ~isempty(varargin)
    fun = varargin{1};
    MPlot_object = fun(MM, varargin{2:end});
end

end