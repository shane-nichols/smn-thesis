function M = mmPartialWave(layerArray, wavelengths, aoi, bReflect, bNorm)

%Required Inputs:
    % layerArray: cell array describing each layer in the model. See
        % input_layerArray.m for details.
    % Lam is the array of wavelengths, in nm.
    % AOI is the angle of incidence, in degrees.
    % b_reflect: boolean value, true for reflection calculation.
% Optional Inputs:
    % varargin: do not put for no plotting. For plotting, all optional
        % parameters valid for MMplot.m can go here. Required plotting
        % parameters are supplied by the calculation.

aoi = aoi.*pi./180;
Nlayers = size(layerArray,2);
g = 1;
while layerArray{g}{4} == 0
    g = g + 1;
end

if strcmp(layerArray{1}{1}, 'air')
    t1 = cos(aoi);
    Psi0 = [t1,0,-t1,0;1,0,1,0;0,1,0,1;0,t1,0,-t1];
    n0 = ones(1, length(wavelengths));
else
    [Psi0, n0] = psiAmbient(layerArray{1}, aoi, wavelengths);
end

kx = n0 .* sin(aoi); % x-component of the wavevector

if strcmp(layerArray{Nlayers}{1}, 'air')
    t1 = cos(aoi);
    Psi2 = [t1,0,-t1,0;1,0,1,0;0,1,0,1;0,t1,0,-t1];
elseif layerArray{Nlayers}{5} == 0
    Psi2 = psiAniso(layerArray{Nlayers}, wavelengths, kx);
elseif layerArray{Nlayers}{2} == Inf
    Psi2 = psiIso(layerArray{Nlayers}, wavelengths, kx);
else
    Psi2 = psiAmbient(layerArray{Nlayers}, aoi, wavelengths);
end

if g > 2 %if there are thin layers before thick one, compute layer matrices
    for m = 2:(g-1)
        layerArray{m}{2} = -layerArray{m}{2}; %change sign of d to invert layer matrix
        if layerArray{m}{5} == 0  %check if layer is anisotropic
            Psi0 = multiprod(layerBerreman(layerArray{m}, wavelengths, kx), Psi0);
        else
            Psi0 = multiprod(layerBerremanIso(layerArray{m}, wavelengths, kx), Psi0);
        end
    end
end

if Nlayers - g > 1   %if there are thin layers after thick one, compute layer matrices
    for m = (Nlayers-1):-1:(g+1)
        if layerArray{m}{5} == 0  %check if layer is anisotropic
            Psi2 = multiprod(layerBerreman(layerArray{m}, wavelengths, kx), Psi2);
        else
            Psi2 = multiprod(layerBerremanIso(layerArray{m}, wavelengths, kx), Psi2);
        end
    end
end

M = partialWave(Psi0, Psi2, layerArray{g}, wavelengths, kx, bReflect);

if bNorm
    M = M ./ M(1,1,:);
end

end