function MM = mmBerremanMap(layerArray, wavelengths, Npts, maxAOI, bReflect, bNorm, bConoscopic) 

% // begin calclation //
Nlayers = size(layerArray,2);
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0];
Ainv = [0.5,0.5,0,0;0,0,0.5,-0.5*1i;0,0,0.5,0.5*1i;0.5,-0.5,0,0];

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


[Psi0, n0] = psiAmbientInvMap(layerArray{1}, kx_t, wavelengths);

if layerArray{Nlayers}{5} == 0
    Psi2 = psiAnisoMap(layerArray{Nlayers}, wavelengths, kx_t, n0);
elseif layerArray{Nlayers}{2} == Inf
    Psi2 = psiIsoMap(layerArray{Nlayers}, wavelengths, kx_t, n0);
else
    Psi2 = psiAmbientMap(layerArray{Nlayers}, kx_t, wavelengths);
end

if Nlayers > 2
    for k = (Nlayers-1):-1:2
        if layerArray{k}{5} == 0  %check if layer is anisotropic
            Psi2 = multiprod(layerBerremanMap(layerArray{k}, wavelengths, kx_t, n0), Psi2);
        else
            Psi2 = multiprod(layerBerremanIsoMap(layerArray{k}, wavelengths, kx_t, n0), Psi2);
        end
    end
end
Psi2 = multiprod(Psi0, Psi2);
if bReflect
    J = multiprod(Psi2([3 4], [1 2], :, :, :), invert2x2(Psi2([1 2], [1 2], :, :, :)));
else
    J = invert2x2(Psi2([1 2], [1 2], :, :, :));
end
temp =  real(multiprod( multiprod(A, bigKron(J)), Ainv)) ;

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

end