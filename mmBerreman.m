function M = mmBerreman(layerArray, wavelengths, aoi, bReflect, bNormalize)

aoi = aoi*pi/180;
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0];
Ainv = [0.5,0.5,0,0;0,0,0.5,-0.5*1i;0,0,0.5,0.5*1i;0.5,-0.5,0,0];
Nlayers = size(layerArray,2);
Nlam = length(wavelengths);
layerMat = zeros(4,4,Nlam,Nlayers-2);

if strcmp(layerArray{1}{1}, 'air')
    t1 = 1 / (2 * cos(aoi));
    Psi_i = [  t1, 1/2,   0,   0 ;...
                0,   0, 1/2,  t1 ;...
              -t1, 1/2,   0,   0 ;...
                0,   0, 1/2, -t1 ];
    n0 = ones(1, length(wavelengths));
else
    [Psi_i, n0] = psiAmbientInv(layerArray{1}, aoi, wavelengths);
end

kx = n0 .* sin(aoi); % x-component of the wavevector

if strcmp(layerArray{Nlayers}{1}, 'air')
    t1 = cos(aoi);
    Psi_e = [t1,0,-t1,0;1,0,1,0;0,1,0,1;0,t1,0,-t1];
elseif layerArray{Nlayers}{5} == 0
    Psi_e = psiAniso(layerArray{Nlayers}, wavelengths, kx);
elseif layerArray{Nlayers}{2} == Inf
    Psi_e = psiIso(layerArray{Nlayers}, wavelengths, kx);
else
    Psi_e = psiAmbient(layerArray{Nlayers}, aoi, wavelengths);
end

for m = 1:Nlayers-2
    if layerArray{m+1}{5} == 0
        layerMat(:,:,:,m) = layerBerreman(layerArray{m+1}, wavelengths, kx);
    else
        layerMat(:,:,:,m) = layerBerremanIso(layerArray{m+1}, wavelengths, kx);
    end
end

% multiply things efficiently without loops
if Nlayers == 2
    P = multiprod(Psi_i, Psi_e);
elseif Nlayers > 3
    P = layerMat(:,:,:,Nlayers-2);
    for k = 1:(Nlayers-3)
        P = multiprod(layerMat(:,:,:,Nlayers-2-k), P);
    end
    P = multiprod( multiprod(Psi_i, P), Psi_e);
else
    P = layerMat(:,:,:,1);
    P = multiprod( multiprod(Psi_i, P), Psi_e);
end
if bReflect
    J = multiprod(P([3 4], [1 2], :), invert2x2(P([1 2], [1 2], :)));
else
    J = invert2x2(P([1 2], [1 2], :));
end
M =  real(multiprod( multiprod(A, bigKron(J)), Ainv)) ;

if bNormalize
    M = M ./ M(1, 1, :);
end

end
