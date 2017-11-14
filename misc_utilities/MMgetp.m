function out = MMgetp(M,parameter)
% This function contains many parameters that one can compute from a
% Mueller matrix (M). In general, M is assumed to be an
% experimental one. Hence, a Mueller-Jones matrix or even a physical M is 
% not assumed. For most parameters, M is first converted to its closest 
% Mueller-Jones matrix, or its Nearest Jones matrix. 

if ndims(M) > 3   % reshape array into 4,4,N
    sz = size(M);
    M = reshape(M,4,4,[]);
else
    sz = 0;
end

switch lower(parameter)
    case 'opt prop'
        J = nearestJones(M);
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
        out = squeeze([LB;LD;LBp;LDp;CB;CD;A]);
    case 'lm'
        J = nearestJones(M);
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
        out = [A,-LD,-LDp,CD ; -LD,A,CB,LBp ; -LDp,-CB,A,-LB ; CD,-LBp,LB,A];
    case 'logm' %log of Mueller matrix with filtering
        Mfiltered = filterM(M);
        out = zeros(size(M));
        for n=1:size(M,3); out(:,:,n) = logm(Mfiltered(:,:,n)); end
    case 'expm' %log of Mueller matrix with filtering
        out = zeros(size(M));
        for n=1:size(M,3); out(:,:,n) = expm(M(:,:,n)); end
    case 'lb'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
    case 'ld'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
    case 'lbp'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
    case 'ldp'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
    case 'cb'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = real(O.*( J(1,2,:) - J(2,1,:) ));
    case 'cd'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        out = -imag(O.*( J(1,2,:) - J(2,1,:) ));
    case 'a' % total mean extinction
        J = nearestJones(M);
        K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
        out = -2*real(log(1./K));
    case 'a_aniso' % anisotropic part of the mean extinction
        J = nearestJones(M);
        K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
        T = acos( K.*( J(1,1,:) + J(2,2,:) )./2); % 2*T = sqrt(L.^2 + Lp.^2 + C.^2)
        O = (T.*K)./(sin(T));
        LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        CD = -imag(O.*( J(1,2,:) - J(2,1,:) ));
        out = sqrt(LD.^2 + LDp.^2 + CD.^2);  % not same as imag(2*T) !
    case 'a_iso'  % isotropic part of the mean extinction
        J = nearestJones(M);
        K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
        T = acos( K.*( J(1,1,:) + J(2,2,:) )./2); % 2*T = sqrt(L.^2 + Lp.^2 + C.^2)
        O = (T.*K)./(sin(T));
        LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        CD = -imag(O.*( J(1,2,:) - J(2,1,:) ));
        out = -2*real(log(1./K)) - sqrt(LD.^2 + LDp.^2 + CD.^2);
    case 'ldmag'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        LD = imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LDp = imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        out = sqrt(LD.^2 + LDp.^2);
    case 'ldang'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        LD = -imag(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LDp = -imag(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        out = atan2(LDp , LD)./2;
        %out = out + pi*(out < 0);
    case 'lbang'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        LB = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LBp = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        out = atan2(LBp , LB)./2;
        out = out + pi*(out < 0);
    case 'lbmag'
        J = nearestJones(M);
        O = jonesAnisotropy(J);
        LB = real(1i.*O.*( J(1,1,:) - J(2,2,:) ));
        LBp = real(1i.*O.*( J(1,2,:) + J(2,1,:) ));
        out = sqrt(LB.^2 + LBp.^2);
    case 'di' % Depolarization Index
        out = (sqrt(squeeze(sum(sum(M.^2,1),2))./squeeze(M(1,1,:)).^2-1)./sqrt(3)).';
    case 'jones' % Jones matrix of a Mueller-Jones matrix
        out = MJ2J(M);
    case 'nearestjones'
        out = nearestJones(M);  % Jones matrix
            % next line just phases the Jones matrix so that the
            % imaginary part of J(1,1) = 0. i.e., it matches case 'jones'
        for n=1:size(out,3); out(:,:,n) = exp( -1i*angle(out(1,1,n)) ) * out(:,:,n); end
    case 'covar' % Cloude carvariance matrix
        out = M2Cov(M);
    case 'covar2m' % Cloude carvariance matrix
        out = Cov2M(M);
    case 'mfiltered' % closest physical Mueller matrix
        out = filterM(M);
end

if size(out,1) == 1 && size(out,2) == 1 %remove extra singletons
    out = squeeze(out).';
end
if sz ~= 0  %  reshape to match input dimensions
    sz2 = size(out);
    out = reshape(out,[sz2(1:(length(sz2)-1)),sz(3:length(sz))]);
end

end  % end parent function

% \\ LOCAL FUNCTIONS \\

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

function Mfiltered = filterM(M)  % M to nearest physical M
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
Mfiltered = Cov2M(C);
end

function O = jonesAnisotropy(J)
K = ( J(1,1,:).*J(2,2,:) - J(1,2,:).*J(2,1,:)).^(-1/2);
T = acos( K.*( J(1,1,:) + J(2,2,:) )./2);
O = (T.*K)./(sin(T));
end