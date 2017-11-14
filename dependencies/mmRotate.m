function MMout = mmRotate(M,radians)
% M is a Mueller matrix array of any dimension. The first two dimension
% must be the Mueller matrix elements. MMout is a Mueller array with the
% same dimension as the input array. 

sz = size(M);
M = reshape(M, 4, 4, []);
MMout = M;
radians=-2*radians;
C2=reshape(cos(radians), 1, 1, []);
S2=reshape(sin(radians), 1, 1, []);
MMout(1,2,:) = M(1,2,:).*C2 + M(1,3,:).*S2;
MMout(1,3,:) = M(1,3,:).*C2 - M(1,2,:).*S2;
MMout(2,1,:) = M(2,1,:).*C2 + M(3,1,:).*S2;
MMout(3,1,:) = M(3,1,:).*C2 - M(2,1,:).*S2;
MMout(2,4,:) = M(2,4,:).*C2 + M(3,4,:).*S2;
MMout(3,4,:) = M(3,4,:).*C2 - M(2,4,:).*S2;
MMout(4,2,:) = M(4,2,:).*C2 + M(4,3,:).*S2;
MMout(4,3,:) = M(4,3,:).*C2 - M(4,2,:).*S2;
MMout(2,2,:) = C2.*(M(3,2,:).*S2 + M(2,2,:).*C2) + S2.*(M(3,3,:).*S2 + M(2,3,:).*C2);
MMout(2,3,:) = C2.*(M(3,3,:).*S2 + M(2,3,:).*C2) - S2.*(M(3,2,:).*S2 + M(2,2,:).*C2);
MMout(3,2,:) = -C2.*(M(2,2,:).*S2 - M(3,2,:).*C2) - S2.*(M(2,3,:).*S2 - M(3,3,:).*C2);
MMout(3,3,:) = S2.*(M(2,2,:).*S2 - M(3,2,:).*C2) - C2.*(M(2,3,:).*S2 - M(3,3,:).*C2);
MMout = reshape(MMout,sz);
end

%%%%This alternative code is faster for very large arrays when radians is scalar.

% sz = size(M);
% M = reshape(M, 4, 4, []);
% radians = reshape(2*radians, 1, []);
% C=cos(radians);
% S=sin(radians);
% if isscalar(radians)
%     R=[1,0,0,0;0,C,-S,0;0,S,C,0;0,0,0,1];
%     M = multiprod(multiprod(R, M), R.');
% else
%     R = zeros(3, 3, size(M, 3));
%     R(1,1,:) = ones(1, size(M, 3));
%     R(4,4,:) = ones(1, size(M, 3));
%     R(2,2,:) = C;
%     R(3,3,:) = C;
%     R(2,3,:) = -S;
%     R(3,2,:) = S;
%     M = multiprod(multiprod(R, M), multitransp(R));
% end
% M = reshape(M,sz);