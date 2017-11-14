function MMout = mmRotateRefl(M,radians)
% M is a 4x4xn array of reflection Mueller matrices
sz = size(M);
M = reshape(M, 4, 4, []);
MMout = M;
radians = 2 .* radians;
C2 = reshape(cos(radians), 1, 1, []);
S2 = reshape(sin(radians), 1, 1, []);
MMout(1,2,:) = M(1,2,:).*C2 - M(1,3,:).*S2;
MMout(1,3,:) = M(1,3,:).*C2 + M(1,2,:).*S2;
MMout(2,1,:) = M(2,1,:).*C2 + M(3,1,:).*S2;
MMout(3,1,:) = M(3,1,:).*C2 - M(2,1,:).*S2;
MMout(2,4,:) = M(2,4,:).*C2 + M(3,4,:).*S2;
MMout(3,4,:) = M(3,4,:).*C2 - M(2,4,:).*S2;
MMout(4,2,:) = M(4,2,:).*C2 - M(4,3,:).*S2;
MMout(4,3,:) = M(4,3,:).*C2 + M(4,2,:).*S2;
MMout(2,2,:) = C2.*(M(3,2,:).*S2 + M(2,2,:).*C2) - S2.*(M(3,3,:).*S2 + M(2,3,:).*C2);
MMout(2,3,:) = C2.*(M(3,3,:).*S2 + M(2,3,:).*C2) + S2.*(M(3,2,:).*S2 + M(2,2,:).*C2);
MMout(3,2,:) = -C2.*(M(2,2,:).*S2 - M(3,2,:).*C2) + S2.*(M(2,3,:).*S2 - M(3,3,:).*C2);
MMout(3,3,:) = S2.*(-M(2,2,:).*S2 + M(3,2,:).*C2) + -C2.*(M(2,3,:).*S2 - M(3,3,:).*C2);
MMout = reshape(MMout, sz);
end