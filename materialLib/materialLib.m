function [epsilon,alpha,mu] = materialLib(material, wavelengths, varargin)
% small library of optical functions for anisotropic materials
Nwl = length(wavelengths);
epsilon = zeros(3,3,Nwl);
mu = setDiag(ones(3,Nwl));
alpha = 0;

switch material
    
    case 'rubrene'
        data = load('rubreneOptfun.mat');
        data = data.filetosave;
        eV = (1239.8)./wavelengths;
        epsilon = setDiag( ((interp1(data(:,1),data(:,2:4),eV)).^2)' );
        
    case '+EDS'
        alpha = zeros(3,3,Nwl);
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 1.06482*lam2./(lam2 - 0.0103027)...
            +2.3712*lam2./(lam2 - 92.3287) + 1.2728;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 1.07588*lam2./(lam2 - 0.0101915)...
            +2.64847*lam2./(lam2 - 94.8497) + 1.2751;
        lam2 = wavelengths.^2;
        alpha(1,1,:) = wavelengths.^3*(0.0146441)./(lam2 - 100.202^2).^2;
        alpha(3,3,:) = wavelengths.^3*(-0.0301548)./(lam2 - 100^2).^2;
        alpha(2,2,:) = alpha(1,1,:);
        
    case '-EDS'
        alpha = zeros(3,3,Nwl);
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 1.06482*lam2./(lam2 - 0.0103027)...
            +2.3712*lam2./(lam2 - 92.3287) + 1.2728;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 1.07588*lam2./(lam2 - 0.0101915)...
            +2.64847*lam2./(lam2 - 94.8497) + 1.2751;
        lam2 = wavelengths.^2;
        alpha(1,1,:) = wavelengths.^3*(-0.0146441)./(lam2 - 100.202^2).^2;
        alpha(3,3,:) = wavelengths.^3*(0.0301548)./(lam2 - 100^2).^2;
        alpha(2,2,:) = alpha(1,1,:);
        
    case 'SYEDS'
        if nargin > 2
            c = varargin{1};
        else
            c = 3.1547;
        end
        alpha = zeros(3,3,Nwl);
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 1.06482*lam2./(lam2 - 0.0103027)...
            +2.3712*lam2./(lam2 - 92.3287) + 1.2728;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 1.07588*lam2./(lam2 - 0.0101915)...
            +2.64847*lam2./(lam2 - 94.8497) + 1.2751;
        data = load('SYEDS');
        data = data.SYEDS;
        data = ((interp1(real(data(1,:)).',[real(data(2,:));imag(data(2,:))].',wavelengths)))';
        epsilon(1,1,:) = squeeze(c*(data(1,:) + 1i*data(2,:))) + squeeze(epsilon(1,1,:)).';
        lam2 = wavelengths.^2;
        alpha(1,1,:) = -wavelengths.^3*(0.0146441)./(lam2 - 100.202^2).^2;
        alpha(3,3,:) = -wavelengths.^3*(-0.0301548)./(lam2 - 100^2).^2;
        alpha(2,2,:) = alpha(1,1,:);
        
    case '+quartz'
        alpha = zeros(3,3,Nwl);
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 1.07044083*lam2./(lam2 - 0.0100585997)...
            +1.10202242*lam2./(lam2 - 100) + 1.28604141;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 1.09509924*lam2./(lam2 - 0.0102101864)...
            +1.15662475*lam2./(lam2 - 100) + 1.28851804;
        lam2 = wavelengths.^2;
        alpha(1,1,:) = wavelengths.^3*(0.0198)./(lam2 - 93^2).^2;
        alpha(3,3,:) = wavelengths.^3*(-0.0408)./(lam2 - 87^2).^2;
        alpha(2,2,:) = alpha(1,1,:);
        
    case '-quartz'
        alpha = zeros(3,3,Nwl);
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 1.07044083*lam2./(lam2 - 0.0100585997)...
            +1.10202242*lam2./(lam2 - 100) + 1.28604141;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 1.09509924*lam2./(lam2 - 0.0102101864)...
            +1.15662475*lam2./(lam2 - 100) + 1.28851804;
        lam2 = wavelengths.^2;
        alpha(1,1,:) = -wavelengths.^3*(0.0198)./(lam2 - 93^2).^2;
        alpha(3,3,:) = -wavelengths.^3*(-0.0408)./(lam2 - 87^2).^2;
        alpha(2,2,:) = alpha(1,1,:);
        
    case 'sapphire'
        osc1A = [1.4313493,0.65054713,5.3414021];
        osc1E = [0.0726631,0.1193242,18.028251].^2;
        osc2A = [1.5039759,0.55069141,6.5927379];
        osc2E = [0.0740288,0.1216529,20.072248].^2;
        lam2 = (wavelengths/1000).^2;
        for n = 1:Nwl
            epsilon(1,1,n) = sum(lam2(n)*osc1A./(lam2(n) - osc1E))+1;
            epsilon(3,3,n) = sum(lam2(n)*osc2A./(lam2(n) - osc2E))+1;
        end
        epsilon(2,2,:) = epsilon(1,1,:);
        
    case 'aBBO'
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = -lam2*0.0155+0.0184./(lam2 - 0.0179)+2.7405;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = -lam2*0.0044+0.0128./(lam2 - 0.0156)+2.373;
        
    case 'KDPnoG' % potassium acid phthalate without gyration
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = lam2.*13.0052/(lam2 - 400)...
            +0.01008956./(lam2 - 0.0129426) + 2.259276;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = lam2.*3.2279924./(lam2 - 400)...
            +0.008637494./(lam2 - 0.012281) + 2.132668;
        
        
    case 'LiNbO3'
        osc1A = [2.6734,1.229,12.614];
        osc1E = [0.01764,0.05914,474.6];
        osc2A = [2.9804,0.5981,8.9543];
        osc2E = [0.02047,0.0666,416.08];
        lam2 = (wavelengths/1000).^2;
        for n = 1:Nwl
            epsilon(1,1,n) = sum(lam2(n)*osc1A./(lam2(n) - osc1E))+1;
            epsilon(3,3,n) = sum(lam2(n)*osc2A./(lam2(n) - osc2E))+1;
        end
        epsilon(2,2,:) = epsilon(1,1,:);
        
    case 'KDP'
        %This includes epsilon from Zernike (1964) and alpha from
        %Konstantinova (2000)
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = (13.00522*lam2./(lam2 - 400))+(0.01008956./(lam2 - 0.0129426))+2.259276;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = (3.2279924*lam2./(lam2 - 400))+(0.008637494./(lam2 - 0.0122810))+2.132668;
        lam2 = wavelengths.^2;
        alpha = zeros(3,3,Nwl);
        alpha(1,2,:) = wavelengths.^3.*(.023)./(lam2 - 10^2).^2;
        alpha(2,1,:) = alpha(1,2,:);
        
    case 'KAP'
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = (-0.013296131.*lam2)+(0.037318762./(lam2 - (0.175493731^2)))+2.663434655;
        epsilon(2,2,:) = 2.670444937+(0.031617528./(lam2-(0.208293225^2)))-(0.004014395.*lam2);
        epsilon(3,3,:) = 2.196073191+(0.015025722./(lam2-(0.190981743^2)))-(0.006100780.*lam2);
        
    case 'TiO2'
        lam2 = (wavelengths/1000).^2;
        epsilon(1,1,:) = 0.2441./(lam2 - 0.0803)+5.913;
        epsilon(2,2,:) = epsilon(1,1,:);
        epsilon(3,3,:) = 0.3322./(lam2 - 0.0843)+7.197;
        
    case 'test'

        epsilon(1,1,:) = ones(Nwl,1)*(2.2);
        epsilon(2,2,:) = ones(Nwl,1)*(2.2);
        epsilon(3,3,:) = ones(Nwl,1)*(2.2);

        
end
end

function out = setDiag(diag)
out = zeros(size(diag, 1), size(diag, 1), size(diag, 2));
for i=1:size(diag, 1)
    out(i,i,:) = diag(i,:);
end
end