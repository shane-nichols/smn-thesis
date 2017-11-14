function n = materialLibIso(material,Lam)

k = length(Lam);

switch material
    case 'air'
        n = ones(k,1);
        
    case 'BK7'
        oscA = [1.03961212,0.231792344,1.01046945];
        oscE = [6000.69867,20017.9144,103560653];
        lam2 = Lam.^2;
        for m = 1:k
            n(m) = sqrt(sum(lam2(m)*oscA./(lam2(m) - oscE))+1);
        end
    case 'Fused Silica'
        Lam = Lam./1000;
        n = sqrt(1+0.6961663./(1-(0.0684043./Lam).^2)+0.4079426./...
            (1-(0.1162414./Lam).^2)+0.8974794./(1-(9.896161./Lam).^2));
    case 'Al_Aspens'
        load Al_Aspens.txt;
        n = interp1(Al_Aspens(:,1),Al_Aspens(:,2),Lam)...
            +1i*interp1(Al_Aspens(:,1),Al_Aspens(:,3),Lam);
        
    case 'Al_Palak'
        load Al_Palak.txt;
        n = interp1(Al_Palak(:,1),Al_Palak(:,2),Lam)...
            +1i*interp1(Al_Palak(:,1),Al_Palak(:,3),Lam);
        
    case 'Ag'
        load Ag.txt;
        Lam = Lam / 1000;
        n = interp1(Ag(:,1),Ag(:,2),Lam)...
            -1i*interp1(Ag(:,1),Ag(:,3),Lam);
        
    case 'AuFilm'
        load AuFilm.txt;
        Lam = Lam / 1000;
        n = interp1(AuFilm(:,1),AuFilm(:,2),Lam)...
            -1i*interp1(AuFilm(:,1),AuFilm(:,3),Lam);
        
    case 'Pt'
        load Pt.txt;
        Lam = Lam / 1000;
        n = interp1(Pt(:,1),Pt(:,2),Lam)...
            -1i*interp1(Pt(:,1),Pt(:,3),Lam);
    case 'Cr'
        data = load('Cr.txt');
        n = interp1(data(:,1),data(:,2),Lam)...
            +1i*interp1(data(:,1),data(:,3),Lam);
        
    case 'Si100'
        data = load('Si100p-type_0.2ohm.txt');
        n(:) = interp1(data(:,1),data(:,2),Lam)...
            +1i*interp1(data(:,1),data(:,3),Lam);
        
    case 'Si110'
        data = load('Si110p-type_11ohm.txt');
        n = interp1(data(:,1),data(:,2),Lam)...
            +1i*interp1(data(:,1),data(:,3),Lam);
    case 'Si111'
        data = load('Si111n-type_0.5ohm.txt');
        n = interp1(data(:,1),data(:,2),Lam)...
            +1i*interp1(data(:,1),data(:,3),Lam);
    case 'SiO2film'
        Lam = Lam/1000;
        data = load('SiO2film.txt');
        n = interp1(data(:,1),data(:,2),Lam)...
            -1i*interp1(data(:,1),data(:,3),Lam);
        
    case 'ITO'
        data = load('ITO.txt');
        n = interp1(data(:,1),data(:,2),Lam)...
            +1i*interp1(data(:,1),data(:,3),Lam);
        
    case 'SU8'
        Lam=Lam/1000;
        n = 1.566 + 0.00796./Lam.^2 +  0.00014./Lam.^4;
        %n = 1.6 + 0.00796./Lam.^2 +  0.00014./Lam.^4;
        
end
end
