%Plot calculation results from NISE

% Linear abs % Luminescence
if isfile('Absorption.dat')
    datlin = load('Absorption.dat');
    figure; plot(datlin(:,1),datlin(:,2));
    title('Absorption & Luminescence');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Intensity (au)');
end
if isfile('Luminescence.dat')
    hold on;
    datlum = load('Luminescence.dat');
    plot(datlum(:,1),datlum(:,2));
    title('Absorption & Luminescence');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Intensity (au)');
end

% LD
if isfile('LD.dat')
    datLD = load('LD.dat');
    figure;
    plot(1e7./datLD(:,1),datLD(:,2));
    title('Linear dichroism');
    xlabel('Wavelength (nm)');
    ylabel('Intensity');
end

% CD
if isfile('CD.dat')
    datCD = load('CD.dat');
    figure;
    plot(1e7./datCD(:,1),datCD(:,2));
    title('Circular dichroism');
    xlabel('Wavelength (nm)');
    ylabel('Intensity');
end

% 2D (absorptive)
if isfile('2D.per.dat') && isfile('2D.par.dat')
    angle = 54.7; % Pump-probe polarization angle [degree]
    dat2Dper = load('2D.per.dat');
    dat2Dpar = load('2D.par.dat');
    Nw = sqrt(size(dat2Dper,1)); % Number of omega elements
    w1 = dat2Dper(1:Nw:end,1);
    w3 = dat2Dper(1:1:Nw,2);
    dat2Dper = reshape(dat2Dper(:,3),Nw,Nw)';
    dat2Dpar = reshape(dat2Dpar(:,3),Nw,Nw)';
    dat2D = struct;
    dat2D.X = 1e7./w3;
    dat2D.Y = 1e7./w1;
    dat2D.T = Tw;
    dat2D.Abs = dat2Dper*sin(angle*180/pi)^2 + dat2Dpar*cos(angle*180/pi)^2;
    plot2Ds(dat2D,Tw);
end

% 2D response (rephasing and non-rephasing)
if isfile('RperI.dat') && isfile('RparI.dat') && isfile('RperII.dat') && isfile('RparII.dat')
    angle = 54.7; % Pump-probe polarization angle [degree]
    datRperI = load('RperI.dat');
    datRparI = load('RparI.dat');
    datRperII = load('RperII.dat');
    datRparII = load('RparII.dat');
    Nt = sqrt(size(datRperI,1));
    t1 = datRperI(1:Nt:end,1);
    t3 = datRperI(1:Nt,3);
    Tw = datRperI(1,2);
    datRperI = reshape(datRperI(:,4),Nt,Nt)';
    datRparI = reshape(datRparI(:,4),Nt,Nt)';
    datRperII = reshape(datRperII(:,4),Nt,Nt)';
    datRparII = reshape(datRparII(:,4),Nt,Nt)';
    figure; contourf(t3,t1,datRperI*sin(angle*180/pi)^2 + datRparI*cos(angle*180/pi)^2);
    xlabel('t_3 (fs)'); ylabel('t_1 (fs)');
    axis equal;
    title('Resphasing response');
    figure; contourf(t3,t1,datRperII*sin(angle*180/pi)^2 + datRparII*cos(angle*180/pi)^2);
    xlabel('t_3 (fs)'); ylabel('t_1 (fs)');
    axis equal;
    title('Non-resphasing response');
end
% Others
% datDif = load('RMSD.dat');
% datPop = load('Pop.dat');
% datPopF = load('PopF.dat');
% datAnalyse = load('LocalDensityMatrix.dat');
