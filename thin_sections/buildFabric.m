close all
clear all
clc

%% this goes through Don's thin section data and predicts sound velocities for the various fabrics at depth
% It also contains code for making schmidt plots
%output: ThinSecDepths, ThinSecVps 


%set(0,'defaulttextinterpreter','latex')


%% this code gets the list of depths
cd /home/dan/Documents/MATLAB/SonicLog/thinSec/C-axisdatabase
dirnames = dir();
Ldir = length(dirnames);


depths = zeros(Ldir,1);
for i = 1:Ldir
    depths(i) = str2double(dirnames(i,1).name);
end

depths(1:2) = []; % get rid of junk
depths = sort(depths); % sort the depths 
Ldepths = Ldir - 2;  % new size

    % build the circle outline for plotting
r = 1;
theta = 1:.1:360;
xc = r*cosd(theta);
yc = r*sind(theta);

%% now cycle through the depth folders

for depthindex = 1:Ldepths
    count = depthindex %keep track of progress
    workDir = num2str(depths(depthindex));
    cd(workDir)
    
    
    [ns, thetas, phis] = importfile('c-axes.txt');
    N = length(ns);
    
    % clean up the data. Get rid of nans, convert to radians
    gns = []; %(indices of good measurements)
    
    for i = 1:N
            if isfinite(thetas(i))
                gns = [gns; i];
            end
    end
    
    N = length(gns); % size of good measures
    phis = str2double(phis);
    phis = phis(gns)*pi/180; % convert the good data to radians
    thetas = thetas(gns)*pi/180;
    
    
    
    %% predict velocity
    
    % first I need to get conventional phi, which I will call phip
    phip = acos(sin(thetas).*sin(phis));
    thetap = 0*thetas; % set to 0 b/c it doesn't matter
    
    sigma = 0*2*pi/360; % reference direction
    
    sps = sp2(thetap, sigma, phip); % slowness of crystals
    spmean = sum(sps)/N;
    vp = 1/spmean*10^6;
    a = length(sps)
    
    vps(depthindex, 1) = vp; %record the velocity
    
    
    
    %% plot fabric


% these are unconventional because it isn't my data
% Dons format (Xd, Yd, Zd) = (Z, X, Y)
Y = sin(thetas).*cos(phis); % 3-d coord
Z = sin(thetas).*sin(phis);
X = cos(thetas);

%for don's plots
Xd = Y;
Yd = Z;
Zd = X;

x = (2./(1+Zd)).^.5 .*Xd/sqrt(2); % lambert azimuthal equal area projection
y = (2./(1+Zd)).^.5 .*Yd/sqrt(2);

% subplot(2,1,1)
% hold off
% plot(xc, yc)
% daspect([1 1 1])
% hold on
% plot(x,y, 'k.')


% reg coords

% first I have to move things to new octrants to fill the ordinary
% hemisphere
for j = 1:N
    if Z(j) < 0
        X(j) = -X(j);
        Y(j) = -Y(j);
        Z(j) = -Z(j);
    end
end

x = (2./(1+Z)).^.5 .*X/sqrt(2); % lambert azimuthal equal area projection
y = (2./(1+Z)).^.5 .*Y/sqrt(2);

% subplot(2,1,2)
f = figure(1);
set(f, 'Color', [1 1 1]);
hold off
plot(xc, yc, 'b-', 'linewidth', 3)
daspect([1 1 1])
hold on
plot(x,y, 'k.')

title(['Depth: ',  num2str(depths(depthindex)), 'm'], 'FontSize', 16, 'FontWeight', 'bold');

axis off

pause()


    %%
    cd ../ %back to main data directory
    
end


plot(vps, depths, 'm.')

cd ../





% save ThinSecDepths depths
% save ThinSecVps vps


