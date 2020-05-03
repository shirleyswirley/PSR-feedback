close all;
clear all;
setup_figs;

%-------------------------------
% 1.) Generate the normalized particle flux curves (eq 13 from DeVries et al., 2014) 
% 2-D matrix, down columns are increasing depth, to the right in rows are
% with increasing epsilon
%-------------------------------

% - Define the epsilon (aka beta = PSD slope) values we want to use
% to generate flux profiles.
epsilon = 3:0.05:6;

% - Define constants from DeVries et al., 2014 mostly table 1
z_0=0; % [m]
z_f=5000; % [m]
z_res=50; % [m]
z=z_0:z_res:z_f; % water column depth range from bottom of euphotic zone to bottom of ocean

c_w = 2.2E5; % coeff in particle sinking rate equation as a fxn of D
c_r = 0.03; % degradation rate of sinking particles
c_1 = c_r/c_w;

DS = 20E-6; % smallest diameter of particles
DL_0 = 2000E-6; % largest diameter of particles in the euphotic zone

eta = 1.17; % exponent in particle sinking rate equation as a fxn of D
zeta = 2.28; % exponent in mass equation as a fxn of D

%-------------------------------
% 2.) Use DeVries et al., 2014 equations to create POC
% flux profiles from epsilon 
%-------------------------------

Fcurves=zeros(length(z),length(epsilon));

for i=1:length(epsilon)
    for j=1:length(z)
        zj=z(j);
        % equation 11
        DL = (max(DL_0^eta-c_1*eta/zeta*zj,0))^(1/eta);
        % equation 13
        F = @(D)(D.^(zeta+eta-epsilon(i))).*(1+(c_1*eta/zeta).*(D.^-eta)*zj).^((1-epsilon(i))/eta);
        if (j==1); F_0 = quad(F,DS,DL); end; % save the particle flux at surface to use for normalization 
        Fcurves(j,i) = quad(F,DS,DL)/F_0;
    end
end

%-------------------------------
% 3.) Visualize the F curve profiles
% depending on epsilon (or beta) at the surface
%-------------------------------
f=figure;
set(f,'color','white','units','inches',...
    'position', [0.5 0.5 6 7]);
colorcell = {'r','g','b'};
icolor = 1;
for idx=[find(epsilon==3.3) find(epsilon==4.3) find(epsilon==5.3)]
    plot(Fcurves(:,idx),-z,'-','linewidth',2,'color',colorcell{icolor});
    hold on;grid on;
    icolor=icolor+1;
end
%plot(((z+100)./100).^-0.858,-z,'k-','linewidth',2);
%legend('\beta = 3.3','\beta = 4.3','\beta = 5.3','Martin Curve');
legend('\beta = 3.3','\beta = 4.3','\beta = 5.3','Location','Southeast');
xlabel('Normalized particulate organic carbon flux');
ylabel('Depth below euphotic zone [m]');
title('PRiSM-generated remineralization profiles');
set(gca,'fontsize',10,'fontweight','bold');
ylim([-1000 0]);

print(f, [fig_save_path 'suppfig1_PSRfbpaper_final.pdf'], '-dpdf', '-r300');
print(f, [fig_save_path 'suppfig1_PSRfbpaper_final.png'], '-dpng', '-r300');
