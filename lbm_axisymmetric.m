% Lattice Boltzmann with axisymmetric condition from Zhoe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a simple implementation of the LBGK model.
%% By Andrey R. da Silva, August 2010
%%
%% The code does not take into account eny specific boundaru condiition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, clc
close all


% Block 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = 250*2+1;                    % Number of lines   (cells in the y direction)
Mc = 500*2+2;                    % Number of columns (cells in the x direction)
N_c = 9;

% Block 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Physical parameters (macro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_p = 340;                    % Sound velocity on the fluid [m/s]
rho_p = 1.2;                  % physical density [kg/m^3]
rho_p = 1;                    % Fluid density  [kg/m^3]
Lx = .5;                      % Maximun dimenssion in the x direction [m]
Ly = 0.0833;                  % Maximun dimenssion on th y direction  [m]
Dx = Lx/Mc                    % Lattice space (pitch)
Dt = (1/sqrt(3))*Dx/c_p       % lattice time step


% Block 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lattice parameters (micro - lattice unities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 1.9;                                      % Relaxation frequency
tau = 1/omega;                                    % Relaxation time
rho_l = 1;                                        % avereged fluid density (latice density
cs = 1/sqrt(3);                                   % lattice speed of sound
e = 1/sqrt(3);                                   % lattice speed of sound
cs2 = cs^2;                                       % Squared speed of sound cl^2
visc = cs2*(1/omega-0.5);                         % lattice viscosity
visc_phy = visc*(Dx^2)/Dt;                        % physical kinematic viscosity

% Array of distribution and relaxation functions
f=zeros(Nr,Mc,N_c);                                 
feq=zeros(Nr,Mc,N_c);
w0 = 4/9.;
w1 = 1/9.;
w2 = 1/36.;
w_alpha = [w1 w2 w1 w2 w1 w2 w1 w2 w0];
% constants velocities (link and cy, cx)
e_alpha(9, 2) = eps;
for link = 1:8
    if mod(link,2) == 1
        lambda_alpha = 1;
        % in y
        e_alpha(link, 1) = lambda_alpha*e*(sin((link-1)*pi/4));
        % in x
        e_alpha(link, 2) = lambda_alpha*e*(cos((link-1)*pi/4));
    elseif mod(link,2) == 0
        lambda_alpha = sqrt(2);
        % in y
        e_alpha(link, 1) = lambda_alpha*e*(sin((link-1)*pi/4));
        % in x
        e_alpha(link, 2) = lambda_alpha*e*(cos((link-1)*pi/4));
    end
end
e_alpha(1,1) = 0;
e_alpha(5,1) = 0;
e_alpha(9,1) = 0;

e_alpha(3,2) = 0;
e_alpha(7,2) = 0;
e_alpha(9,2) = 0;

% constant K
K = sum(e_alpha(:,1).^2)/cs2;

% radius
radius(Nr,Mc) = eps;
for x = 1:Mc
    radius(:, x) = ([0:(Nr-1)]);
end

% omega alpha
omega_alpha(Nr,Mc,N_c) = eps;
for link = 1:N_c
    omega_alpha(:,:,link) = omega*(1 + (2*tau - 1)*e_alpha(link,1)./(2.*radius));
end
omega_alpha(1,:,:) = omega;
%visc = cs2*(1./omega_alpha-0.5);
visc = ((e^2)/6)*(2*tau - 1);

% Filling the initial distribution function (at t=0) with initial values
f(:,:,:)=rho_l/9;   
ux(Nr, Mc) = eps;
uy(Nr, Mc) = eps;

% set wall
xl = [29 29 30+100*2];%Lt=230 voxel; Ltubo=200 voxel
yl = [1 20*2+1 20*2+1];
%xl = [(150-100) (150-100) (150+100) (150+100) (150-100)];
%yl = [(150+100) (150-100) (150-100) (150+100) (150+100)];
[vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8] = crossing3_axis(Nr,Mc,xl,yl);

% 4.0.1 - Adding conditions anechoic
distance = 30;
% condicao anecoica para cima
growth_delta = 0.5;
[sigma_mat9_cima Ft_cima] = build_anechoic_condition_axis(Mc, ... 
Nr, distance, growth_delta, e, e_alpha, w_alpha);
growth_delta = 0.5;
[sigma_mat9_cima_kill Ft_cima_kill] = build_anechoic_condition_axis(Mc, ... 
Nr, Nr, growth_delta, e, e_alpha, w_alpha);
% condicao anecoica para esquerda
growth_delta = -1;
[sigma_mat9_esquerda Ft_esquerda] = build_anechoic_condition_axis(Mc, ... 
Nr, distance, growth_delta, e, e_alpha, w_alpha);
% condicao anecoica para direita
growth_delta = 1;
[sigma_mat9_direito Ft_direito] = build_anechoic_condition_axis(Mc, ... 
Nr, distance, growth_delta, e, e_alpha, w_alpha);

%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construindo chirp
a=40;
total_time = round(10*Mc*sqrt(3)); % meia hora = 20*Mc*sqrt(3)
times = 0 : total_time - 1;
initial_frequency = 0;
frequency_max_lattice = 4*cs/(2*pi*a);
source_chirp = chirp(times, ... 
initial_frequency, times(end), frequency_max_lattice);
% probe to verify
probe_verify(1:total_time) = 0; 

% vetores de pressao e velocidade de particula
pressure1(1:total_time) = 0;
pressure2(1:total_time) = 0;
pressure3(1:total_time) = 0;
pressure4(1:total_time) = 0;
particle_velocity1(1:total_time) = 0;
particle_velocity2(1:total_time) = 0;
particle_velocity3(1:total_time) = 0;
particle_velocity4(1:total_time) = 0;
for ta = 1 : total_time
    
    % Block 5.1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% propagation (streaming)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    f(:,:,1) = [f(:,1:2,1) f(:,2:Mc-1,1)];
    f(:,:,2) = [f(:,1:2,2) f(:,2:Mc-1,2)];
    f(:,:,2) = [f(1:2,:,2);f(2:Nr-1,:,2)];
    f(:,:,3) = [f(1:2,:,3);f(2:Nr-1,:,3)];
    f(:,:,4) = [f(:,2:Mc-1,4) f(:,Mc-1:Mc,4)];
    f(:,:,4) = [f(1:2,:,4);f(2:Nr-1,:,4)];
    f(:,:,5) = [f(:,2:Mc-1,5) f(:,Mc-1:Mc,5)];
    f(:,:,6) = [f(:,2:Mc-1,6) f(:,Mc-1:Mc,6)];
    f(:,:,6) = [f(2:Nr-1,:,6);f(Nr-1:Nr,:,6)];
    f(:,:,7) = [f(2:Nr-1,:,7);f(Nr-1:Nr,:,7)];
    f(:,:,8) = [f(:,1:2,8) f(:,2:Mc-1,8)];
    f(:,:,8) = [f(2:Nr-1,:,8);f(Nr-1:Nr,:,8)];

    if ta >= 1
        % set bounce backs
        G=f;
        f(vec1)=G(vec5);
        f(vec5)=G(vec1);
        f(vec2)=G(vec6);
        f(vec6)=G(vec2);
        f(vec3)=G(vec7);
        f(vec7)=G(vec3);
        f(vec4)=G(vec8);
        f(vec8)=G(vec4);
    end

    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ta == 2120
       %f(150, 150, 9) = rho_p;
    end
    rho=sum(f,3);

    %% Calculando uma fonte ABC dentro do duto
    % direita = 0.5
    density_source = rho_l + 0.1*source_chirp(ta);
    Ux_t = 0; %(0.0001*source_chirp(ta))/sqrt(3) + 0.1*e;
    Uy_t = 0;%(0.0001*source_chirp(ta))/sqrt(3) + 0.1*e;
    point_y = 1;
    distance_y = 39;
    point_x = 30;
    distance_x = 30;
    direction = 0.5;
    [sigma_source Ft_source] = build_source_anechoic_axis(Nr, Mc, ...
    density_source, Ux_t, Uy_t, point_y, ...
    point_x, distance_x, distance_y, direction, e, e_alpha, w_alpha);

    % Determining the velocities according to Eq.() (see slides)
    a1 = e_alpha(1,2)*f(:,:,1);
    b1 = e_alpha(2,2)*f(:,:,2);
    c1 = e_alpha(3,2)*f(:,:,3);
    d1 = e_alpha(4,2)*f(:,:,4);
    e1 = e_alpha(5,2)*f(:,:,5);
    f1 = e_alpha(6,2)*f(:,:,6);
    g1 = e_alpha(7,2)*f(:,:,7);
    h1 = e_alpha(8,2)*f(:,:,8);
    ux = (a1 + b1 + c1 + d1 + e1 + f1 + g1 + h1)./rho;
    %
    a2 = e_alpha(1,1).*f(:,:,1);
    b2 = e_alpha(2,1).*f(:,:,2);
    c2 = e_alpha(3,1).*f(:,:,3);
    d2 = e_alpha(4,1).*f(:,:,4);
    e2 = e_alpha(5,1).*f(:,:,5);
    f2 = e_alpha(6,1).*f(:,:,6);
    g2 = e_alpha(7,1).*f(:,:,7);
    h2 = e_alpha(8,1).*f(:,:,8);
    uy = (a2 + b2 + c2 + d2 + e2 + f2 + g2 + h2)./rho;

    % Block 5.3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Determining the relaxation functions for each direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for link = 1:9
        c1 = 3/(e^2);
        C1 = e_alpha(link,2)*ux + e_alpha(link,1)*uy;
        c2 = 9/(2*e^4);
        C2 = (e_alpha(link,2)^2)*(ux.^2) + 2*e_alpha(link,1)*e_alpha(link,2)*uy.*ux ...
        + (e_alpha(link,1)^2)*(uy.^2);
        c3 = 3/(2*e^2);
        C3 = ux.^2 + uy.^2;
        feq(:,:,link)= w_alpha(link)*rho .*(1 + c1*C1 + c2*C2 - c3*C3);
        %mean(mean(feq(:,:,link)))
        %link
    end

    % Block 5.4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision (relaxation) step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Be careful with lHopital’s rule
    teta = radius;
    teta(2:end,:) = -(rho(2:end,:).*uy(2:end,:))./radius(2:end,:);
    for link = 1:9
        % Be careful with lHopital’s rule
        term_force = radius;
        term_force(2:end,:) = e_alpha(link, 2).*(-(rho(2:end,:).*ux(2:end,:).*uy(2:end,:))./radius(2:end,:)) ...
         + e_alpha(link, 1).*(-((rho(2:end,:).*uy(2:end,:).^2)./radius(2:end,:)) - 2.*rho(2:end,:).*visc.*uy(2:end,:)./radius(2:end,:).^2);
        
         if ta <= Nr*sqrt(3)
             f(:,:,link) = (1 - omega_alpha(:,:,link)).*f(:,:,link) + omega_alpha(:,:,link).*feq(:,:,link) ...
         + w_alpha(link)*teta + term_force./K*(e^2) ...
         - sigma_mat9_cima_kill(:,:,link).*(feq(:,:,link) - Ft_cima_kill(:,:,link)) ...
         - sigma_mat9_esquerda(:,:,link).*(feq(:,:,link) - Ft_esquerda(:,:,link)) ...
         - sigma_mat9_direito(:,:,link).*(feq(:,:,link) - Ft_direito(:,:,link));
         else
             % collide itself
            f(:,:,link) = (1 - omega_alpha(:,:,link)).*f(:,:,link) + omega_alpha(:,:,link).*feq(:,:,link) ...
             + w_alpha(link)*teta + term_force./K*(e^2) ...
             - sigma_mat9_cima(:,:,link).*(feq(:,:,link) - Ft_cima(:,:,link)) ...
             - sigma_mat9_esquerda(:,:,link).*(feq(:,:,link) - Ft_esquerda(:,:,link)) ...
             - sigma_mat9_direito(:,:,link).*(feq(:,:,link) - Ft_direito(:,:,link)); ... 
             - sigma_source(:,:,link).*(feq(:,:,link) - Ft_source(:,:,link));
         end
        
         %mean(mean(w_alpha(link)*teta))
    end

        
    % Ploting the results in real time   
    %surf(rho-1), view(2), shading flat, axis equal, caxis([-.00001 .00001])
    if mod(ta, 100) == 0
        imagesc(flip(rho-1)); axis equal;
        grid off
        pause(.0001)
        disp('Progresso: ');
        disp((ta/total_time*100));
    end

     %% Acquisition Data
    
    x_probe=170; y_probe=35;
    % Face data
    
    pressure1(ta) = (mean(rho(2:40, x_probe)-1))*cs2;
    particle_velocity1(ta) = mean(ux(2:40, x_probe));
    
    particle_velocity3(ta) = ux(y_probe, 228);
    pressure3(ta) = (rho(y_probe, 228)-1)*cs2;

    particle_velocity4(ta) = ux(y_probe, 230);
    pressure4(ta) = (rho(y_probe, 230)-1)*cs2;
    % point data
    
    particle_velocity2(ta) = ux(y_probe, x_probe);
    pressure2(ta) = (rho(y_probe, x_probe)-1)*cs2;
    

end %  End main time Evolution Loop

% Postprocessing

% Radiation
fft_pressure = fft(pressure1);
fft_particle_velocity = fft(particle_velocity1);
ZL = fft_pressure./fft_particle_velocity;
ka_max = (2*pi*a)/cs;
ka = linspace(0, ka_max, length(fft_pressure));

L=55; 
Zo = rho_p*cs;
Zr = Zo*1i*tan(atan(ZL./(1i*Zo))-(ka/a)*L);
Rr=(Zr-Zo)./(Zr+Zo);
%la=(-1/(2*1i*(ka/a))).*log(Rr/abs(Rr));
%% Radiation impedance
open impe.fig
hold on
plot(ka,real(Zr),'--blue')
hold on
plot(ka,imag(Zr),'--red')
axis([0 3 0 1]);
ylabel('Imped\E2ncia, Zr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Real(Analitico)', 'Imaginaria(Analitico)','Real(LBM)', 'Imaginaria(LBM)');
hold off

%% Refletion Coeff. 
open abs_r.fig
hold on
plot(ka,abs(Rr))
axis([0 3 0 1]);
ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
xlabel('Numero de Helmholtz, ka','FontSize',20);
legend('Analitico','LBM');
hold off

% End duct

%% Radiation impedance end plane
% fft_pressure = fft(pressure3);
% fft_particle_velocity = fft(particle_velocity3);
% Z_end = fft_pressure./fft_particle_velocity;
% R_end=(Z_end+Zo)./(Z_end+Zo);
% 
% open impe.fig
% %max_normalization = max([max(abs(real(Zr))) max(abs(real(Zr)))]);
% hold on
% plot(ka,real(Z_end),'--blue')
% hold on
% plot(ka,imag(Z_end),'--red')
% axis([0 3 0 1]);
% ylabel('Imped\E2ncia, Z_end','FontSize',20);
% xlabel('Numero de Helmholtz, ka','FontSize',20);
% legend('Real(Analitico)', 'Imaginaria(Analitico)','Real(LBM)', 'Imaginaria(LBM)');
% hold off
% 
% %%
% open abs_r.fig
% hold on
% plot(ka,abs(R_end))
% axis([0 3 0 1]);
% ylabel('Coeficiente de Reflex\E3o, Rr','FontSize',20);
% xlabel('Numero de Helmholtz, ka','FontSize',20);
% legend('Analitico','LBM');
% hold off