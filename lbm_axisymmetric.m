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

Nr = 300;                    % Number of lines   (cells in the y direction)
Mc = 300;                    % Number of columns (cells in the x direction)
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

%Block 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the iteractive process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ta = 1 : 10000
    
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

    % Block 5.2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recalculating rho and u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ta == 1
       f(150, 150, 9) = rho_p;
    end
    rho=sum(f,3);

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
        
        % collide itself
        f(:,:,link) = (1 - omega_alpha(:,:,link)).*f(:,:,link) + ...
         omega_alpha(:,:,link).*feq(:,:,link) + w_alpha(link)*teta + term_force./K*(e^2); 
         %mean(mean(w_alpha(link)*teta))
    end

        
    % Ploting the results in real time   
    surf(rho-1), view(2), shading flat, axis equal, caxis([-.00001 .00001])
    %imagesc(rho)
    grid off
    pause(.0001)
    ta

end %  End main time Evolution Loop
