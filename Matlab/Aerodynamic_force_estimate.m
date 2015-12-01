% ------------------
% TMohren 2015/09/03
% Aerodynamic force estimate for a flapping and rotating flat plate, based on:
% Sane, S. P. and Dickinson, M. H. The aerodynamic efects of wing rotation and
% a revised quasi-steady model of flapping fight. Journal of Experimental Biology, 2002.
% ------------------

syms t 
R       = 0.05;
C       = 0.02;
S       = R*C;
W       = 1.27e-4;
rho     = 1.225;
rot     = 3;
A_2     = 0.5;  % (R^2*C^2)/2
A_1     = 1;    % R*(R*C^2)
R_2     = 1/3;  % non dimensional radius of the 2nd moment of wing area
Crot    = 0;    

%% define angles and velcities
phi     = deg2rad(15)*-cos(2*pi*10*t);	% stroke angle 
phi_d1  = diff(phi);
phi_d2  = diff(phi_d1);

Uf      = phi_d1 * R;                   % tip velocity
Uind    = 2*pi*rot*R;

alpha   = abs(atan(Uf/Uind));      
alpha_d1 = diff(alpha);
alpha_d2 = diff(alpha_d1);


%% Compute forces 
Fa      =  rho*pi/4*R^2*C^2*A_2*...
    ( phi_d2.*sin(alpha) + phi_d1.*alpha_d1.*cos(alpha) ) - ...
    alpha_d2 * rho*pi/16*C^3*R * A_1  ;
CLt     = 0.225 + 1.58*sin(2.13*alpha-7.2);
CDt     = 1.92 - 1.55 *cos(2.04*alpha-9.81);
F_trans = rho*S*Uf.^2*R_2/2 .* (CLt.^2 + CDt.^2).^(0.5);
Frot    = Crot * rho * Uf .* abs(phi_d1) * C^2 * R * A_2;

%% plot figure
t       = linspace(0,0.05,100);
f       = t*10;
fig1 = figure('position',[200 200 800 500]);
subplot(211)
    plot(f,eval(alpha)*180/pi)
	hold on
    plot(f,eval(phi)*180/pi)
        legend('\alpha','\phi')
        xlabel('Time [s]')
        ylabel('Angle [deg]')
    ax.XTick = linspace(0,0.5,3);
subplot(212)
    plot(f,eval(Fa),'r')
        hold on
    plot(f,eval(F_trans),'b')
    plot(f,eval(Fa+F_trans+Frot),'k','LineWidth',3)
        legend('F added mass','F translation', 'F total' )
        ylabel('Force [N]')
        xlabel('$\displaystyle\frac{t}{T}$[-]','interpreter','latex')
    ax = gca;
    ax.XTick = linspace(0,0.5,3);
