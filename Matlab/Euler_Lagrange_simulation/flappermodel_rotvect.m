function [x,y,z,strainx,strainy,strainxy,T,flappingangle] = flappermodel_rotvect(rot_vec,frot,periodic,fflap,cycles)
%%
% Function to solve ordinary differential equations for plate model of wing flapping with concurrent body rotation 
% Created by Annika Eberle 
% September 20, 2013
% Modified August 4, 2014
% Modified by Thomas Mohren, 2015-07-02


% Outputs: 
% x = x position
% y = y position 
% z = z position (equal to surface displacement, i.e. w(x,t,t)) 
% T = time; 
% flappingangle = input flapping; 
% rotation angle = input rotation; 
%
% Inputs: 
% flapamp = amplitude of flapping (in deg) 
% rotamp = amplitude of rotation (in deg)
% fflap = frequency of flapping (in Hz)
% frot = frequency of rotation (in Hz)
% cycles = number of rotation cycles desired to be run
% sampfreq = sampling frequency for output strains and displacements (in Hz)
% periodic = toggle for type of angular velocity (if 1, then periodic angular velocity; if 0, then constant angular velocity)

%% Parameterize model 
%Geometry and nodes
    nodes   = 4;
    a       = 1; %half chord in cm
    b       = 2.5; %half span in cm
    xpos    = [-a a a -a]; %x position of plate nodes
    ypos    = [0 0 2*b 2*b];%y position of plate nodes

%Material properties 
    E       = 3e9*10^-2;    %Young's modulus (converting from kg/m/s2 to kg/cm/s2) - currently for acrylic (3 GPa), but for moth:500e9*10^-2
    nu      = 0.35;         %Poisson's ratio for the plate - currently for acrylic. for moth: 0.49 
    G       = E/(2*(1+nu)); %Shear modulus for isotropic material
    h       = 1.27e-2;      %plate height in cm -- currently for acrylic sheet, but for moth:50e-4
    density = 1180*(1e-2)^3;%density of plate in kg/cm^3 (converting from m^3)
%Simulation parameters 
    dampingfactor = 63;    %multiplier for velocity-proportional damping via mass matrix 
    sampfreq = 1000; 
    flapamp = 15;
    rotamp  = 30;
    

    
%Kinematics of input angles for flapping and rotation 
    syms x y t %create symbolic variables 
    
%Specify properties for sigmoidal startup
    sigprop = [1;10;3];
    sigd = sigprop(1);
    sigc = sigprop(2);
    sign = sigprop(3);
    sigmoid = (sigd.*(2*pi*fflap*t).^sign)./(sigc+sigd.*(2*pi*fflap*t).^sign);
%     sigmoid = 1; % wing is sufficiently damped, in principle could
%     discard ramp function
    
    % LOCAL FLAPPING 
    phi     = flapamp/180*pi.*sin(2*pi*fflap*t).* sigmoid;
    theta   = 0;
    gamma   = 0;

    % GLOBAL BODY ROTATION 
    if periodic == 1
        globalangle(1) = rot_vec(1)*rotamp/180*pi.*sin(2*pi*frot*t).*sigmoid;
        globalangle(2) = rot_vec(2)*rotamp/180*pi.*sin(2*pi*frot*t).*sigmoid;
        globalangle(3) = rot_vec(3)*rotamp/180*pi.*sin(2*pi*frot*t).*sigmoid;
    elseif periodic == 0
        globalangle(1) = rot_vec(1)*2*pi*frot*t.*sigmoid;
        globalangle(2) = rot_vec(2)*2*pi*frot*t.*sigmoid;
        globalangle(3) = rot_vec(3)*2*pi*frot*t.*sigmoid;
    else
        disp('value for periodic toggle must be 0 or 1')
    end
%Velocity and acceleration of the body (i.e. center base of plate)
    v0  = [0 0 0];
    dv0 = [0 0 0];

%Shape functions and their spatial derivatives 
    N = shapefunc2(a,b,nodes,xpos,ypos); %generate shape functions 
    N = [N(3,:).';N(4,:).']; %put into matrix form 

    for i = 1:6
        dxi(:,i) = [diff(N(i),x,2);diff(N(i),y,2);2*diff(diff(N(i),x),y)]; %compute second spatial derivative 
    end
    

    

%% Generate function with equations for ODEs
%Delete prior PlateODE.m file     
delete('PlateODE.m')
clear PlateODE
[M K Ma Ia Q] = createODEfile_rotvect(a,b,E,G,nu,h,density,dampingfactor,phi,theta,gamma,globalangle,N,dxi);

pause(2) %make sure file saves before solving the ODE 
iter =1; 
while exist('PlateODE.m', 'file') ~= 2 && iter<5
    pause(2)
    iter = iter+1; 
end 

%% Solve ODE 
%Specify initial conditions and damping matrix
    init = zeros(2*6,1); %zero initial conditions

%Solve ODE
    disp('solving ode')
    pause(1)
    options = odeset('RelTol',1e-5);
    teval = 0:1/sampfreq:3; %time matrix at which solution will be evaluated
    [T,Y] = ode45(@(t,y) PlateODE(t,y,v0,dv0,M,K,Ma,Ia,Q,cycles,frot),teval,init,options);


%% Postprocess results 
disp('postprocessing')
pause(1)
%Specify spatial locations where the solution will be evaluated
    xeval = -a:2*a/10:a;
    yeval = 0:2*b/10:2*b;
    [x,y]=meshgrid(xeval,yeval);

%Evaluate shape functions and their derivatives for strains, and spatial derivatives of strain
    for i = 1:6
        Ndisc(i,:,:) = eval(N(i)); %shape function
        strainxi(i,:,:) = eval(dxi(1,i)); %for normal strain along x axis
        strainyi(i,:,:) = eval(dxi(2,i))'; %for normal strain along y axis
        strainxyi(i,:,:) = eval(dxi(3,i))'; %for shear strain in x-y
    end

%Multiply by solution to ODE and sum over all components to solve for actual strains and displacements
    for j = 1:length(Y(:,1))
        %Multiply over all components
        for i = 1:6
            xxnew1(i,:,:) = Ndisc(i,:,:)*squeeze(Y(j,i));
            strainx1(i,:,:)= strainxi(i,:,:)*squeeze(Y(j,i));
            strainy1(i,:,:)= strainyi(i,:,:)*squeeze(Y(j,i));
            strainxy1(i,:,:)= strainxyi(i,:,:)*squeeze(Y(j,i));        
        end

        %Sum over all components
        z(j,:,:) = squeeze(sum(xxnew1,1));
        strainx(j,:,:) = squeeze(sum(strainx1,1))*-h/2;
        strainy(j,:,:) = squeeze(sum(strainy1,1))*-h/2;
        strainxy(j,:,:) = squeeze(sum(strainxy1,1))*-h/2;
       
    end

%Evaluate flapping input, phi(t)
    t = T;
    if phi ~= 0
        flappingangle = eval(phi);
    else
        flappingangle = zeros(length(t),1);
    end

%Evaluate rotating input, phi(t)
%     t = T;
%     if globalangle ~= 0
%         rotationangle = eval(globalangle);
%     else
%         rotationangle = zeros(length(t),1);
%     end
    
file_loc = ''; % use for example: 'data\';, to create a folder with all the results
filename = sprintf('MATLAB_%d%d%dper%d_flapf%d.mat',rot_vec(3)*frot,rot_vec(2)*frot,rot_vec(1)*frot, periodic,fflap);
save( [file_loc, filename] ,'flappingangle','x','y','z','strainx','strainy','strainxy','T')
