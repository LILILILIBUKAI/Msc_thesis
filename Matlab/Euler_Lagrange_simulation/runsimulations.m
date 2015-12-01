% -------------------------
% Run Euler-Lagrange simulation of flat plate and plot strain at wing base
% TMohren 2015/12/01
% Based on model by:
% Eberle, A., Dickerson, B., Reinhall, P., and Daniel, T. 
% A new twist on gyroscopic sensing: body rotations lead to torsion in flapping, flexing insect wings. 
% Journal of The Royal Society Interface, 12(104):20141088, 2015.
% -------------------------


% % % Constant rotation
[x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([0,0,1],3,0,10,6);
% [x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([0,1,0],3,0,10,6);
% [x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([1,0,0],3,0,10,6);

% % %Periodic rotation
% [x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([0,0,1],3,1,10,6);
% [x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([0,1,0],3,1,10,6);
% [x,y,z,exx,eyy,exy,T,phi] = flappermodel_rotvect([1,0,0],3,1,10,6);


%% plot results of this simulation
figure()

eyy     = eyy(1:2000,:,:);
T       = T(1:2000);

subplot(311)
plot(T,eyy(:,2,2))
    xlabel('time [s]')
    ylabel('\epsilon_y at left side of base')

subplot(312)
d_strain = eyy(:,2,2)-eyy(:,10,2);
plot(T,d_strain)
    xlabel('time [s]')
    ylabel('\Delta \epsilon_{y} at base')

subplot(313)
Fs      = 1000;
L       = length(T);
NFFT    = 2^nextpow2(L); 
Y       = fft(d_strain,NFFT)/L;
f       = Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(Y(1:NFFT/2+1))) 
    axis([0 60 0 1e-5])
    ylabel('fft \Delta \epsilon_y')
    xlabel('frequency [Hz]')