% testing function for Bloch simulator using RK4
%%
% Author: Mingming Wu
% Date created: 25.11.2013
% Part of the Master Thesis work at ETH Zurich
% Supervision: Markus Weiger

close all, clear all

%% initialization of simulation parameters
% nt: number of time intervals for the simulating the block pulse
nt=200; 
% dt: time interval duration in ms
dt=0.01;
% initialization of phase and magnitude for block pulse array
phase=zeros(1,nt);
ha=ones(1,nt);

%% B1 magnitude/RF pulse, here a block pulse with a certain flip angle is defined
flip = pi/2; % flip angle in rad
omega1 = flip*ha/sum(ha)/dt; % flip angle divided by pulse duration 

%% B0 off resonance frequency
froff=linspace(-10,10,100); % off resonance frequency in kHz

%% T1 and T2 relaxation times 
% if set to 0, infinitely long T1/T2 will be considered
t1=0;
t2=0;

%% initial magnetization
mxin = 0;
myin = 0;
mzin = 1;

% final longitudinal magnetization after application of the RF pulse
mzend=zeros(1,length(froff));


for i=1:length(froff)
[mx,my,mz]=blochRK4(mxin,myin,mzin,omega1,phase,dt,froff(i),t1,t2);
mzend(i)=mz(end);
end

plot(froff,mzend);
xlabel('B0 off resonance (kHz)');
ylabel('M_z');