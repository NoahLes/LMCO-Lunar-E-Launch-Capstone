% chatgpt was used to debug the code and for minor formatting help

clear all
clc
format long
%%%%%%%%%%%%%%%%%%%%%%% Multistage Picture and Explanation %%%%%%%%%%%%%%%%%%%%
%   R Coil 1        Coil 2         Coil 3        Coil 4
%   ^  L_c           L_c            L_c          L_c
%   |      mult_dist      mult_dist     mult_dist      mult_dist
%   |------|------|------|-------|------|------|------|------|
%   |______|      |______|       |______|      |______|      |
%   |        L_a                                             | 
%   |      |----|                                            |
%   |      |____|                                            | 
%   |      Projectile                                        |
%   |                                                        |
%   |________________________________________________________| --> z
%
%    This is an axisymmetric view of the coilgun, where the horizontal axis
%    is the axial or z-axis, and the vertical axis is the radial axis. Imagine 
%    this picture rotated into/out of the page 360 degrees. 
%
%    Each coil is L_c long, and spaced a distance of mult_dist apart, once
%    the front of the projectile has passed a trip point you set in if
%    statements in the time loop, the next coil will fire. 
%     
%    After the front of the projectile has reached the end of the domain it
%    will be put back at the left of the domain maintining velocity and
%    tracking overall position. It will travel through the 4 coils again,
%    but overall these would be coils 5, 6, 7, and 8. This process
%    continues until the stoptime is reached. The time loop if statement
%    can also be modified to terminate after a certain number of cools,
%    because total number of coils is also tracked. 
%     
%    This could kind of be described as periodic bounday conditions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Efficiency Factors %%%%%%
eta_v_est = 1;
eta_v = 1;

eta_F_est = 1;
eta_F = 1;

eta_B_est = 1;
eta_B = 1;

% first run all will be 1, and you will calculate eta_v, eta_B, and eta_F
eta_v_calc = 0;
eta_B_calc = 0;
eta_F_calc = 0;

%%%%%%%%%%%%%%%%%%%%% READ INPUT PARAMETERS FROM INPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = readtable('INPUTPARAMETERS.csv'); %Read input file
inputs = inputs{:,2};
%%% Coil Properties %%%
N_c = inputs(1);  % Number of coil windings
d_c = inputs(2); % Coil wire diameter
Lay_c = inputs(3); % Number of layers in the coil

%%% Projectile Properties %%%
R_a = inputs(4); % Armature Inner Radius
L_a = inputs(5); % Armature Length
d_a = inputs(6); % Armature wire diameter
N_a = inputs(7); % Number of Armature Windings 
Res_a = inputs(8); % Resistance of Armature
m = inputs(9); % Projectile Mass
pos_0 = inputs(10); %Initial Position of the front of the projectile relative to the end of the first coil
Lay_a = inputs(11); % Number of layers in the armature
%   |------|------|------|-------|------|------|------|------|
%   |______|      |______|       |______|      |______|      |
%   |        L_a                                             | 
%   |      |----|                                            |
%   |      |____|                                            | 
%   |      Projectile                                        |
%   |-----------|                                            |
%   |  ^                                                     |
%   |  |                                                     |
%   |pos_0 is this distance__________________________________|

%%% Drive Circuit Properties %%%
Res = inputs(12); % Resistance of Drive Circuit
C = inputs(13); % Capacitance of Drive Circuit
V_0 = inputs(14); % Initial Voltage of Drive Circuit

%%% Multi-stage Properties %%%
mult_dist = inputs(15); % Distance between stages
stages = inputs(16); % Number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Simulation Properties %%%%%%%%%%%%%%%%%
StopTime = 0.005; % stop time
dt = 1e-6;   % time step
entries = round(StopTime/dt);
time = linspace(0,StopTime,entries);
simTime = 0;
coil_fire_Time = 0;
loops = 1;
drho = 0.0001;
dz = 0.0001;
%%%% Physical Constants %%%
mu_o = (4.*pi.*10.^-7); % magnetic permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Projectile Initializations %%%%
pos = pos_0;
c_d = 0.42; % coefficient of drag for a hemisphere
A_a = pi.*(R_a).^2; % armature area
A_aa = L_a.*d_a.*Lay_a; % area of current carrying in armature
indu_a = mu_o.*(N_a.^2.*A_a)./L_a * 10; % armature inductance
indu_a = 14.583e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Coil Initializations %%%%
R_c = 0.0179;
L_c = 0.04036; % coil length needs to same as was used in Ansys so can't change
A_c = pi.*(R_c).^2; % coil area
indu_c = mu_o.*(N_c.^2.*A_c)./L_c; % coil inductance
indu_c = 60*10^-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_domain_size = (4*L_c) + (4*mult_dist);
r_domain_size = (R_c);
%%%% Multistage Initializations %%%%
coils = 1; % initializing total number of coils that have fired
coil = 1; % initializing tracking which of the 4 coils is firing
z_pos = pos_0; % initializing z_pos which tracks where in the domain the front of the projectile is
% as opposed to pos, which tracks the total position that the projectile
% has traveled
fourstage_pos = pos_0;
%%%%%%%%%%%%%%%%%%%%%

%%%% Circuit Properties %%%
Res_d = 0.1; % resistance of diode 
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% General Initializing %%%%%%%%%%%%%%%%%
%%% Arrays %%%
velocity = zeros(1,entries);
position = zeros(1,entries);
acceleration = zeros(1,entries);
I_armature = zeros(1, entries);
V_armature = zeros(1, entries);
coil_current = zeros(stages, entries);
cap_voltage = zeros(stages, entries);
r_domain_entries = round(r_domain_size/drho);
z_domain_entries = 3214;
Bz_t = zeros(r_domain_entries, z_domain_entries);
Br_t = zeros(r_domain_entries, z_domain_entries);

Bz_tones = zeros(r_domain_entries, z_domain_entries);
Br_tones = zeros(r_domain_entries, z_domain_entries);
Bz_ttwos = zeros(r_domain_entries, z_domain_entries);
Br_ttwos = zeros(r_domain_entries, z_domain_entries);
Bz_tthrees = zeros(r_domain_entries, z_domain_entries);
Br_tthrees = zeros(r_domain_entries, z_domain_entries);
Bz_tfours = zeros(r_domain_entries, z_domain_entries);
Br_tfours = zeros(r_domain_entries, z_domain_entries);

Flux_t = zeros(1, entries);
Force = zeros(1, entries);
max_Br = zeros(1, entries);
max_Bz = zeros(1, entries);
z_Bs = [dz:dz:z_domain_size];
rho_Bs = [drho:drho:r_domain_size];
cur_density = zeros(1,entries)

%%% Scalars %%%
r_flux = 0;
z_flux = 0;
r_force = 0
z_force = 0;
fourstage_num = 0;
flip = 0;
coil_phase2_time = 0;
Q = 0;
zeta = 0;
omega_d = 0;
track1 = 0;
track2 = 0;
track3 = 0;
track4 = 0;
coil_1_time = 0;
coil_2_time = 0;
coil_3_time = 0;
coil_4_time = 0;

volt_coil1_time = 0;
volt_coil2_time = 0;
volt_coil3_time = 0;
volt_coil4_time = 0;

flux_ct_layer1 = 0;
flux_ct_layer2 = 0;

coil1_phase2_time = 0;
coil2_phase2_time = 0;
coil3_phase2_time = 0;
coil4_phase2_time = 0;


coil1 = 1; %coil 1 starts on all others start off
coil2 = 0;
coil3 = 0;
coil4 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Magnetic Field Load in %%%%%%%%%%%%%%%%%
Br = readtable('SubscaleBrField.csv');
Bz = readtable('SubscaleBzField.csv');
Br = Br(2:end, 2:end);
Bz = Bz(2:end, 2:end);
z = Br{:,1};
Br_vals = Br{:,2:end};
Bz_vals = Bz{:,2:end};
Br_mirrored = [fliplr(-Br_vals), Br_vals]; 
Bz_mirrored = [fliplr(Bz_vals),  Bz_vals];
num_pad = (z_domain_entries/2) - 637 ;
zero_pad = zeros(size(Br_vals, 1), num_pad);
Br_padded = [zero_pad, Br_mirrored, zero_pad];
Bz_padded = [zero_pad, Bz_mirrored, zero_pad];

save('Br_Ansys_subscale.mat', 'Br_padded');
save('Bz_Ansys_subscale.mat', 'Bz_padded');

load('Bz_Ansys_subscale.mat');
load('Br_Ansys_subscale.mat');

Bz_profile = Bz_padded/100;
Br_profile = Br_padded/100;

max_Br_t = 0;
max_Bz_t = 0;

L_c_entries = round(L_c/dz);
mult_dist_entries = round(mult_dist/dz);
initial_shift = round(L_c_entries/2 - (length(Bz_profile(1,:)))/2);

Br_profile1 = circshift(Br_profile, [0, initial_shift]);
Bz_profile1 = circshift(Bz_profile, [0, initial_shift]);

Br_profile2 = circshift(Br_profile, [0, initial_shift + L_c_entries + mult_dist_entries]); % shift columns
Bz_profile2 = circshift(Bz_profile, [0, initial_shift + L_c_entries + mult_dist_entries]);

Br_profile3 = circshift(Br_profile, [0, initial_shift + 2*L_c_entries + 2*mult_dist_entries]);
Bz_profile3 = circshift(Bz_profile, [0, initial_shift + 2*L_c_entries + 2*mult_dist_entries]);

Br_profile4 = circshift(Br_profile, [0, initial_shift + 3*L_c_entries + 3*mult_dist_entries]);
Bz_profile4 = circshift(Bz_profile, [0, initial_shift + 3*L_c_entries + 3*mult_dist_entries]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% TIME LOOP BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:entries
    if fourstage_num > 0 && z_pos > L_c && flip == 1%logic for when coil 1 is fired after projectile has passed coil 4
        coil = 1;
        coils = coils + 1;
        coil1 = 1;
        track1 = 0;
        volt_coil1_time = 0;
        coil1_phase2_time = 0;
        coil_1_time = 0;
        flip = 0;
    end
    if z_pos > (L_c + mult_dist + L_c + L_c/2) && coil == 1 % logic/position for when coil 2 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil2 = 1;
        track2 = 0;
        volt_coil2_time = 0;
        coil2_phase2_time = 0;
        coil_2_time = 0;
    end
    if z_pos > (L_c*3 + mult_dist*2 + L_c/2) && coil == 2 % logic/position for when coil 3 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil3 = 1;
        track3 = 0;
        volt_coil3_time = 0;
        coil3_phase2_time = 0;
        coil_3_time = 0;
    end
    if z_pos > (L_c*4 + mult_dist*3 + L_c/2) && coil == 3 % logic/position for when coil 4 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil4 = 1;
        track4 = 0;
        volt_coil4_time = 0;
        coil4_phase2_time = 0;
        coil_4_time = 0;

    end

    if z_pos > z_domain_size-2*dz    % logic for when projectile is put back to the left of the domain
        fourstage_num = fourstage_num + 1;
        z_pos = 2*dz;
        track1 = 0;
        flip = 1;
        coil4 = 0;  %need to turn off all coils from previous domain, somewhat unrealistic but best I could think of
        coil3 = 0;
        coil2 = 0;
        coil1 = 0;
        Bz_tones = zeros(r_domain_entries, z_domain_entries);
        Br_tones = zeros(r_domain_entries, z_domain_entries);
        Bz_ttwos = zeros(r_domain_entries, z_domain_entries);
        Br_ttwos = zeros(r_domain_entries, z_domain_entries);
        Bz_tthrees = zeros(r_domain_entries, z_domain_entries);
        Br_tthrees = zeros(r_domain_entries, z_domain_entries);
        Bz_tfours = zeros(r_domain_entries, z_domain_entries);
        Br_tfours = zeros(r_domain_entries, z_domain_entries);
    end
    % Circuit calculations
    % know circuit is underdamped
    if coil1 == 1 || coil2 == 1 || coil3 == 1 || coil4 == 1
        alpha = Res/(2*indu_c);
        omega_0 = 1/(sqrt(indu_c*C));
        zeta = (Res/2)*sqrt(C/indu_c);
        omega_d = omega_0*sqrt(1-zeta^2);
        b = [0; (-V_0/indu_c)];
        s_1 = -alpha + sqrt(-(omega_0^2 - alpha^2));
        s_2 = -alpha -sqrt(-(omega_0^2 -alpha^2));
        A = [1, 1,;
        s_1, s_2]\b;
        A_1 = A(2);
        A_2 = A(1);
        compute_current_phase1 = @(t) real(A_1*(exp(s_1*t)) + A_2*(exp(s_2*t)));
        compute_voltage_phase1 = @(t) V_0*exp(-zeta*omega_0*t).*(cos(omega_d*t) + (zeta/sqrt(1 - zeta^2))*sin(omega_d*t));
        if coil1 == 1
           I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase1(volt_coil1_time);
           volt_now = compute_voltage_phase1(coil_1_time);
           coil_current((fourstage_num*4) + 1,loops) = I_coil_fire_now;
           cap_voltage((fourstage_num*4) + 1,loops) = volt_now;
           coil_1_time = coil_1_time + dt;
           volt_coil1_time = volt_coil1_time + dt;
           track1 = track1 + 1;
        end
        if coil2 == 1
            I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase1(volt_coil2_time);
            volt_now = compute_voltage_phase1(coil_2_time);
            coil_current((fourstage_num*4) + 2,loops) = I_coil_fire_now;
            cap_voltage((fourstage_num*4) + 2,loops) = volt_now;
            coil_2_time = coil_2_time + dt;
            volt_coil2_time = volt_coil2_time + dt;
            track2 = track2 + 1;
        end
        if coil3 == 1
            I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase1(volt_coil3_time);
            volt_now = compute_voltage_phase1(coil_3_time);
            coil_current((fourstage_num*4) + 3,loops) = I_coil_fire_now;
            cap_voltage((fourstage_num*4) + 3,loops) = volt_now;
            coil_3_time = coil_3_time + dt;
            volt_coil3_time = volt_coil3_time + dt;
            track3 = track3 + 1;
        end
        if coil4 == 1
            I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase1(volt_coil4_time);
            volt_now = compute_voltage_phase1(coil_3_time);
            coil_current((fourstage_num*4) + 4,loops) = I_coil_fire_now;
            cap_voltage((fourstage_num*4) + 4,loops) = volt_now;
            coil_4_time = coil_4_time + dt;
            volt_coil4_time = volt_coil4_time + dt;
            track4 = track4 + 1;
        end
    end
    % circuit switches to slightly different because of a diode which is
    % modeled by just switching to a slightly different RLC circuit
    if coil1 == 2
        alpha_2 = (Res+Res_d)/(2*indu_c);
        omega_02 = 1/(sqrt(indu_c*C));
        b_2 = [init_current1; init_indu_volt1];
        s_1_2 = -alpha_2 + sqrt(-(omega_02^2 - alpha_2^2));
        s_2_2 = -alpha_2 -sqrt(-(omega_02^2 -alpha_2^2));
        A_2 = [1, 1,;
        s_1_2, s_2_2]\b_2;
        A_1_2 = A_2(2);
        A_2_2 = A_2(1);
        compute_current_phase2 = @(t) real(A_1_2*(exp(s_1_2*t)) + A_2_2*(exp(s_2_2*t)));
        I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase2(coil1_phase2_time);        
        coil_current((fourstage_num*4) + 1,loops) = I_coil_fire_now;
        volt_now = compute_voltage_phase1(volt_coil1_time);
        cap_voltage((fourstage_num*4) + 1,loops) = volt_now;
        volt_coil1_time = volt_coil1_time + dt;
        coil1_phase2_time = coil1_phase2_time + dt;
    end
    if coil2 == 2
        alpha_2 = (Res + Res_d)/(2*indu_c);
        omega_02 = 1/(sqrt(indu_c*C));
        b_2 = [init_current2; init_indu_volt2];
        s_1_2 = -alpha_2 + sqrt(-(omega_02^2 - alpha_2^2));
        s_2_2 = -alpha_2 -sqrt(-(omega_02^2 -alpha_2^2));
        A_2 = [1, 1,;
        s_1_2, s_2_2]\b_2;
        A_1_2 = A_2(2);
        A_2_2 = A_2(1);
        compute_current_phase2 = @(t) real(A_1_2*(exp(s_1_2*t)) + A_2_2*(exp(s_2_2*t)));
        I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase2(coil2_phase2_time);
        volt_now = compute_voltage_phase1(volt_coil2_time);
        coil_current((fourstage_num*4) + 2,loops) = I_coil_fire_now;
        cap_voltage((fourstage_num*4) + 2,loops) = volt_now;
        coil2_phase2_time = coil2_phase2_time + dt;
        volt_coil2_time = volt_coil2_time + dt;
    end
    if coil3 == 2
        alpha_2 = (Res + Res_d)/(2*indu_c);
        omega_02 = 1/(sqrt(indu_c*C));
        b_2 = [init_current3; init_indu_volt3];
        s_1_2 = -alpha_2 + sqrt(-(omega_02^2 - alpha_2^2));
        s_2_2 = -alpha_2 -sqrt(-(omega_02^2 -alpha_2^2));
        A_2 = [1, 1,;
        s_1_2, s_2_2]\b_2;
        A_1_2 = A_2(2);
        A_2_2 = A_2(1);
        compute_current_phase2 = @(t) real(A_1_2*(exp(s_1_2*t)) + A_2_2*(exp(s_2_2*t)));
        I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase2(coil3_phase2_time);
        volt_now = compute_voltage_phase1(volt_coil3_time);
        coil_current((fourstage_num*4) + 3,loops) = I_coil_fire_now;
        cap_voltage((fourstage_num*4) + 3,loops) = volt_now;
        coil3_phase2_time = coil3_phase2_time + dt;
        volt_coil3_time = volt_coil3_time + dt;
    end
    if coil4 == 2
        alpha_2 = (Res + Res_d)/(2*indu_c);
        omega_02 = 1/(sqrt(indu_c*C));
        b_2 = [init_current4; init_indu_volt4];
        s_1_2 = -alpha_2 + sqrt(-(omega_02^2 - alpha_2^2));
        s_2_2 = -alpha_2 -sqrt(-(omega_02^2 -alpha_2^2));
        A_2 = [1, 1,;
        s_1_2, s_2_2]\b_2;
        A_1_2 = A_2(2);
        A_2_2 = A_2(1);
        compute_current_phase2 = @(t) real(A_1_2*(exp(s_1_2*t)) + A_2_2*(exp(s_2_2*t)));
        I_coil_fire_now = (eta_B/eta_B_est)*compute_current_phase2(coil4_phase2_time);
        volt_now = compute_voltage_phase1(volt_coil4_time);
        coil_current((fourstage_num*4) + 4,loops) = I_coil_fire_now;
        cap_voltage((fourstage_num*4) + 4,loops) = volt_now;
        coil4_phase2_time = coil4_phase2_time + dt;
        volt_coil4_time = volt_coil4_time + dt;
    end

    if track1 > 4 % check if switch coil 1 to phase 2
       if (3*coil_current((fourstage_num*4) + 1, loops)- 4*coil_current((fourstage_num*4) + 1, (loops-1)) + coil_current((fourstage_num*4) + 1, (loops-2))) < 0 && coil1 == 1
            coil1 = 2;
            init_cap_volt1 = cap_voltage((fourstage_num*4) + 1, loops);
            init_current1 = coil_current((fourstage_num*4) + 1, loops); 
            init_indu_volt1 = (-(Res)*init_current1 - init_cap_volt1)/indu_c;
            coil_1_time = 0;
        end
    end

    if track2 > 4 % check if switch coil 2 to phase 2
       if (3*coil_current((fourstage_num*4) + 2,loops)- 4*coil_current((fourstage_num*4) + 2,(loops-1)) + coil_current((fourstage_num*4) + 2,(loops-2))) < 0 && coil2 == 1
            coil2 = 2;
            init_cap_volt2 = cap_voltage((fourstage_num*4) + 2, loops);
            init_current2 = coil_current((fourstage_num*4) + 2, loops);      
            init_indu_volt2 = (-(Res)*init_current2 - init_cap_volt2)/indu_c;
            coil_2_time = 0;
        end
    end

    if track3 > 4 % check if switch coil 3 to phase 2
       if (3*coil_current((fourstage_num*4) + 3, loops)- 4*coil_current((fourstage_num*4) + 3, (loops-1)) + coil_current((fourstage_num*4) + 3, (loops-2))) < 0 && coil3 == 1
            coil3 = 2;
            init_cap_volt3 = cap_voltage((fourstage_num*4) + 3, loops);
            init_current3 = coil_current((fourstage_num*4) + 3, loops);     
            init_indu_volt3 = (-(Res)*init_current3 - init_cap_volt3)/indu_c;
            coil_3_time = 0;
        end
    end

    if track4 > 4 % check if switch coil 4 to phase 2
       if (3*coil_current((fourstage_num*4) + 4, loops)- 4*coil_current((fourstage_num*4) + 4, (loops-1)) + coil_current((fourstage_num*4) + 4, (loops-2))) < 0 && coil4 == 1
            coil4 = 2;
            init_cap_volt4 = cap_voltage((fourstage_num*4) + 4, loops);
            init_current4 = coil_current((fourstage_num*4) + 4, loops);     
            init_indu_volt4 = (-(Res)*init_current4 - init_cap_volt4)/indu_c;
            coil_4_time = 0;
        end
    end
    % check if any circuits should be off
    if coil1 == 2
        if coil_current((fourstage_num*4) + 1, loops) < 1
            coil1 = 0;
        end
    end
    if coil2 == 2
        if coil_current((fourstage_num*4) + 2, loops) < 1
            coil2 = 0;
        end
    end
    if coil3 == 2
        if coil_current((fourstage_num*4) + 3, loops) < 1
            coil3 = 0;
        end
    end
    if coil4 == 2
        if coil_current((fourstage_num*4) + 4, loops) < 1
            coil4 = 0;
        end
    end

    % Current Scaling
    % shift B profiles to match the position of each of the 4 coils
    if coil1 == 1 || coil1 == 2
        CB1 = coil_current((fourstage_num*4) + 1, loops);
        Bz_t1 = CB1 .* Bz_profile1;
        Br_t1 = CB1 .* Br_profile1;
        Br_tones = Br_t1;
        Bz_tones = Bz_t1;
    end
    if coil2 == 1 || coil2 == 2
        CB2 = coil_current((fourstage_num*4) + 2, loops);
        Bz_t2 = CB2 .* Bz_profile2;
        Br_t2 = CB2 .* Br_profile2;
        Br_ttwos = Br_t2;
        Bz_ttwos = Bz_t2;
    end
    if coil3 == 1 || coil3 == 2
        CB3 = coil_current((fourstage_num*4) + 3, loops);
        Bz_t3 = CB3 .* Bz_profile3;
        Br_t3 = CB3 .* Br_profile3;
        Br_tthrees = Br_t3;
        Bz_tthrees = Bz_t3;
    end
    if coil4 == 1 || coil4 == 2
        CB4 = coil_current((fourstage_num*4) + 4, loops);
        Bz_t4 = CB4 .* Bz_profile4;
        Br_t4 = CB4 .* Br_profile4;
        Br_tfours = Br_t4;
        Bz_tfours = Bz_t4;
    end
    
    Bz_t = Bz_tones + Bz_ttwos + Bz_tthrees + Bz_tfours;
    Br_t = Br_tones + Br_ttwos + Br_tthrees + Br_tfours;

    % Magnetic Flux through armature calculation
    Flux_ct_layer1 = 0;
    Flux_ct_layer2 = 0;
    r_flux = 0;

    for j = 1:N_a
        r_flux = 0;
        for i = 1:(R_a/drho)
            r_flux = r_flux +drho;
            r_idx = round(r_flux/drho);
            z_idx = round((z_pos - j*d_a)/dz);
                if z_idx == 0
                    z_idx = 1;
                end
                if z_idx < 0 
                    z_idx = round(z_domain_size/dz - abs(z_idx));
                end
            Flux_ct_layer1 = Flux_ct_layer1 + Bz_t(r_idx, z_idx)*r_flux*drho;
        end
    end

    for j = 1:N_a
        r_flux = 0;
        for i = 1:((R_a + d_a)/drho)
            r_flux = r_flux +drho;
            r_idx = round(r_flux/drho);
            z_idx = round((z_pos - j*d_a)/dz);
                if z_idx == 0
                    z_idx = 1;
                end
                if z_idx < 0 
                    z_idx = round(z_domain_size/dz - abs(z_idx));
                end
            Flux_ct_layer2 = Flux_ct_layer2 + Bz_t(r_idx, z_idx)*r_flux*drho;
        end
    end

    Flux(loops) = (Flux_ct_layer1 + Flux_ct_layer2)*2*pi;

    %Voltage induced across armature
    if loops == 1
        emf =  - ((Flux(loops) - 0)/dt);
    end
    if loops == 2
        emf =  -((Flux(loops) - Flux(loops-1))/dt);
    end
    if loops == 3
        emf = - (3*Flux(loops) - 4*Flux(loops-1) + Flux(loops-2)) / (2*dt); 
    end
    if loops == 4
        emf = - (3*Flux(loops) - 4*Flux(loops-1) + Flux(loops-2)) / (2*dt);
    end
    if loops > 4
        emf = - (3*Flux(loops) - 4*Flux(loops-1) + Flux(loops-2)) / (2*dt); 
    end

    V_armature(loops) = (emf);

    %%%%%%%%%%%%% Current in armature calculation (RL Circuit) %%%%%%%%%%%%%%%%%%%

    if loops == 1
        I_armature(loops) = V_armature(loops) / Res_a;  % initilize with Ohm's law
    else
        I_armature(loops) = (indu_a*I_armature(loops-1) + dt*V_armature(loops))/(indu_a + Res_a*dt);
    end

    % Force Calculation
    f = 0;
    r_force = R_a - (d_a*Lay_a);
    for i = 1:(d_a*Lay_a)/drho
        z_force = z_pos - L_a;
        r_force = r_force + drho;
        for j = 1:L_a/dz
            z_force = z_force+dz;
            r_idx = min(max(round(r_force / drho), 1), size(Br_t, 1));
            z_idx = min(max(round(z_force / dz), 1), size(Br_t, 2));
            if z_idx == 0
               z_idx = 1;
            end
            if z_idx < 0 
               z_idx = z_domain_size/dz - abs(z_idx);
            end
            cur_dens = I_armature(loops)/(d_a*Lay_a*L_a);
            cur_density(loops) = cur_dens;
            f = f - (cur_dens*Br_t(r_idx, z_idx))*(drho*dz*r_force);
        end
    end

    f = f*2*pi;
    
    Force(loops) = (eta_F/eta_F_est)*f;
 
    % acceleration
    a = (Force(loops))./m;
    acceleration(loops) = a; 

    %total velocity
    if loops == 1
        velocity(loops) = (eta_v/eta_v_est)*((1/2)*(0 + acceleration(loops)) * dt);
        position(loops) = pos_0 + (1/2)*acceleration(loops)*dt^2;
        z_pos = z_pos + (1/2)*acceleration(loops)*dt^2;
        z_position(loops) = z_pos;
        fourstage_pos = fourstage_pos + (1/2)*acceleration(loops)*dt^2;
    else
        velocity(loops) = (eta_v/eta_v_est)*(velocity(loops-1) + (1/2)*(acceleration(loops-1) + acceleration(loops))*dt);
        position(loops) = position(loops-1) + velocity(loops-1)*dt + (1/2)*acceleration(loops-1)*dt^2;
        z_pos = z_pos + velocity(loops-1)*dt + (1/2)*acceleration(loops-1)*dt^2;
        z_position(loops) = z_pos;
        fourstage_pos = fourstage_pos + velocity(loops-1)*dt + (1/2)*acceleration(loops-1)*dt^2;
    end
    if loops > 1
        dz_now = z_position(loops) - z_position(loops-1); 
        Work(loops) = Work(loops-1) + abs(0.5*(Force(loops) + Force(loops-1))*dz_now);
    else
        Work(loops) = 0;
    end
    loops = loops + 1;
    simTime = simTime + dt
end

% calculate eta_v, eta_F, and eta_B
I_peak1 = max(coil_current(1,:));
eta_B_calc = ((1/2)*(indu_c)*I_peak1^2)/((1/2)*C*V_0^2)
eta_F_calc = Work(round(simTime/dt))/(1/2*indu_c*I_peak1^2)
eta_v_calc = (1/2*m*velocity(round(simTime/dt))^2)/(Work(round(simTime/dt)))

%% 1D Velocity/ Current/ Position/ Plots

% Set default font and font size so plots look good
set(groot, 'DefaultAxesFontName', 'Times New Roman')
set(groot, 'DefaultTextFontName', 'Times New Roman')
set(groot, 'DefaultAxesFontSize', 25)
set(groot, 'DefaultTextFontSize', 25)

% Coil 1 Current
figure;
scatter(time,coil_current(1, :), '.')
title('Coil 1 Current')
xlabel('Time (s)')
ylabel('Current (A)')
grid on

% Velocity
figure;
scatter(time,velocity)
title('Velocity')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

% Position
figure;
scatter(time,position)
title('Axial Position')
xlabel('Time (s)')
ylabel('Axial Position (m)')

%% Supplemental Plots 

% Armature Current
figure;
scatter(time,I_armature, '.')
title('Armature Current')
xlabel('Time (s)')
ylabel('Current (A)')


% Current Density
figure;
scatter(time,cur_density)
title('Current density')
xlabel('Time (s)')
ylabel('Current density A/m^2')

% Acceleration
figure;
scatter(time,acceleration)
title('Acceleration')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

% Force
figure;  
scatter(time, Force, '.');  
title('Axial Force'); 
xlabel('Time (s)')
ylabel('Force (N)')

% Flux
figure;
scatter(time,Flux, '.');
title('Flux layer 1')
xlabel('Time (s)')
ylabel('Flux (Wb or V$\cdot$s)', 'Interpreter', 'latex')

% Armature Voltage
figure;
scatter(time,V_armature, '.');
title('Armature Voltage')
xlabel('Time (s)')
ylabel('Armature Voltage (V)')

% Capacitor Voltage
figure;
scatter(time,cap_voltage(1, :), '.');
title('Capacitor 1 Voltage')
xlabel('Time (s)')
ylabel('Voltage (V)') 
grid on
% 
% figure;
% scatter(time,cap_voltage(2, :), '.');
% title('Capacitor 2 Voltage')
% xlabel('Time (s)')
% ylabel('Voltage (V)') 
% 
% figure;
% scatter(time,cap_voltage(3, :), '.');
% title('Capacitor 3 Voltage')
% xlabel('Time (s)')
% ylabel('Voltage (V)') 
% 
% figure;
% scatter(time,cap_voltage(4, :), '.');
% title('Capacitor 4 Voltage')
% xlabel('Time (s)')
% ylabel('Voltage (V)') 

% figure;
% scatter(time,cap_voltage(5, :), '.');
% title('Capacitor 5 Voltage')
% xlabel('Time (s)')
% ylabel('Voltage (V)')

% % Coil 2 Current
% figure;
% scatter(time,coil_current(2, :), '.')
% title('Coil 2 Current')
% xlabel('Time (s)')
% ylabel('Current (A)')
% 
% % Coil 3 Current
% figure;
% scatter(time,coil_current(3, :), '.')
% title('Coil 3 Current')
% xlabel('Time (s)')
% ylabel('Current (A)')
% 
% % Coil 4 Current
% figure;
% scatter(time,coil_current(4, :), '.')
% title('Coil 4 Current')
% xlabel('Time (s)')\
% ylabel('Current (A)')

% % Coil 5 Current
% figure;
% scatter(time,coil_current(5, :), '.')
% title('Coil 5 Current')
% xlabel('Time (s)')
% ylabel('Current (A)')