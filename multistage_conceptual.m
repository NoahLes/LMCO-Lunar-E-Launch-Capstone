% chatgpt was used to help debug the code and for minor formatting suppport

clear all
clc
format long

%%%%%%%%%%%%%%%%%%%%%%% Multistage Picture and Explanation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   R Coil 1        Coil 2         Coil 3        Coil 4        Coil 5       Coil 6        Coil 7 
%   ^  L_c           L_c            L_c          L_c
%   |      mult_dist      mult_dist     mult_dist      mult_dist    mult_dist      mult_dist     mult_dist
%   |------|------|------|-------|------|------|------|------|------|------|------|------|------|-------|
%   |______|      |______|       |______|      |______|      |______|      |______|      |______|       |
%   |        L_a                                                                                        |                                          
%   |      |----|                                                                                       |                                           
%   |      |____|                                                                                       |                                         
%   |      Projectile                                                                                   |                                         
%   |                                                                                                   |
%   |___________________________________________________________________________________________________| --> z
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
%    tracking overall position. It will travel through the 7 coils again,
%    but overall these would be coils 8, 9, 10, etc. This process
%    continues until the s is reached. The time loop if statement
%    can also be modified to terminate after a certain number of cools,
%    because total number of coils is also tracked. 
%     
%    This could kind of be described as periodic bounday conditions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Simulation Properties %%%%%%%%%%%
StopTime = 2; % stop time
dt = 1e-7;   % time step
entries = round(StopTime/dt);
time = linspace(0,StopTime,entries);
simTime = 0;
coil_fire_Time = 0;
loops = 1;

%%% Conversion Factors %%%%%
alpha_B = 1 % provided by subscale
alpha_F = 1 % provided by subscale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Physical Constants %%%%%%%%%%%%%%%%%
mu_o = (4.*pi.*10.^-7); % magnetic permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Projectile Properties %%%%%%%%%%%%%%%%
Res_a = 0.01; % resistance of wire
m = 6000; % projectile mass
pos_0 = 0; % initial position of the front of the projectile
vel_0 = 0; % initial veloctiy of the projectile
pos = pos_0;
N_a = 3; % number of armature windings
R_a =  0.75; % armature radius
L_a = 1; % armature length
d_a = 0.36; % armature wire diameter
Lay_a = 1; % number of armature layers
A_a = pi.*(R_a).^2; % armature area
A_aa = L_a.*d_a.*Lay_a; % area of current carrying in armature
indu_a = mu_o.*(N_a.^2.*A_a)./L_a; % armature inductance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Coil Properties %%%%%%%%%%%%%%%%%%%%%
N_c = 1; % number of coil windings
R_c = 0.7995; % coil radius
L_c = 0.36; % coil length
d_c = 0.3604; % coil wire guage 
Lay_c = 1; % number of coil layers
A_c = pi.*(R_c).^2; % coil area
indu_c = 2.23e-6; % coil inductance
% indu_c = 60*10^-6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Multistage Properties %%%%%%%%%%%%%%%%%%%%%
mult_dist = 0.12; % distance between stages
coils = 1; % initializing total number of coils that have fired
coil = 1; % initializing tracking which of the 7 coils is firing
z_pos = pos_0; % initializing z_pos which tracks where in the domain the back of the projectile is
% as opposed to pos, which tracks the total position that the projectile
% has traveled
sevenstage_pos = pos_0;
stages = 800; % number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_domain_size = (7*L_c) + (7*mult_dist);
r_domain_size = (R_c);

drho = 0.01;
dz = 0.001;

z_domain_entries = round(z_domain_size/dz);
r_domain_entries = round(r_domain_size/drho);
%%%%%%%%%%%%%%%% Circuit Properties %%%%%%%%%%%%%%%%%%%
Res = 0.010;  % total resistance of drive circuit
Res_d = 0.01; % resistance of diode 
C = 0.0616; % capacitance of drive circuit
V_0 = 49000;  % initial voltage in the capacitors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Initializing %%%%%%%%%%%%%%%%%
% Arrays
velocity = zeros(1,entries);
velocity(1) = vel_0;
position = zeros(1,entries);
acceleration = zeros(1,entries);
I_armature = zeros(1, entries);
V_armature = zeros(1, entries);
coil_current = zeros(1, entries);
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
Bz_tfives = zeros(r_domain_entries, z_domain_entries);
Br_tfives = zeros(r_domain_entries, z_domain_entries);
Bz_tsixs = zeros(r_domain_entries, z_domain_entries);
Br_tsixs = zeros(r_domain_entries, z_domain_entries);
Bz_tsevens = zeros(r_domain_entries, z_domain_entries);
Br_tsevens = zeros(r_domain_entries, z_domain_entries);

Flux_t = zeros(1, entries);
Force = zeros(1, entries);
z_Bs = [dz:dz:z_domain_size];
rho_Bs = [0:drho:r_domain_size];

% Scalars
r_flux = 0;
z_flux = 0;
r_force = 0
z_force = 0;
sevenstage_num = 0;
flip = 0;
coil_phase2_time = 0;
track1 = 0;
track2 = 0;
track3 = 0;
track4 = 0;
track5 = 0;
track6 = 0;
track7 = 0;

coil_1_time = 0;
coil_2_time = 0;
coil_3_time = 0;
coil_4_time = 0;
coil_5_time = 0;
coil_6_time = 0;
coil_7_time = 0;


volt_coil1_time = 0;
volt_coil2_time = 0;
volt_coil3_time = 0;
volt_coil4_time = 0;
volt_coil5_time = 0;
volt_coil6_time = 0;
volt_coil7_time = 0;

coil1_phase2_time = 0;
coil2_phase2_time = 0;
coil3_phase2_time = 0;
coil4_phase2_time = 0;
coil5_phase2_time = 0;
coil6_phase2_time = 0;
coil7_phase2_time = 0;

coil1 = 1; %coil 1 starts on all others start off
coil2 = 0;
coil3 = 0;
coil4 = 0;
coil5 = 0;
coil6 = 0;
coil7 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coil1_midpoint = L_c/2;
coil2_midpoint = L_c + mult_dist + L_c/2 - L_a/3;
coil3_midpoint = L_c*2 + mult_dist*2 + L_c/2 - L_a/3;
coil4_midpoint = L_c*3 + mult_dist*3 + L_c/2 - L_a/3;
coil5_midpoint = L_c*4 + mult_dist*4 + L_c/2 - L_a/3;
coil6_midpoint = L_c*5 + mult_dist*5 + L_c/2 - L_a/3; 
coil7_midpoint = L_c*6 + mult_dist*6 + L_c/2; %rerun now because L_a/3 changed

switch_location1 = coil1_midpoint;
switch_location2 = coil2_midpoint;
switch_location3 = coil3_midpoint;
switch_location4 = coil4_midpoint;
switch_location5 = coil5_midpoint;
switch_location6 = coil6_midpoint;
switch_location7 = coil7_midpoint;

%%%%%%%%%%%%%%%% Magnetic Field Load in %%%%%%%%%%%%%%%%%
Br = readtable('ConceptualBrField.csv');
Bz = readtable('ConceptualBzField.csv');
Br = Br(2:10:800, 2:end);
Bz = Bz(2:10:800, 2:end);
z = Br{:,1};
Br_vals = Br{:,2:end};
Bz_vals = Bz{:,2:end};
Br_good = [fliplr(-Br_vals), Br_vals]; 
Bz_good = [fliplr(Bz_vals),  Bz_vals];

Br_profile = Br_good/(1*10^5);
Bz_profile = Bz_good/(1*10^5);

L_c_entries = round(L_c/dz);
mult_dist_entries = round(mult_dist/dz);
initial_shift = round(L_c_entries/2 - (length(Bz_profile(1,:)))/2);
% Loaded in magnetic field has coil in the center, so have to shift that
% coil to where we are placing coil 1, coil2, ..., coil 7
Br_profile1 = circshift(Br_profile, [0, initial_shift]);
Bz_profile1 = circshift(Bz_profile, [0, initial_shift]);

Br_profile2 = circshift(Br_profile, [0, initial_shift + L_c_entries + mult_dist_entries]); % shift columns
Bz_profile2 = circshift(Bz_profile, [0, initial_shift + L_c_entries + mult_dist_entries]);

Br_profile3 = circshift(Br_profile, [0, initial_shift + 2*L_c_entries + 2*mult_dist_entries]);
Bz_profile3 = circshift(Bz_profile, [0, initial_shift + 2*L_c_entries + 2*mult_dist_entries]);

Br_profile4 = circshift(Br_profile, [0, initial_shift + 3*L_c_entries + 3*mult_dist_entries]);
Bz_profile4 = circshift(Bz_profile, [0, initial_shift + 3*L_c_entries + 3*mult_dist_entries]);

Br_profile5 = circshift(Br_profile, [0, initial_shift + 4*L_c_entries + 4*mult_dist_entries]); 
Bz_profile5 = circshift(Bz_profile, [0, initial_shift + 4*L_c_entries + 4*mult_dist_entries]);

Br_profile6 = circshift(Br_profile, [0, initial_shift + 5*L_c_entries + 5*mult_dist_entries]);
Bz_profile6 = circshift(Bz_profile, [0, initial_shift + 5*L_c_entries + 5*mult_dist_entries]);

Br_profile7 = circshift(Br_profile, [0, initial_shift + 6*L_c_entries + 6*mult_dist_entries]);
Bz_profile7 = circshift(Bz_profile, [0, initial_shift + 6*L_c_entries + 6*mult_dist_entries]);
%%%%% Circuit Calcs and setup Analytical Formulas %%%%%%%%%
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
compute_voltage_phase1 = @(t) V_0 * exp(-zeta * omega_0 * t).*(cos(omega_d * t) + (zeta / sqrt(1 - zeta^2)) * sin(omega_d * t));
I_peak = (V_0)/(omega_d*indu_c)*exp(-alpha*pi/(2*omega_d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:entries
    if sevenstage_num > 0 && z_pos > switch_location1 && flip == 1%logic for when coil 1 is fired after projectile has passed coil 4
        coil = 1;
        coils = coils + 1;
        coil1 = 1;
        track1 = 0;
        coil1_phase2_time = 0;
        coil_1_time = 0;
        flip = 0;
        coil7 = 0;
        Bz_tsevens(:) = 0;
        Br_tsevens(:) = 0; 
    end
    if z_pos > (switch_location2) && coil == 1 % logic/position for when coil 2 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil2 = 1;
        track2 = 0;
        coil2_phase2_time = 0;
        coil_2_time = 0;
        coil1 = 0;
        Bz_tones(:) = 0;
        Br_tones(:) = 0;
    end
    if z_pos > (switch_location3) && coil == 2 % logic/position for when coil 3 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil3 = 1;
        track3 = 0;
        coil3_phase2_time = 0;
        coil_3_time = 0;
        coil2 = 0;
        Bz_ttwos(:) = 0;
        Br_ttwos(:) = 0;
    end
    if z_pos > (switch_location4) && coil == 3 % logic/position for when coil 4 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil4 = 1;
        track4 = 0;
        coil4_phase2_time = 0;
        coil_4_time = 0;
        coil3 = 0;
        Bz_tthrees(:) = 0;
        Br_tthrees(:) = 0;
    end
    if z_pos > (switch_location5) && coil == 4 % logic/position for when coil 5 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil5 = 1;
        track5 = 0;
        coil5_phase2_time = 0;
        coil_5_time = 0;
        coil4 = 0;
        Bz_tfours(:) = 0;
        Br_tfours(:) = 0;
    end
    if z_pos > (switch_location6) && coil == 5 % logic/position for when coil 6 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil6 = 1;
        track6 = 0;
        coil6_phase2_time = 0;
        coil_6_time = 0;
        coil5 = 0;
        Bz_tfives(:) = 0;
        Br_tfives(:) = 0;
    end
    if z_pos > switch_location7 && coil == 6 % logic/position for when coil 7 is fired
        coil = coil + 1;
        coils = coils + 1;
        coil7 = 1;
        track7 = 0;
        coil7_phase2_time = 0;
        coil_7_time = 0;
    end


    if z_pos > z_domain_size    % logic for when projectile is put back to the left of the domain
        sevenstage_num = sevenstage_num + 1;
        z_pos = 2*dz;
        track1 = 0;
        flip = 1;
        coil6 = 0;
        coil5 = 0;
        coil4 = 0; 
        coil3 = 0;
        coil2 = 0;
        coil1 = 0;
        Bz_tsixs(:) = 0;
        Br_tsixs(:) = 0;     
    end

    % Circuit calculations
    % know circuit is underdamped
    if coil1 == 1 || coil2 == 1 || coil3 == 1 || coil4 == 1 || coil5 == 1 || coil6 == 1 || coil7 == 1; % this if statement is probably unnecessary but I'm afraid to remove it in case things mess up
        if coil1 == 1
           I_coil_fire_now = alpha_B*compute_current_phase1(coil_1_time);
           volt_now = compute_voltage_phase1(coil_1_time);
           coil_1_time = coil_1_time + dt;
           volt_coil1_time = volt_coil1_time + dt;
           track1 = track1 + 1;
        end
        if coil2 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_2_time);
            volt_now = compute_voltage_phase1(coil_2_time);
            coil_2_time = coil_2_time + dt;
            track2 = track2 + 1;
        end
        if coil3 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_3_time);
            volt_now = compute_voltage_phase1(coil_3_time);
            coil_3_time = coil_3_time + dt;
            volt_coil3_time = volt_coil3_time + dt;
            track3 = track3 + 1;
        end
        if coil4 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_4_time);
            volt_now = compute_voltage_phase1(coil_4_time);
            coil_4_time = coil_4_time + dt;
            volt_coil4_time = volt_coil4_time + dt;
            track4 = track4 + 1;
        end
        if coil5 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_5_time);
            volt_now = compute_voltage_phase1(coil_5_time);
            coil_5_time = coil_5_time + dt;
            volt_coil5_time = volt_coil5_time + dt;
            track5 = track5 + 1;
        end
        if coil6 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_6_time);
            volt_now = compute_voltage_phase1(coil_6_time);
            coil_6_time = coil_6_time + dt;
            volt_coil6_time = volt_coil6_time + dt;
            track6 = track6 + 1;
        end
        if coil7 == 1
            I_coil_fire_now = alpha_B*compute_current_phase1(coil_7_time);
            volt_now = compute_voltage_phase1(coil_7_time);
            coil_7_time = coil_7_time + dt;
            volt_coil7_time = volt_coil7_time + dt;
            track7 = track7 + 1;
        end        
    end
    % circuit switches to slightly different because of a diode which is
    % modeled by just switching to a slightly different RLC circuit
    if coil1 == 2
        alpha_2 = (Res + Res_d)/(2*indu_c);
        omega_02 = 1/(sqrt(indu_c*C));
        b_2 = [init_current1; init_indu_volt1];
        s_1_2 = -alpha_2 + sqrt(-(omega_02^2 - alpha_2^2));
        s_2_2 = -alpha_2 -sqrt(-(omega_02^2 -alpha_2^2));
        A_2 = [1, 1,;
        s_1_2, s_2_2]\b_2;
        A_1_2 = A_2(2);
        A_2_2 = A_2(1);
        compute_current_phase2 = @(t) real(A_1_2*(exp(s_1_2*t)) + A_2_2*(exp(s_2_2*t)));
        I_coil_fire_now = alpha_B*compute_current_phase2(coil1_phase2_time);
        volt_now = compute_voltage_phase1(volt_coil1_time);
        volt_coil1_time = volt_coil1_time + dt;
        coil1_phase2_time = coil1_phase2_time + dt;
    end
    if coil2 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil2_phase2_time);
        coil2_phase2_time = coil2_phase2_time + dt;
    end
    if coil3 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil3_phase2_time);
        coil3_phase2_time = coil3_phase2_time + dt;
    end
    if coil4 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil4_phase2_time);
        coil4_phase2_time = coil4_phase2_time + dt;
    end
    if coil5 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil5_phase2_time);
        coil5_phase2_time = coil5_phase2_time + dt;
    end
    if coil6 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil6_phase2_time);
        coil6_phase2_time = coil6_phase2_time + dt;
    end
    if coil7 == 2
        I_coil_fire_now = alpha_B*compute_current_phase2(coil7_phase2_time);
        coil7_phase2_time = coil7_phase2_time + dt;
    end 
  
    if track1 > 4 % check if switch coil 1 to phase 2
       if I_coil_fire_now > I_peak && coil1 == 1
            coil1 = 2;
            init_cap_volt1 = 3.3356e4;
            init_current1 = 3.3732e6; 
            init_indu_volt1 = -3.80084e10;
            coil_1_time = 0;
        end
    end

    if track2 > 4 % check if switch coil 2 to phase 2
       if I_coil_fire_now > I_peak  && coil2 == 1
            coil2 = 2;
            coil_2_time = 0;
        end
    end

    if track3 > 4 % check if switch coil 3 to phase 2
       if I_coil_fire_now > I_peak && coil3 == 1
            coil3 = 2;
            coil_3_time = 0;
        end
    end

    if track4 > 4 % check if switch coil 4 to phase 2
       if I_coil_fire_now > I_peak && coil4 == 1
            coil4 = 2;
            coil_4_time = 0;
        end
    end

    if track5 > 4 % check if switch coil 5 to phase 2
       if I_coil_fire_now > I_peak && coil5 == 1
            coil5 = 2;
            coil_5_time = 0;
        end
    end 

    if track6 > 4 % check if switch coil 7 to phase 2
       if I_coil_fire_now > I_peak && coil6 == 1
            coil6 = 2;
            coil_6_time = 0;
        end
    end 
    if track7 > 4 % check if switch coil 7 to phase 2
       if I_coil_fire_now > I_peak && coil7 == 1
            coil7 = 2;
            coil_7_time = 0;
        end
    end
    % check if any circuits should be off
    if coil1 == 2
        if I_coil_fire_now < 1
            coil1 = 0;
        end
    end
    if coil2 == 2
        if I_coil_fire_now < 1
            coil2 = 0;
        end
    end
    if coil3 == 2
        if I_coil_fire_now < 1
            coil3 = 0;
        end
    end
    if coil4 == 2
        if I_coil_fire_now < 1
            coil4 = 0;
        end
    end
    if coil5 == 2
        if I_coil_fire_now < 1
            coil5 = 0;
        end
    end
    if coil6 == 2
        if I_coil_fire_now < 1
            coil6 = 0;
        end
    end
    if coil7 == 2
        if I_coil_fire_now < 1
            coil7 = 0;
        end
    end
    % Current Scaling
    % scale the magnetic profile that corresponds to the coil that is firing
    if coil1 == 1 || coil1 == 2
        CB1 = I_coil_fire_now;
        Bz_t1 = CB1 .* Bz_profile1;
        Br_t1 = CB1 .* Br_profile1;
        Br_tones = Br_t1;
        Bz_tones = Bz_t1;
    end
    if coil2 == 1 || coil2 == 2
        CB2 = I_coil_fire_now;
        Bz_t2 = CB2 .* Bz_profile2;
        Br_t2 = CB2 .* Br_profile2;
        Br_ttwos = Br_t2;
        Bz_ttwos = Bz_t2;
    end
    if coil3 == 1 || coil3 == 2
        CB3 = I_coil_fire_now;
        Bz_t3 = CB3 .* Bz_profile3;
        Br_t3 = CB3 .* Br_profile3;
        Br_tthrees = Br_t3;
        Bz_tthrees = Bz_t3;
    end
    if coil4 == 1 || coil4 == 2
        CB4 = I_coil_fire_now;
        Bz_t4 = CB4 .* Bz_profile4;
        Br_t4 = CB4 .* Br_profile4;
        Br_tfours = Br_t4;
        Bz_tfours = Bz_t4;
    end
    if coil5 == 1 || coil5 == 2
        CB5 = I_coil_fire_now;
        Bz_t5 = CB5 .* Bz_profile5;
        Br_t5 = CB5 .* Br_profile5;
        Br_tfives = Br_t5;
        Bz_tfives = Bz_t5;
    end  
    if coil6 == 1 || coil6 == 2
        CB6 = I_coil_fire_now;
        Bz_t6 = CB6 .* Bz_profile6;
        Br_t6 = CB6 .* Br_profile6;
        Br_tsixs = Br_t6;
        Bz_tsixs = Bz_t6;
    end
    if coil7 == 1 || coil7 == 2
        CB7 = I_coil_fire_now;
        Bz_t7 = CB7 .* Bz_profile7;
        Br_t7 = CB7 .* Br_profile7;
        Br_tsevens = Br_t7;
        Bz_tsevens = Bz_t7;
    end   
    Bz_t = Bz_tones + Bz_ttwos + Bz_tthrees + Bz_tfours + Bz_tfives + Bz_tsixs + Bz_tsevens;
    Br_t = Br_tones + Br_ttwos + Br_tthrees + Br_tfours + Br_tfives + Br_tsixs + Br_tsevens;
    % Magnetic Flux through armature calculation
    Flux_ct = 0;
    Flux_cz = 0;
    r_flux = 0;
    coil_current(loops) = I_coil_fire_now;

    for j = 1:N_a
        r_flux = 0;
        for i = 1:(R_a/drho)
            r_flux = r_flux + drho;
            r_idx = round(r_flux/drho);
            z_idx = round((z_pos + j*d_a)/dz);
            if z_idx == 0
               z_idx = 1;
            end
            if z_idx < 0 
               z_idx = z_domain_size/dz - abs(z_idx);
            end
            if z_idx > z_domain_entries
                z_idx = z_idx - z_domain_entries ;
            end
            Flux_ct = Flux_ct + Bz_t(r_idx, z_idx)*r_flux*drho;
        end
    end

    Flux_t(loops) = Flux_ct*2*pi;
    %%%%%%%%%%%%% Current in armature calculation %%%%%%%%%%%%%%%%%%%

    I_armature(loops) = - Flux_t(loops)/indu_a;  % superconduncting armature

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Force Calculation %%%%%%%%%%%%%%%%%%%%%%%%%
    f = 0;
    r_force = R_a - (d_a*Lay_a);
    for i = 1:(d_a*Lay_a)/drho
        z_force = z_pos;
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
            if z_idx > z_domain_entries
                z_idx = z_idx - z_domain_entries ;
            end
            cur_dens = I_armature(loops)/(d_a*Lay_a*L_a);
            f = f - (cur_dens*Br_t(r_idx, z_idx))*drho*dz*r_force;
        end
    end

    f = f*2*pi;
    
    Force(loops) = alpha_F*f;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Drag calculation
    % Drag = -0.5*1.225*v^2*c_d*(3*pi*R_a^2);
    % Calculate inductance of coil at this time step
  
    % acceleration
    a = (Force(loops))./m;
    acceleration(loops) = a; 
    if acceleration(loops) > 7135 
        acceleration(loops) = 7135;
    end

    %%%%%%%%%%%%%% Velocity and Position Calculations %%%%%%%%%%% 
    if loops == 1
        velocity(loops) = (1/2)*(0 + acceleration(loops)) * dt;
        position(loops) = pos_0 + (1/2)*acceleration(loops)*dt^2;
        z_pos = z_pos + (1/2)*acceleration(loops)*dt^2;
    else
        velocity(loops) = velocity(loops-1) + (1/2)*(acceleration(loops-1) + acceleration(loops))*dt;
        position(loops) = position(loops-1) + velocity(loops-1)*dt + (1/2)*acceleration(loops-1)*dt^2;
        z_pos = z_pos + velocity(loops-1)*dt + (1/2)*acceleration(loops-1)*dt^2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update time step and iteration number
    if velocity(loops) > 2353
        break
    end
    loops = loops + 1;
    simTime = simTime + dt;

    %%%% only display some sim-times, keeps track of where in the
    %%%% simulation you are at but saves time by not displaying every time
    %%%% step
    tolerance = 1e-6;
     if abs(mod(simTime, 1e-3)) < tolerance
         simTime
     end 
end


%% 1D Velocity/ Acceleration Plotting

% Set default font and font size so plots look good
set(groot, 'DefaultAxesFontName', 'Times New Roman')
set(groot, 'DefaultTextFontName', 'Times New Roman')
set(groot, 'DefaultAxesFontSize', 31)
set(groot, 'DefaultTextFontSize', 31)
% tip_linex = [0.0002177, 0.0585161, 0.108526, 0.143809, 0.170994, 0.193081, 0.224401, 0.239315, 0.24333, ...
%     0.254913, 0.26542, 0.274967, 0.283761, 0.305497, 0.307858, 0.314906, 0.321489, 0.327681, 0.333523, 0.348302, ...
%     0.349974, 0.354964, 0.359714, 0.364302, 0.368626, 0.379737, 0.381052, 0.384916, 0.388638, 0.395638, 0.404581, ...
%     0.408784, 0.4148, 0.41766, 0.425067, 0.428625, 0.43373, 0.442537, 0.445621, 0.452208, 0.457721, 0.464407, 0.471899, ...
%     0.477242, 0.484, 0.488862, 0.49043, 0.495035, 0.49802, 0.507911, 0.510602, 0.514517, 0.519582, 0.524446, 0.531424, ...
%     0.534787, 0.539135, 0.543354, 0.546432, 0.554335, 0.557192, 0.562735, 0.567206, 0.56894, 0.571014, 0.572386, 0.575732, ...
%     0.576926, 0.578191, 0.581407, 0.582491, 0.583764, 0.586864, 0.588018, 0.589133, 0.592123, 0.593359, 0.594314, 0.597199...
%     , 0.598371, 0.599321, 0.602117, 0.603218, 0.604168, 0.606881, 0.607902, 0.608876, 0.611509, 0.612667, 0.617104, ...
%     0.613446, 0.616021, 0.620408, 0.621376, 0.629598, 0.693942, 0.722333, 0.785046, 0.981492];
% 
% tip_liney = [7009.85, 4850.47, 4872.93, 4895.28, 4917.53, 4923.33, 7135, 6400.16, 5017.26, ... 
%     5039.33, 5045.31, 5083.46, 5105.5, 6334.89, 5134.4, 5173.81, 5197.41, 5220.6, 5243.24, 6302.87, ...
%     5283.06, 5312.23, 5315.53, 5340.4, 5365.43, 6215.61, 5415.45, 5450.13, 5474.03, 5515.58, 6177.46, ...
%     5578.12, 5633.46, 5640.55, 6113.91, 5718.37, 5764.4, 6028.27, 5840.62, 5903.08, 5978.07, 6006.18, 6053.12, ...
%     6115.89, 6163.29, 6217.95, 6234.41, 6265.14, 6298.34, 6388.4, 6415.85, 6442.23, 6497.08, 6528.32, 6589.61, ...
%     6623.9, 6653.41, 6695.59, 6711.23, 6771.87, 6800.64, 6833.12, 6853.37, 6865.85, 7135, 6894.92, 6926.17, ...
%     7135, 6931, 6964.66, 7135, 6964.64, 6999.5, 7135, 6994.31, 7028.64, 7135, 7020.45, 7055.3, 7135 ...
%     7043.38, 7078.52, 7135, 7062.91, 7100.23, 7135, 7080.56, 7116.43, 7135, 7093.36, 7134.27, 7135,  7108, 7133.37,...
%     7135, 7135, 7073.39, 6895.4, 6267.29];

skim_tip_linex = [0.0002177, 0.0585161, 0.108526, 0.143809, 0.170994, 0.193081, 0.24333, ...
    0.254913, 0.26542, 0.274967, 0.283761, 0.307858, 0.314906, 0.321489, 0.327681, 0.333523, ...
    0.349974, 0.354964, 0.359714, 0.364302, 0.368626, 0.381052, 0.384916, 0.388638, 0.395638, ...
    0.408784, 0.4148, 0.41766, 0.428625, 0.43373, 0.445621, 0.452208, 0.457721, 0.464407, 0.471899, ...
    0.477242, 0.484, 0.488862, 0.49043, 0.495035, 0.49802, 0.507911, 0.510602, 0.514517, 0.519582, 0.524446, 0.531424, ...
    0.534787, 0.539135, 0.543354, 0.546432, 0.554335, 0.557192, 0.562735, 0.567206, 0.56894, 0.572386, 0.575732, ...
    0.578191, 0.581407, 0.583764, 0.586864, 0.589133, 0.592123,  0.594314, 0.597199...
    , 0.599321, 0.602117,  0.604168, 0.606881, 0.608876, 0.611509, ...
    0.613446, 0.616021, 0.617897, 0.619769, 0.621376, 0.629598, 0.693942, 0.722333, 0.785046, 0.981492];

skim_tip_liney = [7009.85, 4850.47, 4872.93, 4895.28, 4917.53, 4923.33, 5017.26, ... 
    5039.33, 5045.31, 5083.46, 5105.5, 5134.4, 5173.81, 5197.41, 5220.6, 5243.24, ...
    5283.06, 5312.23, 5315.53, 5340.4, 5365.43, 5415.45, 5450.13, 5474.03, 5515.58, ...
    5578.12, 5633.46, 5640.55, 5718.37, 5764.4, 5840.62, 5903.08, 5978.07, 6006.18, 6053.12, ...
    6115.89, 6163.29, 6217.95, 6234.41, 6265.14, 6298.34, 6388.4, 6415.85, 6442.23, 6497.08, 6528.32, 6589.61, ...
    6623.9, 6653.41, 6695.59, 6711.23, 6771.87, 6800.64, 6833.12, 6853.37, 6865.85,  6894.92, 6926.17, ...
     6931, 6964.66, 6964.64, 6999.5, 6994.31, 7028.64,  7020.45, 7055.3, ...
    7043.38, 7078.52, 7062.91, 7100.23, 7080.56, 7116.43, 7095.37, 7135,  7107.88 ,7113.51, 7133.37,...
    7135, 7135, 7073.39, 6895.4, 6267.29];

% Acceleration
figure;
plot(time(1:9831003),acceleration(1:9831003))
hold on
yline(7500, 'r-', '7500', 'linewidth',4, 'FontSize', 31);     
plot(skim_tip_linex, skim_tip_liney, 'g-', 'linewidth', 3)
title('Acceleration (m/s^2)')
xlabel('Time')                
ylabel('Acceleration') 

% Velocity
figure;
plot(time(1:9831003),velocity(1:9831003), 'linewidth', 4)
hold on
yline(2353, 'r-', '2353', 'linewidth', 4, 'FontSize', 31);  
title('Velocity')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

%% Supplemental Plotting that could be useful

% Coil Current
% figure;
% scatter(time,coil_current, '.')
% title('Coil 1 Current')
% xlabel('Time (s)')
% ylabel('Current (A)')

% Armature Current
figure;
scatter(time,I_armature, '.')
title('Armature Current')
xlabel('Time (s)')
ylabel('Current (A)')

% Position
figure;
scatter(time,position)
title('Axial Position')
xlabel('Time (s)')
ylabel('Axial Position (m)')

% Force
figure;  
scatter(time, Force, '.');  
title('Axial Force'); 
xlabel('Time (s)')
ylabel('Force (N)')

% Flux
figure;
scatter(time,Flux_t, '.');
title('Flux')
xlabel('Time (s)')
ylabel('Flux (Wb or V$\cdot$s)', 'Interpreter', 'latex')