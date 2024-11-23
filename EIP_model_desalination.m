clear all; clc;
% %%
% % Constants definition
F         = 96485;                  % Faraday constant
RT        = 8.3144 * 298;           % gas constant*temperature
V_T       = RT / F;                 % thermal voltage
E         = 300;                    % the pore ion correlation energy
C_St_vol  = 1.45e8;                 % F m-3, Zero-charge capacitance
Alfa      = 30;                     % Factor for charge dependence of Stern capacitance
L_elec    = 450e-6;                 % Electrode thickness
p_pore    = 0.3;                    % pore volume fraction
p_IEP     = 0.6;                    % polymer volume fraction
X_poly    = 2500;                   % charge density
Eta       = 0.05;                   % diffusion reduction factor
Di_Na     = Eta* 1.33 * 10^-9;      % m2 s-1, Na+ effective diffusion coefficient
Di_Cl     = Eta* 2.03 * 10^-9;      % m2 s-1, Cl- effective diffusion coefficient
z_Na      = 1;                      % Na+ valance
z_Cl      = -1;                     % Cl- valance
% operational conditions
c_di_0    = 100;                    % diluate concentration
c_con_0   = 100;                    % concentrate concentration
c_feed_di = c_di_0;                 % initial concentration in diluate tank
c_feed_con= c_con_0;                % initial concentration in concentrate tank
c_sp_di   = c_feed_di;              % initial concentration in diluate spacer channel
c_sp_con  = c_feed_con;             % initial concentration in concentrate spacer channel
n         = 1;                      % cycle number 
TT        = 400;                    % half cycle time
np        = 100;                    % position grid number
mp        = 100;                    % time grid number
dx        = L_elec / np;            % position grid
dTT       = TT / mp;                % time grid
I_max_0   = 20 * 0.9;               % A m-2, current density
L_spacer  = 10e-3;                  % m, thickness of spacer channel
Area      = 7e-4;                   % m2, electrode area
flowrate  = 5e-6/60;                % m3 s-1, flowrate
V_sp      = L_spacer*Area;          % m3, spacer channel volume
V_feed    = 30e-6-V_sp;             % m3, feed tank volume
tau_sp    = V_sp/flowrate;          % s, hydraulic retention time in the spacer
tau_feed  = V_feed /flowrate;       % s, hydraulic retention time in the tank
EER       = 0.0035;                 % external resistance
%%
% create list to store results
c_IEP_Na_list            = zeros(2*n*mp,np+1);
c_IEP_Cl_list            = zeros(2*n*mp,np+1);
c_pore_Na_list           = zeros(2*n*mp,np+1);
c_pore_Cl_list           = zeros(2*n*mp,np+1);
c_v_list                 = zeros(2*n*mp,np+1);
phi_IEP_list             = zeros(2*n*mp,np+1);
phi_interfacial_list     = zeros(2*n*mp,np+1);
phi_Donnan_list          = zeros(2*n*mp,np+1);
phi_Stern_list           = zeros(2*n*mp,np+1);
phi_sp_half_list         = zeros(2*n*mp,1);
phi_EER_list             = zeros(2*n*mp,1);
phi_carbon_list          = zeros(2*n*mp,1);
j_ion_flux_list          = zeros(2*n*mp,1);
Current_list             = zeros(2*n*mp,1);
Aver_sigma_list          = zeros(2*n*mp,1);
j_Na_left_flux_list      = zeros(2*n*mp,1);
j_Na_right_flux_list     = zeros(2*n*mp,1);
j_Cl_left_flux_list      = zeros(2*n*mp,1);
j_Cl_right_flux_list     = zeros(2*n*mp,1);
c_feed_di_list           = zeros(2*n*mp,1);
c_feed_con_list          = zeros(2*n*mp,1);
c_sp_di_list             = zeros(2*n*mp,1);
c_sp_con_list            = zeros(2*n*mp,1);
%%
% Initialize
c_pore_Na_ini = (0:np)*(c_con_0-c_di_0)/np+c_di_0;
c_pore_Cl_ini = (0:np)*(c_con_0-c_di_0)/np+c_di_0;
c_IEP_Na_ini = c_pore_Na_ini.*exp(-z_Na.*asinh(-X_poly./(2.*c_pore_Na_ini)));
c_IEP_Cl_ini = c_pore_Cl_ini.*exp(-z_Cl*asinh(-X_poly./(2.*c_pore_Cl_ini)));
c_v_ini =  (0:np)*(c_con_0-c_di_0)/np+c_di_0;
phi_IEP_ini = asinh(-X_poly./(2.*c_pore_Na_ini));
phi_interfacial_ini = asinh(X_poly./(2.*c_v_ini));
phi_Donnan_ini = asinh(-(c_pore_Na_ini-c_pore_Cl_ini)./(2.*c_pore_Na_ini));
phi_Stern_ini = -(c_pore_Na_ini-c_IEP_Cl_ini)./(C_St_vol.*V_T./F+Alfa.*(c_pore_Na_ini-c_IEP_Cl_ini).^2.*V_T./F);

x0 = [c_IEP_Na_ini , c_IEP_Cl_ini, phi_IEP_ini,   0,0];
options = optimoptions('fsolve', 'MaxFunEvals', 100000, 'Maxiter', 100000, 'Display', 'off','Algorithm', 'trust-region','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12);
[x, fval, exitflag] = fsolve(@(x) Donnandialysis(x, x0, np, dx, X_poly, Di_Na, Di_Cl, c_di_0,c_con_0,  z_Na, z_Cl), x0, options);
disp(exitflag)

c_IEP_Na_ini = x(1:np+1);
c_IEP_Cl_ini = x(np+2:2*np+2);
phi_IEP_ini =x(2*np+3:3*np+3);
Na_ion_flux_donnan_dialysis = x(3*(np+1)+2);
c_v_ini = sqrt( (c_IEP_Na_ini+c_IEP_Cl_ini).^2 -X_poly^2)./2;
phi_interfacial_ini = asinh(X_poly./2./c_v_ini);

x0 = [c_pore_Na_ini, c_pore_Cl_ini,  phi_Donnan_ini, phi_Stern_ini,0];% Initial conditions vector
options = optimoptions('fsolve', 'MaxFunEvals', 100000, 'Maxiter', 100000, 'Display', 'off','Algorithm', 'trust-region','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12);

I_max =0;
charging_discharge = 1;
[x, fval, exitflag] = fsolve(@(x) deal_equations_ini(x, x0, np, V_T, F, C_St_vol, Alfa, z_Na, z_Cl,E,phi_IEP_ini,c_v_ini, phi_interfacial_ini), x0, options);
c_pore_Na_ini = x(1:np+1);
c_pore_Cl_ini = x(np+2:2*np+2);
phi_Donnan_ini =  x(2*np+3:3*np+3);
phi_Stern_ini = x(3*np+4:4*np+4);
phi_carbon_ini = x(4*(np+1)+1);

j_ion_flux_ini = 0;
Aver_sigma_ini = 0;

%%
% main 
x0 = [c_IEP_Na_ini , c_IEP_Cl_ini, c_pore_Na_ini, c_pore_Cl_ini, c_v_ini, phi_IEP_ini, phi_interfacial_ini, phi_Donnan_ini, phi_Stern_ini, phi_carbon_ini, j_ion_flux_ini,Aver_sigma_ini,0,0,0,0 ];% Initial conditions vector

for  ii = 1:n

    charging_discharge = 1;
   
    for j = 1:mp
      
        I_max =I_max_0 *(1-exp(-j*30/mp)-exp(-(100-j)*30/mp));% current density, goes back to zero
        [x, fval, exitflag] = fsolve(@(x) deal_equations(x, x0, np, dx, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na, Di_Cl, I_max, L_elec, charging_discharge, c_sp_di,c_sp_con,  Alfa, z_Na, z_Cl,E), x0, options);
        % disp(exitflag)
        x0 = x;
        j_ion_flux_list((ii-1)*2*mp+j) = x(9*(np+1)+2);
        phi_carbon_list((ii-1)*2*mp+j) = V_T*x(9*(np+1)+1);
        Aver_sigma_list((ii-1)*2*mp+j) = x(9*(np+1)+3);
        Current_list((ii-1)*2*mp+j) = I_max;
    
        c_IEP_Na_list((ii-1)*2*mp+j,1:np+1)= x(1:np+1);
        c_IEP_Cl_list((ii-1)*2*mp+j,1:np+1)= x(np+2:2*np+2);
        c_pore_Na_list((ii-1)*2*mp+j,1:np+1)= x(2*np+3:3*np+3);
        c_pore_Cl_list((ii-1)*2*mp+j,1:np+1)= x(3*np+4:4*np+4);
        c_v_list((ii-1)*2*mp+j,1:np+1) = x(4*np+5:5*np+5);
        phi_IEP_list((ii-1)*2*mp+j,1:np+1) = x(5*np+6:6*np+6);
        phi_interfacial_list((ii-1)*2*mp+j,1:np+1) =x(6*np+7:7*np+7);
        phi_Donnan_list((ii-1)*2*mp+j,1:np+1) = x(7*np+8:8*np+8);
        phi_Stern_list((ii-1)*2*mp+j,1:np+1) = x(8*np+9:9*np+9);
        phi_sp_half_list((ii-1)*2*mp+j) = -I_max_0*L_spacer/4/20/sqrt(Di_Cl*Di_Na)/c_sp_di/F;
        phi_EER_list((ii-1)*2*mp+j) = -I_max_0*EER./(RT./F);
        j_Na_left_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+4);
        j_Na_right_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+5);
        j_Cl_left_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+6);
        j_Cl_right_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+7);

        c_sp_di=c_sp_di+(1/tau_sp*(c_feed_di-c_sp_di)-x(9*(np+1)+4)/L_spacer)*dTT;
        c_sp_con=c_sp_con+(1/tau_sp*(c_feed_con-c_sp_con)+x(9*(np+1)+5)/L_spacer)*dTT;
        c_feed_di=c_feed_di+1/tau_feed*(c_sp_di-c_feed_di)*dTT;
        c_feed_con=c_feed_con+1/tau_feed*(c_sp_con-c_feed_con)*dTT;
        c_feed_di_list((ii-1)*2*mp+j) = c_feed_di;
        c_feed_con_list((ii-1)*2*mp+j) = c_feed_con;
        c_sp_di_list((ii-1)*2*mp+j) = c_sp_di;
        c_sp_con_list((ii-1)*2*mp+j) = c_sp_con;
    end

    charging_discharge = 2;
  
    for j = mp+1:2*mp
        I_max =I_max_0 *(exp(-(j-100)*30/mp)+exp(-(200-j)*30/mp)-1);% current density, goes back to zero
        [x, fval, exitflag] = fsolve(@(x) deal_equations(x, x0, np, dx, dTT,X_poly, V_T, F, C_St_vol, p_pore, p_IEP, Di_Na, Di_Cl, I_max, L_elec, charging_discharge,c_sp_di,c_sp_con,  Alfa, z_Na, z_Cl,E), x0, options);
        % disp(exitflag)
        x0 = x;
        j_ion_flux_list((ii-1)*2*mp+j) = x(9*(np+1)+2);
        phi_carbon_list((ii-1)*2*mp+j) = V_T*x(9*(np+1)+1);
        Aver_sigma_list((ii-1)*2*mp+j) = x(9*(np+1)+3);
        Current_list((ii-1)*2*mp+j) = I_max;

        c_IEP_Na_list((ii-1)*2*mp+j,1:np+1)= x(1:np+1);
        c_IEP_Cl_list((ii-1)*2*mp+j,1:np+1)= x(np+2:2*np+2);
        c_pore_Na_list((ii-1)*2*mp+j,1:np+1)= x(2*np+3:3*np+3);
        c_pore_Cl_list((ii-1)*2*mp+j,1:np+1)= x(3*np+4:4*np+4);
        c_v_list((ii-1)*2*mp+j,1:np+1) = x(4*np+5:5*np+5);
        phi_IEP_list((ii-1)*2*mp+j,1:np+1) = x(5*np+6:6*np+6);
        phi_interfacial_list((ii-1)*2*mp+j,1:np+1) =x(6*np+7:7*np+7);
        phi_Donnan_list((ii-1)*2*mp+j,1:np+1) = x(7*np+8:8*np+8);
        phi_Stern_list((ii-1)*2*mp+j,1:np+1) = x(8*np+9:9*np+9);
        phi_sp_half_list((ii-1)*2*mp+j) = I_max_0*L_spacer/4/20/sqrt(Di_Cl*Di_Na)/c_sp_di/F;
        phi_EER_list((ii-1)*2*mp+j) = I_max_0*EER./(RT./F);
        j_Na_left_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+4);
        j_Na_right_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+5);
        j_Cl_left_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+6);
        j_Cl_right_flux_list((ii-1)*2*mp+j) =  x(9*(np+1)+7);
 
        c_sp_di=c_sp_di+(1/tau_sp*(c_feed_di-c_sp_di)-x(9*(np+1)+4)/L_spacer)*dTT;
        c_sp_con=c_sp_con+(1/tau_sp*(c_feed_con-c_sp_con)+x(9*(np+1)+5)/L_spacer)*dTT;
        c_feed_di=c_feed_di+1/tau_feed*(c_sp_di-c_feed_di)*dTT;
        c_feed_con=c_feed_con+1/tau_feed*(c_sp_con-c_feed_con)*dTT;
        c_feed_di_list((ii-1)*2*mp+j) = c_feed_di;
        c_feed_con_list((ii-1)*2*mp+j) = c_feed_con;
        c_sp_di_list((ii-1)*2*mp+j) = c_sp_di;
        c_sp_con_list((ii-1)*2*mp+j) = c_sp_con;

    end
end

phi_carbon_list = phi_IEP_list + phi_interfacial_list + phi_Donnan_list +phi_Stern_list;
phi_CC_list = 0.36 +(phi_carbon_list(:,1) + phi_sp_half_list + phi_EER_list).*RT./F;

