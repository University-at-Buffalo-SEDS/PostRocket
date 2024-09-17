function [s,x,o,t] = HRAP_Function( MotorConfiguration,RunTime,BurnTime,Timestep, Display,...
                                    GrainOD,GrainID,GrainLength,CEfficiency,AmbientPressure, ...
                                    TankVolume,OxidizerFill,TankTemp, ...
                                    ThroatDiameter,NozzleExit,NozzleCd,NozzleEfficiency, ...
                                    NumberofInjectors,InjectorDiameter,InjectorCd)
    addpath(genpath("core/"))
    addpath(genpath("HRAP/"))
    addpath(genpath("HRAP"))
    addpath(genpath("core"))
    addpath(genpath("core"))

    s = struct(); % input structure
    o = struct(); % output structure
    x = struct(); % state structure
    
    

    % initialize variables for HRAP
    
    % Loading Paraffin as our current fuel source
    paraffin = load('Paraffin.mat').s;
    s.prop_nm = paraffin.prop_nm;
    s.prop_k = paraffin.prop_k;
    s.prop_M = paraffin.prop_M;
    s.prop_OF = paraffin.prop_OF;
    s.prop_Pc = paraffin.prop_Pc;
    s.prop_Reg = paraffin.prop_Reg;
    s.prop_Rho = paraffin.prop_Rho;
    s.prop_T = paraffin.prop_T;
    s.opt_OF = paraffin.opt_OF;

    s.grn_OD = GrainOD;
    s.grn_ID = GrainID;
    s.grn_L = GrainLength;
    s.cstar_eff = CEfficiency;

    s.mtr_nm = MotorConfiguration;
    s.tmax = RunTime;
    s.tburn = BurnTime;
    s.dt = Timestep;
    s.Pa = AmbientPressure;

    s.tnk_V = TankVolume;
    s.cmbr_V = GrainLength*0.25*pi*(GrainOD)^2;
    
    s.regression_model = @(s,x) shift_OF(s,x);

    s.noz_Cd = NozzleCd;
    s.noz_thrt = ThroatDiameter;
    s.noz_ER = (NozzleExit)^2/(s.noz_thrt)^2;
    s.noz_eff = NozzleEfficiency;
    s.inj_CdA = 0.25*pi*(InjectorDiameter)^2*InjectorCd;
    s.inj_N = NumberofInjectors;
    
    s.vnt_S = 0;
    % For now, ignoring vent calculations
    % if VentState == "None"
    %     s.vnt_S = 0;
    % elseif VentState == "External"
    %     s.vnt_S = 1;
    %     s.vnt_CdA = 0.25*pi*(VentDiameter*u.vnt_D)^2*VentCd;
    % else
    %     s.vnt_S = 2;
    %     s.vnt_CdA = 0.25*pi*(VentDiameter*u.vnt_D)^2*VentCd;
    % end
    
    s.mp_calc = 0;
    % For now, ignoring the built-in CG calculations
    %     if MassProperties == 1
    %     s.mtr_m = MotorMass*u.mtr_m;
    %     s.mtr_cg = MotorCG*u.mtr_cg;
    %     s.tnk_X = TankLocation*u.tnk_X;
    %     s.cmbr_X = GrainLocation*u.cmbr_X;
    %     s.mp_calc = 1;
    % else
    %     s.mp_calc = 0;
    % end
    
    x.T_tnk = TankTemp;
    x.ox_props                  = NOX(x.T_tnk);
    x.m_o                       = (OxidizerFill)*s.tnk_V*x.ox_props.rho_l + (1-OxidizerFill)*s.tnk_V*x.ox_props.rho_v;
    
    x.P_tnk                     = x.ox_props.Pv;
    x.P_cmbr                    = s.Pa; % Atmospheric starting chamber pressure
    x.mdot_o                    = 0;
    x.mLiq_new                  = (s.tnk_V - (x.m_o/x.ox_props.rho_v))/((1/x.ox_props.rho_l)-(1/x.ox_props.rho_v));
    x.mLiq_old                  = x.mLiq_new + 1;
    x.m_f                       = 0.25*pi*(s.grn_OD^2 - s.grn_ID^2)*s.prop_Rho*s.grn_L;
    x.m_g                       = 1.225*(s.cmbr_V - 0.25*pi*(s.grn_OD^2 - s.grn_ID^2)*s.grn_L);

    x.OF                        = 0;

    x.mdot_f                    = 0;
    x.mdot_n                    = 0;
    x.rdot                      = 0;
    x.grn_ID                    = s.grn_ID;
    
    t                           = 0;
    
    o.t                         = zeros(1,s.tmax/s.dt + 1);
    o.m_o                       = zeros(1,s.tmax/s.dt + 1);
    o.P_tnk                     = zeros(1,s.tmax/s.dt + 1);
    o.P_cmbr                    = zeros(1,s.tmax/s.dt + 1);
    o.mdot_o                    = zeros(1,s.tmax/s.dt + 1);
    o.mdot_f                    = zeros(1,s.tmax/s.dt + 1);
    o.OF                        = zeros(1,s.tmax/s.dt + 1);
    o.grn_ID                    = zeros(1,s.tmax/s.dt + 1);
    o.mdot_n                    = zeros(1,s.tmax/s.dt + 1);
    o.rdot                      = zeros(1,s.tmax/s.dt + 1);
    o.m_f                       = zeros(1,s.tmax/s.dt + 1);
    o.F_thr                     = zeros(1,s.tmax/s.dt + 1);
    o.dP                        = zeros(1,s.tmax/s.dt + 1);
    
    o.m_o(1)                    = x.m_o;
    o.P_tnk(1)                  = x.P_tnk;
    o.P_cmbr(1)                 = x.P_cmbr;
    o.mdot_o(1)                 = x.mdot_o;
    o.mdot_f(1)                 = x.mdot_f;
    o.OF(1)                     = x.OF;
    o.grn_ID(1)                 = x.grn_ID;
    o.mdot_n(1)                 = x.mdot_n;
    o.rdot(1)                   = x.rdot;
    o.m_f(1)                    = x.m_f;
    
    if s.mp_calc == 1
        o.m_t                   = zeros(1,s.tmax/s.dt + 1);
        o.cg                    = zeros(1,s.tmax/s.dt + 1);
    
        mp                          = mass(s,x);
    
        o.m_t(1)                = mp(1);
        o.cg(1)                 = mp(2);
    end
    
    % RUN HRAP SIMULATION
    
    [s,x,o,t] = sim_loop(s,x,o,t);
    
    % Sim Results
    
    totalImpulse = trapz(o.t(o.F_thr>0),o.F_thr(o.F_thr>0));

    [motorClass,percent] = impulse(totalImpulse);
    %Display Motor Performance
    
    if Display
        burnTime = o.t(sum(o.F_thr>0));
        peakF_thr = max(o.F_thr);
        averageF_thr = mean(o.F_thr(o.F_thr>0));
        peakPressure = max(o.P_cmbr)/10^5;
        averagePressure = mean(o.P_cmbr(o.P_cmbr>0))/10^5;
        fuelConsumed = (o.m_f(1)-o.m_f(end));
        oxidizerConsumed = (o.m_o(1)-o.m_o(end));
        averageOF = oxidizerConsumed/fuelConsumed;
        intPressure = trapz(o.t,o.P_cmbr);
        Cstar = intPressure*0.25*pi*s.noz_thrt^2/(fuelConsumed+oxidizerConsumed);
        specificImpulse = totalImpulse/(((oxidizerConsumed+fuelConsumed))*9.81);
        fuelLeft = x.grn_ID*100;
    
        fprintf(['Motor Name: %s\n    ' ...
                 'Propellant: %s\n    ' ...
                 'Oxidizer Tank Volume: %0.0f cc\n    ' ...
                 'Burn Time: %0.3f s\n    ' ...
                 'Peak Thrust: %0.1f N\n    ' ...
                 'Average Thrust: %0.1f N\n    ' ...
                 'Total Impulse: %0.2f N-s\n    ' ...
                 'Peak Chamber Pressure: %0.3f bar\n    ' ...
                 'Average Chamber Pressure: %0.3f bar\n    ' ...
                 'Port Diameter at Burnout: %0.3f cm\n    ' ...
                 'Fuel Consumed: %0.3f kg\n    ' ...
                 'Oxidizer Consumed: %0.3f kg\n    ' ...
                 'Average OF Ratio: %0.3f\n    ' ...
                 'Characteristic Velocity: %0.1f m/s\n    ' ...
                 'Specific Impulse: %0.1f s\n    ' ...
                 'Motor Classification: %3.0f%% %s%0.0f\n'], s.mtr_nm,s.prop_nm,s.tnk_V*10^6,burnTime,peakF_thr,averageF_thr,totalImpulse,peakPressure,averagePressure,fuelLeft,fuelConsumed,oxidizerConsumed,averageOF,Cstar,specificImpulse,percent,motorClass,averageF_thr);
    end