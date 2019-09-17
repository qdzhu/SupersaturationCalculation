classdef Cal_Saturation
    properties
        a0;
        a1;
        Bw;
        Aw;
        c;
        bw;
    end
    
    methods
        function obj= Cal_Saturation(T_initial)
            g = 9.8; % gravity coefficient m/s2
            Lw = 2260*1e3; % latent heat of water vapor J/kg
            Ra = 286; % gas constant of air  J/kg*k
            Rv = 461; % gas constant of water vapor J/kg*k
            Cp = lookup_Cp(); % specific heat capacitance of moist air J/Kg*K
            qv = lookup_Qv(T_initial); %mixing ratio of water vapor kg/Kg
            qw = 0; %mixing ratio of water liquid
            phow = 997 ; %density of liquid water kg/m3
            phoa = 1.225; %density of air kg/m3
            T = T_initial+273;
            Ew = lookup_Ew(T_initial);%saturation water vapor Nm-2 
            k = lookup_k(T_initial);%coefficient of air heat conductivity Jm-1s-1K-1
            D = lookup_D(T_initial);%coefficient of water vapor diffusiont in the air m2/s
            obj.Aw = (phow*Lw^2/(k*Rv*T^2)+phow*Rv*T/(Ew*D))^(-1);
            obj.a0 = g*(Lw*Ra/(Cp*Rv*T)-1)/(Ra*T);
            obj.a1 = 1/qv+Lw^2/(Cp*Rv*T^2);
            obj.Bw = 4*pi*phow*obj.Aw/phoa;
            obj.c = 1; %ice particle shape factor characterizing capacitance, c=1 for spheres
            obj.bw = obj.a1*obj.Bw;
        end
        function [t,Sw,Sw_appro] = Main_Cal_Saturation(obj,Uv,Sw_initial,rw_initial,Nw)
            dt = 0.001;%s
            t_total = 200;%s
            n = t_total/dt;
            t = dt:dt:t_total;
            Sw = zeros(n,1);
            Sw(1) = Sw_initial;
            rw2 = zeros(n,1);
            rw2(1) = rw_initial^2;
            for i =1:n-1
                Sw(i+1) = Sw(i)+dt*(Sw(i)+1)*(obj.a0*Uv-obj.a1*obj.Bw*Nw*Sw(i)*sqrt(rw2(i)));
                rw2(i+1) = rw2(i)+2*obj.c*obj.Aw*Sw(i+1)*dt;
            end
            Sqsw = obj.a0*Uv/(obj.bw*Nw*rw_initial);
            thu = (obj.a0*Uv+obj.bw*Nw*rw_initial)^-1;
            C0 = (Sqsw-Sw_initial)/(1+Sw_initial);
            Sw_appro = (Sqsw - C0*exp(-t/thu))./(1+C0*exp(-t/thu));
        end
        function Sqsw = Only_Quasi_Supersaturation(obj,Uv,rw_initial,Nw)
            Sqsw = obj.a0*Uv/(obj.bw*Nw*rw_initial);
        end
    end
    
end

function Qv = lookup_Qv(T_initial)
    if T_initial == 0 % celsius degress
        Qv = 3.84*1e-3; %kg/kg
    end
    if T_initial == 5
    end
end

function Cp = lookup_Cp()
    Cp = 1.005; %heat capacity for dry air kJ/Kg*K
    %cs = 1.005 + 1.82H where 1.005 kJ/kg°C is the heat capacity of dry air, 
    %1.82 kJ/kg°C the heat capacity of water vapor, and H water kg/kg air,
    %rather small
    Cp = Cp*1e3; % conversion from kJ/Kg*K to J/Kg*K
end

function Ew = lookup_Ew(T)
    % August-Roche-Magnus equation
    Ew = 0.61094*exp(17.625*T/(T+243.04)); % T in celsius degree
    Ew = Ew*1e3; %conversion from kPa to Pa
end

function k = lookup_k(T)
    if T == 0
        k = 24.36;%mW/(m K)
    end
    k = k*1e-3; % conversion from mW/mk to J/m s k
end

function D = lookup_D(T)
    if T ==0
        D = 0.219*1e-4; %diffusion coefficient m2/s
    end
end