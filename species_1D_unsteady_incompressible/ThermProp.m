%-------------------------------------------------------
%    T   = Temperature (vector optionally)
%    Sp  = struct
%        Sp.pol  : polynomial factors, size: [1:14]
%        Sp.Ts   : 'Switching' temperature [K], size: [1]
%        Sp.Mass : Molar Mass [kg/mole], size: [1]
%-------------------------------------------------------
function [output] = ThermProp(T,Sp,Prop,Dim)
global Runiv

output = zeros(1,length(T));
switch Prop    
    case 'Cv'
        for i = 1:length(T)
            if (T(i) <= Sp.Ts)
                a = Sp.pol(1,:);
            else
                a = Sp.pol(2,:);
            end
            coeff = a(5);
            for j = 4:-1:1
                coeff = coeff*T(i)+a(j);
            end
            coeff = coeff-1;
            output(i) = coeff;
        end
        if (strcmp(Dim,'Mole') == 1)
            output = output*Runiv;
        elseif (strcmp(Dim,'Mass') == 1)
            output = (output*Runiv)./Sp.Mass;
        end
        
    case 'Cp'
        for i = 1:length(T)
            if (T(i) <= Sp.Ts)
                a = Sp.pol(1,:);
            else
                a = Sp.pol(2,:);
            end
            coeff = a(5);
            for j = 4:-1:1
                coeff = coeff*T(i)+a(j);
            end
            output(i) = coeff;
        end
        if (strcmp(Dim,'Mole') == 1)
            output = output*Runiv;
        elseif (strcmp(Dim,'Mass') == 1)
            output = (output*Runiv)./Sp.Mass;
        end
        
    case 'Gamma'
        Cv = ThermProp(T,Sp,'Cv',Dim);
        Cp = ThermProp(T,Sp,'Cp',Dim);
        output = Cp./Cv; 
    
    case 'Enthalpy'
        for i = 1:length(T)
            if (T(i) <= Sp.Ts)
                a = Sp.pol(1,:);
            else
                a = Sp.pol(2,:);
            end
            coeff = a(5)/5;
            for j = 5:-1:2
                coeff = coeff*T(i)+a(j-1)/(j-1);
            end
            coeff = coeff*T(i) + a(6);
            output(i) = coeff;
        end
        if (strcmp(Dim,'Mole') == 1)
            output = output*Runiv;
        elseif (strcmp(Dim,'Mass') == 1)
            output = (output*Runiv)./Sp.Mass;
        end
        
    case 'Entropy'
        for i = 1:length(T)
            if (T(i) <= Sp.Ts)
                a = Sp.pol(1,:);
            else
                a = Sp.pol(2,:);
            end
            coeff = a(5)/4;
            for j = 4:-1:2
                coeff = coeff*T(i)+a(j)/(j-1);
            end
            coeff = coeff*T(i) + a(1)*log(T(i)) + a(7);
            output(i) = coeff;
        end
        if (strcmp(Dim,'Mole') == 1)
            output = output*Runiv;
        elseif (strcmp(Dim,'Mass') == 1)
            output = (output*Runiv)./Sp.Mass;
        end
        
    case 'Gibbs'
        H = ThermProp(T,Sp,'Enthalpy',Dim);
        S = ThermProp(T,Sp,'Entropy',Dim);
        output = H-T.*S;
        
    case 'Conductivity'
        for i = 1:length(T)
            coeff = Sp.cond(Sp.condord);
            for j = Sp.condord-1:-1:1
                coeff = coeff*log(T)+Sp.cond(j);
            end;
            output(i) = exp(coeff);
        end
        
    case 'Viscosity'   
        for i = 1:length(T)
            coeff = Sp.visc(Sp.viscord);
            for j = Sp.viscord-1:-1:1
                coeff = coeff*log(T)+Sp.visc(j);
            end;
            output(i) = exp(coeff);
        end        
        
    otherwise
        output = [ ];
end

return