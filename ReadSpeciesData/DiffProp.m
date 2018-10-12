%-------------------------------------------------------
%    T   = Temperature (vector optionally)
%    Sp  = struct
%        Sp.pol  : polynomial factors, size: [1:14]
%        Sp.Ts   : 'Switching' temperature [K], size: [1]
%        Sp.Mass : Molar Mass [kg/mole], size: [1]
%-------------------------------------------------------
function [output] = DiffProp(T,Sp, iSp,Prop,Dim)
global Runiv
output = zeros(1,length(T));
switch Prop    

    case 'Diffusivity'   
        for i = 1:length(T)
            coeff = Sp.diff(iSp,Sp.difford(iSp));
            for j = Sp.difford(iSp)-1:-1:1
                coeff = coeff*log(T)+Sp.diff(iSp,j);
            end;
            output(i) = exp(coeff);
        end  
        
    otherwise
        output = [ ];        
end
return