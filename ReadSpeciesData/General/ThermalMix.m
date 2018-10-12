function [h,u,s,Cp, Cv ] = ThermalMix(T,Y,P,Sp)
% ThermalMix  computes specific (if Y is fraction) or total (if Y is mass)
%             species properties [h,u,s,Cp,Cv].
%   quantities on  input (T,Y,P,Sp), T maybe an array
%                 output (h,u,s,Cp,Cv) where s is the thermal part of the entropy. Nb. Integral Cp/T 
global Pref Runiv
% hif=zeros(size(Y,1),size(T,2));
n=size(Y,1);m=size(Y,2);
if m > n
    Y=Y';
end
hif=zeros(size(Y,1),size(T,2));
uif=hif;Cpi=hif;Cvi=hif;
M = [Sp.Mass]';
Mave = 1/sum(Y./M);
X = (Y./M)*Mave
for i=1:length(Sp)
    Pi = P*X(i)/sum(X);
    if (Pi <=0 ) 
        Pi=Pref;
    end
    Rgi = Runiv/M(i);
    sif(i,:)= SNasa(T,Sp(i));%-Rgi*log(Pi/Pref);
    hif(i,:)= HNasa(T,Sp(i));
    uif(i,:)= UNasa(T,Sp(i));
    Cpi(i,:)= CpNasa(T,Sp(i));
    Cvi(i,:)= CvNasa(T,Sp(i));
end
h=Y'*hif;u=Y'*uif;Cv=Y'*Cvi;Cp=Y'*Cpi;s=Y'*sif;

end

