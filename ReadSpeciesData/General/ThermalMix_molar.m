function [h,u,s,Cp, Cv ] = ThermalMix_molar(T,X,P,Sp)
%ThermalMix_molar computes molar (if X is fraction) or total (if X is moles) thermal part of entropy (nb Integral Cp/T)
%             species properties [h,u,s,Cp,Cv].
%   quantities on  input (T,X,P,Sp), T maybe an array
%                 output (h,u,s,Cp,Cv) where s is the entropy. Nb. Integral Cp/T - Ru*log(pi/pref)
global Pref Runiv
hif=zeros(size(X,1),size(T,2));
uif=hif;Cpi=hif;Cvi=hif;sif=hif;
for i=1:length(Sp)
    Pi = P*X(i)/sum(X);
    Mi = Sp(i).Mass;
    sif(i,:)= SNasa(T,Sp(i))-Runiv*log(Pi/Pref);
    hif(i,:)= HNasa(T,Sp(i))*Mi;
    uif(i,:)= UNasa(T,Sp(i))*Mi;
    Cpi(i,:)=CpNasa(T,Sp(i))*Mi;
    Cvi(i,:)=CvNasa(T,Sp(i))*Mi;
end
h=X*hif;u=X*uif;Cv=X*Cvi;Cp=X*Cpi;s=X*sif;
end

