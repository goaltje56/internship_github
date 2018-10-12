% ChemEquil function to compute chemical equilibrium composition
% Giel Ramaekers, version 17-09-2007

function [Xeq,Teq] = ChemEquil(Xin,Tn,p,Xiss,Tss)

if (nargin == 3)
    Xiss = Xin; Tss = Tn;  % Start solutions
end

warning off all

global Runiv El Sp
pref = 1.01325e5;  % Reference pressure
Nit_max = 100; Nit_max_Xi = 20; Nit_max_T = 20;

% ---------------------------------------------------------------------
Nsp = length(Sp); Ne = length(El); 
Xit0 = sum(Xin);

b0 = zeros(1,Ne);
for i = 1:Ne
    a = zeros(1,Nsp);
    for j = 1:Nsp
        a(j) = Sp(j).elcomp(i);
    end
    b0(i) = dot(Xin,a);
end

h = zeros(1,Nsp);
for j = 1:Nsp
    h(j) = ThermProp(Tn,Sp(j),'Enthalpy','Mole');
end
h0 = dot(Xin,h);

% ---------------------------------------------------------------------
Xin(:) = max(Xiss(:),1e-9); Tn = Tss;
it = 0; Conv = 1;
while ((Conv > 1e-12) && (it < Nit_max))
    Xi_conv = 1; it_sub = 0;
    while ((Xi_conv > 1e-6) && (it_sub < Nit_max_Xi))
        Xi = Xin;
        A = zeros(Nsp+Ne+1); A(1:Nsp,1:Nsp) = eye(Nsp);
        for i = 1:Nsp
            for j = 1:Ne
                A(i,Nsp+j) = -Sp(i).elcomp(j);
                A(Nsp+j,i) = Sp(i).elcomp(j)*Xi(i);
            end
        end
        A(1:Nsp,end) = -1;
        A(end,1:Nsp) = Xi(1:Nsp);
        A(end,end) = -sum(Xi);

        b = zeros(Nsp+Ne+1,1);
        for i = 1:Nsp
            b(i) = -ThermProp(Tn,Sp(i),'Gibbs','Mole')/(Runiv*Tn)-log(p*Xi(i)/(pref*sum(Xi)));
        end
        for i = 1:Ne
            a = zeros(1,Nsp);
            for j = 1:Nsp
                a(j) = Sp(j).elcomp(i);
            end
            b(Nsp+i) = b0(i) - dot(Xi,a);
        end
        b(end) = Xit0 - sum(Xi);

        x = A\b;
        update = min(1,2/max(5*abs(x(end)),max(abs(x(1:Nsp))))); % orig: 5,1
        DeltaXi = zeros(1,Nsp);
        for i = 1:Nsp
            Xin(i) = exp(log(Xi(i))+update*x(i));
            if (isinf(Xin(i)) == 1)||(isnan(Xin(i)) == 1)
                Xeq = Xi; Teq = Tn;
                return
            end            
            DeltaXi(i) = abs((Xin(i)-Xi(i))/Xi(i));
        end
        Xi_conv = max(DeltaXi);
        it_sub = it_sub+1;
    end

    T_conv = 1; it_sub = 0;
    while ((T_conv > 1e-12) && (it_sub < Nit_max_T))  % Iteratie voor temperatuur
        cp = zeros(1,Nsp); h = zeros(1,Nsp);
        for j = 1:Nsp
            cp(j) = ThermProp(Tn,Sp(j),'Cp','Mole');
            h(j) = ThermProp(Tn,Sp(j),'Enthalpy','Mole');
        end
        cpn = dot(Xin,cp); hn = dot(Xin,h); 
        T = Tn+((h0-hn)/cpn);
        T_conv = abs((h0-hn)/(cpn*T));
        Tn = T;
        it_sub = it_sub+1;
    end

%     Conv = max(0,abs(1-update));
    Conv = 0; % Modification by JvO
    for j = 1:Nsp
        Conv = max(Conv,abs((Xin(j)-Xi(j))/max(Xi(j),1e-9)));
    end
    it = it+1;   
end

% fprintf('\nNumber of iterations for chemical equilibrium solver: %3i\n',it);
Xeq = Xin; Teq = Tn;
return