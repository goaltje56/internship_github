function fi = solve_eq(NPI, aE, aW, aP, b, fi, Istart)
    % equivalence with variables in eq 7.1-7.6
    % BETA = aW(I)
    % D    = aP(I)
    % ALFA = Ari(I)
    % C    = Cri
    % C'   = Cmri(I)
    % b    = b(I)

    Ari(Istart-1) = 0;
    Cmri(Istart-1) = fi(Istart-1);
    
    % forward substitution
    for I = Istart:NPI+1
        Ari(I) = aE(I)/(aP(I)-aW(I)*Ari(I-1));
        Cri    = b(I);
        Cmri(I)= (aW(I)*Cmri(I-1)+Cri)/(aP(I)-aW(I)*Ari(I-1));
    end
    
    % back substitution
    for I = NPI+1:-1:Istart
        if I == NPI +1
            fi(I) = Cmri(I);
        else
            fi(I) = Ari(I)*fi(I+1) + Cmri(I);
        end
    end
end