%% To correct the pressure and the volocity see eq 6.24, 6.25
  
function u  = velcorr(NPI, pc, u, d_u)
    for I=2:NPI+1
        i = I;
        
%         p(I) = p(I) + relax_pc*pc(I);
        
        if I~= 2
            u(i) = u(i) + d_u(i)*(pc(I-1) - pc(I));
        end
               
    end
end