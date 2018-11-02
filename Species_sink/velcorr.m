%% To correct the pressure and the volocity see eq 6.24, 6.25
  
function [p u pc] = velcorr(NPI, pc, p, u, relax_pc, d_u)
    for I=2:NPI+1
        i = I;
        
        p(I) = p(I) + relax_pc*pc(I);
        
        if I~= 2
            u(i) = u(i) + d_u(i)*(pc(I-1) - pc(I));
        end
               
    end
end