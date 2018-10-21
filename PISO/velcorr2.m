%% To correct the pressure and the volocity see eq 6.24, 6.25
  
function [p, u] = velcorr2(NPI, pc,pc2, p, u, b_p, d_u)
    for I=2:NPI+1
        i = I;
        
        p(I) = p(I) + pc(I) + pc2(I);
        
        if I~= 2
            u(i) = u(i) + d_u(i)*(pc(I-1) - pc(I)) + b_p(i) + d_u(i)*(pc2(I-1)- pc2(I));
        end
               
    end
end