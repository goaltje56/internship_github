function u_pseudo  = pseudo(n, x, aW, aE, aP, b, u_pseudo) 
T1 = 0;
T2 = 0;
r  = 1;
    for i = 3:n
        for j = 1:n
            if i == j+1
                T1 = (aW(i)/aP(i))*x(i-1);
            end
            if i == j-1
                T2 = (aE(i)/aP(i))*x(i+1);
            end
        end
        
        u_pseudo(i) = T1 + T2 + b(i)/aP(i);
        T1 = 0;
        T2 = 0;
    end
    
% for i = 1:n
%     if i==1
%         r(i) = b(i) - aP(i)*x(i) + aE(i)*x(i+1);
%     elseif i==n
%         r(i) = b(i) - aP(i)*x(i) + aW(i)*x(i-1);
%     else 
%         r(i) = b(i) - aP(i)*x(i) + aE(i)*x(i+1) + aW(i)*x(i-1);
%     end
% end

end