function [result] = Newton_angle(init_guess,Nt, e, tol)
%NEWTON'S METHOD to sole next value of E with Newton's method.
%       

    E1 = init_guess;
    E2 = E1 - (E1 - e*sin(E1)- Nt)/(1 - e*cos(E1));
    
    i = 0;
    while(abs(E1-E2) > tol)
        E1 = E2;
        E2 = E1 - (E1 - e*sin(E1)- Nt)/(1 - e*cos(E1));
        i = i+1;
        if(i == 1e5)
            error("Newton's method did not manage to converge in 100 000 steps.\n")
        end
    end
    result = E2;
end

