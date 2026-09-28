function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit) 
 
    % Convert two-sided linear constraints 
    % bl <= A*x <= bu 
    % into quadprog form Aineq*x <= bineq 


    Aineq = [ A; 
             -A ]; 
 
    bineq = [ bu; 
             -bl ]; 
 
    % no equality constraints 
    Aeq = []; 
    beq = []; 
 
    options = optimoptions('quadprog','Display','off'); 
 
    [x,fval,exitflag,output,lambda] = ... 
        quadprog(H,g,Aineq,bineq,Aeq,beq,l,u,xinit,options); 
 
    % solver information 
    info.fval     = fval; 
    info.exitflag = exitflag; 
    info.output   = output; 
    info.lambda   = lambda; 
 
end