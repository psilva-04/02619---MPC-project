function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)

    % Convert two-sided constraints into Aineq*x <= bineq
    Aineq = [A; -A];
    bineq = [bu; -bl];

    % Remove inequalities with an infinite upper limit
    keep = (bineq ~= inf);
    Aineq = Aineq(keep,:);
    bineq = bineq(keep);

    Aeq = [];
    beq = [];

    options = optimoptions('quadprog','Display','off');

    [x,fval,exitflag,output,lambda] = quadprog( ...
        H,g,Aineq,bineq,Aeq,beq,l,u,xinit,options);

    info.fval = fval;
    info.exitflag = exitflag;
    info.output = output;
    info.lambda = lambda;

end