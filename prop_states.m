function [y]=prop_states(x0,t)
    % This function integrates numerically using ode45 the desired process for
    % the filter.
    % Inputs:
    % x0: initial state before propagation;
    % dt: time for propagation [s];
    % Outputs:
    % y : new state after propagation;

    % setting tolerance for numerical integration:
        options=odeset('AbsTol',1e-6,'RelTol',1e-12);
    % ode integration:
        [~,Y]=ode113(@(t,y) process(t,y),[0,t],x0,options);

    % extraction of final state:
        y=Y(end,:)';
end