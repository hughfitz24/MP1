%% NewtonRaphson.m
% M-file creating the function that implements the 
% Newton-Raphson algorithm for the solution of systems of equations.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MP1. 
% NOTE: Original authorship was for the completion of MS1, this M-file has
% been re-used. 

function [x, iterations] = NewtonRaphson(f, J, x0, tol, maxIter)
    %% Inputs:
    %   f       - Function handle for system of equations
    %   J       - Function handle for Jacobian matrix
    %   x0      - Initial guess
    %   tol     - Tolerance (1e-9)
    %   maxIter - Maximum iterations (20)
    % Outputs:
    %   x       - Solution vector

    %% Input validation
    if tol < 1e-11 
        error("ERR: Tolerance is excessively small. Consider range of 1e-9.")
    elseif tol > 1e-6
        error("ERR: Tolerance is too large. Consider range of 1e-9.")
    end

    if maxIter > 20
        error("ERR: Maximum iterations should not exceed 20.")
    elseif maxIter < 0
        error("ERR: Maximum iterations should be a positive integer.")
    end

    % Initial guess
    x = x0;

    %% Newton-Raphson iteration
    % Iterative over range of iterations provided by user
    for iter = 1:maxIter
        % Evaluate function and Jacobian
        F = f(x);
        Jacobian = J(x);

        % Solve Jacobian * delta_x = -F using Gaussian elimination, as
        % detailed in Equation 7.
        deltaX = GaussianElimination(Jacobian, -F);

        % Update solution
        x = x + deltaX;

        %% Check for convergence
        % Check if deltaX is less than the tolerance 
        if norm(deltaX) < tol
            fprintf('Newton Raphson algorithm converged after %d iterations.\n', iter);
            iterations = iter;
            return;
        end
    end
    %% Error Messaging
    % If algorithm did not successfully converge, display error msg
    error('Newton-Raphson did not converge after %d iterations.', maxIter);
end

