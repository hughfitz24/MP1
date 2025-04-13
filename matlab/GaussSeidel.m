%% GaussSeidel.m
% M-file creating the function that implements
% Gauss-Siedel for the iterative solution of systems of equations.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MP1.

function [x, iterations] = GaussSeidel(A, b, tolerance, maxIterations, x0)

    % Check for square matrix
    [m, n] = size(A);

    if m ~= n
        error('Matrix A must be square.');
    end

    % Get transpose of A manually

    A_T  = zeros(n, n);

    for i = 1:n
        for j = 1:n
            A_T(i, j) = A(j, i);
        end
    end

    A = A_T*A;

    b = A_T*b;

    % Decompose matrix A

    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);

    % Initial guess

    x = x0;
    iterations = 0;

    % Precompute D inverse

    D_inv = diag(1 ./ diag(D));

    % % Display D_inv, L and U for debugging
    % disp('D_inv matrix:')
    % disp(D_inv)
    % disp('L matrix:')
    % disp(L)
    % disp('U matrix:')
    % disp(U)


    % Start iterations

    for k = 1:maxIterations
        x_old = x;
        % Update x using the Gauss-Seidel formula
        for i = 1:n
            % Compute each term of the formula for element i
            sumL = L(i, :) * x; % uses updated x
            sumU = U(i, :) * x_old; % uses previous x
            x(i) = D_inv(i, i) * (b(i) - sumL - sumU);
        end
        iterations = iterations + 1;
        % Check for convergence
        if norm(x - x_old, inf) < tolerance
            break;
        end
    end
end