%% GaussianElimination.m
% M-file creating the function that implements 
% Gaussian Elimination for the solution of systems of equations.
% Written by Hugh Fitzpatrick, S.N. 22341351 for the completion of MP1. 

function x = GaussianElimination(A, b)

    %% Initial definitions
    Ra = size(A); %% Find size of matrix A
    n = Ra(1); % Number of rows (and columns) of A, assuming A is square
    Id = eye(n); % n x n identity matrix
    Q_tot = Id; % Accumulates the permutation matrix during iteration
 
    % Copy matrices for pivoting, so that originals can be used for final
    % solution
    A1 = A; 
    b1 = b;
 
    for m = 1:n-1 % Loop processes each column in the matrix
        % Index of current row
        row_index = m:n; % Considers rows from m to n 
        % Yc contains max values, Ir contains indices of the values 
        [Yc, Ir] = max(abs(A1(row_index, row_index)));
        [~, Ic1] = max(Yc);
 
        %% Resolving ties
        % Find all occurences of the max value
        tie_indices = find(Yc == max(Yc));
        Ic1 = tie_indices(1);  % Select the first occurence of the max value
 
        % Store row and col number of the pivot point
        row_num = Ir(Ic1); 
        col_num = Ic1;
 
        %% Swap matrix so pivot point is primary element
        % Swap mth row with row containing the pivot 
        p = 1:n;
        p(m) = m-1+row_num;
        p(m-1+row_num) = m;
 
        % Swap mth col with col containing the pivot
        q = 1:n;
        q(m) = m-1+col_num;
        q(m-1+col_num) = m;
 
        % Create permutation matrices corresponding to row and col swaps
        P = Id(p, :);
        Q = Id(:, q);
 
        A1 = P*A1*Q;
        b1 = P*b1;
        Q_tot = Q_tot*Q;
 
        %% Eliminiation
 
        % Build vector that is used to eliminate elements below pivot
        % All values in I_vec are 0 except for the 1st entry in col m
        I_vec = zeros(1,n);
        I_vec(1,m) = 1;
 
        % Eliminate elements under pivot point
        if A1(m,m) ~= 0
         L1 = Id - ([zeros(m,1); A1(m+1:n,m)]*I_vec)/A1(m,m);
        else
            error("ERROR Pivot is equal to 0")
        end 
 
        % Perform elimination 
        A1 = L1*A1;
        b1 = L1*b1;
    end
 
    %% Back substitution 
 
    % Solution vector 
    x = zeros(n,1);
 
    % Iterates in reverse, from last row to first
    for m = 0:n-1
        if A1(n-m,n-m) ~= 0
             x(n-m,1) = (b1(n-m,1) - (A1(n-m,:)*x))/A1(n-m,n-m);
        else 
            error("ERR: Diagonal element is equal to 0. A1 is not an upper triangular matrix.")
        end
        
    end
 
    % Reverse column swaps
     x = Q_tot*x;
 end 
 