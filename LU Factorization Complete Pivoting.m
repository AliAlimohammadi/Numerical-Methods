%  AUTHOR:      Ali Alimohammadi
%    DATE:      6 Apr 2020
%  GitHub:      @AliAlimohammadi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare A as an empty matrix.
A = [];


function [L, U, P, Q] = lu_cp(A)
% This function performs LU factorization for a matrix A. The function
% returns the lower and upper triangular matrices as separate matrices
% to make checking easier. It also returns matrices P and Q indicating
% the row and column exchanges, respectively. The function uses complete
% pivoting, such that PAQ = LU. The matrix must be nonsingular and square.
    
    % We will need to know the dimension of the matrix A.
    dim = length(A);

    % Initialize P and Q to the identity matrices.
    P = eye(dim);
    Q = eye(dim);

    % For each row in A,
    for i = 1 : (dim - 1),

        % Find the element with largest magnitude in each
        % submatrix which will be the new pivot.
        pivot = max(max(abs(A([i : dim], [i : dim]))));

        % Find the indeces of the new pivot
        [X, Y] = find(abs(A([i : dim], [i : dim])) == pivot);
        if i ~= 1,
            X(1) = X(1) + (i - 1);
            Y(1) = Y(1) + (i - 1);
        end;

        % Interchange the rows and columns of the new pivot
        % with the old one
        A([i, X(1)], :) = A([X(1), i], :);
        A(:, [i, Y(1)]) = A(:, [Y(1), i]);

        % Store the changes in the matrices P and Q
        P([i, X(1)], :) = P([X(1), i], :);
        Q(:, [i, Y(1)]) = Q(:, [Y(1), i]);

        % Compute the factor.
        A([(i + 1) : dim], i) = A([(i + 1) : dim], i) / A(i, i);

        % Multiply the "nonzero" elements of row i by the
        % factor. Subtract this result from the "nonzero"
        % elements of row j.
        A([(i + 1) : dim], [(i + 1) : dim]) = A([(i + 1) : dim], [(i + 1) : dim]) -  A([(i + 1) : dim], i) * A(i, [(i + 1) : dim]);
        
    end;
    
    % For each row under row i,
    for j = (i + 1) : dim,
    
        % Check to see if we will encounter division by zero.
        if abs(A(i, i)) <= 1e-12,
            disp('The matrix is singular.');
            U = NaN;
            L = NaN;
            return;
        end;
        
    end;

    % The U factor is the upper triangle of A, with zeros
    % in the lower triangle.
    U = triu(A);

    % The L factor is the lower triangle of A, with zeros
    % in the lower triangle and ones along the main diagonal.
    L = tril(A, -1) + eye(dim);

end;


function [V] = vandermonde(C)
% This function initializes the Vandermonde matrix for
% a given set of constants c<1>, c<2>, c<3>, ..., c<n>.
% The function returns the corresponding Vandermonde matrix.

    % Build the Vandermonde matrix using MATLAB's in-build functions.
    % vander(C) returns the Vandermonde matrix whose columns are
    % powers of the vector C. Then the function finds the alternate
    % form of the Vandermonde matrix using fliplr().
    V = fliplr(vander(C));
    
end;


function [W] = willkinson(n)
% This function creates an n*n Willkinson matrix for the
% corresponding value of n. Then returns the result.

    % A lower triangular matrix filled with -1
    % and ones along the main diagonal.
    W = tril(-ones(n), -1) + eye(n);
    
    % Last column filled with ones.
    W(:, n) = 1; 

end;


function [x] = solve_linear_system(A, b)
% This function performs the required operations to solve
% the system of linear algebraic equations Ax = b, for the
% corresponding square matrix A (coefficient matrix)
% and b as a column vector with n entries.

    % LU decomposition for corresponding coefficient matrix A.
    [L, U, P, Q] = lu_cp(A);
    
    % We want to solve the equation Ax = b for x, given A and b.
    % We have already obtained the LU decomposition of matrix A 
    % such that PAQ = LU. Now we shall find the unknowns.
    
    % In this case the solution is done in two logical steps:
    % First, we solve the equation Ly = Pb for y.
    y = L \ (P * b);
    
    % Second, we solve the equation Ux = y for x.
    x = Q * (U \ y);
    
end;


% Driver Code
fprintf('Enter "n" for your n*n matrix:\n');
n = input('')

fprintf('\nEnter your matrices using standard MATLAB notation.\n');

fprintf('For example:\n\n\t[1, 2, 3; 4, 5, 6; 7, 8, 9]\n\nis interpreted as:\n\n');
disp([1, 2, 3; 4, 5, 6; 7, 8, 9]);

fprintf('\nEnter your constants c<1>, c<2>, c<3>, ..., c<n> for Vandermonde matrix: (use "comma" to seperate entries)\n');
C = input('')

fprintf('\nEnter b (at Ax = b) as a single column vector: (use "semicolon" to seperate entries)\n');
b = input('')

fprintf('\nSolution for Vandermonde matrix:\n');
x = solve_linear_system(vandermonde(C), b)

fprintf('\nSolution for Willkinson matrix:\n');
x = solve_linear_system(willkinson(n), b)