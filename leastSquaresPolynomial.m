
function [polyValues, polyCoeffs, Rsqu] = leastSquaresPolynomial(xValues,yValues,order)

    % Determine coefficients of k-th order polynomial using least squares
    % algorithm
    
    % length(x) = length(y) = n
    n = length(xValues);
    
    % Construct Y matrix
    Y = yValues';
    
    % Construct X matrix
    X = constructPowerMatrix(xValues,n,order);
    
    % Peform matrix multiplication
    % a = (X^T * X)^-1 * X^T * Y
    polyCoeffs = (X'*X)\X'*Y;
    
    % Calculate new curve
    polyValues = X*polyCoeffs;
    
    % Compute R squared
    yMean = sum(yValues)/length(yValues);
    SSres = sum((yValues - polyValues).^2);
    SStot = sum((yValues - yMean).^2);
    Rsqu = 1 - SSres/SStot;
    
end

% 
function X = constructPowerMatrix(xValues,n,k)
    
    % Construct X matrix
    X = zeros(n,k+1);
    
    % Add zero-th powers of x elements
    X(:,1) = ones(n,1);
    
    % Fill up X matrix
    % Fill columns and then rows
    for i = 1:n
        for j = 1:k
            X(i,j+1) = xValues(i)^j;
        end
    end
end