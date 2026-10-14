clear                                         % Clear variables
clc                                           % Clear command window


%% --------- Interpolation Methods for Bessel Function J_0 --------- %%
% Predefined function values (Bessel function J0) at discrete points
% Reference: Abramowitz and Stegun, Table 9.1 (pp. 390-396)
J_0 = [1.00000, 0.99750, 0.99002, 0.97762, 0.96039...
           0.93846, 0.91200, 0.88120, 0.84628, 0.80752...
           0.76519, 0.71962, 0.67113, 0.62008, 0.56685...
           0.51182, 0.45540, 0.39798, 0.33998, 0.28181...
           0.22389, 0.16660, 0.11036, 0.05553, 0.00250...
          -0.04838,-0.09680,-0.14244,-0.18503,-0.22431...
          -0.26005,-0.29206,-0.32018,-0.34429,-0.36429...
          -0.38012,-0.39176,-0.39923,-0.40255,-0.40182...
          -0.39714,-0.38866,-0.37655,-0.36101,-0.34225...
          -0.32054,-0.29613,-0.26933,-0.24042,-0.20973...
          -0.17759];

xData = linspace(0,5,length(J_0));           % Generate evenly spaced x-values corresponding to J_0
xEval = 0:0.01:5;                            % Points where interpolations will be evaluated

% Vandermonde Direct Inversion
yInterp = IntPoly(xData, J_0, xEval);        % Compute interpolation using matrix inversion

% Lagrange Interpolation (Equispaced)
yLagrange = LagrangePoly(xData, J_0, xEval); % Compute interpolation using Lagrange method

% Chebyshev Node Resampling & Lagrange Interpolation
N = 20;                                      % Set number of Chebyshev nodes to 20
i = 0:N;                                     % Generate index array from 0 to N for node generation
xChebyshev = - cos((2*i + 1)/(2*N + 2) * pi); % Chebyshev nodes in [-1,1]
xChebyshev = 5/2 * (xChebyshev + 1);          % Rescale nodes to [0,5]
xRounded = round(xChebyshev,1);               % Round to nearest 0.1 to match available sample grid
indices = round(xRounded * 10) + 1;           % Map rounded nodes to indices in J0 array
yData = J_0(indices);                         % J0 values nearest to Chebyshev nodes
yChebyshev = LagrangePoly(xRounded, yData, xEval); % Compute Lagrange interpolation at Chebyshev nodes

% MATLAB Standard Polyfit/Polyval
N_deg = length(xData) - 1;                    % Max degree of polynomial = number of points - 1
p = polyfit(xData, J_0, N_deg);              % Calculate polynomial coefficients using polyfit
yPoly = polyval(p, xEval);                    % Perform interpolation using polyval on xEval array

% Custom Natural Cubic Spline
ySpline = IntSpline(xData, J_0, xEval);       % Evaluate cubic spline using custom IntSpline implementation

% Built-in MATLAB Spline
yMatlabSpline = spline(xData, J_0, xEval);    % Evaluate cubic spline using MATLAB builtin function (spline)

%% Plot Comparisons for Polynomial Interpolation Methods
figure
hold on
plot(xData, J_0,'ko','MarkerSize',5)               % Plot original data points from textbook
plot(xEval, yInterp,'y--','LineWidth',1.5)         % Plot Vandermonde interpolation (matrix inversion)
plot(xEval, yLagrange,'c.-','LineWidth',1.5)       % Plot equispaced Lagrange interpolation
plot(xEval, yChebyshev,'r:','LineWidth',1.8)       % Plot Chebyshev Lagrange interpolation
plot(xEval, yPoly,'b-','LineWidth',1.5)            % Plot Polyfit interpolation
axis([0 5 -0.5 1])                                 % Set axis limits
legend("Data", "Inversion Interpolation", ...
    "Lagrange Interpolation", "Chebyshev " + ...
    "Interpolation","Polyfit Interpolation")       % Populate legend
xlabel("x")                                        % Formatting and naming
ylabel("J_0(x)")
title("Bessel Function obtained from Polynomial Interpolation Methods")
hold off

%% Plot Comparisons for Spline Interpolation Methods
figure
hold on
plot(xData, J_0,'ko','MarkerSize',5)               % Plot original data points from textbook
plot(xEval, ySpline,'r--','LineWidth',1.5)         % Plot custom spline interpolation
plot(xEval, yMatlabSpline,'b.-','LineWidth',1.5)   % Plot builtin spline interpolation
legend("Data", "Cubic Spline Interpolation", ...
    "Built-in Spline Interpolation")               % Populate legend
xlabel("x")                                        % Formatting and naming
ylabel("J_0(x)")
title("Bessel Function obtained from Cubic Spline Interpolation Methods")
hold off

%% Plot Absolute Differences Between Polynomial Interpolation Methods (Log Scale)
figure
semilogy(xEval, abs(yLagrange - yInterp),'y--','LineWidth',1.5)  % Difference betw. IntPoly vs Lagrange
hold on
semilogy(xEval, abs(yLagrange - yPoly),'c.-','LineWidth',1.5)    % Difference betw. yPoly vs Lagrange
semilogy(xEval, abs(yPoly - yInterp),'r:','LineWidth',1.5)        % Difference betw. IntPoly vs yPoly
legend("Lagrange - IntPoly", "Polyfit - Lagrange", "Polyfit - IntPoly")         % Populate legend
xlabel("x")                                                       % Formatting and naming
ylabel("Absolute Difference (log scale)")
title("Absolute Differences Between Polynomial Interpolation Methods")
hold off

%% Plot Absolute Differences Between Spline Interpolation Methods
figure
plot(xEval, abs(ySpline - yMatlabSpline),'y--','LineWidth',1.5)
legend("|Custom Spline - Built-in Spline|")                            % Populate Legend
xlabel("x")                                                            % Formatting and naming
ylabel("Absolute Difference (linear scale)")
title("Absolute Differences Between Spline Interpolation Methods")

figure
semilogy(xEval, abs(ySpline - yMatlabSpline),'y--','LineWidth',1.5)   % Difference betw. splines
legend("|Custom Spline - Built-in Spline|")                            % Populate Legend
xlabel("x")                                                            % Formatting and naming
ylabel("Absolute Difference (log scale)")
title("Logarithmic Absolute Differences Between Spline Interpolation Methods")

%% --------- Numerical Quadrature Benchmarks for Bessel Function Integration --------- %%
n = 0;                              % Evaluating J0 Bessel function (order n = 0)
N_vec = [5, 11, 17, 19, 35, 37, 51]; % Array of integration node counts
xEval = 0:0.1:6;                    % Evaluation grid for x in [0, 6]
                
% Reference Numerical Integration via MATLAB integral()
yIntegral = zeros(1, length(xEval));                        % Preallocate array for reference J0 values
for j = 1:length(xEval)                                     % Loop through each x in xEval array 
    integrand = @(tau) cos(n * tau - xEval(j) * sin(tau));  % Integrand for J0
    yIntegral(j) = (1/pi) * integral(integrand, 0, pi);     % Compute J0 via inbuilt integral function
end

for i = 1:length(N_vec)
    % Trapezoidal Method
    yTrap = trapzJ(n, N_vec(i), xEval);               % Evaluate integral using trapezoid method
    
    % Simpson's 1/3 Rule
    ySimp = simpsonJ(n, N_vec(i)-1, xEval);           % Evaluate integral using Simpson's 1/3 rule
    
    % Gauss-Legendre Quadrature
    yGauss = gaussJ(n, N_vec(i), xEval);              % Evaluate integral using Gaussian quadrature
    
    % Plot results for current N
    figure('Position',[200 200 1000 400])
    tiledlayout(1,2)                                 % Tiled layout to organize graphs better
    sgtitle(['Node count N = ', num2str(N_vec(i))])  % Add title to each tile group w/ node count
    
    % Plot J0 approximations for each method
    nexttile                                         % Left plot
    hold on
    plot(xEval,yTrap,'--r',xEval,ySimp,'--b',xEval,yGauss,'--g',xEval,yIntegral,'--y')
    xlabel('x')                                      % Formatting and naming
    ylabel('J_0(x)')
    legend("Trapezium", "Simpson", "Gauss", "MATLAB Reference")
    grid on
    hold off
    
    % Calculate absolute errors between method and inbuilt integral function
    nexttile                                         % Right plot
    errTrap = abs(yTrap - yIntegral);
    errSimp = abs(ySimp - yIntegral);
    errGauss = abs(yGauss - yIntegral);
    
    % Plot absolute errors
    semilogy(xEval, errTrap, '--r')                  % Plot trapezoid method approximation on log scale
    hold on
    semilogy(xEval, errSimp, '--b')                  % Plot Simpson's rule approximation on log scale
    semilogy(xEval, errGauss, '--g')                 % Plot Gaussian Quadrature approximation on log scale
    xlabel('x')                                      % Formatting and naming
    ylabel('Absolute Error')
    legend("Trapezium", "Simpson", "Gauss")
    hold off
end

% Numerical Quadrature Convergence Analysis:
% The Trapezoidal Method settles to machine precision (~10^-16) by N = 17 due to exponential 
% convergence for smooth periodic integrands over full periods. Simpson's Rule converges 
% around N = 35. High-order Gaussian Quadrature suffers from ill-conditioning in polynomial 
% root-finding for high N, leading to floating-point instability.
% Note: Zero-error values produce log(0) undefined behavior on logarithmic plots; bounding with
% eps (e.g., max(err, eps)) preserves plot continuity.


%% --------- Root-Finding Algorithms & Application to Membrane Vibrations --------- %%

format long                                 % Display numbers with maximum precision
x = 0:0.1:6;                                % Create x array from 0 to 6
n = 0;                                      % Bessel order = 0
J0 = trapzJ(n, 51, x);                      % Compute approximate J0 using trapzJ function

% Bisection Method
tic                                         % Start timer for bisection method
rootsB = bisectionJ(x, J0, n);              % Compute roots using bisection
tBisection = toc;                           % Record execution time

% Newton-Raphson Method
tic                                         % Start timer for Newton method
rootsN = newtonJ(x, J0, n);                 % Compute roots using Newton
tNewton = toc;                              % Record execution time

% Secant Method
tic                                         % Start timer for secant method
rootsS = secantJ(x, J0, n);                 % Compute roots using secant method
tSecant = toc;                              % Record execution time

% MATLAB Built-in fzero Method
tic                                         % Start timer for fzero method
rootsZ = fzeroJ(x, J0, n);                  % Compute roots using fzero
tZero = toc;                                % Record execution time

% Reference Roots (Abramowitz & Stegun Table 9.5)
rootsActual = [2.4048255577; 5.5200781103]; 

% Calculate Mean Relative Percentage Error
errBisection = 100 * mean(abs(rootsActual - rootsB)./rootsActual); % Compute % error for bisection
errNewton    = 100 * mean(abs(rootsActual - rootsN)./rootsActual); % Compute % error for Newton
errSecant    = 100 * mean(abs(rootsActual - rootsS)./rootsActual); % Compute % error for secant
errZero      = 100 * mean(abs(rootsActual - rootsZ)./rootsActual); % Compute % error for fzero

%% Root-Finding Performance Visualization
methods = {'Bisection','Newton','Secant','fzero', 'Actual'};         % List of method names for plotting
rootsAll = {rootsB, rootsN, rootsS, rootsZ};                         % Store roots
timesAll = [tBisection, tNewton, tSecant, tZero];                    % Store execution times
errsAll  = [errBisection, errNewton, errSecant, errZero];            % Store errors

for k = 1:4                                    % Loop over each method for plotting
    figure;                                    % Create new figure for each method
    plot(x, J0, 'b', 'LineWidth', 1.5);        
    hold on;
    yline(0,'k--');                            % Horizontal line at y=0 for reference
    plot(rootsAll{k}, zeros(size(rootsAll{k})), 'rx', 'MarkerSize', 8, 'LineWidth',2);     % Plot roots
    for i = 1:length(rootsAll{k})
        r = rootsAll{k}(i);
        text(r, 0, sprintf('%.8g', r), ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','center', ...
            'FontSize',8);
    end
    % Add label showing execution time and mean error
    text(0.5, 0.9*max(J0), sprintf('Time: %.4f s\nMean Error: %.2e %%', timesAll(k), errsAll(k)), 'Units','normalized');
    xlabel('x'); ylabel('J_0(x)'); 
    title(sprintf('%s Method: Roots vs J_0(x)', methods{k}));
    legend('J_0(x)', 'y = 0', 'Roots', 'Location','best');
    hold off;
end

%% Tabulated Performance Results
values = zeros(length(rootsActual), 5);        % Initialise values array
for i = 1:length(rootsActual)                  % Loop through each root
    values(i, :) = [rootsB(i), rootsN(i), rootsS(i), rootsZ(i), rootsActual(i)];
end                                            % Put roots obtained from each method into same row

% Print Roots Table (10 decimal places)
fprintf('\nCalculated Roots Across Methods:\n');    
fprintf('%-12s %-12s %-12s %-12s %-12s\n', 'Bisection', ...
    'Newton', 'Secant', 'fzero', 'Actual');          
for i = 1:size(values, 1)                            
    fprintf('%-12.10f %-12.10f %-12.10f %-12.10f %-12.10f\n', ...
            values(i,1), values(i,2), values(i,3), values(i,4), values(i,5)); 
end

% Print Execution Time & Relative Error Table
fprintf('\n%-12s | %-18s | %-18s\n', 'Method', ...
    'Execution Time (s)', 'Relative Error (%)');         
for i = 1:length(timesAll)                           
    fprintf('%-12s | %-18.2e | %-18.4e\n', ...
        methods{i}, timesAll(i), errsAll(i));        
end

%% Physical Application: Circular Membrane Mode Frequencies
% Dynamic physical parameters for a clamped circular membrane
h = 0.0001;                 % Membrane thickness (m)
rho = 1000;                 % Density (kg/m^3)
D = 0.1;                    % Diameter (m)
N_tension = 12;             % Membrane tension resultant (N/m)
a = 0.05;                   % Radius (m), derived from diameter
c = sqrt(N_tension/(rho*h)); % Wave propagation speed

% Evaluation points for Bessel function approximation
xEval = 0.1:0.1:15;

% Compute Bessel functions of order 0, 1, and 2 using 17-point Trapezoidal Quadrature
J0 = trapzJ(0, 17, xEval); % Evaluate J0(xEval)
J1 = trapzJ(1, 17, xEval); % Evaluate J1(xEval)
J2 = trapzJ(2, 17, xEval); % Evaluate J2(xEval)

% Find fundamental roots using Newton-Raphson Method
r0 = newtonJ(xEval, J0, 0); % Find J0 roots
r1 = newtonJ(xEval, J1, 1); % Find J1 roots
r2 = newtonJ(xEval, J2, 2); % Find J2 roots

% Combine roots corresponding to modes (0,1), (1,1), (2,1), (0,2)
roots = [r0(1); r1(1); r2(1); r0(2)]; 
roots = sort(roots);                   % Sort roots in ascending order to align with increasing frequency

% Calculate Modal Frequencies (Hz): f_mn = (alpha_mn * c) / (2 * pi * a)
freq = roots * c / (2*pi*a);

% Reference Benchmark Frequencies
real_freq = [82.2, 158, 217, 227];

% Display Computed Frequencies
fprintf('\nCalculated Circular Membrane Frequencies (Hz):\n');
for i = 1:length(freq)
    fprintf('Mode %d: %.1f Hz (Reference: %.1f Hz)\n', i, freq(i), real_freq(i));
end

%% Functions
function y = IntPoly(xData, yData, xEval)       % Initialize function w/ inputs x, y and xEval and return y array

    n = length(yData);          % Number of data points

    VDB = ones(n);              % Initialize the Vandermonde matrix with ones instead of zeros so multiplications does smth (include this point)
    for row = 1:n               % Loop through data points
        for col = 1:n           % Loop through polynomial powers
            VDB(row,col) = xData(row)^(col - 1);    % Fill matrix: x^0 to x^(n-1) from class 4 slides
        end
    end

    coeffs = VDB \ yData';      % Solve linear system V*a = y for poly coefficients using matrix inversion

    a = flipud(coeffs);         % Since ascending flip so that descending (for polyval function)
    y = polyval(a, xEval);      % Evaluate p(x) at xEval points

end


function y = LagrangePoly(xData, yData, xEval)  % Initialize Lagrange function w/ inputs x, y and xEval and return y array
    
    n = length(xData);                  % Number of data points

    y = zeros(size(xEval));             % Initialise y array
    for i = 1:n                         % Loop over each basis polynomial 

        product = ones(size(xEval));    % Initialise Li(x) (ones used to not multiply by 0)
                                     
        for j = 1:n                     % Loop over all data points
            if i ~= j                   % When i = j, dy = 0 so nothing added
                product = product .* (xEval - xData(j)) / (xData(i) - xData(j));    % Populate Li(x) w/ formula
            end
        end

        y = y + yData(i) * product;     % Add ith term to final y

    end
end


function y = IntSpline(xData, yData, xEval)     % Initialize function w/ inputs x, y and xEval and return y array

    n = length(xData);              % Number of data points
    h = xData(2) - xData(1);        % Step size (assumes uniform spacing)
    m = n - 2;                      % Number of interior points

    A = zeros(m,m);                 % System matrix for second derivatives
    RHS = zeros(m,1);               % Right-hand side vector

    for i = 1:m                     % Loop thru interior points
        A(i,i) = 2/3;               % Assign center weight for the second-derivative continuity
        
        if i < m                    % If not at the rightmost value
            A(i, i+1) = 1/6;        % Link node to the curvature of the neighbor on the right
        end
        
        if i > 1                    % If not at the leftmost value
            A(i, i-1) = 1/6;        % Link node to the curvature of the neighbor on the left
        end
  
        RHS(i) = (yData(i+2) - 2*yData(i+1) + yData(i))/h^2; % Calculate RHS using formula from class 4 slides
    end

    y_dd_raw = A \ RHS;             % Solve for 2nd derivatives using inversion

    y_dd = zeros(n,1);              % Initialise second derivative array
    y_dd(2:n-1) = y_dd_raw;         % Fill w/ 2nd der values and apply free runout boundary conditions (ends = 0)

    y = zeros(size(xEval));         % Initialise second derivative array

    for k = 1:length(xEval)         % Loop thru xEval array

        x = xEval(k);               % Current evaluation point

        for i = 1:n - 1                          % Loop thru intervals
            if x >= xData(i) && x <= xData(i+1)  % Check if interval contains x
                
                a = xData(i);       % Left endpoint of interval
                b = xData(i+1);     % Right endpoint of interval

                % Apply cubic spline formula
                y(k) =   y_dd(i)/6   * ((b - x)^3/(b - a) - (b - a) * (b - x)) + ...
                       + y_dd(i+1)/6 * ((x - a)^3/(b - a) - (b - a) * (x - a)) + ...
                       + yData(i)   * (b - x)/(b - a) + ...
                         yData(i+1) * (x - a)/(b - a);

                break
            end
        end
    end

end


function I = trapzJ(n, N, x)             % Initialize function w/ inputs n, N, x, and return I array.
    tau = linspace(0, pi, N+1);          % Generate N+1 equispaced points in range 0 to pi
    
    F = cos(n * tau - x.' * sin(tau));   % Evaluate integrand for all combinations of x and tau, 
                                         % transpose x to allow for matrix multiplication (x --> col vector, sin(tau) --> row vector)
    
    I = trapz(tau, F, 2);                % Integrate over each row of tau using trapezium rule function
                                         % to get integral array w/ a J0 estimate for each x in x vector

    I = I / pi;                          % Divide by pi (formula)

    I = I';                              % Turn into column vector (like the input x)
end


function I = simpsonJ(n, N, x)              % Initialize function w/ inputs n, N, x, and return I array.
    
    h = pi / N;                             % Step size
    I = zeros(size(x));                     % Initialize integral variable
    f = @(tau) cos(n * tau - x * sin(tau)); % Create function handle for integrand

    for i = 1:N-1                           % Loop through all indices except 0 and N (dealt w/ later)
        tau = i * h;                        % Assign current tau value
        if mod(i, 2) == 0                   % If even index
            I = I + 2 * f(tau);             % Add f(tau) w/ weighting of 2
        else                                % If odd index
            I = I + 4 * f(tau);             % Add f(tau) w/ weighting of 4
        end
    end

    I = h/3 * (f(0) + I + f(pi));           % Add endpoints and multiply by h/3 (formula)

    I = I / pi;                             % Divide by pi (formula)
end


function L = LegendrePoly(N)                                    % Function returns L polynomial array and uses N, the order of Legendre Polynomial, as input
    P0 = 1;                                                     % Set P0 = 1
    P1 = [1 0];                                                 % Set P1 = x

    if N == 0                                                   % If n = 0
        L = P0;                                                 % Set L to P0=1, as L = 1
        return                                                  % Exit the function
    elseif N == 1                                               % If n = 1
        L = P1;                                                 % Set L to P1=[1 0], as L = x
        return                                                  % Exit the function
    end

    for i = 1:N-1                                               % Loop from 1 to N-1 (as per formula)

        xPn = conv(P1, [1 0]);                                  % Polynomial multiplication by x
        
        P0_padded = [zeros(1, length(xPn) - length(P0)) P0];    % Ensure P0 has same length by adding zeros to the left until matching size w/ xPn

        P2 = ((2*i + 1) * xPn - P0_padded * i) / (i + 1);       % Carry out recurrence relation

        P0 = P1;                                                % Shift forward in the recurrence by reassigning P0
        P1 = P2;                                                % Shift forward in the recurrence by reassigning P1
    end

    L = P1;                                                     % Assign the current polynomial P1 to L array
end


function I = gaussJ(n, N, x)                                    % Initialize function w/ inputs n, N, x, and return I array.
    L = LegendrePoly(N);                                        % Extract Legendre polynomial coefficients from LegendrePoly() function
    quad_points = sort(roots(L));                               % Sort roots in ascending order in case unsorted 

    w = zeros(1,N);                                             % Initialize weight array for memory efficiency
    der = polyder(L);                                           % Derivative array uses polyder on Lagrange polynomial to obtain dL_N/dx
                                                                % to calculate weights.
    
    for i = 1:N                                                 % Compute weights for all values of i from 1 to N
        der_eval = polyval(der, quad_points(i));                % Calculate derivative at x = eta_i as per formula
        w(i) = 2 / ((1 - quad_points(i)^2) * (der_eval)^(2));   % Calculate weight(i) as per formula on slide 3.46
    end

    f = @(tau) cos(n * tau - x * sin(tau));                     % Create function handle for integrand
    tau_map = @(u) (pi/2) * (u + 1);                            % Transform tau values such that bounds 0 to pi become -1 to 1

    I = zeros(size(x));                                         % Initialize Integral variable
    for i = 1:N                                                 % Add to integral from i = 1 to N                  
        I = I + w(i) * f(tau_map(quad_points(i)));              % Iteratively add to the integral by adding weight(i) * f(eta_i)
    end
    I = 1/2 * I;                                                % Normalise by dividing by 2 (1/pi (formula) * pi/2 (jacobian from transformation) = 1/2)
end


function roots = fzeroJ(x, J, n) % Find roots using MATLAB's fzero for each interval with a sign change

    indices = find(J(1:end-1) .* J(2:end) < 0);  % Find indices where J changes sign
    tStart = tic;                                % Start timer to enforce time cap
    roots = zeros(length(indices), 1);           % Initialise array for root results

    for i = 1:length(indices)                    % Loop over each interval with sign change
        if toc(tStart) > 600, warning('Time cap reached'); break; end % Stop if time > 10min

        idx = indices(i);                        % Current index of sign change
        xInterval = [x(idx), x(idx+1)];          % Extract x values for the current interval

        f = @(x) trapzJ(n, 51, x);               % Function handle for integral evaluation
        roots(i) = fzero(f, xInterval);          % Find root using fzero w/ input interval and function handle
    end
end


function roots = bisectionJ(x, J, n)              % Find roots using the bisection method

    tol = 1e-10;                                  % Convergence tolerance for iterations
    tStart = tic;                                 % Start timer
    indices = find(J(1:end-1) .* J(2:end) < 0);   % Find sign-change interval indices
    roots = zeros(length(indices), 1);            % Preallocate roots array

    for i = 1:length(indices)                     % Loop over each sign-change interval
        if toc(tStart) > 600, warning('Time cap reached'); break; end % Time limit check

        a = x(indices(i));                        % Left endpoint of interval
        b = x(indices(i)+1);                      % Right endpoint of interval
        f = @(x) trapzJ(n, 51, x);                % Function handle for integral

        err = Inf;                                % Initialize error to infinity (to ensure no early exit)
        c = (a+b)/2;                              % Initial midpoint of the interval

        while err > tol                           % Iterate until convergence
            c0 = c;                               % Store previous midpoint
            if f(c) * f(a) > 0                    % Check which half contains the root
                a = c;                            % Root in [c,b], update lower bound to c
            else
                b = c;                            % Root in [a,c], update upper bound to c
            end
            c = (a+b)/2;                          % Compute new midpoint
            err = abs(c - c0);
        end
        roots(i) = c;                             % Store computed root in roots array
    end
end


function roots = newtonJ(x, J, n)               % Find roots using Newton-Raphson method

    tol = 1e-10;                                % Convergence tolerance
    tStart = tic;                               % Start timer for time cap

    indices = find(J(1:end-1) .* J(2:end) < 0); % Find intervals where J changes sign
    roots = zeros(length(indices), 1);          % Initialise array for root results

    for i = 1:length(indices)                   % Loop thru each sign-change interval
        if toc(tStart) > 600, warning('Time cap reached'); break; end % Stop if time >10 min

        x0 = (x(indices(i)) + x(indices(i)+1))/2; % Initial guess: midpoint of interval
        err = Inf;                              % Initialize error to infinity (to ensure no early exit)

        while err > tol                         % Continue until convergence (when err < tol)
            if toc(tStart) > 600, warning('Time cap reached'); break; end % Check timer

            f  = trapzJ(n, 51, x0);             % Evaluate func at guess
            df = 0.5*(trapzJ(n-1,51,x0) - trapzJ(n+1,51,x0)); % Approximate derivative using formula

            if abs(df) < 1e-12, df = eps; end   % Deal w/ dividing by 0 if derivative is very small

            x1 = x0 - f/df;                     % Newton Raphson formula
            err = abs(x1 - x0);                 % Calculate difference to check convergence
            x0 = x1;                            % Update guess for next iteration
        end
        roots(i) = x0;                          % Store converged root for this interval
    end
end


function roots = secantJ(x, J, n)               % Find roots using the secant method

    tol = 1e-10;                                % Convergence tolerance
    tStart = tic;                               % Start timerfor time cap

    indices = find(J(1:end-1) .* J(2:end) < 0); % Find intervals where J changes sign
    roots = zeros(length(indices), 1);          % Initialise array to store roots

    f = @(a) trapzJ(n, 51, a);                  % Function handle for evaluation using integration

    for i = 1:length(indices)                   % Loop over each sign-change interval
        if toc(tStart) > 600, warning('Time cap reached'); break; end % Stop if t>10min

        x0 = x(indices(i));                     % First initial guess (left endpoint)
        x1 = x(indices(i)+1);                   % Second initial guess (right endpoint)

        err = Inf;                              % Initialize error to infinity (to ensure no early exit)
        while err > tol                         % Continue until convergence (when err < tol)
            if toc(tStart) > 600, warning('Time cap reached'); break; end % Check timer

            % Secant method formula: compute new estimate based on previous two
            x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0)); 

            err = abs(x2 - x1);                 % Calculate error (diff in consecutive roots)
            x0 = x1;                            % Update previous guess
            x1 = x2;                            % Update current guess
        end
        roots(i) = x2;                          % Store computed root in roots array
    end
end