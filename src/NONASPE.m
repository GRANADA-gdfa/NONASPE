%==========================================================================
%                                                                         
% NONASPE                                                                 
% NON-newtonian Accurate Solver for Planetary Estuarines                        
%                                                                         
% Copyright © 2023 Llorente Lázaro, Víctor Javier                         
% Last Update: May 24, 2023                                               
% Website: https://sites.google.com/view/vjllorente                       
% Contact: victor.javier.llorente@gmail.com                               
%                                                                         
%--------------------------------------------------------------------------
%
% ------------
% Description:
% ------------
% This is an in-house CFD code which solves the boundary-value problem:
%
%   -                                                          o z = 0
%  | d  /      / du \ ^ n \                                    |
%  | -- | K0 * | -- |     | = -g * etax , z in ( 0, -H )       |
%  | dz \      \ dz /     /                                    |
%  |                                                           |
% -|            / du \ ^ n                                     |
%  | rho * K0 * | -- |     = taux       , at z = 0             |
%  |            \ dz /                                         |
%  |                                                           o z = -H
%  | u = 0                              , at z = -H
%   -
%
% where K0 is an eddy viscosity (assume constant), rho is the density of
% the estuarine, etax is a barotropic pressure gradient (assume constant), 
% and tauw is the wind stress. g is the gravity constant.
%
% An effective eddy viscosity can be formulated, 
%
%                | du | ^ n - 1
%    Keff = K0 * | -- |         
%                | dz |
%
% Therefore, the differential equation is re-written as follows:
%
%  d  /        du \                                   
%  -- | Keff * -- | = -g * etax    
%  dz \        dz /    
%
% -------
% INPUTS:
% -------
% Double                :: H            - Depth [m]
% Double                :: K0           - Eddy viscosity [m^2/s]
% Double                :: rho          - Water density [kg/m^3] 
% Double                :: n            - Power exponent [-]
% Double                :: etax         - Barotropic pressure gradient [-]
% Double                :: tauw         - Wind stress [kg/ms^(n+1)]
% Double                :: g            - Gravity [m/s^2]
% String                :: induced      - Solve wind-driven for n = 1 [wind-driven / no wind-driven]
% String                :: viscosity    - Choose a viscosity profile [visco 0 / visco 1]
% Double                :: zm           - Position of the maximum for visco 1
% Double                :: zh           - Position of leyer separation for visco 1
% Double                :: m            - Decay exponent for visco 1
% String                :: solver       - Choose the solver [TDMA / no TDMA]
% String                :: trick        - Bound the values of the matrix [bounded matrix / no bounded matrix]
% Double                :: Amax         - Upper bound value
% Double                :: Amin         - Lower bound value
% Double                :: N            - Number of nodes [-]
% Double                :: IterMax      - Maximum of iterations [-]
% Double                :: Tol          - Tolerance [-]
% 
% --------
% OUTPUTS:
% --------
% Double, Nx1 array     :: U            - Velocity [m/s]
% Double                :: Q            - Flow rate [m^3/s]
% Double, Nx1 array     :: dUdz         - Velocity gradient [1/s]
% Double, Nx1 array     :: d2Udz2       - Velocity laplacian [1/ms]
% Double, Nx1 array     :: vK           - Physical eddy viscosity vector [m^2/s]
% Double, Nx1 array     :: vKeff        - Effective eddy viscosity vector [m^2/s^n]
% Double, iterx1 array  :: RES          - Residuals vector [-]
%
% ---------
% COMMENTS:
% ---------
% The solver works apparently for a range of 0.3 < n < 2 (H = 10) with the 
% solver of Matlab and 0.3 <= n <= 2 using TDMA. The solution presents a 
% inflection point at some point in the domain. 
% 1) For n <= 0.3, the numerical Keff at the discontinuity is so large 
%    that make the matrix singular (bad condition). 
% 2) For n >= 2, the numerical Keff at the discontinuity tends to zero 
%    making a zero-row in the matrix (det(A) = 0). 
% 
% A patch could be bound the values of the matrix. It works for n <= 0.3  
% but lower than 0.2 the solution is not the correct one. Is still does
% not work for n > 2.
% * Use preconditioning?
% * Use Finite Volume Method and a flux limiter for the diffusive flux?
% * Use the exact solution to build a discrete version?
%
% If the solution has an inflection point (flex), then d^2u/dz^2 = 0 at the
% flex point. Therefore, the differential equation writes as follows:
%
%  dKeff |       du |                                   
%  ----- |     * -- |      = -g * etax    
%   dz   |_flex  dz |_flex    
% 
% It looks like at the flex point the behaviour of the fluid is like a pure 
% convective flow. Determine how the information is propagated with the 
% sign of dKeff/dz and discretize the above differential equation.
%
%==========================================================================

clear all
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain
H = 50; %30
% Parameters
K0 = 1e-3;
rho = 1000;
n = 0.5;
etax = -1e-6;
tauw = 1e-2;
g = 9.807;
% Choose
induced = 'no wind-driven';
viscosity = 'visco 1';
zm = -17; %-5
zh = -32; %-10
m = 14; %5
solver = 'no TDMA';
powerlaw = 'no modified';
ShearRateMin = 1e-6;
ShearRateMax = 1e+5;
% Discretization
N = 200;
IterMax = 1000;
Tol = 1e-5;
IterMaxQ = 1000;
TolQ = 1e-3;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case of a pure wind-driven flow (only for n = 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( induced, 'wind-driven' )
% For this case, Q = 0 -> int(u,-H,0) = 0    
    n = 1;
    etax = - ( 3 * tauw ) / ( 2 * rho * g * H );
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesher
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = linspace( 0, -H, N )';
Dz = ( H - 0 ) / ( N - 1 );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective eddy viscosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( viscosity, 'visco 0' )
    K =@( z ) K0;
elseif strcmp( viscosity, 'visco 1' )
    alpha = 1 / ( 2 * zh * ( zm - zh ) / m - zh * ( zh - 2 * zm ) );
    beta  = 2 * alpha * zh * ( zm - zh ) / m;
    % Danger! K( z( 1 ) ) = NaN
    K =@( z ) ( K0 .* ( 1 - 2 .* alpha .* zm .* z + alpha .* z .^2 ) ) .* ( z > zh ) + ...
              ( K0 .* beta .* abs( z ./ zh ) .^ ( - m ) ) .* ( z <= zh );
end
Keff =@( dudz, z ) K( z ) .* abs( dudz ) .^ ( n - 1 );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uini =@( z ) 1 - ( z ./ H ) .^ 2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Allocation
b  = zeros( N, 1 );
Ac = zeros( N, 1 );
Aw = zeros( N, 1 );
Ae = zeros( N, 1 );
% Source term at boundary
b ( 1 ) =  Dz * ( tauw / ( rho * K0 ) ) ^ ( 1 / n );
% Coefficients of the matrix for the boundaries
Ac( 1 ) =  1;
Ae( 1 ) = -1;
Ac( N ) =  1;
% Coefficient 2: ||A||_1 * |Q| <= k_1(A) * [ C0 + C1 * etax ]
C1 = ( N - 2 ) * Dz^2 * g;
for iterQ = 1 : IterMaxQ
    fprintf( '========== ITERATION Q %d ==========\n', iterQ )
    % Source term
    b( 2 : N - 1 ) = -Dz ^ 2 * g * etax;
    % Inicialization
    Uold = uini( z );
    % Residual vector
    RES = 0;
    % Bounding the value of the matrix
    if strcmp( powerlaw, 'modified' )
        gamma2 = ShearRateMax;
        gamma1 = ShearRateMin;
    else
        gamma2 = realmax;
        gamma1 = realmin;
    end
    % Main loop   
    for iter = 1 : IterMax
        % Computing the modified gradient of the velocity field
    %     gammaDot(1) -> \dot{\gamma}_{1/2}, gammaDot(2) -> \dot{\gamma}_{3/2}, ...
        gammaDot = ( Uold( 2 : N ) - Uold( 1 : N - 1 ) ) ./ Dz;
        gradU = min( max( abs( gammaDot ),  gamma1 ), gamma2 );
        gradU( gradU == gamma1 ) = 1.0;
        % Calculation of the coefficients of the matrix
        Aw( 2 : N - 1 ) = Keff( gradU( 1 : N - 2 ), ( z( 1 : N - 2 ) + z( 2 : N - 1 ) ) / 2 );
        Ae( 2 : N - 1 ) = Keff( gradU( 2 : N - 1 ), ( z( 2 : N - 1 ) + z( 3 : N     ) ) / 2 );
        Ac( 2 : N - 1 ) = - ( Aw( 2 : N - 1 ) +  Ae( 2 : N - 1 ) );
        % Build the matrix in matricial form
        A = diag( Ac, 0 ) + diag( Aw( 2 : N     ), -1 ) + ...
                            diag( Ae( 1 : N - 1 ) , 1 );
        % Calculation of the velocity field 
        if strcmp( solver, 'TDMA' )
            % TDMA solver
            U = TDMA( Aw, Ac, Ae, b );
        else
            % Matlab direct solver
            U = A \ b;
        end
        % Calculation of the maximum velocity (in absolute value)
        Umax = max( abs( U ) );
        % Calculation of the residual
        error = ( U - Uold ) ./ U;
        res = sum( abs( error( 1 : end - 1 ) ) ); 
        RES = [ res; RES ];
        fprintf( 'iter = %d res = %d\n', iter, res ) 
        % Stop the iterative procedure
        if res <= Tol
            break
        end
        % old <- new
        Uold = U;
        % Plotting the solution for each iteration
        plot( U ./ Umax, z ./ H, 'b-' )
        grid on
        set( gcf, 'color', 'w' )
        xlabel( '$u/|u_{\mathrm{max.}}|$', 'interpreter', 'latex' )
        ylabel( '$z/H$', 'interpreter', 'latex' )
        %pause(0.1)
        drawnow
    end
% Volumetric flow rate
Q = trapz( z, U );
fprintf( 'Q = %d m^3/s\n', Q ) 
% Stop the iterative procedure
if abs(Q) <= TolQ
    break
end
% Increment of etax
Detax = norm( A, 1 ) * abs( Q ) / ( cond( A, 1 ) * C1 );
% Update etax
etax = etax + sign( Q ) * Detax;
fprintf( ' \n ' )
end
toc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Calculation of the flow rate
% Q = trapz( z, U );
% Calculation and plotting the gradient and laplacian of the velocity field
dUdz = gradient( U, Dz );
dUdzmax = max( abs( dUdz ) );
d2Udz2 = del2( U, Dz );
d2Udz2max = max( abs( d2Udz2 ) );
figure
subplot(1,2,1)
plot( dUdz ./ dUdzmax, z ./ H, 'b-' )
grid on
set(gcf,'color','w')
xlabel('d$u/$d$z$ / $|($d$u/$d$z)_{\mathrm{max.}}|$','interpreter','latex')
ylabel('$z/H$','interpreter','latex')
subplot(1,2,2)
plot( d2Udz2 ./ d2Udz2max, z ./ H, 'b-' )
xlim( [ 0 1 ] )
grid on
set(gcf,'color','w')
xlabel('d$^2u/$d$z^2$ / $|($d$^2u/$d$z^2)_{\mathrm{max.}}|$','interpreter','latex')
ylabel('$z/H$','interpreter','latex')
% Plotting the eddy viscosity and the effective eddy viscosity
vK( 1 : N ) = K( z( 1 : N ) );
vKmax = max( vK );
vKeff = K( z ) .* abs( dUdz ) .^ ( n - 1 );
vKeffmax = max( vKeff );
figure
subplot(1,2,1)
plot ( vK ./ vKmax, z ./ H, 'b-' )
grid on
set(gcf,'color','w')
xlabel('$K/K_{\mathrm{max.}}$','interpreter','latex')
ylabel('$z/H$','interpreter','latex')
subplot(1,2,2)
plot ( vKeff ./ vKeffmax, z ./ H, 'b-' )
grid on
set(gcf,'color','w')
xlabel('$K_{\mathrm{eff.}}/K_{\mathrm{eff.,max.}}$','interpreter','latex')
ylabel('$z/H$','interpreter','latex')
% Plotting the residual
RES = flip( RES );
figure
iterations = [ 1 : 1 : length( RES ) ];
plot( iterations, RES )
grid on
set(gcf,'color','w')
xlabel('Iteration [-]','interpreter','latex')
ylabel('Residual [-]','interpreter','latex')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
