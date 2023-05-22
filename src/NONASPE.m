%==========================================================================
%                                                                         
% NONASPE                                                                 
% NON-newtonian Accurate Solver for Planetary Estuarines                        
%                                                                         
% Copyright © 2023 Llorente Lázaro, Víctor Javier                         
% Last Update: May 22, 2023                                               
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
%    Keff = K0 * | -- |         .
%                | dz |
%
% Therefore, the differential equation is re-written as follows:
%
%  d  /        du \                                   
%  -- | Keff * -- | = -g * etax    
%  dz \        dz /    
%
%
% ------
% INPUT:
% ------
% Double            :: H            - Depth [m]
% Double            :: K0           - Eddy viscosity [m^2/s]
% Double            :: rho          - Water density [kg/m^3] 
% Double            :: n            - Power exponent [-]
% Double            :: etax         - Barotropic pressure gradient [-]
% Double            :: g            - Gravity [m/s^2]
% String            :: induced      - Solve wind-driven for n = 1 [wind-driven / no wind-driven]
% String            :: solver       - Choose the solver [TDMA / no TDMA]
% String            :: trick        - Bound the values of the matrix [bounded matrix / no bounded matrix]
% Double            :: N            - Number of nodes [-]
% Double            :: IterMax      - Maximum of iterations [-]
% Double            :: Tol          - Tolerance [-]
% Double            :: Amax         - Upper bound value
% Double            :: Amin         - Lower bound value
% 
%
% -------
% OUTPUT:
% -------
% Double, Nx1 array :: U            - Velocity [m/s]
% Double, Nx1 array :: dUdz         - Velocity gradient [1/s]
% Double            :: Q            - Flow rate [m^3/s]
%
% ---------
% COMMENTS:
% ---------
% The solver works apparently for a range of 0.3 < n < 2 with the solver 
% of Matlab and 0.3 <= n <= 2 using TDMA. The gradient of u presents a 
% discontinuity at some point in the domain. 
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
%==========================================================================

clear all
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain
H = 10;
% Parameters
K0 = 1e-3;
rho = 1000;
n = 1;
etax = -1e-6;
tauw = 1e-2;
g = 9.807;
% Choose
induced = 'wind-driven';
solver = 'no TDMA';
trick = 'no bounded matrix';
% Discretization
N = 500;
IterMax = 1000;
Tol = 1e-5;
Amax = 1e+7;
Amin = 1e-5;

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
zi = linspace( 0, -H, N )';
Dz = ( H - 0 ) / ( N - 1 );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective eddy viscosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Keff =@( a, b ) K0 .* abs( ( a - b ) ./ Dz ) .^ ( n - 1 );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guess solution
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
% Source term
b( 1         ) =  Dz * ( tauw / ( rho * K0 ) ) ^ ( 1 / n );
b( 2 : N - 1 ) = -Dz ^ 2 * g * etax;
% Inicialization
Uold = uini( zi );
% Coefficients of the matrix for the boundaries
Ac( 1 ) =  1; 
Ae( 1 ) = -1;
Ac( N ) =  1;
% Residual vector
RES = 0;
% Bounding the value of the matrix
if strcmp( trick, 'bounded matrix' )
    UpperBound = Amax;
    LowerBound = Amin;
else
    UpperBound = realmax;
    LowerBound = realmin;
end
% Main loop   
for iter = 1 : IterMax
    % Calculation of the coefficients of the matrix
    Aw( 2 : N - 1 ) = max( min( Keff( Uold( 1 : N - 2 ), Uold( 2 : N - 1 ) ), UpperBound ), LowerBound ); 
    Ae( 2 : N - 1 ) = max( min( Keff( Uold( 2 : N - 1 ), Uold( 3 : N     ) ), UpperBound ), LowerBound );
    Ac( 2 : N - 1 ) = - ( Aw( 2 : N - 1 ) +  Ae( 2 : N - 1 ) );
    % Calculation of the velocity field 
    if strcmp( solver, 'TDMA' )
        % TDMA solver
        U = TDMA( Aw, Ac, Ae, b );
    else
        % Build the matrix in matricial form
        A = diag( Ac, 0 ) + diag( Aw( 2 : N     ), -1 ) + ...
                            diag( Ae( 1 : N - 1 ) , 1 );
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
    % new -> old
    Uold = U;
    % Plotting the solution for each iteration
    plot( U ./ Umax, zi ./ H, 'b-')
    grid on
    set(gcf,'color','w')
    xlabel('$u/|u_{max.}|$ [m/s]','interpreter','latex')
    ylabel('$z/H$ [m]','interpreter','latex')
    pause(0.1)%drawnow
end
toc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the flow rate
Q = trapz( zi, U );
% Calculation and plotting of the gradient of the velocity field
dUdz = gradient( U );
dUdzmax = max( abs( dUdz ) );
figure
plot( dUdz ./ dUdzmax, zi ./ H, 'b-' )
grid on
set(gcf,'color','w')
xlabel('d$u/$d$z$ / $|($d$u/$d$z)_{max.}|$ [1/s]','interpreter','latex')
ylabel('$z/H$ [m]','interpreter','latex')
% Plotting the effective eddy viscosity
vKeff = K0 .* abs( dUdz ) .^ ( n - 1 );
vKeffmax = max ( vKeff );
figure
plot ( vKeff ./ vKeffmax, zi ./ H, 'b-' )
grid on
set(gcf,'color','w')
xlabel('$K_{eff.}/K_{eff.,max.}$ [m$^2$/s]','interpreter','latex')
ylabel('$z/H$ [m]','interpreter','latex')
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
