%==========================================================================
% NONASPEv21_Bathymetry_ZetaC
%
% ------------
% Description:
% ------------
%
% This script computes and plot the value of the parameter ζ_c along the 
% horizontal domain for a given rheological index n. The local bathymetry 
% of Ligeia Mare is incorporated to determine the variation of depth across 
% the basin.
%       
%                                                                        
% NONASPE                                                                 
% NON-newtonian Accurate Solver for Planetary Estuarines                        
%                                                                         
% Copyright © 2025 Llorente Lázaro, Víctor Javier                         
% Last Update: May 8, 2025                                               
% Website: https://sites.google.com/view/vjllorente                       
% Contact: vjllorente@unizar.com                               
%                                                                         
%--------------------------------------------------------------------------
%
% ------------
% Description:
% ------------
% This is an in-house CFD code which solves the boundary-value problem:
%
%   -                                                          o z = 0
%  | d  /      / du \ ^ n \       dη                           |
%  | -- | K0 * | -- |     | = g * -- , z in ( 0, -H )          |
%  | dz \      \ dz /     /       dx                           |
%  |                                                           |
% -|          / du \ ^ n                                       |
%  | ρ * K0 * | -- |     = 𝜏w        , at z = 0                |
%  |          \ dz /                                           |
%  |                                                           o z = -H
%  | u = 0                           , at z = -H
%   -
%
% where K0 is an eddy viscosity (assume constant), ρ is the density of
% the estuarine, dη/dx is a barotropic pressure gradient (assume constant), 
% and 𝜏w is the wind stress. g is the gravity constant.
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
%  -- | Keff * -- | = g * etax    
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
%==========================================================================

clear all
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
K0 = 1e-3;
rho = 420;
tauw = 1e-2;
g = 1.352;
I = 0.05;
Cmu = 0.09;


% Choose
boundscheme = 'order2';
induced = 'no wind-driven';
viscosity = 'visco 0';
flow = 'turbulent';
zm = -17;
zh = -32;
m = 10;
% Discretization
N = 400;
IterMax = 1000;
Tol = 1e-5;
IterMaxQ = 10000;
TolQ = 1e-6;

% Initialization of bathymetry-related values
% Rheological index
n_val = 1;
% Bathymetry profile (Mastrogiuseppe et al. (2014))
x_dom = linspace(0, 300, 300);
bathy_profile = zeros(1, length(x_dom));
for i = 1:length(x_dom)
    x = x_dom(i);
    bathy_profile(i) = -0.000018 * x.^3 + 0.012322 * x.^2 - 1.955119 * x - 60.654995;
end

% Store ζ_c results
zeta_c = zeros(size(x_dom));
z_c = zeros(size(x_dom));

for j = 1:length(x_dom)
    H = -bathy_profile(j);
    n = n_val;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Case of a pure wind-driven flow (only for n = 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp( induced, 'wind-driven' )
    % For this case, Q = 0 -> int(u,-H,0) = 0    
        n = 1;
        etax = ( 3 * tauw ) / ( 2 * rho * g * H );
    end
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mesher
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z = linspace( 0, -H, N )';
    Dz = ( H - 0 ) / ( N - 1 );
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Effective eddy viscosity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp( flow, 'laminar' )
        kavg =@( dUdz, z ) K0 .* abs( dUdz ) .^ ( n - 1 );
        kt =@(z) 0;
        kn =@(uprom, dudz, z) 0;
    
    elseif strcmp( flow, 'turbulent' )
        kavg =@( dUdz, z ) K0 .* abs( dUdz ) .^ ( n - 1 );
        
        alpha = 1 / ( 2 * zh * ( zm - zh ) / m - zh * ( zh - 2 * zm ) );
        beta  = 2 * alpha * zh * ( zm - zh ) / m;
        % Danger! K( z( 1 ) ) = NaN
        kt =@( z ) ( K0 .* ( 1 - 2 .* alpha .* zm .* z + alpha .* z .^2 ) ) .* ( z > zh ) + ...
                  ( K0 .* beta .* abs( z ./ zh ) .^ ( - m ) ) .* ( z <= zh );
        
        eps =@( uprom, z ) Cmu * I.^4 ./ kt( z ) * 1/4 .* uprom.^4;
        kn =@( uprom, dUdz, z ) ( ((n-1) * eps( uprom, z )) ./ ((dUdz).^2 + eps( uprom, z ) ./ kavg( dUdz, z )) );
    end
    Keff =@( uprom, dUdz, z ) kavg( dUdz, z ) + kt( z ) + kn( uprom, dUdz, z );
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial guess
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zetac = 0.67;
    etax = tauw / ( zetac * rho * g * H );
    al = ( g * etax ) / K0;
    be = tauw / ( rho * K0 );
    ga = n / ( al * ( n + 1 ) );
    p  = 1 + 1 / n; 
    uini =@( z ) ga .* ( abs( be + al .* z ) .^ p - abs( be - al .* H ) .^ p ); %1 - ( z ./ H ) .^ 2;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    % Allocation
    b   = zeros( N, 1 );
    Ac  = zeros( N, 1 );
    Aw  = zeros( N, 1 );
    Ae  = zeros( N, 1 );
    Aee = zeros( N, 1 );
    % Source term at boundary
    b( 1 ) =  Dz * abs( tauw / ( rho * K0 ) ) ^ ( 1 / n );
    % Coefficients of the matrix for the boundaries
    if strcmp( boundscheme, 'order1' )
        Ac( 1 ) =  1;
        Ae( 1 ) = -1;
    else
        Ac ( 1 ) =  3/2;
        Ae ( 1 ) = -2;
        Aee( 1 ) =  1/2;
    end
    Ac( N ) =  1;
    % Coefficient 2: ||A||_1 * |Q| <= k_1(A) * [ C0 + C1 * etax ]
    C1 = ( N - 2 ) * Dz ^ 2 * g;
    % Inicialization
    Uold = uini( z );
    %
    for iterQ = 1 : IterMaxQ
        fprintf( '========== ITERATION Q %d ==========\n', iterQ )
        % Source term
        b( 2 : N - 1 ) = Dz ^ 2 * g * etax;
        % Residual vector
        RES = 0;
        % Main loop   
        for iter = 1 : IterMax
            % Computing the modified gradient of the velocity field
        %   gradU(1) -> \frac{dU}{dz}_{1/2}, gradU(2) -> \frac{dU}{dz}_{3/2}, ...
            gradU = abs( ( Uold( 2 : N ) - Uold( 1 : N - 1 ) ) ./ Dz );
            % Calculation of the coefficients of the matrix
            uprom = [ ( Uold( 1 : N - 1 ) + Uold( 2 : N ) ) / 2; Uold( N ) ]; 
    
            Aw( 2 : N - 1 ) = Keff( uprom( 1 : N - 2 ), gradU( 1 : N - 2 ), ( z( 1 : N - 2 ) + z( 2 : N - 1 ) ) / 2 );
            Ae( 2 : N - 1 ) = Keff( uprom( 2 : N - 1 ), gradU( 2 : N - 1 ), ( z( 2 : N - 1 ) + z( 3 : N     ) ) / 2 );
            Ac( 2 : N - 1 ) = - ( Aw( 2 : N - 1 ) +  Ae( 2 : N - 1 ) );
            % Build the matrix in matricial form
            A = diag( Ac, 0 ) + diag( Aw( 2 : N     ), -1 ) + ...
                                diag( Ae( 1 : N - 1 ),  1 );
            if strcmp( boundscheme, 'order2' )
                A = A + diag( Aee( 1 : N - 2 ), 2 );
            end
            % Calculation of the velocity field 
            U = A \ b;
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
            %{
            % Plotting the solution for each iteration
            plot( U ./ Umax, z ./ H, 'b-' )
            grid on
            set( gcf, 'color', 'w' )
            xlabel( '$u/|u_{\mathrm{max.}}|$', 'interpreter', 'latex' )
            ylabel( '$z/H$', 'interpreter', 'latex' )
            %pause(0.1)
            drawnow
            %}
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
    etax = etax - sign( Q ) * Detax;
    fprintf( ' \n ' )
    end
    toc


    zeta_c(j) = tauw ./ (rho * g * etax * H);
    z_c(j) = -zeta_c(j) * H;
end


   

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
vK( 1 : N ) = kavg( z( 1 : N ) );
vKmax = max( vK );
vKeff = kavg( z ) + K0 .* abs( dUdz ) .^ ( n - 1 );
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


% Bathymetry mesh (X–Z grid)
[X, Z] = meshgrid(x_dom, linspace(0, max(bathy_profile), N));
for i = 1:length(x_dom)
    Z(:, i) = linspace(0, bathy_profile(i), N);
end

% Plotting
figure
plot( x_dom, bathy_profile, 'b-' )
grid on
set( gcf, 'color', 'w' )
xlabel('Distancia horizontal [km]', 'Interpreter', 'latex')
ylabel('Profundidad [m]', 'Interpreter', 'latex')
hold on
    
plot( x_dom, z_c, 'r-' )
grid on
set( gcf, 'color', 'w' )
hold off


    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
