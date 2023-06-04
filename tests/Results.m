clear all
close all
clc
%%
%
Sm2n1   = load( 'Kvariable_m2_n1.mat'   );
Sm2n04  = load( 'Kvariable_m2_n04.mat'  );
Sm2n06  = load( 'Kvariable_m2_n06.mat'  );
Sm2n08  = load( 'Kvariable_m2_n08.mat'  );
Sm2n12  = load( 'Kvariable_m2_n12.mat'  );
Sm2n14  = load( 'Kvariable_m2_n14.mat'  );
Sm2n16  = load( 'Kvariable_m2_n16.mat'  );
Sm2n18  = load( 'Kvariable_m2_n18.mat'  );
%
Sm8n1   = load( 'Kvariable_m8_n1.mat'   );
Sm8n04  = load( 'Kvariable_m8_n04.mat'  );
Sm8n06  = load( 'Kvariable_m8_n06.mat'  );
Sm8n08  = load( 'Kvariable_m8_n08.mat'  );
Sm8n12  = load( 'Kvariable_m8_n12.mat'  );
Sm8n14  = load( 'Kvariable_m8_n14.mat'  );
Sm8n16  = load( 'Kvariable_m8_n16.mat'  );
Sm8n18  = load( 'Kvariable_m8_n18.mat'  );
%
Sm14n1  = load( 'Kvariable_m14_n1.mat'  );
Sm14n04 = load( 'Kvariable_m14_n04.mat' );
Sm14n06 = load( 'Kvariable_m14_n06.mat' );
Sm14n08 = load( 'Kvariable_m14_n08.mat' );
Sm14n12 = load( 'Kvariable_m14_n12.mat' );
Sm14n14 = load( 'Kvariable_m14_n14.mat' );
Sm14n16 = load( 'Kvariable_m14_n16.mat' );
Sm14n18 = load( 'Kvariable_m14_n18.mat' );
%%
Un1max   = max( [ max( abs( Sm2n1.U    ) ), ...
                  max( abs( Sm8n1.U    ) ), ...
                  max( abs( Sm14n1.U   ) ) ] );
Unle1max = max( [ max( abs( Sm2n04.U   ) ), ...
                  max( abs( Sm2n06.U   ) ), ...
                  max( abs( Sm2n08.U   ) ), ...
                  max( abs( Sm8n04.U   ) ), ...
                  max( abs( Sm8n06.U   ) ), ...
                  max( abs( Sm8n08.U   ) ), ...
                  max( abs( Sm14n04.U  ) ), ...
                  max( abs( Sm14n06.U  ) ), ...
                  max( abs( Sm14n08.U  ) ) ] );
Unge1max = max( [ max( abs( Sm2n12.U   ) ), ...
                  max( abs( Sm2n14.U   ) ), ...
                  max( abs( Sm2n16.U   ) ), ...
                  max( abs( Sm2n18.U   ) ), ...
                  max( abs( Sm8n12.U   ) ), ...
                  max( abs( Sm8n14.U   ) ), ...
                  max( abs( Sm8n16.U   ) ), ...
                  max( abs( Sm8n18.U   ) ), ...
                  max( abs( Sm14n12.U  ) ), ...
                  max( abs( Sm14n14.U  ) ), ...
                  max( abs( Sm14n16.U  ) ), ...
                  max( abs( Sm14n18.U  ) )] );
H = max( abs( Sm2n1.z ) );
%%
figure
% Newtonian
subplot( 3, 3, 1 )
plot( Sm2n1.U  / Un1max, Sm2n1.z / H, 'b-' )
grid on
title( 'Newtonian' )
lgd = legend( { '$n = 1$' }, 'Location', 'southeast', ...
                             'Interpreter', 'latex' );
lgd.FontSize = 9;
ylabel( '$\zeta$', 'Interpreter', 'latex' )
text( -0.4, -0.1, '(a.1) $p = 2$', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
subplot( 3, 3, 4 )
plot( Sm8n1.U  / Un1max, Sm2n1.z / H, 'b-' )
grid on
ylabel( '$\zeta$', 'Interpreter', 'latex' )
text( -0.4, -0.1, '(a.2) $p = 8$', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
subplot( 3, 3, 7 )
plot( Sm14n1.U / Un1max, Sm2n1.z / H, 'b-' )
grid on
ylabel( '$\zeta$', 'Interpreter', 'latex' )
xlabel( '$u/|u_{\mathrm{max.}(n = 1)}|$', 'Interpreter', 'latex' )
text( -0.4, -0.1, '(a.3) $p = 14$', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
% Shear-thinning
subplot( 3, 3, 2 )
plot( Sm2n04.U / Unle1max, Sm2n04.z / H, ...
      Sm2n06.U / Unle1max, Sm2n06.z / H, ...
      Sm2n08.U / Unle1max, Sm2n08.z / H )
grid on
title( 'Shear-thinning' )
lgd = legend( { '$n = 0.4$', '$n = 0.6$', '$n = 0.8$' }, 'Location', 'southeast', ...
                                                         'Interpreter', 'latex' );
lgd.FontSize = 9;
text( -0.5, -0.1, '(b.1)', 'Interpreter', 'latex' )
xlim( [ -0.6 1 ] )
subplot( 3, 3, 5 )
plot( Sm8n04.U  / Unle1max, Sm8n04.z   / H, ...
      Sm8n06.U  / Unle1max, Sm8n06.z   / H, ...
      Sm8n08.U  / Unle1max, Sm8n08.z   / H )
grid on
text( -0.5, -0.1, '(b.2)', 'Interpreter', 'latex' )
xlim( [ -0.6 1 ] )
subplot( 3, 3, 8 )
plot( Sm14n04.U / Unle1max, Sm14n04.z  / H, ...
      Sm14n06.U / Unle1max, Sm14n06.z  / H, ...
      Sm14n08.U / Unle1max, Sm14n08.z  / H )
grid on
xlabel( '$u/|u_{\mathrm{max.}(n \leq 1)}|$', 'Interpreter', 'latex' )
text( -0.5, -0.1, '(b.3)', 'Interpreter', 'latex' )
xlim( [ -0.6 1 ] )
% Shear-thickening
subplot( 3, 3, 3 )
plot( Sm2n12.U  / Unge1max, Sm2n12.z  / H, ...
      Sm2n14.U  / Unge1max, Sm2n14.z  / H, ...
      Sm2n16.U  / Unge1max, Sm2n16.z  / H, ...
      Sm2n18.U  / Unge1max, Sm2n18.z  / H )
grid on
title( 'Shear-thickening' )
lgd = legend( { '$n = 1.2$', '$n = 1.4$', '$n = 1.6$', '$n = 1.8$' }, ...
              'Location', 'southeast', 'Interpreter', 'latex' );
lgd.FontSize = 9;
text( -0.4, -0.1, '(c.1)', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
subplot( 3, 3, 6 )
plot( Sm8n12.U  / Unge1max, Sm8n12.z  / H, ...
      Sm8n14.U  / Unge1max, Sm8n14.z  / H, ...
      Sm8n16.U  / Unge1max, Sm8n16.z  / H, ...
      Sm8n18.U  / Unge1max, Sm8n18.z  / H )
grid on
text( -0.4, -0.1, '(c.2)', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
subplot( 3, 3, 9 )
plot( Sm14n12.U / Unge1max, Sm14n12.z / H, ...
      Sm14n14.U / Unge1max, Sm14n14.z / H, ...
      Sm14n16.U / Unge1max, Sm14n16.z / H, ...
      Sm14n18.U / Unge1max, Sm14n18.z / H )
grid on
xlabel( '$u/|u_{\mathrm{max.}(n \geq 1)}|$', 'Interpreter', 'latex' )
text( -0.4, -0.1, '(c.3)', 'Interpreter', 'latex' )
xlim( [ -0.5 1 ] )
%
set(gcf,'color','w')
adjustpdfpage(gcf,20)
print(gcf,'figura2','-dpdf')