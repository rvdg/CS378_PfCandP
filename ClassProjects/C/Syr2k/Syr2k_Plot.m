output
close all

% gflops = number of flops (in billions )
gflops = 2 * data_ref( :,1 ) .* data_ref( :, 1 ) .* data_ref( :, 1 ) ...
    * 1.0e-9;  
plot( data_ref( :,1 ), gflops ./ data_ref( :, 2 ), 'k-x' ); 

hold on    % Plot additional data in same figure

plot( data_unb_var1( :,1 ), gflops ./ data_unb_var1( :, 2 ), 'b-.o' ); 

hold off   % Stop plotting additional data in same figure% Plot additional data in same figure

% Add x- and y-axis labels and legend
xlabel( 'matrix dimension sn' );
ylabel( 'GFLOPS' );
legend( 'FLA_Syr2k', ...
        'Syr2k_xx_unb_var1', ...
        'Location', 'NorthWest' );
 print('Plot', '-dpdf')
 


