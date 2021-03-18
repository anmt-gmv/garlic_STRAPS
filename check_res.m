clear all
close all
clc

%% Plot results
set (0, "defaultaxesfontsize", 16)
set (0, "defaulttextfontsize", 16)
set (0, "defaultlinelinewidth", 2)


state_a = load('.\Session604\Session604\lsq_state.txt');
state_b = load('.\Session604_noglo\Session604_noglo\lsq_state.txt');

n_epochs = size(state_a(:,1),1);
dt = 0.5;
time = (1 : n_epochs ) * dt;

med_a = median(state_a(:,end));
med_b = median(state_b(:,end));

nsat_a = state_a( : , 4 : 7 );
mean_nsa_a = mean( nsat_a );
tot_sat_a = sum( nsat_a , 2 );

nsat_b = state_b( : , 4 : 7 );
mean_nsa_b = mean( nsat_b );
tot_sat_b = sum( nsat_b , 2 );

n_sat_diff = tot_sat_a - tot_sat_b;

figure(1)
plot( time, state_a(:,end) )
hold on
plot( time, state_b( : , end ) )
axis([time(1),time(end),0,100])
legend(sprintf('GPS-BEI-GLO, median=%d',med_a),sprintf('GPS-BEI, median=%d',med_b))
title('Ublox - Session 604 - TIR 10^{-7}')
ylabel('IBPL [m]')
xlabel('time [s]')

figure( 2 )
plot( time , state_a( : , 4 : 7 ) )
legend( sprintf('GPS, mean=%d',mean_nsa_a(1)),sprintf('GLO, mean=%d',mean_nsa_a(2)),...
sprintf('GAL, mean=%d',mean_nsa_a(3)),sprintf('BEI, mean=%d',mean_nsa_a(4))  )
title('Ublox - Session 604')
ylabel('# sats')
xlabel('time [s]')

figure( 3 )
plot( time , tot_sat_a ) 
hold on
plot( time , tot_sat_b ) 
ylabel('# sats')
xlabel('time [s]')
title('Ublox - Session 604')
legend('GPS-BEI-GLO','GPS-BEI')

figure( 4 )
plot( time , n_sat_diff )
ylabel('# sats')
xlabel('time [s]')
title('Ublox - Session 604')
legend('# GLO sats ')
