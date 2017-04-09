#/bin/sh
# Execute experiments for paper 2

echo Executing tests 1 to 4
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_01.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_02.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_03.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_04.csv &
wait

echo Executing tests 5 to 8
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_05.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_06.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_07.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_08.csv &
wait


echo Executing tests 9 to 12
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_09.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_10.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_11.csv &
./simulate -f parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_12.csv &
wait
