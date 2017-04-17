#/bin/sh
# Execute experiments for paper 2

echo Executing tests 1 to 4
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_01.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_02.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_03.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_04.csv &
wait

echo Executing tests 5 to 8
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_05.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_06.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_07.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_08.csv &
wait


echo Executing tests 9 to 12
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_09.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_10.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_11.csv &
./faststi -f ../parms/parms_40m_3yr_1sim.txt -s $RANDOM > output_40m_3yr_1sim_12.csv &
wait
