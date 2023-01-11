#!/bin/bash

## FDIST island model coalescent simulations | individuals grouped by sample location in data set

# simulate 3 demes 
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes3.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes3_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes3_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res


# simulate 4 demes
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes4.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes4_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes4_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res


# simulate 9 demes
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes9.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes9_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes9_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res


# simulate 25 demes
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes25.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes25_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes25_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res


# simulate 50 demes
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes50.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes50_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes50_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res


# simulate 100 demes
./scr/arlecore3522_64bit data/POPGEN/BMA_by_pop.arp scr/fdist2_nsim20kndemes100.ars

# move and rename results files
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_ObsOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes100_ObsOut.txt
cp /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res/fdist2_simOut.txt /home/soleary/GAFFTOPS/BMA_POPGEN/results/fdist2_nsim20kndemes100_simOut.txt
rm -r /home/soleary/GAFFTOPS/BMA_POPGEN/data/POPGEN/BMA_by_pop.res

