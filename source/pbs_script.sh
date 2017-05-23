# This is just a sample PBS script. You should adapt the
# following such that you can run all your experiments.
# For more information about Torque and PBS, check the following sources:
#    https://support.cc.gatech.edu/facilities/instructional-labs/how-to-run-jobs-on-the-instructional-hpc-clusters
#    http://www.democritos.it/activities/IT-MC/documentation/newinterface/pages/runningcodes.html

#PBS -q class
#PBS -l nodes=jinx10+jinx11+jinx12+jinx13+jinx15+jinx17+jinx18
#PBS -l walltime=02:00:00
#PBS -N cse6220-prog2

# allocate 4 of the (24 total) sixcore nodes for up to 5 minutes

# TODO: change this to your project directory relative to your home directory
#       (= $HOME)
# TODO:
# try various different configurations and different input files
# it could be a good idea to create loops over different input files
# (or randomly generated on the fly using `generate_input` using different
# input sizes)
# and trying different computation depths for the master
# (again, you might want to add this to a loop and increase the total all
#  time accordingly)

#---------------------------------------------------
# Set by adsouza31
# 
export PBS_O_WORKDIR=$HOME/cse6220-prog3
EXE=$PBS_O_WORKDIR/sort
#export PBS_NODEFILE=$PBS_O_WORKDIR/host-mpich.txt
#
#
#---------------------------------------------------

# loop over number of processors (just an example, uncomment to use the loop)
for n in 1000 10000 15000 35000 50000 75000 100000 250000 500000 1000000
 do
  for p in 3 4 5 8 10 14 16 18 20 22 24 25 30
  do
    OMPI_MCA_mpi_yield_when_idle=0 mpirun -np $p $EXE -o $PBS_O_WORKDIR/output/output-$n.txt $PBS_O_WORKDIR/input/input-$n.txt > $PBS_O_WORKDIR/output/p-$p-n-$n.txt 2>&1
  done
done

