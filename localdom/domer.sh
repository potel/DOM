module load GSL/1.15
#make irreducible
#./irreducible
make dpdom
#./dpdom parametros48Ca.txt

qsub dpdom.sub
#cp misc1.txt /mnt/home/potel/Documents/gromacs_new/gromacs_install/bin
#cp misc2.txt /mnt/home/potel/Documents/gromacs_new/gromacs_install/bin
#cp misc3.txt /mnt/home/potel/Documents/gromacs_new/gromacs_install/bin
#cp misc4.txt /mnt/home/potel/Documents/gromacs_new/gromacs_install/bin
