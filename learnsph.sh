#!/usr/local_rwth/bin/zsh

###SBATCH --partition=c23ms

### Job name
#SBATCH --job-name=learnSPH

### File / path where STDOUT & STDERR will be written
###  %J is the job ID, %I is the array ID
#SBATCH --output=learnSPH.%J.txt

### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]:minute,
### that means for 80 minutes you could also use this: 1:20
#SBATCH --time=0:10:00

### Request memory you need for your job per PROCESS in MB


### Request the number of compute slots you want to use
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

### Change to the work directory
$PWD

### Execute your application
### module load cmake
module load intel

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
cd build/app
./learnSPH_app

