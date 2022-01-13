#!/bin/bash
# 
# Help
Help()
{
   # Display Help
   echo "DIPEND is a Python-based pipeline to build random structures for intrinsically disordered protein segments, making use of neighbor-dependent backbone conformation preferences described in Ting et al. PLoS Comput Biol 6(4): e1000763. "
   echo
   echo "It uses ChimeraX, Scwrl4 and GROMACS as external tools."
   echo
   echo "Input parameteres can be supplied as options or listed in the Data/parameters.in file."
   echo "-t --threadnum: the number of the threads (parallel running processes), 1 by default"
}
angletoadd=-99
base=""
cycle=-99
dataset=""
gmxcheck=-99
keep=-99
mode=""
numofstructures=-99
proline=-99
sequence=""
random=-99
threads=-99
unknotmax=-99
paramfile="/home/zita/Scripts-Research/DIPEND/Data/parameters.in"
# Get the options
while getopts "a:b:c:d:g:k:m:n:p:r:s:t:u:h" option; do
   case $option in
      a) angletoadd=$OPTARG;;
      b) base=$OPTARG;;
      c) cycle=$OPTARG;;
      d) dataset=$OPTARG;;
      g) gmxcheck=$OPTARG;;
      k) keep=$OPTARG;;
      m) mode=$OPTARG;;
      n) numofstructures=$OPTARG;;
      p) proline=$OPTARG;;
      r) random=$OPTARG;;
      s) sequence=$OPTARG;;
      t) threads=$OPTARG;;
      u) unknotmax=$OPTARG;;
      h) # display Help
         Help
         exit;;
   esac
done
# read paramfile
while read -r line; do 
set -- "$line" 
IFS=" " && arrline=($*) 
case ${arrline[0]} in
      \#) true;;
      a) 
        if [ $angletoadd -eq -99 ];  then 
            angletoadd=${arrline[1]} 
        fi 
        ;;
      b) 
        if [ -z $base ]; then
            base=${arrline[1]}
        fi
        ;;
      c) if [ $cycle -eq -99 ]; then
            cycle=${arrline[1]}
        fi
        ;;
      d) if [ -z $dataset ]; then
            dataset=${arrline[1]}
        fi
        ;;
      g) if [ $gmxcheck -eq -99 ]; then
            gmxcheck=${arrline[1]}
        fi
        ;;
      k) if [ $keep -eq -99 ]; then
            keep=${arrline[1]}
        fi
        ;;
      m) if [ -z $mode ]; then
            mode=${arrline[1]}
        fi
        ;;
      n) if [ $numofstructures -eq -99 ]; then
            numofstructures=${arrline[1]}
        fi
        ;;
      p) if [ $proline -eq -99 ]; then
            proline=${arrline[1]}
        fi
        ;;
      r) if [ $random -eq -99 ]; then
            random=${arrline[1]}
        fi
        ;;
      s) if [ -z $sequence ]; then
            sequence=${arrline[1]}
        fi
        ;;
      t) if [ $threads -eq -99 ]; then
            threads=${arrline[1]}
        fi
        ;;
      u) if [ $unknotmax -eq -99 ]; then
            unknotmax=${arrline[1]}
        fi
        ;;
esac
done < $paramfile
# calculate distribution of structures into thread groups
if [ $threads -gt $numofstructures ]
then
    threads=$numofstructures
fi
if [ $threads -lt 1 ]
then
    threads=1
fi
groupsize=1
groupsize=$((numofstructures/threads))
# beginning execution
first=1
x=1
y=$((groupsize+1))
if [ $y -gt $numofstructures ]
then
  y=$numofstructures
fi
# first time running
    nohup python3 ~/Scripts-Research/DIPEND/Dipend.py -a $angletoadd -b $base -c $cycle -d $dataset -f 1 -g $gmxcheck -k $keep -m $mode -n $numofstructures -p $proline -r $random -s $sequence -t $threads -u $unknotmax -x $x -y $y > out1 2>&1 &
sleep 10 # (in seconds) we have to wait a bit for the first thread to make the initial structure for the other threads
x=$((x+groupsize))
y=$((y+groupsize))
# from the second thread (if any)
if [ $threads>=2 ] 
then
for (( i=2; i<$threads; i++ ))
do
    nohup python3 ~/Scripts-Research/DIPEND/Dipend.py -a $angletoadd -b $base -c $cycle -d $dataset -f 0 -g $gmxcheck -k $keep -m $mode -n $numofstructures -p $proline -r $random -s $sequence -t $threads -u $unknotmax -x $x -y $y > out$i 2>&1 &
x=$((x+groupsize))
y=$((y+groupsize))
if [ $y -gt $numofstructures ]
then
  y=$numofstructures
fi
done
# last time running
    nohup python3 ~/Scripts-Research/DIPEND/Dipend.py -a $angletoadd -b $base -c $cycle -d $dataset -f 2 -g $gmxcheck -k $keep -m $mode -n $numofstructures -p $proline -r $random -s $sequence -t $threads -u $unknotmax -x $x -y $y > out$threads 2>&1 &
fi