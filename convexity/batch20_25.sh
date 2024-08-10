#!/usr/local/bin/zsh
zmodload zsh/mathfunc
for ((k=0.0; k <= 0.05; k+=0.005)); do
  let n=k+0.2
  let l=0.005-k/50
  #echo $l
  let a=n+0.005 
  let i=$(( int(rint(n*1000)) ))
  ./convex $n $a 2.1 2.100001 $l 0.01 > convex_n_$i 
done
