#!/usr/local/bin/zsh
zmodload zsh/mathfunc
for ((k=0.025; k <= 0.050; k+=0.0025)); do
  let n=k+0.4
  let l=0.0015-k/200
  #echo $l
  let a=n+0.0025 
  let i=$(( int(rint(n*1000)) ))
  ./convex $n $a 2.1 2.100001 $l 0.01 > convex_n_$i 
done
