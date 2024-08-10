#!/usr/local/bin/zsh
zmodload zsh/mathfunc
for ((n=0.00; n <= 0.1; n+=0.005)); do
  let a=n+0.005 
  let i=$(( int(rint(n*1000)) ))
  ./convex $n $a 2.1 2.100001 0.005 0.01 > convex_$i 
done
