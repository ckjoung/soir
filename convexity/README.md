# Convexity

## Compilation
```
    g++ -O2 main.cpp -o convex `~/lib/capd-5.3.0/build/bin/capd-config --cflags --libs`
```

## Test cases
```
    ./convex 0.0 0.005 2.1 2.100001 0.005 0.01 
    
```

## Producing the results from the paper:
To reproduce the results of the paper, run the zsh shell scripts:
zsh batch00_10.sh
zsh batch10_20.sh
zsh batch20_25.sh
zsh batch25_30.sh
zsh batch30_35.sh
zsh batch35_40r.sh
zsh batch400_425n.sh
zsh batch425_450n.sh
zsh batch450_475n.sh
zsh batch475_50n.sh

These are for the ranges
mu=0.0 to 0.1
mu=0.1 to 0.2
mu=0.2 to 0.25
mu=0.25 to 0.3
mu=0.3 to 0.35
mu=0.35 to 0.4
mu=0.4 to 0.425
mu=0.425 to 0.45
mu=0.45 to 0.475
mu=0.475 to 0.5


