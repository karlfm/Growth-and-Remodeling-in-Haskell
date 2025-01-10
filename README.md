# G&R Solver

A Haskell solver for growth and remodeling. Currently, I am trying to reproduce *On the definition and modeling of incremental, cumulative,
and continuous growth laws in morphoelasticity* by Alain Goriely and Martine Ben Amar (https://pubmed.ncbi.nlm.nih.gov/17123061/). 

## Use

Open GHCI and write
```Haskell
ghci> :l incremental_growth.hs 
```
and then simply write
```Haskell
ghci> main
```
to run the code. 

## Issues

Currently I am just setting a constant growth rate, such that growth is happening at a faster rate each step. 
In the paper they change the growth rate $\mu$ such that growth rate will be 1% each step. 
The only way I can think of to do this is by running two root finding methods at the same time, which resluts in very slow performance. 
