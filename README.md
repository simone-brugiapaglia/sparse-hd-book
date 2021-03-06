# sparse-hd-book

This is a repository associated with the book:

*Sparse Polynomial Approximation of High-Dimensional Functions*
by Ben Adcock, Simone Brugiapaglia and Clayton G. Webster

to be published by SIAM in late 2021. It contains Matlab scripts to generate all figures in the book. 

If you have questions or comments about the code, or would like to view a draft copy of the manuscript, please contact ben_adcock@sfu.ca or simone.brugiapaglia@concordia.ca.

## Code organization

Files are organized into four main directories:

### src

Contains the main Matlab files used to create figures in each chapter

Organized in chapters

### utils

Contains various Matlab functions needed across main scripts

### data

Contains .mat files generated by the scripts in src

Organized in chapters

### figs

Contains .eps files generated by the scripts in src

Organized in chapters


### Filenames

These generically take the form 

fig\_[chpt][number]\_[col] or figs\_[chpt][number1]\_[chpt][number2]\_[row]\_[col] 

where [chpt] is the chapter number, [number], [number1] and [number2] are the figure numbers and [row] and [col] are the row number and column number is multi-panel figures.



