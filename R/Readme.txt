This folder is R code to achieve the results for the project Correspondence analysis: handling cell-wise outliers via the reconstitution algorithm

Section 3 How outliers originate:

smallmadeupexample.R

Section 6 Empirical studies:

The CA for the original matrix: origbrands.R, origop.R. 
The CA for the reconstitution of order h: reconbrands.R, reconop.R. 
The CA for the MacroPCA: macropcabrands.R, macropcaop.R.
The CA for the supplementary points method: supplebrands.R, suppleop.R.

In addition, reconca.R is a function to implement reconstitution of order h; the code is written by adjusting imputeCA function in the package missMDA version 1.19.