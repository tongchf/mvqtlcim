# Introduction
MVQTLCIM is a software package for composite inverval mapping QTLs with multivariate phenotype data in an outbred full-sib family. The software utilizes a genetic linkage map constructed with different segregation molecular markers such as 1:1, 1:2:1 and 1:1:1:1, and assumes that QTL may segregate in the five different segregation patterns on a specific position of the genetic map. It allows users to select the best QTL segregation pattern with AIC, BIC and TIC for a significant QTL. It also provides command line parameters to be chosen for alternative analyses, including the number of background markers, window size ([Basten, et al., 1994](http://statgen.ncsu.edu/qtlcart/WinQTLCart.pdf)), QTL segregation type, genetic map function and number of permutations. Specifically, when performing permutations to determine the empirical threshold of significant QTLs, `mvqtlcim` permits to use multithreads to accelerate computing speed. When an analysis completes, the software will generate two files for each QTL model, of which one contains the parameter estimates and the corresponding statistic values at every 1 cM on the genome, and the other saves the maximum LR value of each permutation. With these result files, we wrote an R script, `lrPlot.r`, to summarize the significant QTL information and generate scatter plots of LR against genome position. These plots can be optionally saved in pdf, jpg, png, tif or bmp format. 
# Usage
With the embedded example input file `example.txt`, users can get started with the following commands:  
`./mvqtlcim -i example.txt -b 5 -m 1`  
`./mvqtlcim -i example.txt -b 5 -w 15.0 -p 1000 -t 20`  
`Rscript lrPlot.r -i example_CIM_M1Rst.txt -p example_CIMPermuM1Rst.txt -w 1`  
  
One can directly run the commands, `./mvqtlcim` and `Rscript lrPlot.r`, to show the usages of `mvqtlcim` and `lrPlot.r`, respectivley, as follows  

    Usage: mvqtlcim <-i inputfile.txt> [options]  
    Options:  
            -b INT  	number of background markers [15]  
            -w FLOAT	window size [10.0]  
            -m STR		qtl models chosen from 1 to 5,such as '134' and '235' [12345]  
            -f INT		genetic map function, 1 for Haldane and 2 for Kosambi [1]  
            -p INT		number of permutaions [0]  
            -t INT		number of threads [1]  
            -a FLOAT	threshold of p-value for a marker to be added in stepwise regression [0.05]  
            -d FLOAT	threshold of p-value for a marker to be deleted in stepwise regression [0.10]  
--- - - -     
    Usage: Rscript lrPlot.r -i qtlrstfile [options]
    Options:
        	-p	str	result file for permutations
	        -w	int	just taking the value of 0, 1 or 2 [0]
			        0: manually giving the threshold value
			        1: genome-wide threshold determined by permutations
			        2: linkage-group-wide threshold determined by permutations
	        -s	float	significant level in(0.0001,0.2) [0.05]
	        -t	str	the plot format of 'pdf','png','jpg','tif' and 'bmp' [pdf]
	        -v	float	giving the threshold value manually

# Input data
