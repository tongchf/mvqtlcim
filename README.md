# Introduction
MVQTLCIM is a software package for composite inverval mapping QTLs with multivariate phenotype data in an outbred full-sib family. The software utilizes a genetic linkage map constructed with different segregation molecular markers such as 1:1, 1:2:1 and 1:1:1:1, and assumes that QTL may segregate in the five different segregation patterns on a specific position of the genetic map. It allows users to select the best QTL segregation pattern with AIC, BIC and TIC for a significant QTL. It also provides command line parameters to be chosen for alternative analyses, including the number of background markers, window size ([Basten, et al., 1994](http://statgen.ncsu.edu/qtlcart/WinQTLCart.pdf)), QTL segregation type, genetic map function and number of permutations. Specifically, when performing permutations to determine the empirical threshold of significant QTLs, `mvqtlcim` permits to use multithreads to accelerate computing speed. When an analysis completes, the software will generate two files for each QTL model, of which one contains the parameter estimates and the corresponding statistic values at every 1 cM on the genome, and the other saves the maximum LR value of each permutation. With these result files, we wrote an R script, `lrPlot.r`, to summarize the significant QTL information and generate scatter plots of LR against genome position. These plots can be optionally saved in pdf, jpg, png, tif or bmp format. 
# Usage
With the embedded example input file `example.txt`, users can get started with the following commands:  
`./mvqtlcim -i example.txt -b 5 -m 1`  
`./mvqtlcim -i example.txt -b 5 -w 15.0 -p 1000 -t 20`  
`Rscript lrPlot.r -i example_CIM_M1Rst.txt -p example_CIMPermuM1Rst.txt -w 1`  
  
One can directly run the commands, `./mvqtlcim` and `Rscript lrPlot.r`, to show the usages of `mvqtlcim` and `lrPlot.r`, respectivley, as follows  

<table><tr><td bgcolor=#7FFFD4>   abcd </td></tr></table> 
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
The input data is formated in a plain text as the following example file, where 

	Pop_Type:	F1  
	Sample_Size:  300  
	Linkage_Group_Number:  5  
	Marker_Number_Within_Groups:  
	LG1	11  
	LG2	11  
	LG3	11  
	LG4	11  
	LG5	11  

	Marker_Data:
	LG	Marker	Distance	P1_Phase	P2_Phase	#1	#2	#3	#4	#5	#6	#7	#8	#9	#10	#11	#12	#13	#14	#15	#16	#17	#18	#19	#20	#21	#22	#23	#24	#25	#26	#27	#28	#29	#30	#31	#32	#33	#34	#35	#36	#37	#38	#39	#40	#41	#42	#43	#44	#45	#46	#47	#48	#49	#50	#51	#52	#53	#54	#55	#56	#57	#58	#59	#60	#61	#62	#63	#64	#65	#66	#67	#68	#69	#70	#71	#72	#73	#74	#75	#76	#77	#78	#79	#80	#81	#82	#83	#84	#85	#86	#87	#88	#89	#90	#91	#92	#93	#94	#95	#96	#97	#98	#99	#100	#101	#102	#103	#104	#105	#106	#107	#108	#109	#110	#111	#112	#113	#114	#115	#116	#117	#118	#119	#120	#121	#122	#123	#124	#125	#126	#127	#128	#129	#130	#131	#132	#133	#134	#135	#136	#137	#138	#139	#140	#141	#142	#143	#144	#145	#146	#147	#148	#149	#150	#151	#152	#153	#154	#155	#156	#157	#158	#159	#160	#161	#162	#163	#164	#165	#166	#167	#168	#169	#170	#171	#172	#173	#174	#175	#176	#177	#178	#179	#180	#181	#182	#183	#184	#185	#186	#187	#188	#189	#190	#191	#192	#193	#194	#195	#196	#197	#198	#199	#200	#201	#202	#203	#204	#205	#206	#207	#208	#209	#210	#211	#212	#213	#214	#215	#216	#217	#218	#219	#220	#221	#222	#223	#224	#225	#226	#227	#228	#229	#230	#231	#232	#233	#234	#235	#236	#237	#238	#239	#240	#241	#242	#243	#244	#245	#246	#247	#248	#249	#250	#251	#252	#253	#254	#255	#256	#257	#258	#259	#260	#261	#262	#263	#264	#265	#266	#267	#268	#269	#270	#271	#272	#273	#274	#275	#276	#277	#278	#279	#280	#281	#282	#283	#284	#285	#286	#287	#288	#289	#290	#291	#292	#293	#294	#295	#296	#297	#298	#299	#300
	1	1	0	b|a	c|d	ac	bd	ad	bd	ad	bd	bd	ad	ac	ac	bc	bc	bc	ac	ac	bc	bc	bc	ac	bc	ac	bc	bd	ac	bc	bd	ac	ac	bc	ac	bc	bc	bd	ac	ad	ac	ac	ad	ad	ac	ad	ac	bd	ad	ac	bc	bd	bd	ac	ad	ac	bd	ac	ad	ad	ad	bc	ac	ac	ad	ac	ad	ac	bd	bd	bd	bd	bc	ad	ac	bc	ac	ad	bd	ad	bd	ad	ad	ad	ac	ad	ad	ac	bc	ad	bd	bd	bc	bc	ac	bc	ac	ad	ad	ac	bc	ac	bc	ac	ac	bd	bd	bd	ac	ad	bd	ac	ad	bd	bc	ad	ad	bc	ad	bd	ac	bc	bd	ac	ad	ac	ad	bc	ac	bc	bc	ac	ad	bc	bc	bc	bc	ad	ad	bd	bc	bc	bc	ac	ad	bd	bc	bd	ac	bc	bc	bd	bc	bc	bc	ac	ad	ac	bc	bc	bd	ad	ac	bd	ac	bc	bd	ad	bc	ac	bc	ac	bd	ad	ac	bd	ac	bd	ad	ad	bd	bc	bd	ac	bc	ac	bd	ad	bc	ad	ac	bd	bc	ac	ac	bc	bd	ad	ad	bd	bc	bd	ad	ac	ad	ad	ad	bc	ad	bd	ac	ac	ad	bd	bc	ac	ad	bc	bc	ad	ad	ad	ad	bc	ad	bc	bc	ac	bd	ac	ac	ac	ad	bc	bd	bc	bc	ac	bd	ac	ad	ac	ad	ac	bc	bc	bd	ad	bd	bd	ac	bd	ad	bd	bc	ac	bd	ac	ac	ac	bd	ac	bc	ad	ac	bc	ac	ad	ac	bd	bd	ac	bc	ac	bd	ad	bc	ac	ad	bc	bd	ad	ad	ac	bc	ad	bc	ad	ac	ac	ac	bd	ac	ac	bd	bc	bc	bd	bc	ac	bc	ad	ad	bd	ad

	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	aaaa bbbbb ccccc
	
