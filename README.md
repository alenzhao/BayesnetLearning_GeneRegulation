GBNet: (G)ibbs sampling enhanced (B)ayesian (Net)work
=====

Reference:
Shen, L., Liu, J. and Wang, W. (2008) GBNet: Deciphering regulatory rules in the co-regulated genes using a Gibbs sampler enhanced Bayesian network approach, BMC Bioinformatics, 9, 395.


Bayesian network learning of gene regulatory rules. This is a toolset for inferring transcriptional rules from a group of co-regulated genes using Bayesian networks. 
--
Li Shen, UCSD
Feb. 11, 2008

This toolset is separated into four programs: func, bayescor, bbnet AND gbnet.

bbnet and gbnet use bayescor's OUTPUT as INPUT.

bayescor, bbnet and gbnet all depend on the functional depth files from func.

******************************************************************************
* func: Prepare the functional depth files for a motif list on all sequences.*
******************************************************************************

INPUT: motif list, PWM folder, genomic sequences file, functional depth folder 
and normalization constant file
OUTPUT: all motifs functional depth files in a folder
Arguments:
-m  motif list
-w  PWM folder where you should put all matrix files
-g  a single file contains all promoter sequences (TAB delimited)
-f  where all functional depth files will go
-n  [optional] a file contains all motifs' normalization constants

Example: func -m motif.list -w pwm -g genome -f func_folder -n norm.txt

************************************************************************************
* bayescor: Calculate the score of each single TF's presence on a cluster of genes.*
************************************************************************************

INPUT: motif list, node gene list, background gene list and folder to store binding information
OUTPUT: single motif score list
Arguments:
-m  motif list
-n  node gene list
-b  background gene list
-f  folder to store binding information
-o  motif score list

Example: bayescor -m motif.list -n cluster.list -b bkg.list -f folder -o scores.list


*************************************************************************
* bbnet: Learn transcriptional regulatory rules from a cluster of genes.*
*************************************************************************

INPUT: motif score list, node gene list, background gene list, folder to store binding information and logK value
OUTPUT: results and parameter settings of Bayesian network running
Arguments:
-s  motif score list
-n  node gene list
-b  background gene list
-f  folder to store binding information
-o  results output file
Optional:
-k  logK value (network complexity penalization)(default = 5.0)
-c  number of candidate motifs (optional, default = 50)
-d  positive negative (for prediction, default = NULL) 
This parameter supplies two files containing known positive and negative cases. bayesnet will
use the regulatory rules learnt from BN to predict these genes' categories and output TP, FP, TN and FN.
-l  all training genes' information of satisfying rules (optional, default = NO output)
-t  all genes' translational/transcriptional start sites locations in TAB delimited format (optional, default = right end)
-rb     bit-string to determine which rules to include.(Default = 111110)
This parameter can be used to "mask" out certain rules that do not make sense in your situation
Rule order: TSS Orientation Second copy Spacing Order   Loop
Default=     1      1            1         1      1  0
-i      Use mutual information instead of Bayesian score.(Default = off)

Example: bbnet -s scores.list -n node.list -b bkg.list -f func -k 6.5 -o results_6.5.txt -c 50


*******************************************
* Explanation of motif binding file format*
*******************************************

1. file name: motif.func, Eg. YY1.func

Motif functional depth varies from 0.05 to 0.95 in a step of 0.05. All motifs' binding 
information must be stored in files as named above so that the programs can find them.

2. For each binding information file, the file format must follow this: 

Gene name\t number of binding sites\t binding site 1\t binding site 2\t ... \t binding site n

Separated by TAB, Example:

Hs.106529   2   F,0.054,772 R,0.97,230

3. For each binding site, the format must follow this:

orientation,binding score,distance to TSS

Separated by comma, Example:
F,0.054,772
That means: the motif is binding in Forward orientation with matrix score 0.054 at 772 bps upstream from TSS.


***************************************************
* gbnet: Gibbs sampler enhanced Bayesian networks.*
***************************************************

GBNet combines Gibbs sampling and Simulated annealing to search in a sequence constraints space 
trying to find a transcriptional module that is best supported by the data.

Pros:
    - Less prone to local minima than greedy search
    - Search exaustively in stochastic fashion
    - Increase Bayesian score and find rules that are more meaningful
Cons:
    - Large computational cost

The usage of GBNet is almost the same with BBNet. The only difference is that
GBNet uses simulated annealing so you'll need to specify the parameters for it
or GBNet will use the default values.

To specify the parameters for SA:

Use: -sa repeats iterations max_changes alpha initial_temperature

E.g. -sa 30 10 300 0.9 10.0 will set SA to start from initial temperature=10.0 
and goes down in exponential fashion by alpha = 0.9. This process will repeat
for at most 30 times. During each repeat, SA will run 10 iterations or make
300 changes to Bayesian network structure, whichever comes first.

Some default parameters:

Starting Temp:  5.0 Temperature at the initial point
Repeat:     20  Number of times that temperature changes
Iteration:  20  Number of iterations at each temperature
Changes:    500 Number of required changes to Bayesian network before
            moving to the next temperature
Alpha:      0.9 Temperature change rate

If Bayesian network doesn't make any change under a certain temperature after enough iterations, 
the process stops assuming the ground zero status is achieved.












