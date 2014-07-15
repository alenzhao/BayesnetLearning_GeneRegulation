** A toy example for BBNet and GBNet **

This shows an example of combinatorial regulation of PAC and RRPE in genes that are 
related with rRNA transcription in yeast. It is taken from Beer's Cell paper: Predicting gene 
expression from sequence. The problem here is small-scale: only 4 candidate motifs are 
involved. The two AlignACE motifs: M600 and M602 are learned from computational method 
and are highly similar with PAC and RRPE, respectively. The PAC and RRPE motifs are 
from literature. I included them to add some "noise" in the dataset. The purpose is to 
illustrate the ability of Bayesian networks to pick up the motifs that can best explain 
the data.

The functional depth files have been already calculated so you can run BBNet and GBNet 
directly. The cluster contains 114 genes while the background contains 1789
genes. What Bayesian network does is to learn the regulatory rules that can
best distinguish the cluster from the background. 

> To calculate Bayesian scores for the 4 motifs, type:

./bscor -m motif.list -n node.list -b bkg.list -f func -o scores.list

in command line. This will calculate the scores for the four motifs and output
them in "scores.list".

> To run BBNet based on the scores and func depth, type:

./bbnet -s scores.list -n node.list -b bkg.list -f func -k 6.5 -o
bb_results_6.5.txt > bb_log_6.5.txt

in command line. This will output results in "bb_results_6.5.txt" and also
record the running status of BBNet in "bb_log_6.5.txt". 

> Similarly, to run GBNet, type:

./gbnet -s scores.list -n node.list -b bkg.list -f func -k 6.5 -o
gb_results_6.5.txt -sa 40 20 500 0.9 10.0 > gb_log_6.5.txt

in command line. The simulated annealing parameters are tuned according to
experience. There is no golden standard to do this! However, the above
parameters will usually give you a good result for problems of similar size.
GBNet usually takes a long time to finish (4-6 hours) depending on the size of
a problem. It also yields a huge number of text lines telling you what he(she?) did
in changing the Bayesian network structure. The number of text lines can
easily exceed 30K so I would usually put them in a text file as record. 

You can build a non-verbose version of BBNet or GBNet. Please refer to README
files under their folders.

You don't need to worry about the number of candidate motifs. We tried with
"-c 50" or even "-c 666". It can still produce good results.

If lower logK values, like 4.0 or 3.0, are used, you can lower the initial temperature to 
some value like 5.0. This can avoid GBNet waste too much time in high
temperature. A good practice is to always set the initial temperature to be
higher than logK value, which allows SA to wander around the search space when
temperature is high.

Note: Beer used logK = 15 in his work but he used natural logrithm. And ln(15) = log10(6.5), 
that's why we set the logK equal 6.5 in this example. We want to compare GBNet and BBNet 
in exactly the same condition as what Beer performed BBNet in his Cell paper.

As you will see in the results of GBNet, it usually gives you two rules like
this:

"""""""""""""""""""""""""""
Bayesian score of the network: -96.1711		(named RSLT1)

Constraints:
1. Distance to TSS of M600:140, 0.6
2. Distance to TSS of M602:240, 0.7

Conditional probability table:
        k = 0   k = 1
00      1723    23
01      5       27
10      61      18
11      0       46
"""""""""""""""""""""""""""

In Beer's paper, he reported two rules like this:

""""""""""""""""""""""""""""""""""""
Bayesian score of the network: -96.2621		(named RSLT2)

Constraints:
1. Distance to TSS of M600:140, 0.6
2. Distance to TSS of M602:240, 0.65

Conditional probability table:
        k = 0   k = 1
00      1707    23
01      5       13
10      77      18
11      0       60
"""""""""""""""""""""""""""""""""""

So RSLT1 is slightly better than RSLT2 but with less genes satisfying both
rules. However, if you use BBNet, you won't be able to learn the two distance
rules. Instead, you will learn a distance rule of M600 and a presence rule of
M602 because when the distance rule of M602 is added, the Bayesian score
becomes worse, so for a Greedy search algorithm, it just stops there. But SA
algorithm has the chance to overcome this obstacle. Cool.


BTW, I haven't decided the gender of GBNet...

--
Li Shen, UCSD, Sep. 12, 2007



