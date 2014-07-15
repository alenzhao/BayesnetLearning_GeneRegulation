GBNet: (G)ibbs sampling enhanced (B)ayesian (Net)work
=====
Li Shen, UCSD
Feb. 11, 2008

Reference:
Shen, L., Liu, J. and Wang, W. (2008) GBNet: Deciphering regulatory rules in the co-regulated genes using a Gibbs sampler enhanced Bayesian network approach, BMC Bioinformatics, 9, 395.


Bayesian network learning of gene regulatory rules. This is a toolset for inferring transcriptional rules from a group of co-regulated genes using Bayesian networks. This toolset is separated into four programs: func, bayescor, bbnet AND gbnet. bbnet and gbnet use bayescor's OUTPUT as INPUT. bayescor, bbnet and gbnet all depend on the functional depth files from func.

See Readme.txt for details. There is also a test run example under the BN_example folder.


* func: Prepare the functional depth files for a motif list on all sequences.

Example: func -m motif.list -w pwm -g genome -f func_folder -n norm.txt

* bayescor: Calculate the score of each single TF's presence on a cluster of genes.

Example: bayescor -m motif.list -n cluster.list -b bkg.list -f folder -o scores.list

* bbnet: Learn transcriptional regulatory rules from a cluster of genes.

Example: bbnet -s scores.list -n node.list -b bkg.list -f func -k 6.5 -o results_6.5.txt -c 50

* gbnet: Gibbs sampler enhanced Bayesian networks.

GBNet combines Gibbs sampling and Simulated annealing to search in a sequence constraints space trying to find a transcriptional module that is best supported by the data.

The usage of GBNet is almost the same with BBNet. The only difference is that GBNet uses simulated annealing so you'll need to specify the parameters for it or GBNet will use the default values.

To specify the parameters for SA:

Use: -sa repeats iterations max_changes alpha initial_temperature

E.g. -sa 30 10 300 0.9 10.0 will set SA to start from initial temperature=10.0 and goes down in exponential fashion by alpha = 0.9. This process will repeat for at most 30 times. During each repeat, SA will run 10 iterations or make 300 changes to Bayesian network structure, whichever comes first.

If Bayesian network doesn't make any change under a certain temperature after enough iterations, 
the process stops assuming the ground zero status is achieved.












