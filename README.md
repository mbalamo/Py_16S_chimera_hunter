# Py_16S_chimera_hunter
Python script to search FASTA long sequences of 16S for chimeric ends



Part 1: mock FASTA creator, chimera_16S_maker_v1.py

The purpose of this script is to make a test FASTA file of 16S long contigs to check for evidence of chimeric ends. The bacteria chosen are based on Zymo's 16S microbial control standards, which we can run as part of a process QC to look for evidence of chimeric molecules created as an unintended side-effect (by too many PCR cycles, bad reagents, etc.). The 16S gene in bacteria is appox 1500 bp composed of 9 variable regions that are unique to a genus or species. What we produce here is a mock file of assembled short into long reads or native long reads for the full 16S gene. The mock file can be used to test the accuracy of the "hunter script", and to see how a process faliure can lead to changes in truncated or chimeric molcueles. In this mock file, we control how many contigs we check, how the diferent bacteria are mixed (at what ratios), what fraction will be chimeric, and which fraction will be truncated. Contigs that are too short (truncated) will be missing the end(s) so we can't test for chimeric end if those ends are missing. Likewise, if a contig was not a perfect match to the "left end" and "right end" anchor sequences, we will report that contig as "undefined" even if full length. All of these "modes" are possible from real samples in a real libraray prep and sequencing event.

Considerations: user parameters for fraction chimeric, truncated, and undefined, must be numbers from 0 up to 100. Any value above or below will be auto-adjusted. The sum of those values must be <= 100. If they are above, each value will be scaled uniformly to become sum == 100. The minumum number of contigs = 100.
