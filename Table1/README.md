Here we detail the steps for reproducing Figure 2: Genome-wide Manhattan plots of pleiotropy signals and Table of Method Comparision. 

1. Open each of the .R scripts and change the top few lines where necessary, such as to point to the correct output folder.

2. Make necessary changes to .lsf files, such as using the correct queue names, email address, and runtimes.

4. Run the submit_summarize_ukb.lsf and submit_threshold_table.lsf jobs.

5. After step 4 is completed, run the run_pleio_a1.sh and run_pleio_a2.sh jobs.

6. Run the table1.R script to create Table of Method Comparision.



 
