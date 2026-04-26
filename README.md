# csmGmm_CRC_IBD_reproduce

This repository includes 5 main parts:

1. Data cleaning

2. Regional pleiotropy analysis for CRC vs IBD
   - Run `Regional Analysis Plot.R`
   - Change the package path and file path, for example:
     ```r
     .libPaths("/home/swang25/R/ubuntu/4.4.1")
     ```

3. Figure 1 -- Manhattan plot
   - First, correct the path in `make_fig4.R`
   - Then submit `run_pleio.lsf`
   - Do not forget to modify the file path
   - Next, change the path in `plot_manhattan.R`
   - Then submit `submit_manhattan.lsf`

4. Table

5. Output files
