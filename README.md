## Fast, threaded scripts for comparing proteomes with Needleman\-Wunsch alignment (all vs all) 

These scripts exploit existing fast implementations of the Needleman-Wunsch global alignment algorithm (opal and parasail) to quickly find all closely related proteins between multiple proteomes and compile these into a small database. This supports an upcoming publication.

### Dependencies
You will need to have installed:
* [opal](https://github.com/Martinsos/opal) (opal\_aligner)
* [parasail](https://github.com/jeffdaily/parasail) (with python3 API)
* Common python3 libraries - see top of build\_homology\_database.py and unpack\_pickle.py for these

### How to use
1. Clone this repo to a folder. The subfolder "./genomes/" is required.
1. After compiling opal put a softlink/shortcut to the program "opal\_aligner" into *this* folder (i.e. into "./" - it is expected by build\_homology\_database.py).
1. Download proteomes and place them into the "./genomes/" subfolder. The expected format is multiline fasta (\*.fa) with standard amino acid code (see Identity\_score\_matrix.txt for expected alphabet).

1. (Optional) Parameters (including number of threads used, percent identity cutoffs, and output file names) can be modified at the top of build\_homology\_database.py and unpack\_pickle.py. The default number of threads is 8 as this was empirically found to run the quickest on a quad-core hyperthreaded CPU. Note: do not use more threads than the square of the number of genomes (ie if using 2 genomes then reduce threads to four or less)
1. Run ./build\_homology\_database.py. This will output a pickle (mhdb\_\*.pickle). With 8 mycobacterial genomes this ran for about 7 hours on an Intel i7-7700HQ CPU @ 2.80GHz.
1. Run ./unpack\_pickle.py. This will output an excel file (mhdb\_\*.xlsx) containing data from the all vs all comparisons.

If you have issues or feedback, please feel free to contact me:

Will Matern - maternwill@gmail.com
