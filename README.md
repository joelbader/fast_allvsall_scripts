## Fast, threaded scripts for comparing proteomes with Needleman\-Wunsch alignment (all vs all) 

These scripts exploit existing fast implementations of the Needleman-Wunsch global alignment algorithm (opal and parasail) to quickly find all closely related proteins between multiple proteomes and compile these into a small database. This supports our upcoming publication:

Matern WM, Bader, JS, Karakousis PC. Genome data for Mycobacterium avium subsp. hominissuis strain 109. Scientific Data. 2018.


### Dependencies
You will need to have installed:
* [opal](https://github.com/Martinsos/opal) (opal\_aligner)
* [parasail-python](https://github.com/jeffdaily/parasail-python) (python3 API for parasail)
* Common python3 libraries - see top of build\_homology\_database.py and unpack\_pickle.py for these

### How to use
1. Clone this repo to a folder. The subfolder "./genomes/" is required.
1. After compiling opal put a softlink/shortcut to the program "opal\_aligner" into *this* folder (i.e. into "./" - it is expected by build\_homology\_database.py).
1. Download proteomes and place them into the "./genomes/" subfolder. The expected format is multiline fasta (\*.fa) with standard amino acid code (see Identity\_score\_matrix.txt for expected alphabet).

1. (Optional) Parameters (including number of threads used, percent identity cutoffs, and output file names) can be modified at the top of build\_homology\_database.py and unpack\_pickle.py. The default number of threads is 8 as this was empirically found to run the quickest on a quad-core hyperthreaded CPU. Note: do not use more threads than the square of the number of genomes (ie if using 2 genomes then reduce threads to four or less)
1. Run ./build\_homology\_database.py. This will output a pickle (mhdb\_\*.pickle). With 8 mycobacterial genomes this ran for about 7 hours on an Intel i7-7700HQ CPU @ 2.80GHz.
1. Run ./unpack\_pickle.py. This will output an excel file (mhdb\_\*.xlsx) containing data from the all vs all comparisons.

### Notes
* If you modify the proteomes contained in "./genomes/" then build\_homology\_database.py will only compute proteome comparisons which are needed in the new database (and does not recompute all comparisons). For example, if you run build\_homology\_database.py and then remove a proteome from ./genomes/ then NW alignment does not actually need to be called. However, the entries of the removed proteome will be deleted from the database.

If you have issues or feedback, please feel free to contact me:

Will Matern - maternwill@gmail.com
