# [Under construction!!]

# Basic principle

Each assembled genome was screened for contaminating scaffolds, i.e. scaffolds that are likely to be of bacterial (or human) origin. For this, we compiled different databases containing 1908 complete bacterial genome sequences, 43 complete insect genome sequences, or XX complete insect genome sequences excluding ants, as well as three databases containing corresponding bacterial or insect (with or without ants) CDS sequences. We next created 2000 bp sliding windows (overlap 500 bp) for each target ant genome and searched the different insect and bacterial databases using `mmseqs search` (release_12-113e3), with `--start-sens 1 --sens-steps 2 -s 7 --search-type 3` followed by `mmseqs convertalis`. Similarly, we searched all sliding windows against a database containing the human genome sequence (release **XX**). We filtered the single best blast hit (sorting by evalue (-k7,7g), then by bitscore (-k8,8gr)) for each sliding window using `sort -k1,1 -k7,7g -k8,8gr | sort -u -k1,1 --merge`. Highly conserved rRNA genes in the ant genomes could spuriously be identified as of bacterial origin. We thus annotated  bacterial and eukaryotic rRNAs with `barrnap` and identified overlapping sliding windows with `bedtools intersect`. We then used `infoseq  -nocolumn -delimiter "\t" -auto -only -name -length -pgc` to calculate GC content of each sliding window and of complete scaffolds of the target genome. We calculated coverage from long or short read genomic data across each sliding window using `minimap2`, `samtools` and `bedtools coverage`. This coverage data was used as further lines of evidence to identify  contaminants (where coverage often deviates from the mean coverage across the rest of the genome) as well as high-quality LGTs (where LGT boundaries should be well-covered by genomic reads). 

Scaffolds were flagged as bacterial contaminants if one of the following was true: (1) all sliding windows of a scaffold have better blast hits against the bacterial than the eukaryotic database (based on bitscores), (2) >50% of sliding windows have a better prokaryotic hit and the difference between total bacterial and insect blastn bitscores is >200, (3) >50% of sliding windows have a better prokaryotic hit, the difference between total bacterial and insect blastn bitscores is >0, the difference between total bacterial and insect blastx bitscores is >50, the scaffold-wide GC content and sequencing coverage do not lie within 95% confidence interval of the respective distribution of other eukaryotic scaffolds of the assembly. 
Scaffolds flagged as bacterial were removed from the assembly before further processing.

To identify putative horizontally transferred genes ("lateral gene transfers", LGTs), we identified all individual sliding windows with a HSP against the bacterial database that did not overlap with a higher-scoring HSP against the insect database. We then filtered out any HSP, where the bitscore difference bacterial and eukaryotic hits were <50 and where the bacterial HSP was shorter than 50 bp (cutoffs were 50 for the `loose` LGTs and 100 for the not-`loose`). Remaining candidate LGT HSPs less than 500 bp apart were merged into LGT loci for further analyses. 

# GAGA-Metagenome-LGT

This repository hosts a series of scripts to screen GAGA (Global Ant Genomics Alliance) assembled genomes for bacterial contigs/contaminants and LGT (lateral gene transfers) candidates based on blast results (mmseqs blastn and blastx).

## Before running

The only tested working environment is on the HPC (High Performance Computing) clusters Computerome (computerome.dk) which is based on TORQUE queuing system (e.g. qsub) and Moab Workload Manager (e.g. msub). The dependencies (often requeired) will need to be loaded (e.g. module load; modules.sourceforge.net) within the HPC environment.

## The pipeline

In order to obtain the desired results, scripts need work as the following sequence.

1. Run batch submission script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh)
   * To run this script in the TORQUE queuing system with  necessary information:

   ```bash
   qsub -A <account_name> -l nodes=<node_count>:ppn=<core_count>:thinnode,mem=<mem_amount>,walltime=<time_requested> -N <Job Name> -e <your_run>.err -o <your_run>.log -v "id=<GAGA_id>,pacbio=<1 or ...>" runLGT.v4.sh
   ```

   * If want to receive email warning only when error occurs, add an option:

   ```bash
    -m n
   ```

   The script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh) does the following for a given genome assembly:
   (A) blast sliding windows against different genome databases.  
   * generate 2.5 kb windows with 500 bp walking steps for all scaffolds
   * search each window against curated `mmseqs` genome databases, a *bacterial db*, an *insects db* with selected insect genome assemblies, an *no-ants insect db* with selected insect genome assemblies excluding ants, and a *human genome db*
   * keep the best hit for each window, sorted first by evalue (-k7,7g), then by bitscore (-k8,8gr)
   (B) Annotate rRNAs (with `barrnap`) in the assembly (to filter false-positive hits against the bacterial database)
   (C) Calculate gc-content and length (relevant only for windows at the end of a scaffold) for each window
   (D) Analyze sequencing raw read coverage (usually PacBio) for each window 
   
   </br>

2. The script invokes a separate R script ([`GAGA_metagenome_pipeline.R`](GAGA_metagenome_pipeline.R)) that processes the raw output (blast results, etc.) and identifies scaffolds to be flagged as bacterial or human contaminations.

   ```bash
   Rscript /home/people/dinghe/github/GAGA-Metagenome-LGT/GAGA_metagenome_pipeline.R ${id}
   ```

   Scaffolds are flagged according to the following criteria:

   ```
   # classify a scaffold as prokaryotic if
   #   1. 100 % of hitted windows are identified as prokaryotic OR
   #   2. more than 50% of windows are identified as prokaryotic AND bsratio > 200
   #   3. more than 50% of windows are identified as prokaryotic AND 0 < bsratio < 200 AND their
   #     3.1 bsratio.x > 50 AND
   #     3.2 GC content does not lie within 95% confidence interval of euk scaffolds GC content distribution AND
   #     3.3 coverage does not lie within 95% confidence interval of euk scaffolds coverage distribution
   # note that it is possible to have ratio 1 for only small amount of proWindow but the rest are uncertain (no hit)
   # for downstream genome annotation, recommending only the scaffolds tagged "pro" should be excluded
   # for downstram microbiome analyses, recommending to include all "pro" and "uncertain" scaffolds
   # Checking for human contamination has added: kingdom will be assigned human if more than half of
   # the windows hitted better against human DB
   ```

3. To identify regions in the genomes that have the signature of being horizontally transferred from a bacterium, the script [`runLGT.v4.sh`](qsub_batch/runLGT.v4.sh) does the following:
   ```
   # -- Intersect euk and pro blast hits -- #
   # get overlapping HSPs between pro (windows.fa.DB.pro.bed) and euk (windows.fa.DB.euk.bed).
   # for each pro hit, all overlapping euk hits are returned (bedtools intersect -a windows.fa.DB.pro.bed -b windows.fa.DB.euk.bed -wao)
   # sort overlap by euk bitscore (sort -t$'\t' -k 1,4 -k 10,10gr)
   # keep only the best (-u) scoring euk HSP (as sorted before) (-k10,10gr) for each pro HSP ( sort -t$'\t' -u -k1,4)
   ```
4. The R-script `analyseLGTs.cr2.Rmd` processes files produced in step 3 and visualizes LGT candidate regions and produces one overview plot per processed genome. 

   The following info is loaded and processed:

   ```
   # genome wide info
   genome.file <- genome.file                           # basic bedtools tsv genome file with the length of each scaffold.
   globalCoverage <- genome.overlappingwindows.cov.tsv  # read coverage along sliding-windows

   # genome-wide besthits
   ebh <- query.GAGAid.DB.euk_blastn.bh                 # best blast hits of sliding-window against the eukaryotic database
   pbh <- query.GAGAid.DB.pro_blastn.bh                 # best blast hits of sliding-window against the prokaryotic database
   nAbh <- query.GAGAid.DB.euk_noAnts_blastn.bh         # best blast hits of sliding-window against the noAnts database

   # lgt candidate info
   lgt.candidates <- LGTs.candidateloci.loose.bed          # bed file with LGT candidates (loose filtering thresholds)  
   lgt.coverage <- LGTs.candidateloci.loose.coverage.bed   # coverage along LGT candidates
   lgt.proteins <- LGTs.candidateloci.loose.proteins.bed   # blastX hits intersecting with LGT candidates
   lgt.candidate.fasta <- LGTs.candidateloci.loose.fa      # LGT candidate fasta sequence
   lgt.complexity <- LGTs.candidateloci.loose.complex      # LGT candidate sequence complexity scores 
   lgt.candidate.bam <- LGTs.candidateloci.loose.PacBio.bam   # PacBio or stLFR read mapping along LGT canidates
   ```
   The script eventually creates a single overview plot per genome that summarizes all identified LGT candidates, such as [this one](data/GAGA-0515/GAGA-0515.euk.lgt.candidates.pdf). In addition, a single detailed plot for every LGT candidate is created, e.g. [this one (, a true LGT)](data/GAGA-0515/GAGA-0515.Scaffold8-1760778-1761241.pdf) or [this one (, which is a misassembly)](data/GAGA-0515/GAGA-0515.Scaffold96-1-22253.pdf). These plots, along with tsv files for all candidates (e.g. [GAGA-0515.euk.lgt.all.candidates.tsv](data/GAGA-0515/GAGA-0515.euk.lgt.all.candidates.tsv)) and for the most promising candidates (e.g. [GAGA-0515.euk.lgt.good.candidates.tsv](data/GAGA-0515/GAGA-0515.euk.lgt.good.candidates.tsv)) can be used for manual and further automatic screening of LGT candidates. 


# Manual assessment of LGT candidate quality
# First steps to do for each promising LGT

Here's a brief step-by-step guide on how we analyze LGT candidates.

Here, we discuss analyzing the LGT candidate on `Scaffold8` from position `1760778-1761241` in the GAGA genome of *C. obscurior*, **GAGA-0515**.
## 1. Is the LGT candidate worth studying in detail?
We first check, how convincing the evidence is that the LGT candidate is a proper piece of bacterial DNA. A lot of the candidates we identified are in fact highly repetitive regions in the ant genome, that just by coincidence align better to a bacterial genome in our database.

For each GAGA genome we screened, we have generated a plot (`GAGA-0515.euk.lgt.candidates.pdf` available [here](data/GAGA-0515/GAGA-0515.euk.lgt.candidates.pdf)), showing relevant info for all identified candidates. The candidates labelled in the PCA  (here the two candidates "`Scaffold8:1760778-1761241`" and "`Scaffold96:1-22253`") are automatically labelled, if they look like promising candidates based on the following criteria:

|  Metric |  condition  |
|-|-|
|GC-content| 10-90 %|
|*ct6* complexity*| > 0.001 |
|Entropy |> 1.5 |
| `bitscore prokaryotic database` - `bitscore eukaryotic database` | > 100  |
|C+G skew   | 0-1   |

(*: ct6 is "*Trifnov's complexity with order 6*", see https://en.wikipedia.org/wiki/Linguistic_sequence_complexity)

However, this automatic characterization of "good candidates" is not perfect, so it is important to have a quick look at the different values, and maybe even the sequence of the LGT.

The four plots in the bottom of [the Figure](data/GAGA-0515/GAGA-0515.euk.lgt.candidates.pdf) show the following:

- **top left:** How much better are the blast hits against the prokaryotic database than against the eukaryotic databse? We want high values here for a good LGT.
- **top right:**: This just shows you the (logartihmically scaled) length of the different candidates. This is not really relevant for assessing the quality of an LGT.
- **bottom left:** This gives the GC content of the different candidates. You see that a lot of the candidates here have a GC content close to 1 (all the green-ish candidates for example). These are very likely bad candidates, as they only contain repetitive DNA that can occur in both eukaryotic and prokaryotic genomes.
- **bottom right:** This shows the entropy. Proper genetic sequence is usually highly entropic, so we expect high entropy values for good candidates.

So, good candidates will have high values in the top left and bottom right plots and intermediate values in the bottom left plot. The candidate `Scaffold8:1760778-1761241` that we focus on here fits this perfectly!

When you look at the actual sequence of the **good candidate**, you can see no clear pattern and no clear bias towards anything. This is a good sign!
```
>Scaffold8:1760778-1761241
GGCCATCGGCGGTACAAACGCAACGCCCGTGTCCCAGGGTAGATCAATCCAGGTGTCCTGCGGAACGTCGATCACGTAATCGTTGACCAGTGGACGGCCGGCAGGTTTGGCGAATAT
GGTAATAAAGTGGGCTTTGGGATACAAGTCGCGGATCGCGGCGGCGGTAACGCCGGTGTCCACTAGGTCGTCCACCACGATGAAGCCTTCGCCGTCCCCCTCCGCTCGCTTAATGAC
GCTTAGGTCACGTTGAACGTCGTGGTCATAACTCGAGATGCAGACGGTGTCCACGTGGCGGATGTCCAGCTCACGCGCCAGTATCGCGGCCGGTACCAGGCCGCCGCGACTGACCGC
GATAATTCCTTTCCACCGTTCGGTAGATGATAGGAGGCGTTTAGCCAGCGTGCGCGTATGAGTTTGCAACATTTCCCAAGTGACGATATATTTTTTGCTCGACAATTCAGTC
```
In contrast, a **bad candidate** usually looks very suspicious and is easy to spot. For example `Scaffold10:95463-95666`:

```
>Scaffold10:95463-95666
GGGTGGGGGTGGGGGTGGGAGGGGGGGGGTTGGGGGGGGGGGGGTGGTTGGGGTGGGGGGGGGTGGGGGGGGGGGGTGGGGGGGGGGGGTGGGGTGGTGGGGGGGGGGGTGGGGTGG
GGGGGGGGGGGTGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGATTGGTGGGGGGGTGGGTGGGGGGGGG
```

The four plots in the detailed figures [here for `Scaffold8:1760778-1761241`](data/GAGA-0515/GAGA-0515.Scaffold8-1760778-1761241.pdf) and [here for `Scaffold96:1-22253`](data/GAGA-0515/GAGA-0515.Scaffold96-1-22253.pdf) show the following info:


## 2. A first detailed look at the good candidates

Once you have selected a good candidate, you can take a closer look. I generated figures for each candidate in each genome. You can find these figures in the `lgt.candidates` folders.

Take a look at [the figure `GAGA-0515.Scaffold8-1760778-1761241.pdf`](data/GAGA-0515/GAGA-0515.Scaffold8-1760778-1761241.pdf).


**At the very top**, you see the full length of the scaffold (here `Scaffold8`) and in red with a black circle the position of the LGT.

Then on the **top left**, you see a closer look at the candidate region and surroundings. The blue boxes show you where we got significant hits against a database containing prokaryotic proteins. So, here you see that 9 prokaryotic proteins in our database have a significant similarity to this one particular region. **In the second plot on the left**, you see the (log2-scaled) relative coverage (from long read PacBio data) for the region. A value of 0 here means that the coverage at a given position is just as the average of the genome. A value of -1 would mean half as much coverage, a value of 1 means twice as much coverage, etc. Good LGTs should not deviate very strongly from 0 with regard to the coverage, because values > 1 and < -1 might hint that the genome assembly was bad at this particular region.

**The third panel on the left** again gives an overview of the blast-scores against a prokaryotic (orange) and eukaryotic (blue) database. A good LGT candidate has a pattern like we see here. In a background of high blast hits against the eukaryotic database (the constantly high blue line), we see a sudden spike of a prokaryotic hit (the orange line). The small red box in this plot is in fact the precise range of the LGT candidate (i.e. position 1760778-1761241).

**The last plot on the left side** shows the long-read PacBio data that were mapped to this region. This can be used to differentiate between a misassembly and a good LGT candidate. Misassemblies should have no or very few reads overlapping the boundaries of the LGT candidate. A good LGT candidate will have many reads of the type `antDNA--LGT--antDNA`.

**The right part of the figure** shows at the top a focused view of the PacBio reads. Below that, you see an overview of the best blast hits. The blue one with the tag `A4W6X0` is the best protein hit against the `SwissProt` protein database. You can see what this is here: https://www.uniprot.org/uniprot/A4W6X0

The red box is the hit against our database of bacterial genomes. According to our blastn searches, the LGT candidate has a very high similarity (e-value 5.50e-93) against a piece of the genome of *Cedecea neteri*. The NCBI accession of that is `(CP009451.1)`: https://www.ncbi.nlm.nih.gov/nuccore/CP009451.1/

**On the very bottom of the right side**,  you see the first 75 bp of the LGT candidate.
