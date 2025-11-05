
<img width="8803" height="2602" alt="TargetMosaic_3" src="https://github.com/user-attachments/assets/1c5bdfb6-6ffc-4450-8344-3970cc70cdf2" />

**TargetMosaic** is a workflow for downloading bacterial whole-genome annotations, identifying target proteins from a reference FASTA, and analyzing their sequence coverage across the pathogen species of interest.  
This helps evaluate the potential of a protein as a **drug target** or **binder target**, facilitating scientific contributions to precision medicine and targeted diagnostics.  

---

## Features  
- Parallel download of genome annotation files from NCBI using SLURM.  
- Automatic flattening of download directories.  
- Target gene search against large genome collections.  
- Strain coverage statistics (absent, present, duplicated).  
- Sequence deduplication and multiple sequence alignment.  
- Domain-level analysis of conservation and loss.  
- Final coverage summary (fully targetable, partially targetable, not targetable).  

---

## Installation  

Make sure you have access to a SLURM cluster and the following software available:  

- `anaconda3`  
- `bakrep` (for downloading annotation files)  
- `rsync`  
- `clustalo` (Clustal Omega for alignments)  
- Python ≥ 3.10  

In your working directory:
```bash
git clone https://github.com/arpita-sahoo/TargetMosaic.git
```
## Environment setup  

This repository requires a conda environment named `myenv`.  
You can create it directly from the provided `environment.yml`:  

```bash
module load anaconda3   # if on a cluster with modules
conda env create -f environment.yml
conda activate myenv
```
If the environment already exists, update it with:
```bash
conda env update -f environment.yml --prune
```
### Alternative setup with pip  

If you prefer a lightweight installation, a `requirements.txt` is also provided:  

```bash
conda create -n myenv python=3.10
conda activate myenv
pip install -r requirements.txt
```
---

## Workflow  

TargetMosaic can be applied to **any target gene** in **any bacterial pathogen**.  

1. **Download annotations**  
   - Retrieve protein annotation files (`.faa`) for all strains of interest.  

2. **Search for target gene**  
   - Provide a FASTA sequence of the target gene (e.g., *spa* in *Staphylococcus aureus*).  
   - Identify homologous sequences across all genomes.  

3. **Count strain coverage**  
   - Report how many strains lack the gene, have a single copy, or have multiple copies.  

4. **Extract unique variants**  
   - Collapse redundant sequences into unique variants.  
   - Headers are annotated as:  
     ```
     >Var<variant_number>_<occurrence_count>
     ```

5. **Multiple sequence alignment**  
   - Perform alignment with Clustal Omega.  

6. **Domain analysis**  
   - Extract and analyze functional domains from aligned sequences.  
   - Assess whether domain loss affects druggability.  

7. **Coverage summary**  
   - Report the percentage of strains fully targetable, partially targetable, or not targetable.  

---

## Step 1. Download genome annotation files  

### Input files (in current working directory)  
- `SAMID.tsv` – TSV file with all SAM IDs (one per line)   

This Workflow uses the bakrep-cli tool; hence, the SAM IDs must be compatible with the bakrep database. To generate the TSV file:
- Visit https://bakrep.computational.bio/browse
- Enter appropriate filters. Important:
     - Genome size: according to the species of interest
     - Contigs: 1-200
     - Completeness: 90-100%
     - Contamination: 1-10%
     - Species: full scientific name of the species of interest
- Click `Export as TSV` and upload the exported metadata TSV file(bakrep-export.tsv) into your working directory
- Extract the SAM IDs using:
```bash
module load anaconda3
source activate myenv
python TargetMosaic/extract_SAMid.py Pseudomonas aeruginosa
```
Optionally with custom input/output:
```bash
module load anaconda3
source activate myenv
python TargetMosaic/extract_SAMid.py Escherichia coli -i bakrep-export.tsv -o SAMID.tsv
```

### Run download job  

```bash
mkdir batch_download_faa
cd batch_download_faa
sbatch ../TargetMosaic/download_single.slurm
```
## Step 2. Flatten folder hierarchy  

After downloads finish, flatten the directory structure.  

```bash
cd faa_files
sbatch ../../TargetMosaic/flatten_copy.slurm
```
## Step 3. Target gene identification and coverage analysis  

With all annotation files downloaded and flattened, the next step is to search for your **target gene** across strains. 
We can come back to the main working directory:
```bash
cd ../..
```

### Input files required in the main working directory  
- `WT.fasta` – FASTA file containing the reference sequence of your target gene.  
- `faa_files/` – directory with all `.faa` annotation files from Step 2.  

### Run sequence search  

```bash
sbatch TargetMosaic/extract_kmer_counter.slurm
```
This will produce:
- `matches.fasta` – all sequences homologous to the target gene.
- `counts.tsv` – table with per-strain presence/absence/duplication counts.

Summarize species coverage
# Count number of strains without the target gene
```bash
awk -F'\t' '$2 == 0 {count++} END {print count}' counts.tsv
```
This provides the number of strains without the target.  
`$2 == 0` can be changed to `$2 == 1` to identify number of strains with one copy of the gene, and so on.

---

### Step 3.1 Extract unique variants  

Collapse identical sequences into unique variants, annotated with their occurrence counts:  

```bash
module load anaconda3
source activate myenv
python TargetMosaic/unique_seq.py matches.fasta
```

Output:  
- `unique_seq.fasta` – FASTA with headers in the form:  
  ```
  >Var<variant_number>_<occurrence_count>
  ```  
- Report of total unique variants.  

---

### Step 3.2 Multiple sequence alignment  

Run multiple sequence alignments on unique variants:  

```bash
sbatch TargetMosaic/run_clustalo.sbatch
```

Output:  
- `unique_seq_aligned.fasta` – aligned FASTA file of all variants.  

---

### Step 3.3 Domain analysis  

Extract and analyze specific functional domains from the aligned sequences:  

```bash
module load anaconda3
source activate myenv
python TargetMosaic/domain_extraction.py --aligned unique_seq_aligned.fasta --domains domains.fasta --reference Var1_x
```
Replace `Var1_x` with whichever variant header corresponds to the reference sequence
Replace `domains.fasta` with a fasta file containing sequence(s) of the domain(s) of interest


Output:  
- Domain-specific FASTA file  
- Summary of unique domain sequence variants with corresponding counts.  

---

### Step 3.4 Coverage summary  

Combine gene presence/absence and domain analysis to classify strains:  
- **Fully targetable** – gene present, with intact functional domains.  
- **Partially targetable** – gene present, but altered domains without loss of function.  
- **Not targetable** – gene absent, or gene present with missing non-redundant domains.  

---

## Example: *spa* (Protein A) in *Staphylococcus aureus*  

We applied TargetMosaic to **SpA**, a virulence factor of *S. aureus*:  

- Total strains analyzed: **102,510**  
- Strains lacking *spa*: ~**10%**  
- Unique sequence variants: **8,440**  
- Domain analysis showed:  
  - ~**78%** of strains retain full target capacity  
  - ~**12%** retain partial capacity  
  - ~**10%** cannot be targeted  

This demonstrates how TargetMosaic can guide **drug target prioritization** and species **diagnostic coverage** analysis.  

---

## Citation  

If you use **TargetMosaic** in your work, please cite this repository.  

---

## License  

This project is licensed under the MIT License.  
