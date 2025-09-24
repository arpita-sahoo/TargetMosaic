# TargetMosaic  

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
   - Report percentage of strains fully targetable, partially targetable, or not targetable.  

---

## Step 1. Download genome annotation files  

### Input files (in current working directory)  
- `SAMID.tsv` – TSV file with all SAM IDs (one per line)   

### Run download job  

```bash
module load anaconda3
source activate myenv
sbatch TargetMosaic/download_single.slurm
```
## Step 2. Flatten folder hierarchy  

After downloads finish, flatten the directory structure.  

From inside `batch_downloads/faa_files/`:  

```bash
sbatch flatten_copy.slurm
```
## Step 3. Target gene identification and coverage analysis  

With all annotation files downloaded and flattened, the next step is to search for your **target gene** across strains.  

### Input files  
- `target.fasta` – FASTA file containing the reference sequence(s) of your target gene.  
- `faa_files/` – directory with all `.faa` annotation files from Step 2.  

### Run sequence search  

Example command (adapt to your environment and tools):  

```bash
sbatch extract_kmer_counter.slurm
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
python unique_seq.py matches.fasta
```

Output:  
- `unique_seq.fasta` – FASTA with headers in the form:  
  ```
  >Var<variant_number>_<occurrence_count>
  ```  
- Report of total unique variants.  

---

### Step 3.2 Multiple sequence alignment  

Run multiple sequence alignment on unique variants:  

```bash
sbatch run_clustalo.sbatch
```

Output:  
- `unique_seq_aligned.fasta` – aligned FASTA file of all variants.  

---

### Step 3.3 Domain analysis  

Extract and analyze specific functional domains from the aligned sequences:  

```bash
python domain_extraction.py --aligned unique_seq_aligned.fasta                             --domains domains.fasta                             --reference Var1_8302
```

Output:  
- Domain-specific FASTA files (e.g., `signal_sequences.fasta`, `z_domains.fasta`).  
- Summary of domain conservation across variants.  

---

### Step 3.4 Coverage summary  

Combine gene presence/absence and domain analysis to classify strains:  
- **Fully targetable** – gene present, with intact functional domains.  
- **Partially targetable** – gene present, but missing or altered domains.  
- **Not targetable** – gene absent.  

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

This demonstrates how TargetMosaic can guide **drug target prioritization** and **species coverage analysis**.  

---

## Citation  

If you use **TargetMosaic** in your work, please cite this repository.  

---

## License  

This project is licensed under the MIT License.  
