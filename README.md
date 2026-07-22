# BIP — Barcode Inference Pipeline

BIP is a containerized DNA barcoding pipeline for use with any genetic marker, primer set, and second-generation or third-generation sequencing platform, with extra features for the COI barcode region (COI-5P). BIP intakes raw sequence data and performs the following: paired-end merging (Illumina only), demultiplexing, primer/UMI trimming, per‑sample OTU clustering, reference‑frame sequence correction (COI-5P only), chimera detection and removal, probabilistic taxonomic assignment, barcode inference, and BIN matching (COI-5P only). All tools (R/Python packages, reference databases) are included in the container. Supports **Illumina paired‑end**, **PacBio**, **Oxford Nanopore (ONT)**, and generic single‑end platforms.

---

## Contents
- [What BIP does](#what-bip-does)
- [How BIP works](#how-bip-works)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Run on your own data](#run-on-your-own-data)
- [Parameters](#parameters)
- [Reference library](#reference-library)
- [Outputs](#outputs)
- [Build the image yourself](#build-the-image-yourself)
- [Interactive container](#interactive-container)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## What BIP does

BIP must be provided with a non-demultiplexed `.fastq.gz` file (or `R1.fastq.gz` + `R2.fastq.gz` if Illumina paired-end data) and a `parameters.xlsx` spreadsheet (describing wells, UMIs, primer-specific parameters, and expected taxonomy) and produces per-well OTU table containing consensus sequences with taxonomy, target barcode versus non-target sequence status, various sequence quality metrics, and BOLD BIN assignments (COI-5P only). A SINTAX-formatted reference library (`.fasta`) is also required — a COI BOLDistilled library is bundled into the Docker image at build time, so it's already available with no extra setup, but you can supply your own for any genetic marker instead (see [Reference library](#reference-library)).

---

## How BIP works

Stages, in order (these are the exact stage names written to `duration_log.txt`):

1. **Curating input data** — validate that at least one `.fastq.gz`/`.fastq` and exactly one `.xlsx` are present in the working directory.
2. **FASTQ file merge** — for Illumina paired-end, merge R1/R2 with `PEAR`; otherwise concatenate all `.fastq.gz`/`.fastq` into a single file.
3. **Primer recognition & read orientation** — orient and split reads by primer pair (`cutadapt`).
4. **Size filtering** — length-filter by the `Min/Max amplicon length` set for each primer set.
5. **Demultiplexing** — assign reads to wells by UMI (`cutadapt`), sanitize/rename sample headers.
6. **OTU clustering** — per-well primary + secondary clustering (`vsearch`), consensus building via multiple-sequence alignment (`muscle`/`DECIPHER`/`msa` in R).
7. **Indel correction** (COI-5P only, full-length amplicons) — reference-frame homopolymer/indel correction against a curated reference set.
8. **Chimera removal** — `vsearch --uchime_denovo`.
9. **Tax assignment** — `vsearch --sintax` against the user-specified (or bundled BOLDistilled) reference library.
10. **Barcode inference** — classify OTUs as on-target barcodes vs. non-target sequences using a decision tree.
11. **BIN matching** — (COI-5P only) — match OTUs to BOLD BINs (`vsearch --usearch_global`).
12. **Report generation** — an R script renders a PDF report (`<RunName>_Report.pdf`) and writes `publication_metrics.txt`.

---

## Requirements
- Supported architectures: **linux/amd64** and **linux/arm64** (Apple Silicon).
- Disk space and memory — the reference library alone is several GB, and large runs generate substantial intermediate files.
- [Docker](https://docs.docker.com/get-docker/) (Desktop or Engine).

### Instructions for easy Docker setup

**Linux:**
```bash
curl -fsSL https://get.docker.com | sh      # installs Docker Engine + compose plugin
sudo usermod -aG docker "$USER"             # so you don't need sudo
newgrp docker                               # apply group now (or just log out/in)
docker run hello-world                      # verify, no sudo
```

**Mac/Windows (with WSL):** install Docker Desktop from https://www.docker.com/products/docker-desktop/

---

## Quick start

Run the bundled demo first to confirm your setup works, before touching real data. Copy the three files from the BIP repo's `workdir/` directory (`BAI109_demo.fastq.gz` + `parameters.xlsx` + `compose.yaml`) into a working directory of your own.

**Pull the image:**
```bash
docker pull ghcr.io/cbg-innov/bip:latest
```

**Run the demo:**
```bash
mkdir -p ~/bip_workdir && cd ~/bip_workdir
cp /path/to/workdir/BAI109_demo.fastq.gz /path/to/workdir/parameters.xlsx /path/to/workdir/compose.yaml .

BIP_DATA="$(pwd)" docker compose -f compose.yaml run --rm \
  -e wkdir=/data \
  bip bash /BIP/SCRIPTS/BIP.sh --platform nanopore --minreads 2 DemoRun
```

Results land in `./DemoRun_results/` inside your working directory. 

---

## Run on your own data

Once that's worked, switch to your own data **in the same directory** — replace the files and rerun the identical command. You may also override advanced parameters here with a flag (e.g., `--chimera_mindiv`), for instance:
```bash
rm BAI109_demo.fastq.gz
cp /path/to/your.fastq.gz /path/to/your_parameters.xlsx .

BIP_DATA="$(pwd)" docker compose -f compose.yaml run --rm \
  -e wkdir=/data \
  bip bash /BIP/SCRIPTS/BIP.sh --platform nanopore --minreads 2 <RunName> --chimera_mindiv 0.001
```

No new flags needed for this. BIP just rescans the working directory for at least one `.fastq.gz` (or `.fastq`) and exactly one `.xlsx`. Your `.fastq.gz`, `.xlsx`, and a persistent `dictionary_BIP.tsv` are never moved or archived, so you can rerun without re-copying anything. If it finds no `.fastq.gz`/`.fastq` or more than one `.xlsx`, it exits with an error naming the problem (this is also why a stray Microsoft Excel lock file, `~$parameters.xlsx`, will break the run if the params file is open elsewhere).

- `--platform` is required: `illumina_pe`, `pacbio`, `nanopore`, or `other`.
- `--minreads N` (default `5`) sets a floor for per-well OTU retention.
- The last argument is the run name you wish to provide, which is essential for naming outputs.
- `wkdir` (working directory), `scripts_directory`, and `reference_library_directory` are environment variables, not flags. Set via `-e` on `docker compose run`, as shown above. Their in-container defaults are `/BIP/Barcoding`, `/BIP/SCRIPTS`, and `/BIP/REFS`.
- For Illumina paired-end data, set `--platform illumina_pe`; BIP merges R1/R2 with `PEAR` before proceeding.


To use your own reference library instead of the bundled one, see [Reference library](#reference-library).

---

## Parameters

The parameters spreadsheet (`.xlsx`) has three tabs: `Instructions`, `UMIs and Primers`, and `Dictionary Update`. `Instructions` explains what each column in the next two tabs means. `UMIs and Primers` must be minimally completed as it provides well-based details; `Dictionary Update` only needs to be completed for primer pairs you're using for the first time or if you're updating primer-specific parameters.

### UMIs and Primers (one row per well)

| Field | Meaning |
|---|---|
| **Plate** | Plate name/number (e.g., `Plate01`) |
| **Well** | Well ID (e.g., `A01`) |
| **Sample** | Sample name for this well |
| **Forward / Reverse UMI** | UMI sequence for this well |
| **Forward / Reverse Primer Name** | Must match a name defined in the persistent `dictionary_BIP.tsv` file or specified in the `Dictionary Update` tab |
| **Negative Control** | `yes`, `y`, `Y`, `1`, etc if this well is a negative control, otherwise leave blank |
| **Kingdom … Species** | Optional expected taxonomy, used to greatly aid barcode inference |

### Dictionary Update (one row per primer pair)

| Field | Meaning |
|---|---|
| **Forward / Reverse Primer Name** | Name referenced in the `UMIs and Primers` tab |
| **Forward / Reverse Primer Sequence** | The actual primer sequence (5' to 3')|
| **Marker** | Locus name (e.g., `COI-5P` for the COI barcode region, anything else for all other markers) |
| **Reference Library** | Prefix of a SINTAX-formatted `.fasta` in `REFS/` (e.g., `BOLDistilled_COI_Apr2026` is fine for referring to `BOLDistilled_COI_Apr2026_SEQUENCES_sintax.fasta`). Must be an unambiguous prefix. If another file in `REFS/` shares it, the run fails with a "multiple files match" error. |
| **Min / Max amplicon length** | Length filter, applied with UMIs and primers still attached |
| **Target amplicon length** | Expected clean amplicon length (no primers/UMIs) |
| **Primary / Secondary OTU clustering threshold** | vsearch clustering identity thresholds (within-well, then across-run) |
| **Tax assign probability threshold** | SINTAX confidence cutoff |
| **BIN assign threshold** | Minimum identity to accept a BOLD BIN match (used for COI-5P only) |

Entries here are merged into a persistent `dictionary_BIP.tsv` in your working directory each run, so previously-used primer pairs accumulate rather than needing to be re-entered every time.

### Advanced tuning

A handful of parameters (UMI/primer overlap minimums, per-base error tolerances, the satellite-OTU read-count floor, core count) are editable as flags. 
| Flag | Default | Description |
|---|---|---|
| `--platform` | *(required)* | `illumina_pe`, `pacbio`, `nanopore`, or `other`. Determines whether R1/R2 are merged with `PEAR` and whether chimera screening runs early (per-sample) or late (post-clustering). |
| `--minreads` | `5` | OTUs below this per-sample read count are discarded during clustering. |
| `--wd` | `/BIP/Barcoding` | Working directory — where BIP looks for the input `.fastq.gz`/`.xlsx` and writes `<RunName>_results/`. |
| `--scripts` | `/BIP/SCRIPTS` | Directory containing `BIP.sh`'s companion `bip1`-`bip4` R/Python scripts. |
| `--refs` | `/BIP/REFS` | Reference library directory - searched for both the SINTAX taxonomy reference and the error-correction reference (`reference_seqs_*.fasta`). |
| `--ref_seq_corr` | *(auto-detect)* | Explicit path to the error-correction reference FASTA. Leave unset to auto-detect the single `reference_seqs_*.fasta` in `--refs`. |
| `--cores_to_leave` | `3` | How many CPU cores to leave free. BIP will use the rest. |
| `--umi_overlap_min` | `0.75` | Multiplier applied to UMI lengths during demultiplexing. Recommended: 0.75 for UMIs >=12 nt, 1.0 for UMIs <12 nt. |
| `--primer_overlap_min` | `0.75` | Multiplier applied to primer lengths during demultiplexing. |
| `--error_umi1` / `--error_umi2` | `0.125` / `0.125` | Max mismatch rate (Cutadapt) for the forward/reverse UMI. Consider lower for UMIs <8 bp. |
| `--error_primer1` / `--error_primer2` | `0.2` / `0.2` | Max mismatch rate (Cutadapt) for the forward/reverse primer. |
| `--minreadssatellite` | `0.00000625` | Controls how "satellite OTUs" are removed, as a fraction of total input reads (default yields a min read count of ~5 on a typical run). |
| `--minoverlap` | `0.75` | Minimum overlap (as a fraction of expected amplicon length) required for a BIN match. |
| `--chimera_abskew` | `10` | `vsearch --uchime_denovo` `abskew` parameter. |
| `--chimera_mindiv` | `0.0005` | `vsearch --uchime_denovo` `mindiv` parameter. |
| `--bin_maxhits` | `1` | `vsearch --usearch_global` `maxhits` parameter (BIN matching). |
| `--bin_maxaccepts` | `1` | `vsearch --usearch_global` `maxaccepts` parameter (BIN matching). |

---

## Reference library

The BOLDistilled COI SINTAX reference set and the error-correction reference set are both downloaded and unpacked into `/BIP/REFS` **at image build time** (not first run). With `docker compose`, `REFS` is mounted as a named volume (`reflib`) so it persists across container recreation.

To use a different reference library, set the `Reference Library` field in your `parameters.xlsx`'s `Dictionary Update` tab to an unambiguous prefix of your reference file's name, then make the file itself available one of two ways:

**From a folder on your local device** (no rebuild needed) — mount it and point `--refs` at it:
```bash
BIP_DATA="$(pwd)" docker compose -f /path/to/workdir/compose.yaml run --rm \
  -e wkdir=/data \
  -v ~/Desktop/REFS:/myrefs \
  bip bash /BIP/SCRIPTS/BIP.sh --platform nanopore --minreads 2 --refs /myrefs <RunName>
```
Since `--refs` replaces the search location entirely rather than adding to it, `~/Desktop/REFS` needs to contain both your SINTAX-formatted `.fasta` *and* the error-correction reference (`reference_seqs_*.fasta`) — or use `--ref_seq_corr /BIP/REFS/reference_seqs_327K.fasta` to keep pulling that one specific file from the built-in volume while everything else comes from your local folder.

**Permanently, in the built-in `reflib` volume** — copy your file directly into `REFS/` before building your own image, so it's baked in alongside the bundled BOLDistilled library.

---

## Outputs

Written to `<wkdir>/<RunName>_results/`:

- **`<RunName>_Report.pdf`** — summary report.
- **`publication_metrics.txt`** — run parameters, software versions, and system info often required for publication.
- **`duration_log.txt`** — per-stage timing.
- **`<PrimerPair>/`** — per-primer-pair OTU details (`<RunName>__<PrimerPair>__OTUDetails.tsv`), consensus FASTA files, and error-correction intermediates.
- **`Miscelleneous_Files/`** — read-count summaries, the merged parameters TSV, and other intermediate files.
- **`Individual_Raw_Fastq_Files.tar.gz`** — per-well demultiplexed reads, often required for SRA archiving for publication.

Your original `.fastq.gz`/`.fastq`, `.xlsx`, `dictionary_BIP.tsv`, and `compose.yaml` (if present) stay at the top level of `wkdir`. They are never moved or archived into the results folder, so nothing is lost if you delete `<RunName>_results/` and rerun.

---

## Build the image yourself

**Locally (single-arch, matches your own machine):**
```bash
cd Docker
docker build -t bip .
```
This installs the environment via `micromamba`, lays out `SCRIPTS/` and `Barcoding/`, and downloads + unpacks the SINTAX reference library and error-correction reference set.

Run a locally-built image by using `bip:latest` instead of `ghcr.io/cbg-innov/bip:latest` in the commands above.

---

## Interactive container
```bash
docker compose -f workdir/compose.yaml up -d      # start a long-running 'bip' container
docker compose exec bip bash                      # drop into a shell (the 'bip' env auto-activates)

# inside the container:
bash /BIP/SCRIPTS/BIP.sh --platform nanopore --minreads 2 <RunName>

docker compose down     # stop & remove the container (the reflib volume persists)
```

To make your own data visible inside the interactive container without `docker compose run`, set `BIP_DATA` before `up -d`, or add a bind mount directly under the `bip` service in `compose.yaml`:
```yaml
    volumes:
      - reflib:/BIP/REFS
      - ~/Desktop/my_run:/data

```

---

## Troubleshooting

- **No output on the host?** Make sure `wkdir=/data` is set (via `-e wkdir=/data`) and that `BIP_DATA` points at a real folder. Otherwise BIP writes to its default in-container `/BIP/Barcoding`, which is lost when the container exits.
- **"Multiple .xlsx files found"?** Check for a Microsoft Excel lock file (`~$<name>.xlsx`) in the same folder and ensure that the file is not open on your device. Close the file in Excel or delete the lock file.
- **"Multiple files match" for the reference library?** Another file in `REFS/` shares the same prefix as your `Reference Library` value (e.g., a metadata PDF). Use a longer, more specific prefix.
- **Disk space.** Large runs can generate many intermediate files; ensure adequate free disk on the Docker host.

---

## Citation

> Sean WJ Prosser, Ken A Thompson, Nicholas W Bard, Robin M Floyd, Emine Ozsahin, and Paul DN Hebert. The Barcode Inference Pipeline (BIP): From sequencer output to DNA barcodes. <i> In prep. </i>
