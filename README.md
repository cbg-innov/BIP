# BIP — Barcode Identification Pipeline

BIP is a containerized DNA barcoding pipeline for use with all barcode sequences, with particular optimisation for the COI barcode. The pipeline covers the following: demultiplexing, primer/UMI trimming, per‑sample OTU clustering, reference‑frame sequence correction, chimera screening, SINTAX taxonomic assignment, and BOLD BIN matching. All tools (R/Python packages, reference databases) are included in the container. Supports **Illumina paired‑end**, **PacBio**, **Oxford Nanopore (ONT)**, and generic single‑end platforms.

---

## Contents
- [What BIP does](#what-bip-does)
- [How BIP works](#how-bip-works)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Run on your own data](#run-on-your-own-data)
- [Interactive container](#interactive-container)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Build the image yourself](#build-the-image-yourself)
- [Reference library](#reference-library)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## What BIP does

BIP can be provided with a non-demultiplexed `.fastq.gz` file and a `parameters.xlsx` spreadsheet (describing wells, UMIs, and primers) and produces per-well OTU consensus sequences with taxonomy and BOLD BIN assignments. 

---

## How BIP works

Stages, in order (these are the exact stage names written to `duration_log.txt`):

1. **Curating input data** — validate that exactly one `.fastq.gz`/`.fastq` and one `.xlsx` are present in the working directory.
2. **FASTQ file merge** — for Illumina paired-end, merge R1/R2 with `PEAR`; otherwise concatenate.
3. **Primer recognition & read orientation** — orient and split reads by primer pair (`cutadapt`).
4. **Size filtering** — length-filter by the `Min/Max amplicon length` set in the dictionary.
5. **Demultiplexing** — assign reads to wells by UMI (`cutadapt`), sanitize/rename sample headers.
6. **OTU clustering** — per-well primary + secondary clustering (`vsearch`), consensus building via multiple-sequence alignment (`muscle`/`DECIPHER`/`msa` in R).
7. **Indel correction** (COI only, full-length amplicons) — reference-frame homopolymer/indel correction against a curated reference set.
8. **Chimera removal** — `vsearch --uchime_denovo`.
9. **Tax assignment** — SINTAX against the bundled BOLDistilled reference library.
10. **Barcode inference** — classify OTUs as on-target barcodes vs. satellite/off-target, using a read-count floor scaled to run depth.
11. **BIN matching** — match OTUs to BOLD BINs (`vsearch --usearch_global`) for COI markers.
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

Change directory to a working directory containing your `.fastq.gz`, parameters `.xlsx`, and `compose.yaml`. Start by downloading the three files small demo dataset included in the BIP repo's `workdir/` directory (`BAI109_demo.fastq.gz` + `parameters.xlsx` + `compose.yaml`), before using real data. Navigate to a new working directory with the files copied and edited as needed, but replace the `fastq.gz` with your own data. Now you are ready to run BIP!

**Pull the image:**
```bash
docker pull ghcr.io/cbg-innov/bip:latest
```

**Set up a working directory and run:**
```bash
mkdir -p ~/bip_workdir && cd ~/bip_workdir
cp /path/to/your.fastq.gz /path/to/parameters.xlsx .

BIP_DATA="$(pwd)" docker compose -f /path/to/Docker/compose.yaml run --rm \
  -e wkdir=/data \
  bip bash /BIP/SCRIPTS/BIP.sh --platform nanopore --minreads 2 <RunName>
```

Results land in `./<RunName>_results/` inside your working directory. Your original `.fastq.gz`, `.xlsx`, and a persistent `dictionary_BIP.tsv` are left in place at the top level — they are not moved or archived, so you can rerun without re-copying anything.

---

## Run on your own data

Same command as above — there are no `--fastq`/`--params` flags. BIP finds its inputs by scanning the working directory for exactly one `.fastq.gz` (or `.fastq`) and exactly one `.xlsx`. If it finds zero or more than one of either, it exits with an error naming the problem (this is also why a stray Microsoft Excel lock file, `~$parameters.xlsx`, will break the run if the params file is open elsewhere — close it first).

- `--platform` is required: `illumina_pe`, `pacbio`, `nanopore`, or `other`.
- `--minreads N` (default `5`) sets a floor for per-well OTU retention.
- The last argument is the run name you wish to provide, which is essential for naming outputs.
- `wkdir` (working directory), `scripts_directory`, and `reference_library_directory` are environment variables, not flags — set via `-e` on `docker compose run`, as shown above. Their in-container defaults are `/BIP/Barcoding`, `/BIP/SCRIPTS`, and `/BIP/REFS`.

For Illumina paired-end data, set `--platform illumina_pe`; BIP merges R1/R2 with `PEAR` before proceeding.

---

## Interactive container

```bash
docker compose -f Docker/compose.yaml up -d      # start a long-running 'bip' container
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

## Parameters

The parameters spreadsheet (`.xlsx`) has three tabs: `Instructions`, `UMIs and Primers`, and `Dictionary Update`. The first two must be filled out; only add rows to `Dictionary Update` for primer pairs you're actually using.

### UMIs and Primers (one row per well)

| Field | Meaning |
|---|---|
| **Plate** | Plate name/number (e.g., `Plate01`) |
| **Well** | Well ID (e.g., `A01`) |
| **Sample** | Sample name for this well |
| **Forward / Reverse UMI** | UMI sequence for this well |
| **Forward / Reverse Primer Name** | Must match a name defined in the `Dictionary Update` tab |
| **Negative Control** | `yes` if this well is a negative control, otherwise leave blank |
| **Kingdom … Species** | Optional expected taxonomy, used for reporting only |

### Dictionary Update (one row per primer pair)

| Field | Meaning |
|---|---|
| **Forward / Reverse Primer Name** | Name referenced from the `UMIs and Primers` tab |
| **Forward / Reverse Primer Sequence** | The actual primer sequence |
| **Marker** | Locus name (e.g., `COI-5P`) |
| **Reference Library** | Prefix of a SINTAX-formatted `.fasta` in `REFS/` (e.g., `BOLDistilled_COI_Apr2026` matching `BOLDistilled_COI_Apr2026_SEQUENCES_sintax.fasta`). Must be an unambiguous prefix — if another file in `REFS/` shares it, the run fails with a "multiple files match" error. |
| **Min / Max amplicon length** | Length filter, applied with UMIs and primers still attached |
| **Target amplicon length** | Expected clean amplicon length (no primers/UMIs) |
| **Primary / Secondary OTU clustering threshold** | vsearch clustering identity thresholds (within-well, then across-run) |
| **Tax assign probability threshold** | SINTAX confidence cutoff |
| **BIN assign threshold** | Minimum identity to accept a BOLD BIN match |

Entries here are merged into a persistent `dictionary_BIP.tsv` in your working directory each run, so previously-used primer pairs accumulate rather than needing to be re-entered every time.

### Advanced tuning

A handful of parameters (UMI/primer overlap minimums, per-base error tolerances, the satellite-OTU read-count floor, core count) are set as constants near the top of `BIP.sh` rather than exposed as flags. Edit them directly in the script if you need to change them.

---

## Outputs

Written to `<wkdir>/<RunName>_results/`:

- **`<RunName>_Report.pdf`** — summary report.
- **`publication_metrics.txt`** — run parameters, software versions, and system info.
- **`duration_log.txt`** — per-stage timing.
- **`<PrimerPair>/`** — per-primer-pair OTU details (`<RunName>__<PrimerPair>__OTUDetails.tsv`), consensus FASTA files, and error-correction intermediates.
- **`Miscelleneous_Files/`** — read-count summaries, the merged parameters TSV, and other intermediate files.
- **`Individual_Raw_Fastq_Files.tar.gz`** — per-well demultiplexed reads.

Your original `.fastq.gz`, `.xlsx`, `dictionary_BIP.tsv`, and `compose.yaml` (if present) stay at the top level of `wkdir` rather than being archived into the results folder.

---

## Build the image yourself

**Locally (single-arch, matches your own machine):**
```bash
cd Docker
docker build -t bip .
```
This installs the environment via `micromamba`, lays out `SCRIPTS/` and `Barcoding/`, and downloads + unpacks the SINTAX reference library and error-correction reference set. `

Run a locally-built image by using `bip:latest` instead of `ghcr.io/cbg-innov/bip:latest` in the commands above.

---

## Reference library

The BOLDistilled COI SINTAX reference set and the error-correction reference set are both downloaded and unpacked into `/BIP/REFS` **at image build time** (not first run), so no internet connection is needed at run time beyond pulling the image. With `docker compose`, `REFS` is mounted as a named volume (`reflib`) so it persists across container recreation.

To use a different reference library, copy your SINTAX-formatted `.fasta` into `REFS/` and set the `Reference Library` field in the `Dictionary Update` tab to an unambiguous prefix of its filename.

---

## Troubleshooting

- **No output on the host?** Make sure `wkdir=/data` is set (via `-e wkdir=/data`) and that `BIP_DATA` points at a real folder — otherwise BIP writes to its default in-container `/BIP/Barcoding`, which is lost when the container exits.
- **"Multiple .xlsx files found"?** Check for a Microsoft Excel lock file (`~$<name>.xlsx`) in the same folder — it matches the same glob as your real parameters file. Close the file in Excel or delete the lock file.
- **"Multiple files match" for the reference library?** Another file in `REFS/` shares the same prefix as your `Reference Library` value (e.g., a metadata PDF). Use a longer, more specific prefix.
- **Disk space.** Large runs can generate many intermediate files; ensure adequate free disk on the Docker host.

---

## Citation

> Sean WJ Prosser, Ken A Thompson,  Nicholas W Bard, Robin M Floyd, Emine Ozsahin, and Paul DN Hebert. The Barcoding Inference Pipeline (BIP): From sequencer output to DNA barcodes. <i> In prep. </i>

BIP_DATA="$(pwd)" docker compose -f compose.yaml   run --rm -e wkdir=/data bip   bash /BIP/SCRIPTS/BIP.sh   --platform nanopore --minreads 2 BAI