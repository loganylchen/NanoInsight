# Snakemake workflow: `NanoInsight`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/loganylchen/NanoInsight/workflows/Tests/badge.svg?branch=main)](https://github.com/loganylchen/NanoInsight/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `Comprehensive Nanopore DRS Data Analysis Workflow`

## Workflow Structure

### 1. Fusion Transcript Analysis
NanoInsight employs multiple fusion detection algorithms optimized for long-read sequencing data:

- **Jaffal**  
  A long-read-specific fusion detection tool designed to identify chimeric transcripts from Nanopore data. Jaffal leverages both alignment and assembly-based approaches to maximize sensitivity while minimizing false positives.

- **Genion**  
  An alignment-free method that detects fusion events by directly analyzing k-mer signatures in long reads. Genion excels at identifying novel fusion junctions with complex structural rearrangements.

- **LongGF**  
  A graph-based fusion detector that constructs splicing graphs from long reads to identify both canonical and non-canonical fusion transcripts. LongGF is particularly effective for detecting fusion events involving alternative splicing.

- **Aeron** (Currently Disabled)  
  Although initially integrated, Aeron has been temporarily disabled due to compatibility issues. Future updates will re-evaluate its integration based on tool maintenance and performance improvements.


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=loganylchen%2FNanoInsight).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) NanoInsightsitory and its DOI (see above).

# TODO

* Replace `loganylchen` and `NanoInsight` everywhere in the template (also under .github/workflows) with the correct `NanoInsight` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `NanoInsight`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `loganylchen` and `NanoInsight` were correctly set.