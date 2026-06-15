"""AMR Gene Cargo in Commercial Probiotic Strains — analysis pipeline.

Implements the seven-phase research plan:
  Phase 2  genome_fetch     download probiotic + pathogen genomes (NCBI)
  Phase 3  amr_screen       detect acquired ARGs in every genome
  Phase 4  mge_analysis     locate plasmid / IS / transposon context of each ARG
  Phase 4  risk_classify    assign the four-tier transfer-risk level
  Phase 5  pathogen_compare build the probiotic<->pathogen gene-sharing network
  Phase 6  figures          heatmap, risk chart, gene context, sharing network
           pipeline         orchestrate everything from a single config file
"""

__version__ = "1.0.0"
