# 2018_sorghum_nematodes
Single-marker analysis pipeline for Harris-Shultz et al, A Novel QTL for Root-Knot Nematode Resistance is identified from a South African Sweet Sorghum Line

This is a short pipeline to do single-marker analysis for sorghum nematode resistance, BRIX, and root weight. The entire pipeline is run by calling 0_RunSorghumQtl.sh in bash. This pipeline depends on the following software (note: specific version numbers are probably less crucial but given just in case):
* TASSEL version 5.2.43
* Python 3.5.2
     * argparse 1.1
     * gzip
     * math
     * matplotlib 2.2.2
     * numpy 1.14.5
     * pandas 0.23.1
     * re 2.2.1
* R 3.4.4
     * argparse 1.1.1
     * parallel 3.4.4
     * qtl 1.42-8
     
For the publication, this pipeline was run on a desktop workstation (4-core (8 virutal cores) Intel Xeon W-2123 CPU with 64 GB RAM) running Linux Mint 18.3.
