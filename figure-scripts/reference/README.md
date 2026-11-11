# Reference files

`uniprotkb_human_reviewed_2025_10_28.fasta.gz` is the UniProtKB/Swiss-Prot human reviewed proteome
(20,420 entries), downloaded from UniProt on 2025-10-28 as
`uniprotkb_Human_AND_reviewed_true_AND_m_2025_10_28.fasta` (human, reviewed entries). It is the FASTA that the sequence-coverage and
PTM-site analyses (`03_ptm_coverage/`, `04_reproducibility/`) use for residue numbering. It is included
here, gzipped, because that exact snapshot is hard to re-download later. The scripts read it as-is; R
reads `.gz` transparently.

Source: The UniProt Consortium, UniProt: the Universal Protein Knowledgebase in 2025, *Nucleic Acids
Res.* 2025. Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see
<https://www.uniprot.org/help/license>.

The DIA-NN searches themselves used the FASTA files deposited with the MassIVE dataset (MSV000103457).
