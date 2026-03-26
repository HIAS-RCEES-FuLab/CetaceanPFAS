This repository provides R-based tools for the Non-Targeted Analysis of PFAS.
It integrates two primary screening strategies: Homologue Series Analysis and Diagnostic Fragment Filtering.
Homologue Series Analysis module identifies PFAS clusters based on repeating structural units, such as the -CF2- unit.
Diagnostic Fragment Filtering module identifies PFAS candidates by tracking characteristic Fragment ions. It extracts EICs for PFAS diagnostic fragments (e.g., C3F7), and calculates the Pearson Correlation Coefficient or Cosine Similarity between the peak shapes of fragments and potential precursor ions. Candidates with a high correlation score are flagged as origin-fragment pairs, helping to reconstruct the parent molecule structure even when the molecular ion is labile.
