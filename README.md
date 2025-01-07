# TCRconvertR

> **Warning**: This project is in **alpha stage**. It is under active development and may be unstable.

**Rename T-cell receptor genes between 10X, Adaptive, and IMGT formats**

TCRconvertR takes T-cell receptor (TCR) data containing V, D, J, and/or C genes from 10X, Adaptive, or other sequencing platforms and renames them from any of these formats to any other one:

* **10X**: TRAV1-2
* **Adaptive**: TCRAV01-02*01
* **IMGT**: TRAV1-2*01

TCRconvertR will work with human, mouse, and rhesus macaque data out-of-the-box, but users can also add their own species.

TCRconverRt helps researchers unify TCR datasets by converting them to a standard naming convention. It is fast, reliable, and prevents errors from manual conversions. Unlike other tools that require custom objects, TCRconvertR works directly with data frames and CSV/TSV files.
