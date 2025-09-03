# fish_eDNA

Gerald fish eDNA - Notes

-most common marker genes used in fish: 12S (MiFish-U, highly specific, 163-185 bp), 16S, cytB, COI (large ref library but off target amplification, 185 bp)
-regionally curated databases perform better (Gold et al., 2021)
-occupancy models need at least two eDNA samples per site
-multi-species occupancy models: mitigate high error rates of meatbarcoding studies
-process-based models
-geographic curation: remove geographically improbable species


Previous study (MSc thesis)
Impact of disturbance of eDNA detection
Compare eDNA to conventional methods
-Experimental design:
--9 sites sampled
--two reps before disturbance, two after
--pre: 1.2 micrometers, post: 5 micrometers filters, target volume 2L
--Fish pos control with 10 marine species for 12S
--Tax with RDP classifier (COI trained on eukaryotic ref seqs from GenBank and BOLD; 12S on vertebrate reference sequences from NCBI nucleotide and MitoFish database)
--Filtered for class Actinopteri A >=90%
--Detection threshold: 20 reads for 12S, 60 reads for COI
--This thesis contains a comprehensive list of Canadian freshwater species!!!

Gerald eDNA dataset
-objective: compare eDNA work to conventional methods (electro fishing)
-questions: 18 samples in raw data (12S: S3-S22, COI: S43-S60), but 22 comercially analyzed
-Not sure if this is the same data, ask Gerald for meta data
-Ask Gerald for traditional tax data

-Previous approach/analyses
--Run dada2 pipeline for 12S and COI dataset
---COI data appears to be bad quality (many short reads, are we sure we have the right primers?): Only few reads survive filtering etc.
--12S: Assign taxonomy with two approaches: 
---rdp classifier (12S fish classifier training set from terrimporter github): many actinopteri, but low higher tax support
---assign tax by blasting against mito tax: questionable, mostly catfish, but low id match
--COI: Assign taxonomy with two approaches
---rdp classifier: Only very few ASVs assigned, tropical hits!
---assign tax by blasting against midori unique database: many tropical hits (short but perfect hits)
---another problem: they don't pass the pseudogene filter
---There is also a midori db directly implemented in dada2 format! Yet to be tried


--Overall project goal: Curate custom databases for all our fish eDNA work

-Next steps
1. Re-run dada2 pipeline to familiarize yourself with data, quality, etc., BARQUE
2. Literature research: what are the newest versions of freshwater fish databases for 12S, COI (and maybe also other marker genes)
3. Assign taxonomy using newest database versions and most up-to-date algorithms
4. Also trial commercial approach that has been used previously
5. Compare outputs (e.g., TAXXI framework)
6. Using the comprehensive list of Canadian species, can we curate our own dbs from public data (Genbank, BOLD, midori?)
7. Compare output to electro fishing id results
8. Community ecology of fish species (after we obtain metadata)
9. What is the effect of a general Canadian database versus a regional (ecozone? ecoprovince?) database?
10. What are our tax blindspots? Can we fill those using stored id'ed fish tissue?

Notes: also look at barque workflow, contact Vincent Fugere, Martin Laporte, Michelle in Rene's lab, meet with Gerald: What has been previously done? Where are projects at?
