
![DICAST](docs/source/img/logo.png )

Alternative splicing is a major contributor to transcriptome and proteome diversity in health and disease. A plethora of tools have been developed for studying alternative splicing in RNA-seq data. Previous benchmarking efforts focused on isoform quantification and mapping, neglecting event detection tools which arguably provide the most detailed insights into the alternative splicing process. 

DICAST closes this gap by offering a modular and extensible alternative splicing framework integrating eleven splice-aware mapping and nine event detection tools, which we benchmark extensively on simulated as well as whole blood RNA-seq data. We further propose the first uniform reporting standard to unify existing formats and to guide future tool development. The performance of event detection tools varies widely with no tool outperforming all others. DICAST allows researchers to employ a consensus approach to consider the most successful tools jointly in robust event detection. 

### Documentation

https://dicast.readthedocs.io/en/latest/index.html

#### Mapping tools
- bb- map
- contextmap2
- crac
- dart
- gsnap
- hisat2
- mapsplice2
- minimap2
- segemehl
- star
- subjunc

#### Alternative Splicing Event Detection Tools
- asgal
- aspli
- eventpointer
- irfinder
- majiq
- sgseq
- spladder
- whippet