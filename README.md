# HPRC year 1 pggb analysis pipeline

## data

Get the assemblies.

```
cd assemblies
time cat ../Year1_F1_assemblies.index.txt| cut -f 3 | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%' | parallel -j 4 'wget -q {} && echo got {}'
```

## preprocessing

Prefix them and adjust their naming:

```
ls *.paternal.f1_assembly_v1.fa.gz | sort -V | while read f; do echo '#!/bin/bash\n'f=$f'\no=$(basename $f .f1_assembly_v1.fa.gz)\nzcat $f | sed s%#1#%.paternal#% | bgzip >$o.fa.gz && samtools faidx $o.fa.gz' | sbatch ; done
ls *.maternal.f1_assembly_v1.fa.gz | sort -V | while read f; do echo '#!/bin/bash\n'f=$f'\no=$(basename $f .f1_assembly_v1.fa.gz)\nzcat $f | sed s%#2#%.maternal#% | bgzip >$o.fa.gz && samtools faidx $o.fa.gz' | sbatch ; done
```

Apply name prefixes to the reference sequences.

```
( fastix -p 'grch38#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) | bgzip >grch38.fna.gz && samtools faidx grch38.fna.gz ) &
( fastix -p 'chm13#' <(zcat chm13.draft_v1.0.fasta.gz) | bgzip >chm13.fa.gz && samtools faidx chm13.fa.gz ) &
```

Combine them into a single reference for competitive assignment of sample contigs to chromosome bins.

```
zcat chm13.fa.gz grch38.fna.gz >chm13+grch38.pan.fa && samtools faidx chm13+grch38.pan.fa
```

Partition the assembly contigs by chromosome by mapping each assembly against the scaffolded references, and then subsetting the graph.

```
cd ..
ls assemblies | grep 'maternal.fa.gz$\|paternal.fa.gz$' | cut -f 1-2 -d . | sort -V | uniq >haps.list
ref=assemblies/chm13+grch38.pan.fa
for hap in $(cat haps.list);
do
    in=assemblies/$hap.fa.gz
    out=alignments/$hap.vs.ref.paf
    sbatch -p lowmem -c 16 --wrap "wfmash -t 16 -m -N -p 90 $ref $in >$out" >>slurm.jobids
done
```

Subset by chromosome.

```
( seq 22; echo X; echo Y; echo M ) | while read i; do awk '$6 ~ "chr'$i'$"' alignments/*.vs.ref.paf | cut -f 1 | sort -V >parts/chr$i.contigs
( seq 22; echo X; echo Y; echo M ) | while read i; do sbatch -p lowmem -c 24 --wrap 'cat parts/chr'$i'.contigs | parallel -j 24 "s=\$(echo {} | cut -f 1 -d\#); samtools faidx assemblies/\$s.fa.gz {}" >parts/chr'$i'.pan.fa' >> slurm.jobids; done
```


## graph generation

Todo.

We apply pggb.