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

Partition the assembly contigs by chromosome by mapping each assembly against the scaffolded references, and then subsetting the graph. Here we use [wfmash](https://github.com/ekg/wfmash) for the mapping.

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

Then we can merge both the reference scaffolds and HPRCy1 contigs:

```
( seq 22; echo X; echo Y; echo M ) | while read i; do sbatch -p lowmem -c 8 --wrap '( samtools faidx assemblies/chm13.fa.gz chm13#chr'$i'; samtools faidx assemblies/grch38.fna.gz grch38#chr'$i'; cat parts/chr'$i'.pan.fa ) >parts/chr'$i'.pan+refs.fa && samtools faidx parts/chr'$i'.pan+refs.fa' >> slurm.jobids; done
```

We will use these files directly in [pggb](https://github.com/pangenome/pggb).


## graph generation

We apply [pggb](https://github.com/pangenome/pggb).

```
( seq 22 ; echo X; echo Y; echo M ) | while read i; do sbatch -c 48 --wrap 'pggb -i parts/chr'$i'.pan+refs.fa -s 15000 -p 98 -l 100000 -w 500000 -j 15000 -e 15000 -n 7 -t 48 -v -Y "#" -k 27 -B 30000000 -I 0.6 -R 0.2 -C 100,1000,10000,10000:parts/refs/chm13+grch38.chr'$i'.txt:y,10000:parts/refs/chm13+grch38.chr'$i'.txt:n -o graphs/chr'$i'.pan+refs' >>wg.jobids; done
```

Todo.

## graph evaluation

We apply [pgge](https://github.com/pangenome/pgge).

Todo.
