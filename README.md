# mutspy
```
Mutation counts in bam files

Usage:
  mutspy.R count [-srDP] (-b bam | -l bamlist) [-o out -m mapq] <mut>...
  mutspy.R -h



Options:
  -h                       Show this screen.
  -s                       Count fw and bw strand separately.
  -r                       Count both bases of overlapping read pairs.
  -D                       Include duplicate reads (bitflag 0x400) [default: False]
  -P                       Print reads instead of counting [default: False]
  -b --bam=bam             Path to an indexed bam file.
  -l --bamlist=bamlist     Path to file listing indexed bam files one per line.
  -o --out=out             Path to file were output are appended [default: stdout].
  -m --mapq=mapq           Only include data with minimum this mapq [default: 10].

Arguments:
  <mut>                   Mutations is in the form <seqname>:<pos>_<ref>/<alt>, 
                          e.g. chr10:114925316_G/A (SNP), chr10:114925316_GGT/G (deletion),
                          chr10:114925316_G/GAT (insertion), and  chr10:114925316_GGT/GAT 
                          (MNP). Several mutations can be given as 'mut1 mut2 mutn'.
                          Alternatively a file with lines of sitemuts can be given.
                          The mutation N expands into all four bases.
                          Count outputs are 'alt total' or, if stranded,
                          'alt.fw alt.bw total.fw total.bw'
                          Note: mutspy is reference agnostic, only knowns the read sequence.

Requires:
  Rsamtools, docopt.

Version:
  0.3.1 
```


## Requires
- R and R packages docopt and Rsamtools

## Example
Clone the repo:
```bash
git clone https://github.com/lindbjerg-group/mutspy mutspy
```

Install required dependencies:
```bash
conda create -n mutspy-env bioconda::r-docopt bioconda::bioconductor-rsamtools
```

Run some mutspy commands:
```bash
conda activate mutspy-env
mutspy/src/mutspy.R --help
```

Make a stranded count three specific positions, one using
the ambiguous base N (expands into all four bases):
```bash
mutspy/src/mutspy.R count -s -b sample.bam 'chr2:233840772_TC/GA chr12:25245347_T/N chrX:108175747_G/GT'
```

## Output

|Sample Name |Chr  |Start     |Ref seq |Alt seq |Fwd alt |Rev alt |Fwd total |Rev total |
|------------|-----|----------|--------|--------|--------|--------|----------|----------|
|sample.bam  |chr2 |233840772 |TC      |GA      |0       |0       |74        |69        |   
|sample.bam  |chr12|25245347  |C       |A       |0       |0       |90        |75        |
|sample.bam  |chrX |108175747 |G       |GT      |0       |0       |34        |43        |
|sample.bam  |chr12|25245347  |C       |T       |14      |14      |90        |75        |
|sample.bam  |chr12|25245347  |C       |C       |76      |61      |90        |75        |
|sample.bam  |chr12|25245347  |C       |G       |0       |0       |90        |75        |
 


## Contact
- mads.heilskov@clin.au.dk

