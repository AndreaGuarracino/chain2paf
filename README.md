# chain2paf
convert CHAIN format to PAF format

## installation

`chain2paf` is built with rust, and so we install using `cargo`:

```shell
https://github.com/AndreaGuarracino/chain2paf
cd chain2paf
cargo install --force --path .
```

## usage

With alignments in `data/hg18ToHg19.over.chain`, we would convert it into a PAF format file using this call:

```shell
chain2paf -i data/hg18ToHg19.over.chain > hg18ToHg19.over.chain.paf
```

If full CIGAR strings (with `=`/`X` operators) are required, we can specify two FASTA files (uncompressed or bgzipped),
the 1-st for the targets and the 2-nd for the queries. For example, if the sequences of the human genome versions 
hg18 and hg19 are available in the `hg18.fa.gz` and `hg19.fa.gz` files, you can execute

```shell
chain2paf -i data/hg18ToHg19.over.chain -f hg18.fa.gz hg19.fa.gz > hg18ToHg19.over.chain.paf
```

If the CHAIN file is the result of a pairwise alignment, you can specify the same FASTA file for both targets and queries:

```shell
chain2paf -i input.chain -f input.fa.gz input.fa.gz > input.chain.paf
```

## info

`chain2paf` performs the reverse operation of [paf2chain](https://github.com/AndreaGuarracino/paf2chain).
