# chain2paf
convert CHAIN format to PAF format

## installation

`chain2paf` is built with rust, and so we install using `cargo`:

```
https://github.com/AndreaGuarracino/chain2paf
cd chain2paf
cargo install --force --path .
```

## usage

With alignments in `data/hg18ToHg19.over.chain`, we would convert it into a PAF format file using this call:

```
chain2paf -i data/hg18ToHg19.over.chain > hg18ToHg19.over.chain.paf
```
