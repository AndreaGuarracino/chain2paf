# chain2paf
convert CHAIN format to PAF format

## installation

`paf2chain` is built with rust, and so we install using `cargo`:

```
https://github.com/AndreaGuarracino/chain2paf
cd chain2paf
cargo install --force --path .
```

## usage

With alignments in `aln.chain`, we would convert it into a PAF format file using this call:

```
chain2paf -i aln.chain > aln.paf
```
