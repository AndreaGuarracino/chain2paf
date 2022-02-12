extern crate clap;

use std::borrow::Borrow;
use std::cmp::max;
use clap::{App, Arg}; //, SubCommand};

use std::fs::File;
use std::io::{BufRead, BufReader};
use rust_htslib::faidx;

fn _emit_paf_row(
    q_name: String, q_size: String, q_start: usize, q_end: usize, q_strand: String,
    t_name: String, t_size: String, t_start: usize, t_end: usize,
    num_matches: usize,
    cigar: String) {
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
        q_name, q_size, q_start, q_end, q_strand,
        t_name, t_size, t_start, t_end,
        num_matches,
        max(q_end - q_start, t_end - t_start),
        "255",
        cigar
    );
}

pub fn _fetch_subsequence(reader: &faidx::Reader, name: &str, begin: usize, end: usize) -> Vec<u8> {
    // 0-based
    reader.fetch_seq(name, begin, end).unwrap().to_owned()
}

pub fn complement(a: u8) -> u8 {
    match a {
        b'a' => b't',
        b'c' => b'g',
        b't' => b'a',
        b'g' => b'c',
        b'u' => b'a',
        b'A' => b'T',
        b'C' => b'G',
        b'T' => b'A',
        b'G' => b'C',
        b'U' => b'A',
        _ => b'N'
    }
}

/// Calculate reverse complement of given text (IUPAC alphabet supported).
///
/// Casing of characters is preserved, e.g. `b"NaCgT"` → `b"aCgTN"`.
/// All `N`s remain as they are.
///
/// ```
/// use bio::alphabets::dna;
///
/// assert_eq!(dna::revcomp(b"ACGTN"), b"NACGT");
/// assert_eq!(dna::revcomp(b"GaTtaCA"), b"TGtaAtC");
/// assert_eq!(dna::revcomp(b"AGCTYRWSKMDVHBN"), b"NVDBHKMSWYRAGCT");
/// ```
pub fn revcomp<C, T>(text: T) -> Vec<u8>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item=C>,
        T::IntoIter: DoubleEndedIterator,
{
    text.into_iter()
        .rev()
        .map(|a| complement(*a.borrow()))
        .collect()
}

fn main() {
    let matches = App::new("chain2paf")
        .version("0.1.0")
        .author("Andrea Guarracino")
        .about("Generate a PAF format file from a CHAIN format file")
        .arg(
            Arg::with_name("CHAIN")
                .required(true)
                .takes_value(true)
                .short("i")
                .long("input")
                .help("input CHAIN file"),
        ).arg(
        Arg::with_name("FASTA")
            .required(false)
            .takes_value(true)
            .short("f")
            .long("fasta")
            .help("input FASTA file to write =/X CIGAR operators (slower)"),
    )
        .get_matches();

    // Open the input CHAIN file
    let path_input_chain = matches.value_of("CHAIN").unwrap();
    let file = File::open(path_input_chain).unwrap();
    let chain_reader = BufReader::new(file);

    // Open the input FASTA file
    let path_input_fasta = matches.value_of("FASTA").unwrap_or("");
    let write_full_cigar = !path_input_fasta.is_empty();
    let fasta_reader = faidx::Reader::from_path(path_input_fasta).unwrap();
    //todo to fix
    // if write_full_cigar {
    //     fasta_reader = faidx::Reader::from_path(path_input_fasta).unwrap();
    // }

    // Initialize PAF information
    let mut t_name = String::from("");
    let mut t_size = String::from("");
    //let mut t_strand;
    let mut t_start: usize = 0;
    let mut t_end: usize = 0;
    let mut q_name = String::from("");
    let mut q_size = String::from("");
    let mut q_strand = String::from("");
    let mut q_start: usize = 0;
    let mut q_end: usize = 0;
    //let id;
    let mut num_matches: usize = 0;
    let mut cigar = String::from("");

    let mut q_seq: Vec<u8> = Vec::new();
    let mut t_seq: Vec<u8> = Vec::new();

    let mut q_offset: usize = 0;
    let mut t_offset: usize = 0;

    for line in chain_reader.lines() {
        let l = &line.unwrap();
        let l_tab = l.replace(" ", "\t"); // There are CHAIN files with space-separated fields
        let l_vec: Vec<&str> = l_tab.split('\t').collect();

        if l_vec[0].eq("chain") {
            //eprintln!("{:?}", l_vec);

            // Emit previous PAF row, if available
            if !cigar.is_empty() {
                _emit_paf_row(
                    q_name, q_size, q_start, q_end, q_strand,
                    t_name, t_size, t_start, t_end,
                    num_matches,
                    cigar,
                );
            }

            // Save query/target information
            // score -- chain score
            // tName -- chromosome (reference sequence)
            // tSize -- chromosome size (reference sequence)
            // tStrand -- strand (reference sequence)
            // tStart -- alignment start position (reference sequence)
            // tEnd -- alignment end position (reference sequence)
            // qName -- chromosome (query sequence)
            // qSize -- chromosome size (query sequence)
            // qStrand -- strand (query sequence)
            // qStart -- alignment start position (query sequence)
            // qEnd -- alignment end position (query sequence)
            // id -- chain ID
            t_name = l_vec[2].to_string();
            t_size = l_vec[3].to_string();
            //t_strand = line_vec[4].to_string();
            t_start = l_vec[5].parse::<usize>().unwrap();
            t_end = l_vec[6].parse::<usize>().unwrap();
            q_name = l_vec[7].to_string();
            q_size = l_vec[8].to_string();
            q_strand = l_vec[9].to_string();
            q_start = l_vec[10].parse::<usize>().unwrap();
            q_end = l_vec[11].parse::<usize>().unwrap();
            //let id = l_vec[12].clone();

            // Initialize PAF fields
            num_matches = 0;
            cigar = String::from("");

            if write_full_cigar {
                //eprintln!("Fetch subquery {} - {}", q_start, q_end);
                if q_strand == "+" {
                    q_seq = _fetch_subsequence(
                        &fasta_reader,
                        &q_name,
                        q_start,
                        q_end,
                    );
                } else {
                    q_seq = _fetch_subsequence(
                        &fasta_reader,
                        &q_name,
                        q_size.parse::<usize>().unwrap() - q_end,
                        q_size.parse::<usize>().unwrap() - q_start - 1,
                    );

                    q_seq = revcomp(q_seq);
                }
                // for x in q_seq.iter() {
                //     eprint!("{}", *x as char);
                // }
                // eprintln!("");

                //eprintln!("Fetch target {} - {}", t_start, t_end);
                t_seq = _fetch_subsequence(
                    &fasta_reader,
                    &t_name,
                    t_start,
                    t_end,
                );
                // for x in t_seq.iter() {
                //     eprint!("{}", *x as char);
                // }
                // eprintln!("");

                q_offset = 0;
                t_offset = 0;
            }
        } else {
            //eprintln!("{:?}", l_vec);

            // size -- the size of the ungapped alignment
            // dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
            // dq -- the difference between the end of this block and the beginning of the next block (query sequence)
            let size_ungapped_alignment = if !l_vec[0].is_empty() { l_vec[0] } else { "0" };
            let diff_in_target = if l_vec.len() > 1 { l_vec[1] } else { "0" };
            let diff_in_query = if l_vec.len() > 2 { l_vec[2] } else { "0" };

            if !size_ungapped_alignment.eq("0") {
                if write_full_cigar {
                    let mut last_op = ' ';
                    let mut len_last_op = 0;

                    for _ in 0..size_ungapped_alignment.parse::<u64>().unwrap() {
                        let q_character = q_seq[q_offset];
                        let t_character = t_seq[t_offset];

                        // eprintln!(
                        //     "{} ({} - {}) == {} ({} = {})? {} -- {} {}",
                        //     q_offset, q_start + q_offset, q_character as char,
                        //     t_offset, t_start + t_offset, t_character as char,
                        //     if q_character == t_character { "YES" } else { "NO"},
                        //     last_op, len_last_op
                        // );

                        if q_character == t_character {
                            if last_op != ' ' && last_op != '=' {
                                //eprintln!("\t ------> {}X", len_last_op);

                                cigar.push_str(&*len_last_op.to_string());
                                cigar.push('X');

                                len_last_op = 0;
                            }

                            last_op = '=';
                            len_last_op += 1;
                        } else {
                            if last_op != ' ' && last_op != 'X' {
                                //eprintln!("\t ------> {}=", len_last_op);

                                num_matches += len_last_op;

                                cigar.push_str(&*len_last_op.to_string());
                                cigar.push('=');

                                len_last_op = 0;
                            }

                            last_op = 'X';
                            len_last_op += 1;
                        }

                        q_offset += 1;
                        t_offset += 1;
                    }

                    if last_op == '=' {
                        num_matches += len_last_op;
                    }

                    //eprintln!("\t ------> {}{}", len_last_op, last_op);
                    cigar.push_str(&*len_last_op.to_string());
                    cigar.push(last_op);
                } else {
                    num_matches += size_ungapped_alignment.parse::<usize>().unwrap();

                    cigar.push_str(size_ungapped_alignment);
                    cigar.push('M');
                }
            }
            if !diff_in_target.eq("0") {
                if write_full_cigar {
                    // Keep the offset updated
                    t_offset += diff_in_target.parse::<usize>().unwrap();
                }

                cigar.push_str(diff_in_target);
                cigar.push('D');
            }
            if !diff_in_query.eq("0") {
                if write_full_cigar {
                    // Keep the offset updated
                    q_offset += diff_in_query.parse::<usize>().unwrap();
                }

                cigar.push_str(diff_in_query);
                cigar.push('I');
            }
        }
    }

    // Emit last PAF row, if available
    if !cigar.is_empty() {
        _emit_paf_row(
            q_name, q_size, q_start, q_end, q_strand,
            t_name, t_size, t_start, t_end,
            num_matches,
            cigar,
        );
    }
}
