extern crate clap;

use std::cmp::max;
use clap::{App, Arg}; //, SubCommand};

use std::fs::File;
use std::io::{BufRead, BufReader};

fn _emit_paf_row(
    q_name: String, q_size: String, q_start: String, q_end: String, q_strand: String,
    t_name: String, t_size: String, t_start: String, t_end: String,
    num_matches: u64,
    cigar: String) {
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
        q_name, q_size, q_start, q_end, q_strand,
        t_name, t_size, t_start, t_end,
        num_matches,
        max(
            q_end.parse::<u64>().unwrap() - q_start.parse::<u64>().unwrap(),
            t_end.parse::<u64>().unwrap() - t_start.parse::<u64>().unwrap(),
        ),
        "255",
        cigar
    );
}

fn main() {
    let matches = App::new("chain2paf")
        .version("0.1.0")
        .author("Andrea Guarracino")
        .about("Generate a PAF format file from a CHAIN format file")
        .arg(
            Arg::with_name("INPUT")
                .required(true)
                .takes_value(true)
                .short("i")
                .long("input")
                .help("input CHAIN file"),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);

    let mut t_name = String::from("");
    let mut t_size = String::from("");
    //let mut t_strand;
    let mut t_start = String::from("");
    let mut t_end = String::from("");
    let mut q_name = String::from("");
    let mut q_size = String::from("");
    let mut q_strand= String::from("");
    let mut q_start = String::from("");
    let mut q_end = String::from("");
    //let id;

    let mut num_matches: u64 = 0;
    let mut cigar = String::from("");

    for line in reader.lines() {
        let l = &line.unwrap();
        let l_vec: Vec<&str> = l.split('\t').collect();

        if l_vec[0].eq("chain") {
            //eprintln!("{:?}", l_vec);

            // Emit previous PAF row, if available
            if !cigar.is_empty() {
                _emit_paf_row(
                    q_name, q_size, q_start, q_end, q_strand,
                    t_name, t_size, t_start, t_end,
                    num_matches,
                    cigar
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
            t_start = l_vec[5].to_string();
            t_end = l_vec[6].to_string();
            q_name = l_vec[7].to_string();
            q_size = l_vec[8].to_string();
            q_strand = l_vec[9].to_string();
            q_start = l_vec[10].to_string();
            q_end = l_vec[11].to_string();
            //let id = l_vec[12].clone();

            // Initialize PAF fields
            num_matches = 0;
            cigar = String::from("");
        } else {
            //eprintln!("{:?}", l_vec);

            // size -- the size of the ungapped alignment
            // dt -- the difference between the end of this block and the beginning of the next block (reference sequence)
            // dq -- the difference between the end of this block and the beginning of the next block (query sequence)
            let size_ungapped_alignment = if !l_vec[0].is_empty() {l_vec[0]} else {"0"};
            let diff_in_target = if l_vec.len() > 1 {l_vec[1]} else {"0"};
            let diff_in_query = if l_vec.len() > 2 {l_vec[2]} else {"0"};

            if !size_ungapped_alignment.eq("0") {
                num_matches += size_ungapped_alignment.parse::<u64>().unwrap();

                cigar.push_str(size_ungapped_alignment);
                cigar.push('M');
            }
            if !diff_in_target.eq("0") {
                cigar.push_str(diff_in_target);
                cigar.push('D');
            }
            if !diff_in_query.eq("0") {
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
            cigar
        );
    }
}
