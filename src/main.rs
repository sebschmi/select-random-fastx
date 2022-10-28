use clap::Parser;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use log::{info, warn, LevelFilter};
use rand::distributions::{Distribution, Uniform, WeightedIndex};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro512PlusPlus;
use seq_io::fastx::dynamic::FastxReader;
use seq_io::fastx::Reader;
use seq_io::{fasta, fastq};
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::{fs, mem};

#[derive(Parser)]
struct Cli {
    /// Fasta or fastq input file (automatically detected).
    /// Pass multiple times for multiple files.
    #[clap(long)]
    input: Vec<PathBuf>,

    /// Output fastx file.
    #[clap(long)]
    output: PathBuf,

    /// The amount of entries to randomly select.
    #[clap(long)]
    amount: Option<u64>,

    /// Instead of selecting random entries, simply concatenate the input files.
    #[clap(long, conflicts_with = "amount", required_unless_present = "amount")]
    concatenate: bool,

    /// If true, allow to randomly select a fastx entry more than once.
    /// If this is not given and there are not enough entries, the program aborts.
    #[clap(long)]
    allow_repetitions: bool,
}

pub fn initialise_logging(log_level: LevelFilter) {
    CombinedLogger::init(vec![TermLogger::new(
        if cfg!(debug_assertions) {
            LevelFilter::Trace
        } else {
            log_level
        },
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )])
    .unwrap();

    info!("Logging initialised successfully");
}

fn progress_bar_style() -> ProgressStyle {
    ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})",
    )
    .unwrap()
    .with_key(
        "eta",
        |state: &ProgressState, w: &mut dyn std::fmt::Write| {
            write!(w, "{:.0}s", state.eta().as_secs_f64()).unwrap()
        },
    )
    .progress_chars("#>-")
}

struct OutputWriter<Output: std::io::Write> {
    writer: Output,
    next_index: u64,
}

fn copy_entries(
    source_file: &Path,
    target_file: &mut OutputWriter<impl std::io::Write>,
    index_iterator: impl Iterator<Item = u64>,
    progress_bar: &ProgressBar,
) -> Result<(), String> {
    let mut fastx_reader =
        Reader::new(BufReader::new(File::open(source_file).map_err(|err| {
            format!("Could not open file {source_file:?} for reading: {err}")
        })?));
    let mut current_index = 0;
    let mut current_record = fastx_reader
        .next()
        .unwrap()
        .map_err(|err| format!("Unable to parse fastx record: {err}"))?;

    for next_index in index_iterator {
        while next_index > current_index {
            current_index += 1;
            current_record = fastx_reader
                .next()
                .unwrap()
                .map_err(|err| format!("Unable to parse fastx record: {err}"))?;
        }

        if let Some(qual_lines) = current_record.opt_qual_lines() {
            fastq::write_iter(
                &mut target_file.writer,
                format!("{}", target_file.next_index).as_bytes(),
                current_record.seq_lines(),
                qual_lines,
            )
            .map_err(|err| format!("Cannot write fastq record: {err}"))?;
        } else {
            fasta::write_iter(
                &mut target_file.writer,
                format!("{}", target_file.next_index).as_bytes(),
                current_record.seq_lines(),
            )
            .map_err(|err| format!("Cannot write fasta record: {err}"))?;
        }

        target_file.next_index += 1;
        progress_bar.inc(1);
    }

    Ok(())
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    initialise_logging(LevelFilter::Info);

    info!("Reading file metadata...");
    let pb = ProgressBar::new(cli.input.len() as u64);
    pb.set_style(progress_bar_style());
    let mut file_size_sum = 0u64;
    for path in &cli.input {
        file_size_sum = file_size_sum
            .checked_add(
                fs::metadata(path)
                    .map_err(|err| format!("Could access metadata of file {path:?}: {err}"))?
                    .len(),
            )
            .ok_or_else(|| {
                "Overflow when summing file sizes, total file size exceeds 2^64".to_string()
            })?;
        pb.inc(1);
    }
    pb.finish_and_clear();

    info!("Counting fastx entries...");
    let pb = ProgressBar::new(file_size_sum);
    pb.set_style(progress_bar_style());
    let mut entry_count_per_file = Vec::new();
    let mut entry_count_sum = 0;
    let mut format = None;
    for path in &cli.input {
        let mut entry_count = 0u64;
        let mut fastx_reader =
            Reader::new(BufReader::new(File::open(path).map_err(|err| {
                format!("Could not open file {path:?} for reading: {err}")
            })?));
        while let Some(record) = fastx_reader.next() {
            let _ = record.map_err(|err| format!("Unable to parse fastx record: {err}"))?;
            entry_count = entry_count.checked_add(1).ok_or_else(|| {
                format!("Overflow (>2^64) when counting records in file {path:?}")
            })?;
        }
        entry_count_per_file.push(entry_count);
        entry_count_sum += entry_count;

        if let Some(format) = format.as_mut() {
            if let Some(reader_format) = fastx_reader.format() {
                if reader_format != *format {
                    return Err(format!("Mismatched sequence formats, the first nonempty file is {format:?}, but {path:?} is {reader_format:?}"));
                }
            }
        } else {
            format = fastx_reader.format()
        }
    }
    pb.finish_and_clear();

    let mut output = OutputWriter {
        writer: BufWriter::new(File::create(&cli.output).map_err(|err| {
            format!(
                "Could not create and open output file {:?} for writing: {err}",
                &cli.output
            )
        })?),
        next_index: 0,
    };
    let mut rng = Xoshiro512PlusPlus::from_entropy();

    if let Some(amount) = cli.amount {
        if cli.allow_repetitions {
            info!("Assigning {amount} random draws with repetition to files...");
            let pb = ProgressBar::new(amount);
            pb.set_style(progress_bar_style());
            let mut amount_per_file = vec![0usize; entry_count_per_file.len()];
            let distribution = WeightedIndex::new(&entry_count_per_file).unwrap();
            for _ in 0..amount {
                let amount_per_file = &mut amount_per_file[distribution.sample(&mut rng)];
                *amount_per_file = (*amount_per_file).checked_add(1).ok_or_else(|| {
                    format!(
                        "Overflow (>2^{}) when assigning random draws to file",
                        mem::size_of::<usize>() * 8
                    )
                })?;
                pb.inc(1);
            }
            pb.finish_and_clear();

            info!("Drawing {amount} records with repetition...");
            let pb = ProgressBar::new(amount);
            pb.set_style(progress_bar_style());

            for (entry_count_per_file, (amount_per_file, source_path)) in entry_count_per_file
                .iter()
                .copied()
                .zip(amount_per_file.iter().copied().zip(cli.input.iter()))
            {
                if amount_per_file == 0 {
                    continue;
                }

                if u128::try_from(amount_per_file).unwrap() > u128::from(entry_count_per_file) {
                    warn!("Using inefficient algorithm for a large amount of repetitious draws");
                }

                let mut indices: Vec<_> = Uniform::new(0, entry_count_per_file)
                    .sample_iter(&mut rng)
                    .take(amount_per_file)
                    .collect();
                indices.sort_unstable();

                copy_entries(source_path, &mut output, indices.iter().copied(), &pb)?;
            }

            pb.finish_and_clear();
        } else {
            info!("Assigning {amount} random draws without repetitions to files...");
            let pb = ProgressBar::new(amount);
            pb.set_style(progress_bar_style());
            let mut remaining_entries_per_file = entry_count_per_file.clone();
            let mut amount_per_file = vec![0usize; entry_count_per_file.len()];
            let mut distribution = WeightedIndex::new(&entry_count_per_file).unwrap();
            for _ in 0..amount {
                let file_index = distribution.sample(&mut rng);

                let amount_per_file = &mut amount_per_file[file_index];
                *amount_per_file = (*amount_per_file).checked_add(1).ok_or_else(|| {
                    format!(
                        "Overflow (>2^{}) when assigning random draws to file",
                        mem::size_of::<usize>() * 8
                    )
                })?;

                remaining_entries_per_file[file_index] -= 1;
                distribution
                    .update_weights(&[(file_index, &remaining_entries_per_file[file_index])])
                    .map_err(|err| format!("Error updating random weights: {err}"))?;

                pb.inc(1);
            }
            pb.finish_and_clear();

            info!("Drawing {amount} records without repetitions...");
            let pb = ProgressBar::new(amount);
            pb.set_style(progress_bar_style());

            for (entry_count_per_file, (amount_per_file, source_path)) in entry_count_per_file
                .iter()
                .copied()
                .zip(amount_per_file.iter().copied().zip(cli.input.iter()))
            {
                if amount_per_file == 0 {
                    continue;
                }

                assert!(
                    u128::from(entry_count_per_file) >= u128::try_from(amount_per_file).unwrap()
                );
                let mut indices: Vec<_> = (0..entry_count_per_file).collect();
                indices.shuffle(&mut rng);
                indices.resize(amount_per_file, 0);
                indices.sort_unstable();

                copy_entries(source_path, &mut output, indices.iter().copied(), &pb)?;
            }

            pb.finish_and_clear();
        }
    } else {
        info!("Concatenating the input files (no random selection)...");
        let pb = ProgressBar::new(entry_count_sum);
        pb.set_style(progress_bar_style());
        for (entry_count, source_path) in entry_count_per_file.iter().copied().zip(cli.input.iter())
        {
            if entry_count == 0 {
                continue;
            }
            copy_entries(source_path, &mut output, 0..entry_count, &pb)?;
        }
        pb.finish_and_clear();
    }

    output
        .writer
        .flush()
        .map_err(|err| format!("Unable to flush output buffer: {err}"))?;

    info!("Done");
    Ok(())
}
