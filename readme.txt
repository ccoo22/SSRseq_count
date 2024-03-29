An str count tool based on NGS SSR sequence

perl str_count.pl
Options: [required]

        --fastq_dir/-i                Paired-end fastq file dir
                                      file name format: sample_name_R1.fastq.gz, sample_name_R2.fastq.gz
        --target_fasta/-t             pcr targeted fasta seq, seq name should be named as our designed.
        --output_dir/-o               output dir
        --parallel                    how many samples will be analyzed in parallel, example : 10
        --software_flash              path/to/flash     flash software directory. please download and install by yourself. http://ccb.jhu.edu/software/FLASH/
        --software_blastn             path/to/blastn      blastn software directory. please download and install by yourself.  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
        --software_fastq_to_fasta     path/to/fastq_to_fasta    fastq_to_fasta software directory. please download and install by yourself. http://hannonlab.cshl.edu/fastx_toolkit/commandline.html

Options: [optional]
        --sample_list/-s          samples that need to be analyzed. If there are multiple samples, separate them with commas. example: sample1,sample2,sample3
                                  we will use all samples in fastq_dir by default.
         --help/-h                help doc


how target_fasta seq named?

example1:
>seq_name1
GTCAATCATGGCACAACAGCACAACTGTGTGAGATAATTTATTTTCCTATCTCATCATCATCATCATCATCGTCACATGACCACACTGCACTCCACCACCATGCAACCACCACATCTGCACCTCCTAACCACCACCCTCCTCCTTCTGCTCCTCTGGGCCGTAACTTTTCCGGCATTCCTCTCCGCCGTCAACCTCCCTGAATTCCGTGAAGCCCCAGCATTCCGAAACGGAAACCAATGCCCCAACCCAACCTCCTCTTCCTCCACCAT
only one motif will be analyzed in seq_name1, and motif seq is TCA, start from 52 to 69 in target seq.
this seq name in fasta file should be named as:  seq_name1|TCA(52-69)


example2:
>seq_name2
GTCAATCATGGCACAACAGCACAACTGTGTGAGATAATTTATTTTCCTATCTCATCATCATCATCATCATCGTCACATGACCACACTGCACTCCACCTCATCATCATCATCATCATCATCATCAAACCACCACCCTCCTCCTTCTGCTCCTCTGGGCCGTAACTTTTCCGGCATTCCTCTCCGCCGTCAACCTCCCTGAATTCCGTGAAGCCCCAGCATTCCGAAACGGAAACCAATGCCCCAACCCAACCTCCTCTTCCTCCACCAT
there are two motif should be analyzed in seq_name2: TCA and TCA
we should named it as:  seq_name2|TCA(52-69);TCA(98-124)


example3:
>seq_name3
GTCAATCATGGCACAACAGCACAACTGTGTGAGATAATTTATTTTCCTATCTCATCATCATCATCATCATCGTCACATGACCACACTGCACTCCACCCCTCCTCCTCCTCCTCCTCCTCCTAACCACCACCCTCCTCCTTCTGCTCCTCTGGGCCGTAACTTTTCCGGCATTCCTCTCCGCCGTCAACCTCCCTGAATTCCGTGAAGCCCCAGCATTCCGAAACGGAAACCAATGCCCCAACCCAACCTCCTCTTCCTCCACCAT
there are two motif should be analyzed in seq_name3: TCA and CCT
we should named it as:  seq_name3|TCA(52-69);CCT(98-121)


Dependencies: 
perl packages
    File::Spec;
    Getopt::Long;
    Bio::SeqIO;
    Bio::SearchIO;
    Parallel::ForkManager;

software
    FLASH: merges reads from paired-end sequencing. http://ccb.jhu.edu/software/FLASH/
    BLASTN: perform BLAST searches. ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    fastq_to_fasta: convert fastq to fasta. http://hannonlab.cshl.edu/fastx_toolkit/commandline.html

example usage: 
perl str_count.pl -i example/fastq_example  -s CC10,CC11,CC12 -t example/target.fa -o ./example/result --parallel 10 --software_flash /path/to/flash --software_blastn /path/to/blastn --software_fastq_to_fasta /path/to/fastq_to_fasta

Citing 
Cui X, Li C, Qin S, Huang Z, Gan B, Jiang Z, Huang X, Yang X, Li Q, Xiang X, Chen J, Zhao Y, Rong J. High-throughput sequencing-based microsatellite genotyping for polyploids to resolve allele dosage uncertainty and improve analyses of genetic diversity, structure and differentiation: A case study of the hexaploid Camellia oleifera. Mol Ecol Resour. 2021 Jul 14. doi: 10.1111/1755-0998.13469. Epub ahead of print. PMID: 34260828.

