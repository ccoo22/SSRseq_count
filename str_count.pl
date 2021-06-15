$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
ssr reads count
Version: v1.0 2020-06-21
Contact: 129 ganb
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Parallel::ForkManager;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($fastq_dir, $sample_list, $target_fasta, $output_dir, $parallel, $SOFT_FLASH, $SOFT_BLASTN, $SOFT_FASTQ_TO_FASTA, $if_help);
GetOptions(
        "fastq_dir|i=s"    => \$fastq_dir,
        "sample_list|s=s"  => \$sample_list,
        "target_seq|t=s"   => \$target_fasta,
        "output_dir|o=s"   => \$output_dir,

        "parallel=s"       => \$parallel,
        "software_flash=s"     => \$SOFT_FLASH,
        "software_blastn=s"    => \$SOFT_BLASTN,
        "software_fastq_to_fasta=s"    => \$SOFT_FASTQ_TO_FASTA,

        "help|h"           => \$if_help,
);
die "
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
        --help/-h                 help doc

\n" if (defined $if_help or not defined $fastq_dir or not defined $target_fasta or not defined $output_dir or not defined $parallel or not defined $SOFT_FLASH or not defined $SOFT_BLASTN);

$sample_list = get_sample_list($fastq_dir) if(not defined $sample_list);
###################################################################### 初始化
make_dir($output_dir);

my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET]  fastq_dir : $fastq_dir
[SET]  sample_list : $sample_list
[SET]  target_fasta : $target_fasta
[SET]  output_dir : $output_dir
[SET]  parallel : $parallel
[SET]  soft_flash : $SOFT_FLASH
[SET]  soft_blastn : $SOFT_BLASTN
[SET]  soft_fastq_to_fasta : $SOFT_FASTQ_TO_FASTA
";
open SAVE, ">>$output_dir/run.info"; 
print SAVE $SCRIPT_INFO.$RUN_INFO; 
close SAVE;
print $RUN_INFO;

die "[Error] lost flash soft\n"          if(is_file_ok($SOFT_FLASH) == 0);
die "[Error] lost blastn soft\n"         if(is_file_ok($SOFT_BLASTN) == 0);
die "[Error] lost fastq_to_fasta soft\n" if(is_file_ok($SOFT_FASTQ_TO_FASTA) == 0);
die "[Error] lost target_fasta\n"        if(is_file_ok($target_fasta) == 0);
die "[Error] lost fastq_dir\n"           if(is_dir_ok($fastq_dir) == 0);

my $blast_plus_dir = `dirname $SOFT_BLASTN`;
   $blast_plus_dir=~s/[\r\n]//g;
my $SOFT_MAKEBLASTDB = "$blast_plus_dir/makeblastdb";
die "[Error] lost makeblastdb soft. please set blastn and makeblastdb at the same directory\n" if(is_file_ok($SOFT_MAKEBLASTDB) == 0);

my @samples = split /,/, $sample_list;

###################################################################### 主程序
check_fastq();  
my $target_motif = check_target_motif();  
my $target_fasta_db = build_blast_index();  

# sample result dir
print "[process] start mapping sample and check fastq str info\n";
my $mapping_dir = "$output_dir/mapping";
make_dir($mapping_dir);

my $pm = Parallel::ForkManager->new($parallel);
foreach my $sample(@samples)
{
    $pm->start() and next;
    sub_run($sample);
    $pm->finish;    
}
$pm->wait_all_children;

print "[process] summary\n";
quality_control();
str_count_summary();


###################################################################### 子程序

sub str_count_summary{
    my %hashSTR;

    # read sample str count
    foreach my $sample(@samples)
    {
        my $str_count_file = "$mapping_dir/$sample/$sample.str_count.txt";
        open STR, $str_count_file;
        my $line1 = <STR>;
           $line1 =~ s/[\r\n]//g;
        my @heads = split /\t/, $line1;

        while(<STR>)
        {
            $_ =~ s/[\r\n]//g;
            next if($_!~/\w/);
            my @datas = split /\t/, $_;
            my %hashTmp = map{ ($heads[$_], $datas[$_])  } (0..$#heads);

            my $target       = $hashTmp{"target"};
            my $aim_str      = $hashTmp{"aim_str"};
            my $str_identify = $hashTmp{"str_identify"};
            my $reads_count  = $hashTmp{"reads_count"};
            $hashSTR{$target}{$aim_str}{$str_identify}{$sample} = $reads_count;            
        }
        close STR;
    }

    # output
    my $str_count_summary = "$output_dir/str_count.txt";
    open OUT, ">$str_count_summary";

    print OUT "Target\tAimSTR\tSTR\t" . (join "\t", @samples) . "\n";
    foreach my $target(sort keys %hashSTR)
    {
        foreach my $aim_str(sort keys %{$hashSTR{$target}})
        {   
            my ($motif, $start, $end) = split /,/, $aim_str;
            my $motif_length = length($motif);
            foreach my $str_identify (sort keys %{$hashSTR{$target}{$aim_str}})
            {
                my $str_motif_count = length($str_identify) / $motif_length;
                my @values = ($target, $aim_str, "$motif($str_motif_count)");
                foreach my $sample(@samples)
                {
                    my $reads_count = exists $hashSTR{$target}{$aim_str}{$str_identify}{$sample} ? $hashSTR{$target}{$aim_str}{$str_identify}{$sample} : "";
                    push @values, $reads_count;
                }
                print OUT (join "\t", @values) . "\n";
            }
        }
    }
    close OUT;

    print "STR Count: $str_count_summary\n";
    
}

sub quality_control{
    my %hashQC;
    my %hashQCT;

    foreach my $sample(@samples)
    {
        my $flash_log         = "$mapping_dir/$sample/flash.log.txt";
        my $statistics_target = "$mapping_dir/$sample/$sample.filter_ok.txt";

        open FLASH, $flash_log;
        while(<FLASH>)
        {
            if($_=~/Total pairs:\s+(\d+)/)
            {
                $hashQC{$sample}{'RawReads'} = $1;
                next;
            }
            if($_=~/Combined pairs:\s+(\d+)/)
            {
                $hashQC{$sample}{'MergeReads'} = $1;
                next;
            }
        }
        close FLASH;

        open TARGET, $statistics_target;
        while(<TARGET>)
        {
            $_=~s/[\r\n]//g;
            my ($target, $count) = split /\t/, $_;
            $hashQCT{$target}{$sample} = $count;
            $hashQC{$sample}{'FilterReads'} += $count;
        }
        close TARGET;

        $hashQC{$sample}{'Sample'} = $sample;
        $hashQC{$sample}{'MergePerc'} = sprintf "%0.2f", $hashQC{$sample}{'MergeReads'} / $hashQC{$sample}{'RawReads'};
        $hashQC{$sample}{'FilterPerc'} = sprintf "%0.2f", $hashQC{$sample}{'FilterReads'} / $hashQC{$sample}{'MergeReads'};

    }

    # output result
    my $qc = "$output_dir/qc.txt";
    my @heads_qc = qw(Sample RawReads MergeReads MergePerc FilterReads FilterPerc);
    open QC, ">$qc";
    print QC (join "\t", @heads_qc) . "\n";
    foreach my $sample(@samples)
    {
        my @values = map{ $hashQC{$sample}{$_} } @heads_qc;
        print QC (join "\t", @values) . "\n";
    }

    my $qct = "$output_dir/qc_target.txt";
    my @targets = sort keys %hashQCT;

    open QCT, ">$qct";
    print QCT "Sample\t" . (join "\t", @targets) . "\n";
    foreach my $sample(@samples)
    {
        my @counts = map{ $hashQCT{$_}{$sample} } @targets;
        print QCT "$sample\t" . (join "\t", @counts) . "\n";
    }
    close QCT;

    print "QC: $qc\n";
    print "QC Target: $qct\n";
}

 

sub sub_run{
    my $sample   = shift @_;
    my $fastq_r1 = "$fastq_dir/$sample\_R1.fastq.gz";
    my $fastq_r2 = "$fastq_dir/$sample\_R2.fastq.gz";

    my $sample_dir = "$mapping_dir/$sample";
    make_dir($sample_dir);

    # Merge
    print "1.Merge Fastq : $sample\n";
    system "$SOFT_FLASH -t 4 --allow-outies -m 15 -M 150 -x 0.1 -d $sample_dir -o $sample $fastq_r1 $fastq_r2 > $sample_dir/flash.log.txt 2>&1"; # M提高，能一定程度上提高短片段的合并成功数量，提高体系富集效率大约1%，但对总体来说几乎无影响

    # filter
    print "2.Filter Fastq : $sample \n";
    my $fastq_merge_mapped    = "$sample_dir/$sample.merge.mapped.fastq";
    my $merge_reads_belong    = "$sample_dir/$sample.merge.reads.belong";
    my %hashTargetFasta  = read_fasta($target_fasta);
    merge_fastq_analysis($target_fasta_db, \%hashTargetFasta, $sample_dir, $sample); # 合并成功reads过滤

    my $final_fastq = "$sample_dir/$sample.fastq";
    system("cat $fastq_merge_mapped > $final_fastq");

    # fastq比对
    print "3.blast final Fastq : $sample \n";
    my $final_blast = "$sample_dir/$sample.blast.gz";
    blast_final($final_fastq, $merge_reads_belong, $target_fasta_db, $final_blast, $sample_dir);

    # str reads count 
    print "4.get str count : $sample \n";
    str_count($merge_reads_belong, $final_blast, $sample_dir, $sample);

}


sub str_count{
    my $reads_belong = shift @_;
    my $blast_out    = shift @_;
    my $sample_dir   = shift @_;
    my $sample       = shift @_;


    my %hashMotif       = read_motif($target_motif);
    my %hashReadsBelong = read_read_belong($reads_belong, ""); # 经过上一步的质控，提取reads属于哪一个参考片段
    
    ### 参照序列STR信息，只需要起始与终止位置
    my %hashResult;   # 记录str信息
    my %hashBlastSeq; # 记录比对序列，用于后期核查
    my %hashFinishReads; # reads 已经处理完毕
    
    open my $fh, "gzip -cd $blast_out |";
    my $blast_in = Bio::SearchIO->new(-fh => $fh, -format => 'blast');
    while( my $r = $blast_in->next_result ) {
        while( my $h = $r->next_hit ) {  
            while( my $hsp = $h->next_hsp ) {  
                my $hit_name   = $h->name;       # 参照序列名称
                my $query_name = $r->query_name ;# reads序列名称
                   $hit_name   =~ s/\s+//g;
                   $query_name =~ s/\s+//g;
                last if(exists $hashFinishReads{$query_name} ); # 该reads已经分析完毕
                last if($hit_name ne $hashReadsBelong{$query_name}); # 该reads质控时对应的参考片段不是当前的参考片段

                my $hit_string   = uc($hsp->hit_string);  # 比对序列,参考序列
                my $qerry_string = uc($hsp->query_string); # 比对序列，reads
                my ($hit_start, $hit_end) = $hsp->range('hit') ; # 参考序列比对范围，永远是从小到大的顺序，即使是Minus
                my $fangxiang = $hsp->strand('hit') ; # 获取参照序列的方向 
                # 如果参照序列是反向匹配的，则需要反转,必须保证参考序列是正向的
                if($fangxiang == -1){ 
                   $qerry_string = reverse $qerry_string;
                   $qerry_string =~ tr/ATCGatcg/TAGCtagc/;
                   $hit_string   = reverse $hit_string; 
                   $hit_string   =~ tr/ATCGatcg/TAGCtagc/;                                        
                }

                # 判定该reads的STR序列
                foreach my $aim_str(keys %{$hashMotif{$hit_name}})# 对需要在该片段上判定的目标STR循环识别
                {   
                    # 获取参考序列上STR情况
                    my $motif         = $hashMotif{$hit_name}{$aim_str}{'motif'};
                    my $motif_length  = $hashMotif{$hit_name}{$aim_str}{'motif_length'};
                    my $str_seq_start = $hashMotif{$hit_name}{$aim_str}{'str_seq_start'};
                    my $str_seq_end   = $hashMotif{$hit_name}{$aim_str}{'str_seq_end'};
                    
                    # 比对区域没有包含目标STR区域
                    next if($str_seq_start < $hit_start or $str_seq_end > $hit_end);
                    
                    # 获取参照STR两端各扩展3bp，因为当indel出现时，可能会造成原始序列STR变长
                    my $extend = 3;
                    $str_seq_start = $str_seq_start - $extend;
                    $str_seq_end   = $str_seq_end + $extend;

                    # (1) 第一步，循环遍历hit_string, 根据参考位置，提取目标STR的序列
                    my $get_seq     = ""; # 实际reads在目标STR区域的序列
                    my $before_seq  = ""; # 目标区域前面的序列
                    my $after_seq   = ""; # 目标区域后面的序列
                    my $target_pos  = $hit_start; # 参考序列上的绝对位置

                    for(my $hit_pos = 0; $hit_pos < length($hit_string); $hit_pos++)
                    {
                        my $hit_base   = substr($hit_string, $hit_pos, 1);
                        my $query_base = substr($qerry_string, $hit_pos, 1);
                        $get_seq    = $get_seq.$query_base    if($target_pos >= $str_seq_start and $target_pos <= $str_seq_end);
                        $before_seq = $before_seq.$query_base if($target_pos < $str_seq_start);
                        $after_seq  = $after_seq.$query_base  if($target_pos > $str_seq_end);
                        $target_pos++ if($hit_base=~/[ATCG]/i);# 只有当hit_base是碱基的时候，绝对位置向后移动一位，因为存在插入缺失，即：字符“-”
                    }

                    # (2) 第二步，判定提取的序列上STR信息
                    # 去掉插入缺失字符
                    $get_seq    =~ s/[-]//g;
                    $before_seq =~ s/[-]//g;
                    $after_seq  =~ s/[-]//g;
                    my ($str_seq, $start, $end, $str_seq_length) = get_seq_str($query_name, $get_seq, $motif, 1, 'Force'); # 获取该片段最大的STR序列,最少重复1次

                    # (3) 第三步，从目标STR区域向两侧扩展，每次左右各扩展time * motifLength 长度，因为有可能前面的extend扩展的不足
                    my $time     = 0; # 左右扩展次数，
                    my $continue = 1; # 控制是否继续循环
                    while($continue){
                       $time++;
                       my $cut_length    = $motif_length * $time; # 本次前后扩展长度
                       my $before_cutseq = (length($before_seq) <= $cut_length) ? $before_seq : substr($before_seq, length($before_seq) - $cut_length, $cut_length); # getSeq前面扩展的序列
                       my $after_cutseq  = (length($after_seq) <= $cut_length)  ? $after_seq : substr($after_seq, 0, $cut_length); # getSeq后面扩展的序列
                       my $combine_seq   = $before_cutseq.$get_seq.$after_cutseq;                  
                       my ($str_seq_tmp, $start_tmp, $end_tmp, $str_seq_length_tmp) = get_seq_str($query_name, $combine_seq, $motif, 1, 'Force'); # 获取扩展后该片段最大的STR序列

                       # 扩展后变长
                       if($str_seq_length < $str_seq_length_tmp){
                          $str_seq_length = $str_seq_length_tmp;
                          $str_seq        = $str_seq_tmp;
                       }else{ # 扩展并不会使STR变长，离开循环
                          $continue = 0;
                       }
                    }
                    # (4) 第四步，结果保存
                    if($str_seq =~ /[agtc]/ig){
                        $hashFinishReads{$query_name}++; # 记住当前reads名称，表明其已分析完毕，避免同源导致多次处理
                        $hashResult{$hit_name}{$aim_str}{$str_seq}++; # 参考片段上记录一个str
                        $hashBlastSeq{$hit_name}{$aim_str}{$str_seq}{$query_name}{'querry'} = $qerry_string; # 详细比对序列，用于后续输出，检查
                        $hashBlastSeq{$hit_name}{$aim_str}{$str_seq}{$query_name}{'hit'}    = $hit_string;
                    }
                }                   
            }
        }
    }
    close $fh;


    # 结果输出
    open STRCOUNT, ">$sample_dir/$sample.str_count.txt"; # 每种STR的reads数量
    open STRCOUNT_DETAIL, "| gzip > $sample_dir/$sample.str_count.detail.xls.gz"; # 每条reads分型的STR详情，用于检查是否分型有问题。压缩输出，减少空间浪费
    open READ_LENGTH, ">$sample_dir/$sample.reads.length.count.txt"; # 每个片段上的reads每种长度reads的数量统计

    print STRCOUNT "target\taim_str\tstr_referrence\tstr_identify\treads_count\n";
    print STRCOUNT_DETAIL "target\taim_str\tstr_identify\treads_count\treads_name\tseq\n";
    print READ_LENGTH "target\taim_str\tlength\tcount\n";
   
    foreach my $target(sort keys %hashResult)
    {
        foreach my $aim_str(sort keys %{$hashResult{$target}})
        {   
            my %hashLength;
            foreach my $str_identify(sort keys %{$hashResult{$target}{$aim_str}})
            {
                print STRCOUNT "$target\t$aim_str\t$hashMotif{$target}{$aim_str}{'str_seq'}\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\n";
                foreach my $querryname(sort keys %{$hashBlastSeq{$target}{$aim_str}{$str_identify}}){
                    print STRCOUNT_DETAIL "$target\t$aim_str\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\t$querryname\_QUERRY\t$hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'querry'}\n";
                    print STRCOUNT_DETAIL "$target\t$aim_str\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\t$querryname\_HIT\t$hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'hit'}\n";            
                    
                    # 记录序列长度
                    my $seq = $hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'querry'};
                       $seq =~ s/-//g;
                    my $length = length($seq);
                    $hashLength{$length}++;
                }
            }
            # 输出长度信息
            my $motif_length = $hashMotif{$target}{$aim_str}{'motif_length'};
            my @lengths      = sort {$a <=> $b} keys %hashLength;
            my $min_length   = $lengths[0] - $motif_length;
            my $max_length   = $lengths[$#lengths] + $motif_length;
            $hashLength{$min_length}++;
            $hashLength{$max_length}++;# 最大、最小长度上下扩展1个motif，使图形更美观
            map{ print READ_LENGTH "$target\t$aim_str\t$_\t$hashLength{$_}\n"; } sort {$a <=> $b} keys %hashLength; # 输出reads长度的数量            
        }
    }
    close STRCOUNT;
    close STRCOUNT_DETAIL;
    close READ_LENGTH;
}



# 最终比对，按片段拆分执行，以节省CPU/存储消耗
sub blast_final{
    my $final_fastq        = shift @_; # 要处理的fastq
    my $reads_belong_merge = shift @_; # 每条reads所属target
    my $target_fasta_db    = shift @_; # target数据库
    my $final_blast        = shift @_; # 最终输出文件
    my $output_dir         = shift @_; 

    # 临时目录，用后即删
    my $tmp_dir = "$output_dir/fastq_split_blast"; 
    make_dir($tmp_dir);


    my %hashTargetFasta = read_fasta($target_fasta_db); # 参考序列名称获取 {'target_name'}= seq
    my %hashReadsBelong = read_read_belong($reads_belong_merge, ""); # reads所属片段信息获取 {'reads_name'} = target_name

    # 创建每个片段fa句柄
    my %hashHandle; 
    foreach my $target_name(keys %hashTargetFasta)
    {
        my ($target_name_prefix) = split /[|]/, $target_name;  # 取片段前缀作为文件名
        my $tmp_fa    = "$tmp_dir/$target_name_prefix.fa";         # 属于当前片段的reads序列
        my $tmp_limit = "$tmp_dir/$target_name_prefix.limit.txt";  # 当前片段的完整名称，用于blastn的输入
        system("echo '$target_name\n' > $tmp_limit");
        open $hashHandle{$target_name}, ">$tmp_fa";
    }

    # final_fastq 拆分到每个target片段句柄里
    open FASTQ, $final_fastq;
    while(my $line1 = <FASTQ>)
    {
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;

           $line1=~ s/[\r\n]//g;
           $line1=(split /\s/,$line1)[0];
           $line1=~ s/^\@//g;
        my $handle = $hashHandle{$hashReadsBelong{$line1}};
        print $handle ">$line1\n$line2";
    }
    close FASTQ;  

    # 关闭句柄
    foreach my $target_name(keys %hashHandle)
    {
        close $hashHandle{$target_name};
    }
    %hashReadsBelong = (); # 清空，释放内存

    # 比对
    my @blast_results; # 所有比对结果文件
    foreach my $target_name(keys %hashTargetFasta)
    {
        my ($target_name_prefix) = split /[\|]/, $target_name;
        my $tmp_fa    = "$tmp_dir/$target_name_prefix.fa";         # 属于当前片段的reads序列
        my $tmp_limit = "$tmp_dir/$target_name_prefix.limit.txt";  # 当前片段的名称，用于blastn的输入
        my $tmp_blast = "$tmp_dir/$target_name_prefix.blast.gz";
        next if(is_file_ok($tmp_fa) == 0); # 没有数据，不必再做
        system("$SOFT_BLASTN -task blastn -query $tmp_fa -evalue 0.00001 -db $target_fasta_db -num_threads 4 -dust no -seqidlist $tmp_limit | gzip > $tmp_blast"); # 压缩，减少存储消耗
        push @blast_results, $tmp_blast;
    }   

    # 合并
    system("cat @blast_results > $final_blast");
    # 清空中间文件
    system("rm -r $tmp_dir");
}

# 提取reads所属信息
sub read_read_belong{
    my $reads_belong = shift @_;
    my $type         = shift @_;

    my %hashReadsBelong;
    open BELONG, $reads_belong;
    while(<BELONG>)
    {
        $_=~s/[\r\n]//g;
        my ($reads_name, $target_name) = split /\t/, $_;        
        if($type eq "FASTQ_TYPE") # 原始fastq命名格式
        {
            $reads_name=~s/:MERGE$//;
            $reads_name=~s/:R\d$//;
            $reads_name="\@$reads_name";
        }
        elsif($type eq 'READS_NAME_TYPE')
        {
            $reads_name=~s/:MERGE$//;
            $reads_name=~s/:R\d$//;            
        }
        $hashReadsBelong{$reads_name} = $target_name;
    }
    close BELONG;
    return %hashReadsBelong;
}


# 合并成功reads过滤
sub merge_fastq_analysis{
    my $target_fasta     = shift @_;
    my $hashTargetFasta  = shift @_;
    my $sample_dir       = shift @_;
    my $sample           = shift @_;

    # fastq 转 fasta
    my $merge_fastq = "$sample_dir/$sample.extendedFrags.fastq";
    my $merge_fasta = "$sample_dir/$sample.extendedFrags.fa";   
    system "$SOFT_FASTQ_TO_FASTA -Q 33 -n -i $merge_fastq -o $merge_fasta";

    # blast
    my $merge_blast = "$sample_dir/$sample.extendedFrags.blast";
    system("$SOFT_BLASTN -task blastn -query $merge_fasta -outfmt 6 -evalue 0.00001 -db $target_fasta -out $merge_blast -num_threads 4 -max_target_seqs 2 -dust no");
    
    # 筛选
    my %hashmergeFastq = read_fastq($merge_fastq);
    my %hashmergeBlast = read_blast(\%hashmergeFastq, $hashTargetFasta, $merge_blast, 'Merge');
    
    my $fastq_mapped    = "$sample_dir/$sample.merge.mapped.fastq";
    my $fastq_un_mapped = "$sample_dir/$sample.merge.unmapped.fastq";
    my $reads_belong    = "$sample_dir/$sample.merge.reads.belong";
    open MAPPED, ">$fastq_mapped";
    open UNMAPPED, ">$fastq_un_mapped";
    open BELONG, ">$reads_belong";
 
    my %hashFilterOK; # qualified reads count in each target
    foreach my $reads_name(keys %hashmergeFastq){
        if(exists $hashmergeBlast{$reads_name}){
            print MAPPED "\@$reads_name:MERGE\n";
            print MAPPED "$hashmergeFastq{$reads_name}{2}\n";
            print MAPPED "$hashmergeFastq{$reads_name}{3}\n";
            print MAPPED "$hashmergeFastq{$reads_name}{4}\n";
            print BELONG "$reads_name:MERGE\t$hashmergeBlast{$reads_name}\n";

            $hashFilterOK{$hashmergeBlast{$reads_name}}++;

        }else{
            print UNMAPPED "\@$reads_name\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{2}\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{3}\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{4}\n";
        }
    }
    close MAPPED;
    close UNMAPPED;
    close BELONG;

    my $statistics_target = "$sample_dir/$sample.filter_ok.txt";
    open STAT, ">$statistics_target";
    foreach my $target_name(sort keys %$hashTargetFasta)
    {
        my $filter_ok = exists $hashFilterOK{$target_name} ? $hashFilterOK{$target_name} : 0;
        print STAT "$target_name\t$filter_ok\n";
    }
    close STAT;
}

# 读取fastq文件
sub read_fastq{
    my $fastq = shift @_;
    my %hashFastq;
    open FASTQ,$fastq;
    while(my $line1 = <FASTQ>){
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;
        $line1=~ s/[\r\n]//g;
        $line2=~ s/[\r\n]//g;
        $line3=~ s/[\r\n]//g;
        $line4=~ s/[\r\n]//g;
        my $name=(split /\s/,$line1)[0];
           $name=~ s/^\@//;
        $hashFastq{$name}{1} = $line1;
        $hashFastq{$name}{2} = $line2;
        $hashFastq{$name}{3} = $line3;
        $hashFastq{$name}{4} = $line4;
    }
    close FASTQ;
    return %hashFastq;
}

# 读取blast结果
sub read_blast{
    my $hashFastq       = shift @_; 
    my $hashTargetFasta = shift @_;
    my $blast           = shift @_;
    my $type            = shift @_;###single表示对merge的reads进行比对。pair表示对没有merge上的R1或R2的reads进行比对。
    
    # 读取blast比对信息，并进行基础过滤
    my %hashBlast;
    open BLAST,$blast;
    while(<BLAST>){
        $_=~ s/[\r\n]//g;
        my ($reads_name, $target_name, $identify_perc, $identify_len, $mismatch, $gap, $reads_start, $reads_end, $target_start, $target_end, $evalue, $score) = split /\t/, $_;
        my @array=split(/\t/,$_);
        my $reads_length  = length( $hashFastq->{$reads_name}{2} );
        my $target_length = length( $hashTargetFasta->{$target_name} );
        ($target_start, $target_end) = ($target_end, $target_start) if($target_end < $target_start);
 
        my $judge_value = $identify_len; # 默认用匹配长度作为后续候选的筛选规则
        #
        # 情况一：Reads覆盖整个目标区域
        # 条件，覆盖超过80%的区域
        # 条件，覆盖起始位置小于开头加5，结束位置大于末尾减5
        #
        if($type eq 'Merge' and $identify_len > $target_length * 0.8 and $target_start <= 5 and $target_end > $target_length - 5){
            $hashBlast{$reads_name}{$judge_value} = $target_name;
        }
    }
    close BLAST;

    # 进一步筛选
    my %hashBest=();

    foreach my $title(keys %hashBlast){
        foreach my $len (sort {$b <=> $a} keys %{$hashBlast{$title}}){
            $hashBest{$title} = $hashBlast{$title}{$len};
            last;
        }
    }       

    return %hashBest;
}



sub check_fastq{
    print "[process] Check fastq \n";
    foreach my $sample(@samples)
    {
        my $fastq_r1 = "$fastq_dir/$sample\_R1.fastq.gz";
        my $fastq_r2 = "$fastq_dir/$sample\_R2.fastq.gz";
        die "[Error] lost fastq: $fastq_r1 $fastq_r2\n" if(is_file_ok($fastq_r1, $fastq_r2) == 0);
    }
}

sub check_target_motif{
    my $hashConfig = shift @_;

    print "[process] Check target motif \n";
    my %hashTarget   = read_fasta($target_fasta);

    # （1）提取motif详细信息
    my %hashResult; # 记录该片段需要查找的所有motif，并判断是否重复
    my $count = 0;
    foreach my $target_full_name(sort keys %hashTarget)
    {
        my $target_seq        = $hashTarget{$target_full_name};
        my $target_seq_length = length($hashTarget{$target_full_name});# 参考序列长度
        my ($target_name, $motif_info_list) = split /\|/, $target_full_name, 2;
        
        foreach my $motif_info(split /;/, $motif_info_list)
        { 
            $count++;
            my $motif = "";
            my ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = ('', '', '', '', '');

            # 自定义区域信息
            if($motif_info =~ /^([ATCG]+)\((\d+)-(\d+)\)$/) 
            {
                $motif     = uc($1);
                my $aim_start = $2;
                my $aim_end   = $3;
                die "region must be in ascending order! $target_full_name\n" if($aim_start > $aim_end);
                ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = get_aim_str($target_full_name, $target_seq, $motif, $aim_start, $aim_end); # 获取指定目标区域motif的STR 
            }
            else
            {
                die "motif format info error, please read readme.txt : $target_full_name -> $motif_info\n";
            }
            
            # 检测是否重复
            my $title = "$motif,$start,$end";
            die "[Error] duplicate STR in $target_full_name: $motif $motif_info\n" if(exists($hashResult{$target_full_name}{$title})); # 重复STR
            
            # 信息保存
            $hashResult{$target_full_name}{$title}{'sort_order'}          = $count;      # 排序
            $hashResult{$target_full_name}{$title}{'original_info'}       = $motif_info; # 当前motif原信息
            $hashResult{$target_full_name}{$title}{'motif'}               = $motif; # motif
            $hashResult{$target_full_name}{$title}{'motif_length'}        = length($motif); # motif
            $hashResult{$target_full_name}{$title}{'str_seq_start'}       = $start; # str序列起始
            $hashResult{$target_full_name}{$title}{'str_seq_end'}         = $end;   # str序列终止
            $hashResult{$target_full_name}{$title}{'str_seq_length'}      = $str_seq_length;      # str序列长度
            $hashResult{$target_full_name}{$title}{'str_seq_motif_count'} = $str_seq_motif_count; # str序列含有的motif数量
            $hashResult{$target_full_name}{$title}{'target_seq_length'}   = $target_seq_length;   # 参考序列长度
            $hashResult{$target_full_name}{$title}{'str_mark'}            = "$motif($str_seq_motif_count)";   # str序列特殊表示形式
            $hashResult{$target_full_name}{$title}{'str_seq'}             = $motif x $str_seq_motif_count;   # str序列
            $hashResult{$target_full_name}{$title}{'target_seq'}          = $target_seq;   # 参考序列
            $hashResult{$target_full_name}{$title}{'target'}              = $target_full_name;   # 片段名称
        }
    }

    # (2) 输出到文件
    my @heads = ('target', 
        'motif', 
        'motif_length', 
        'str_mark',
        'str_seq_start', 
        'str_seq_end', 
        'str_seq_length', 
        'str_seq_motif_count', 
        'target_seq_length', 
        'str_seq',
        'target_seq',
        );
    
    open OUT, ">$output_dir/target.fasta.STR.txt";
    print OUT (join "\t", @heads) . "\n";
    foreach my $target_full_name(sort keys %hashResult)
    {
        foreach my $title(sort {$hashResult{$target_full_name}{$a}{'sort_order'} <=> $hashResult{$target_full_name}{$b}{'sort_order'}} keys %{$hashResult{$target_full_name}})
        {
            my @datas = map{ $hashResult{$target_full_name}{$title}{$_} } @heads;
            print OUT (join "\t", @datas) . "\n";
        }
    }
    close OUT;

    return "$output_dir/target.fasta.STR.txt";
}



# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 提取fasta序列
sub read_fasta{
    my $fasta = shift @_;
    my %hashFasta;
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}= $seq;
    }
    return %hashFasta;
}

# 获取指定目标区域motif的STR
sub get_aim_str{
    my $target_full_name = shift @_;
    my $target_seq = shift @_;
    my $motif      = shift @_;
    my $aim_start  = shift @_;
    my $aim_end    = shift @_;

    my $start          = $aim_start;
    my $end            = $aim_end;
    my $motif_length   = length($motif);
    my $str_seq        = substr($target_seq, $aim_start - 1, $aim_end - $aim_start + 1);
    my $str_seq_length = length($str_seq);

    my $str_seq_motif_count      = $str_seq_length / $motif_length; # STR序列中包含的motif数量
    my ($target_seq_motif_count) = ($target_seq =~ s/$motif/$motif/ig); # 参考序列中包含的motif数量
    my ($str_seq_motif_count_R)  = ($str_seq =~ s/$motif/$motif/ig); # STR序列中包含的motif数量, 根据实际获取的STR序列计算
    
    die "[Error] str length can not be divided by motif len!>$target_full_name\n$target_seq\n" if($str_seq_length % $motif_length!=0);# 非整数倍
    die "[Error] aimed STR Region Error(motif count in region $aim_start-$aim_end not equal regionLength/motifLength):>$target_full_name\n$target_seq\n" if($str_seq_motif_count != $str_seq_motif_count_R);# 指定区域的motif数量跟实际不符
    die "[Error] no motif exists in Seq :>$target_full_name\n$target_seq\n" if($target_seq_motif_count == 0 or $str_seq_motif_count_R == 0);# 参考序列中不含有motif
        
    return ($str_seq, $aim_start, $aim_end, $str_seq_length, $str_seq_motif_count);  
}


# blast+ 建库
sub build_blast_index{
    print "[process] build blast DB\n";

    my $target_fasta_db_dir = "$output_dir/target_fasta_db";
    my $target_fasta_db     = "$target_fasta_db_dir/seq.fa";
    make_dir($target_fasta_db_dir);
    system("cp $target_fasta $target_fasta_db");
    
    system("$SOFT_MAKEBLASTDB -in $target_fasta_db -dbtype nucl  -parse_seqids > $target_fasta_db_dir/log.txt 2>&1");
    return  $target_fasta_db;  
}


# 读取片段motif信息
sub read_motif{
    my $motif_file = shift @_;
    my $model      = (exists $_[0]) ? $_[0] : "full_name";
    
    die "Lost $motif_file\n" if(is_file_ok($motif_file) == 0);
    my %hashMotif;
    open MOTIF, $motif_file;
    my $line1 = <MOTIF>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    while(<MOTIF>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        
        my %hashTmp = map{ ($heads[$_], $datas[$_])} (0..$#heads);

        my $target        = $hashTmp{'target'};
        my $motif         = $hashTmp{'motif'};
        my $str_seq_start = $hashTmp{'str_seq_start'};
        my $str_seq_end   = $hashTmp{'str_seq_end'};
        ($target) = split /[\|]/, $target if($model eq 'clean_name');

        my $aim_str       = "$motif,$str_seq_start,$str_seq_end";
        map{ $hashMotif{$target}{$aim_str}{$_} = $hashTmp{$_} } keys %hashTmp;              
    }
    close MOTIF;

    return %hashMotif;
}


# 获取序列中motif对应的最长STR
sub get_seq_str{
    my $target_full_name = shift @_;
    my $target_seq = shift @_;
    my $motif      = shift @_;
    my $minrep     = shift @_; # 最小重复次数
    my $type       = exists $_[0] ? $_[0] : 'Check';

    my $motif_length = length($motif);
       $minrep       = $minrep - 1;
    my $regexp       = "(($motif)\\2{".$minrep.",})";

    # 匹配最长str
    my $bestlength = 0;
    my $getinfo    = "";
    while($target_seq =~ /$regexp/ig)
    {
        my $str_seq        = $1; 
        my $str_seq_length = length($str_seq);
        my $end            = pos($target_seq);
        my $start          = $end - $str_seq_length + 1;
        my $str_seq_motif_count = $str_seq_length / $motif_length; # STR序列中包含的motif数量
 
        if($str_seq_length > $bestlength)
        {
           $bestlength = $str_seq_length;
           $getinfo    = "$str_seq|$start|$end|$str_seq_length|$str_seq_motif_count";
        }
        pos($target_seq) = $start; # 重新设置检测起点，以单碱基步频移动  
    }
    my ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = ('', '', '', 0, 0);
       ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = split /\|/, $getinfo if($getinfo =~ /\w/);

    my ($target_seq_motif_count) = ($target_seq =~ s/$motif/$motif/ig); # 参考序列中包含的motif数量
    die "no motif '$motif' exists in Seq >$target_full_name\n$target_seq \n" if($type eq 'Check' and ($target_seq_motif_count == 0 or $bestlength == 0) ); # 参考序列中不含有motif
    my $str_seq_motif_perc = ($target_seq_motif_count == 0) ? 0 : sprintf "%0.4f", $str_seq_motif_count / $target_seq_motif_count;# str序列中motif占参考序列中motif数量的比例

    return ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count, $str_seq_motif_perc);
}


sub get_sample_list{
    my $fastq_dir = shift @_;

    my @samples;
    opendir FASTQ_DIR, $fastq_dir;
    while(my $file = readdir FASTQ_DIR)
    {
        my ($sample) = $file =~ /(.*)_R1.fastq.gz/;
        push @samples, $sample if(defined $sample);
    }
    close FASTQ_DIR;
    my $sample_list = join ",", @samples;
    return $sample_list;
}