# SCycDB
# SCycDB: A Curated Functional Gene Database for Metagenomic Profiling of Sulfur Cycling Pathways
The sulfur (S) cycle driven by microorganisms is an important biogeochemical cycling process of the Earth's biosphere, and it is usually coupled with carbon, nitrogen and metal cycling in natural ecosystems. Shotgun metagenome sequencing has opened a new avenue to advance our understanding of S cycling microbial communities. However, accurate metagenomic profiling of S cycling microbial communities remains technically challenging, mainly due to low coverage of S cycling genes/pathways, difficulties in distinguishing homologous genes and a long research time on publicly available orthology databases. It is essential to develop a comprehensive and accurate database for characterizing S cycling microbial communities in metagenomic studies. To solve those problems, we constructed a manually curated sulfur cycling database (SCycDB) for metagenome sequencing data analysis of S cycling microbial communities in the environment. 

The developed SCycDB contains 207 gene families and 585,055 representative sequences affiliated with 52 phyla and 2684 genera of bacteria/archaea, and 20,761 homologous orthology groups were also included to reduce false positive sequence assignments. 

Three files are generated for SCycDB:

<b>1. SCycDB_2020Mar.zip</b>: fasta format representative sequences obtained by clustering curated sequences at 100% sequence identity. This file can be used for "BLAST" searching SCycDB genes in shotgun metagenomes.

<b>2. id2gene.2020Mar.map</b>: a mapping file that maps sequence IDs to gene names, only sequences belonging to SCycDB gene families are included. Sequences for SCycDB homologs are not included. This file is used to generate SCycDB profiles from BLAST-like results against the SCycDB database. 

<b>3. SCycDB_FunctionProfiler.PL</b>: a perl script for functional profiling of S cycling genes.

<b>4. SCycDB_TaxonomyProfiler.PL</b>: a perl script for taxonomical profiling of S cycling microbial communities.

<b>DOWNLOAD/INSTALLATION</b>

git clone https://github.com/qichao1984/SCycDB.git

<b>Dependencies and Tools</b>

<i>Perl modules that can be easily installed via cpan:</i>
<p>List::Util</p>
<p>Getopt::Long</p>
<i>Dependencies for SCycDB_FunctionProfiler.PL, currently supported database searching tools are: </i>
<p>usearch: https://www.drive5.com/usearch/download.html
<p>diamond: https://github.com/bbuchfink/diamond/releases
<p>blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz</p>
<i>Dependencies for SCycDB_TaxonomyProfiler.PL:</i>
<p>seqtk: https://github.com/lh3/seqtk.git
<p>kraken2: https://github.com/DerrickWood/kraken2.git

<b>USAGE</b>

Before getting started, please modify both scripts (SCycDB_FunctionProfiler.PL, SCycDB_TaxonomyProfiler.PL) at lines 6-16 to specify the locations of third party tools and their parameters. If the tools are in the system path, no revision is needed for the path of these tools. By default, basic parameters are used for these tools. Users are encouraged to make revisions in cases of short reads and/or expecting more strict/relaxed results. We also encourage users to develop useful implementations based on SCycDB.


<b>Example for using SCycDB_FunctionProfiler.PL:</b>

perl SCycDB_FunctionProfiler.PL -d \<workdir\> -m \<diamond|usearch|blast\> -f \<filetype\> -s \<seqtype\> -si \<sample size info file\> -rs \<random sampling size\> -o \<outfile\>
  
Detailed explanations: 

-d : specify the directory where your fasta/fastq (or gzipped) files are located. 

-m : specify the database searching program you plan to use, currently diamond, usearch and blast are supported. 

-f : specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s : sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

-o : the output file for N cycle gene profiles.   


<b>Example for using SCycDB_TaxonomyProfiler.PL:</b>

perl SCycDB_TaxonomyProfiler.PL -d \<workdir\> -m \<diamond|usearch|blast\> -f \<filetype\> -s \<seqtype\> -si \<sample size info file\> -rs \<random sampling size\> 
  
Detailed explanations: 

-d : specify the directory where your fasta/fastq (or gzipped) files are located. 

-m : specify the database searching program you plan to use, currently diamond, usearch and blast are supported. 

-f : specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s : sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have, note that no file extensions should be included here.

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

  
