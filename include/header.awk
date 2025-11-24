@namespace "header"

function run() {
help="\
    Update the VCF header with contigs from a new reference. \n\
    Same as doing `bcftools reheader -f file.fa.fai file.vcf.gz`. \n\
    Requires tabix in PATH. The VCF must be bgzipped and tabixed. \n\
    Arguments are supplied positionally (first VCF, then FASTA) and/or using named options. \n\
    EXAMPLE: \n\tulitka header file.vcf.gz file.fa > new_header.vcf"
  arg::add_argument("v", "vcf", 0, "VCF file to be lifted over")
  arg::add_argument("f", "fasta", 0, "new reference FASTA file")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  if ( system("tabix --version > /dev/null") != 0 ) {
       print "Can't run tabix"; exit 1
  }
  if ( "vcf" in arg::args ) {
    vcf = arg::args["vcf"]
  } else {
    vcf = arg::nonargs[1]
  }
  if ( "fasta" in arg::args ) {
    fasta = arg::args["fasta"]
  } else {
    fasta = arg::nonargs[ length(arg::nonargs) ]
  }
  if ( fasta == "" ) {
    print "Warning: FASTA file not provided, header will not be updated" > "/dev/stderr"
  }
  ulitka::assert( vcf != "", "VCF file required" )
  exit main(vcf, fasta)
}

function main(vcf, fasta,    fai, firstcontig, cmd, m) {
  if (fasta != "") {
    fai=fasta".fai"
    if ( system("test -f " fai) != 0 ) {
      faidx::faindex(fasta)
    }
    firstcontig=1
  }

  cmd="tabix -H " vcf
  while( cmd | getline > 0 ) {
    if (match($0, /^#(#contig|CHROM)/, m)) {
      if ( firstcontig ) {
        firstcontig=0
        printcontigs(fai)
      }
      if (fasta != "" && m[0]=="##contig") {
        continue
      }
    } 
    print
  }
  close(cmd)
}

function printcontigs(fai) {
  while( getline < fai > 0 ) {
    print "##contig=<ID="$1",length="$2">"
  }
  close(fai)
}

