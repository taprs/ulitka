@namespace "easy"

function run() {
help="\
    Run all steps of `ulitka` VCF liftover, optionally parallelized. \n\
    Runs `ulitka header`, `ulitka lift`, `sort` and `bgzip`. \n\
    Requires tabix, GNU parallel (for parallel execution) and bgzip (for compression) in PATH. \n\
    The VCF must be bgzipped and tabixed. \n\
    Arguments are supplied positionally (first chain file, then VCF) and/or using named options. \n\
    EXAMPLE: \n\tulitka easy --bgzip --fasta file.fa file.chain file.vcf.gz > file_liftover.vcf.gz"
  arg::add_argument("v", "vcf", 0, "VCF file to be lifted over")
  arg::add_argument("c", "chain", 0, "chain file")
  arg::add_argument("f", "fasta", 0, "new reference FASTA file")
  arg::add_argument("j", "jobs", 0, "number of parallel jobs to run (default 1)")
  arg::add_argument("z", "bgzip", 1, "compress output with bgzip")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  OFS="\t"
  if ( system("tabix --version > /dev/null") != 0 ) {
       print "Can't run tabix" > "/dev/stderr"; exit 1
  }
  if ( system("bgzip --version > /dev/null") != 0 ) {
       print "Can't run bgzip" > "/dev/stderr"; exit 1
  }
  if ( system("parallel --version > /dev/null") != 0 ) {
     print "Warning: parallel not found in PATH, using one job" > "/dev/stderr"
     arg::args["jobs"]=1
  }
  fasta=arg::args["fasta"]
  if ( "vcf" in args ) { 
    vcf=arg::args["vcf"] 
  } else {
    vcf=arg::nonargs[ length(arg::nonargs) ]
    delete arg::nonargs[ length(arg::nonargs) ]
  }
  if ( "chain" in arg::args ) { chains[arg::args["chain"]]=nchains=1 }
  for(i in arg::nonargs) {
    chains[arg::nonargs[i]]=1
    nchains++
  }
  ulitka::assert( vcf != "", "VCF file required" )
  ulitka::assert( nchains != 0, "chain file required" )
  if ( fasta == "" ) {
    print "Warning: FASTA file not provided, header will not be updated" > "/dev/stderr"
  }
  exit main()
}

function main() {
  cmd="header.awk "vcf" "fasta
  if ( arg::args["bgzip"]==1 ) {
    cmd = cmd " | bgzip"
  }
  print "Executing command: "cmd > "/dev/stderr"
  system(cmd)
  if (arg::args["jobs"]<2) {
    cmd = "ulitka"
    if (fasta != "") {
      cmd = cmd" --fasta "fasta
    }
    for (c in chains) {
      cmd = cmd" "c
    }
    cmd = cmd " " vcf
  } else {
    cmd = "parallel"
    for (c in chains) {
      cmd = cmd" -a "c
    }
    cmd = cmd " -j "arg::args["jobs"]" --pipepart --recend '\\n\\n' ulitka"( fasta != "" ? " --fasta "fasta : "" )" - "vcf
  }
    cmd = cmd " | sort -k1,1 -k2,2n"
  if ( arg::args["bgzip"]==1 ) {
    cmd = cmd " | bgzip"
  }
    print "Executing command: "cmd > "/dev/stderr"
    system(cmd)
    print "Done!" > "/dev/stderr"
}
