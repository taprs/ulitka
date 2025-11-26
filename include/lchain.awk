@namespace "lchain"

function run() {
  help="\
  Write a (sorted) linearized chain ('lchain') indexable by target genome coordinates. \n\
  If multiple input files are given, they are concatenated. \n\
  If run on an lchain.gz file, indexes it by old ref coordinates: `tabix -s3 -b6 -e7 file.lchain.gz`. \n\
  EXAMPLES: \n\tulitka lchain --bgzip file.chain > file.lchain.gz; ulitka lchain file.lchain.gz \n\
  \tzcat file.lchain.gz | ulitka lchain -r - > file.chain \n\
  "
  arg::add_argument("r", "rev", 1, "convert lchain to chain instead")
  arg::add_argument("z", "bgzip", 1, "compress output with bgzip (no effect with --rev)")
  arg::parse_args(2, help)
  arg::parse_nonargs()
  if ( system("sort --version > /dev/null") != 0 ) {
    print "Can't run sort" > "/dev/stderr"
    exit 1
  }
  # Everything after options is input files
  for( i in arg::nonargs ) {
    files[arg::nonargs[i]]=1
    nfiles++
  }
  if ( arg::args["rev"] && arg::args["bgzip"] ) {
    print "Warning: option --bgzip has no effect with --rev" > "/dev/stderr"
  }
  ulitka::assert( nfiles > 0, "at least one input file required" )
  exit main()
}

function main() {
  if ( arg::args["rev"] ) {
    FS="\t"
    OFS=" "
    for ( f in files ) {
      while ( getline < f ) {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13
        for (nf=14;nf<NF+2;nf++) { print $nf }
      }
      close(f)
    }
  } else {
    PROCINFO["BUFFERPIPE"]=1
    RS="\n\n"
    FS="\n"
    OFS="\t"
    sort="sort -k3,3 -k6,7n"
    if (arg::args["bgzip"]) { 
      sort = sort " | bgzip"
    }
    for ( f in files ) {
      if ( f ~ /lchain\.gz$/ ) {
        if ( system("tabix -s3 -b6 -e7 "f) != 0 ) {
          print "Error: indexing failed for file "f > "/dev/stderr"
          exit 1
        }
        continue
      }
      while ( getline < f > 0 ) {
        gsub(/ /, "\t", $1) 
        print | sort
      }
      close(f)
    }
    close(sort)
  }
}
