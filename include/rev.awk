@namespace "rev"

function run() {
  help="\
      Reverse chain file, i.e. swap query and target genomes. \n\
      The output represents the opposite liftover operation. \n\
      Arguments are supplied positionally and/or using named options. \n\
      Multiple inputs are concatenated. \n\
      EXAMPLE: \n\tulitka rev file.chain > file_rev.chain"
  arg::add_argument("c", "chain", 0, "chain file")
  arg::parse_args(2, help) 
  arg::parse_nonargs()

  if ( "chain" in arg::args ) { chains[arg::args["chain"]]=nchains=1 }
  for(i in arg::nonargs) {
    chains[arg::nonargs[i]]=1
    nchains++
  }
  ulitka::assert( nchains > 0, "chain file required" )
  exit main()
}

function main() {
  for (chain in chains) {
    while ( getline < chain > 0 ) {
      if ( $1 == "chain" ) {
        print $1,$2,$8,$4,$5,$11,$12,$3,$9,$10,$6,$7,$13
        if ( $5 == $10 ) {
          while ( getline < chain > 0 ) {
            if ( $0 == "" ) { break }
            print $1" "$3" "$2
          }
        } else {
          nr=0
          delete partone
          delete parttwo
          while ( getline < chain > 0 ) {
            if ( $0 == "" ) { break }
            partone[nr]=$1
            parttwo[nr]=$3" "$2
            nr++
          }
          while ( nr-- ) { print partone[nr]"\t"parttwo[nr-1] }
        }
      }
      print ""
    }
    close(chain)
  }
}

