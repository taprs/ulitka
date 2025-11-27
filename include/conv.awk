@namespace "convert"

function run() {
help="\
    Convert genomic coordinates using a chain file. \n\
    Accepts BED format (0-based) and chr:start-end notation (1-based); output is in BED format. \n\
    Arguments are supplied positionally (first chain, then region) and/or using named options. \n\
    Multiple input regions can be concatenated. \n\
    EXAMPLE: \n\tulitka convert file.chain file.bed > file_liftover.bed \n\
    GNU PARALLEL EXAMPLES: \n\t parallel --pipepart -a file.chain --recend '\\n\\n' ulitka convert {} file.bed > file_liftover.bed \n\
    \t parallel --pipepart -a file.bed ulitka convert file.chain {} > file_liftover.bed"
  arg::add_argument("c", "chain", 0, "chain file")
  arg::add_argument("r", "reg", 0, "chr:start-end or BED file")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  if ( "chain" in arg::args ) { 
    chain=arg::args["chain"] 
  } else {
    chain=arg::nonargs[1]
    delete arg::nonargs[1]
  }
  if ( "reg" in arg::args ) { regs[arg::args["reg"]]=nregs=1 }
  for(i in arg::nonargs) {
    regs[arg::nonargs[i]]=1
    nregs++
  }

  ulitka::assert( chain != "", "chain file required" )
  ulitka::assert( nregs >  0, "region or region file required" )
  exit main()
}

function main() {
  for (reg in regs) {
    while ( getline < reg > 0 ) {
      isfile=1
      print convert_coords(chain, $1, $2, $3)
    }
    if ( ! isfile ) {
      faidx::parse_region(reg)
      print convert_coords(chain, faidx::chr, faidx::start, faidx::end)
    }
  }
}

function convert_coords(chainfile, chr, pos1, pos2,
                        pos, npos, gotpos, newpos, id) {

  ulitka::assert( pos1!="", "please provide at least one numeric coordinate: "chr", "pos1", "pos2)
  pos[1]=pos1
  npos=1
  if ( pos2!="" ) { pos[2]=pos2; npos++ }
  gotpos=0

  while ( getline < chainfile > 0 && gotpos < npos ) {
    if ( $1 != "chain" || $3!=chr || ( $6 > pos1 && $6 > pos2 ) || ( $7 < pos1 && $7 < pos2 ) ) { continue }
    if ( chain::parse_next_chain(chainfile) != 0 ) { 
      print "Error: invalid chain "$0 > "/dev/stderr"
      exit 1
    }
    for (id=0;id<=chain::idx;id++) {
      for ( p in pos ) {
        if ( pos[p] >= chain::start[id] && pos[p] < chain::end[id] ) {
          gotpos++
          newpos[p]=pos[p]+chain::delta[id]
          if ( !samestrand ) { 
            newpos[p]-=2*(pos[p]-chain::end[id]+1)
          }
          delete pos[p]
        }
      }
      if ( gotpos >= npos ) { break }
    }
  }
  close(chainfile)
  if ( gotpos < npos ) {
    print "Warning: position "chr":"(1 in pos ? pos[1]+1 : pos[2])" not covered by the chain file" > "/dev/stderr"
    return ""
  }
  if ( npos == 1 ) {
    return chain::qchr"\t"newpos[1]"\t"newpos[1]
  } else if ( newpos[2] >= newpos[1] ) {
    return chain::qchr"\t"newpos[1]"\t"newpos[2]
  } else {
    return chain::qchr"\t"newpos[2]"\t"newpos[1]
  }
}

