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
  arg::add_argument("i", "index", 1, "use tabix index of the regions file")
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
  if ( "index" in arg::args ) {
    if ( system("bgzip --version > /dev/null") != 0 ) {
         print "Can't run bgzip"; exit 1
    }
    if ( system("tabix --version > /dev/null") != 0 ) {
         print "Can't run tabix"; exit 1
    }
    for (reg in regs) {
      convert_coords_index(chain, reg)
    }
  } else {
    for (reg in regs) {
      while ( getline < reg > 0 ) {
        isfile=1
        convert_coords(chain, $1, $2, $3)
      }
      if ( ! isfile ) {
        faidx::parse_region(reg)
        convert_coords(chain, faidx::chr, faidx::start, faidx::end)
      }
    }
  }
}

function slide_chain(p1, p2,    id, p) {
  ulitka::assert( p1!="", "please provide at least one numeric coordinate: "chain::tchr", "pos1", "pos2)
  pos[1]=p1
  npos=1
  if ( p2 != "" ) { pos[2]=p2; npos++ }
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
    if ( gotpos >= npos ) { 
      break 
    }
  }
  return gotpos
}

function return_bed() {
  if ( npos == 1 ) {
    print chain::qchr"\t"newpos[1]"\t"newpos[1]
  } else if ( newpos[2] >= newpos[1] ) {
    print chain::qchr"\t"newpos[1]"\t"newpos[2]
  } else {
    print chain::qchr"\t"newpos[2]"\t"newpos[1]
  }
}

function convert_coords_index(chainfile, regfile) {
  regidx = prepare_index(regfile)
  while ( getline < chainfile > 0 ) {
    if ( chain::parse_next_chain(chainfile) != 0 ) { 
      print "Error: invalid chain "$0 > "/dev/stderr"
      exit 1
    }
    cmd="tabix "regidx" "chain::tchr":"chain::t0+1"-"chain::t1
    while ( cmd | getline > 0 ) {
      gotpos=0
      if ( !slide_chain($2, $3) ) {
        print "Warning: position "chain::tchr":"(1 in pos ? pos[1]+1 : pos[2])" not covered by the chain file" > "/dev/stderr"
        continue
      }
      return_bed()
    }
    close(cmd)
  }
  close(chainfile)
}

function prepare_index(regfile) {
  if ( system("test -f "regfile".tbi") == 0 ) {
    print "Index found: "regfile".tbi" > "/dev/stderr"
    return regfile
  } else {
    print regfile".tbi not found, assuming unsorted input" > "/dev/stderr"
    if ( tmpdir=="" ) {
      tmpcmd="mktemp -d ulitkatmp.XXXXXXXX"
      if ( tmpcmd | getline tmpdir <= 0 ) {
        say("Error: failed to create temporary directory!")
        exit 1
      }
      close(tmpcmd)
      # We can't trap an exec'ed script but we can do this:
      system( "{ while kill -0 "PROCINFO["pid"]" 2>/dev/null; do sleep 1; done; rm -r "tmpdir"; } &" )
    }
    tmpf=tmpdir"/"regfile".gz"
    cmd="gzip -cdfq "regfile" | sort -k1,1 -k2,3n | bgzip > "tmpf" && tabix -p bed "tmpf
    print "Executing command: "cmd > "/dev/stderr"
    system ( cmd ) 
    return tmpf
  }
}

function convert_coords(chainfile, chr, pos1, pos2) {
  gotpos=0
  while ( getline < chainfile > 0 ) {
    if ( chain::parse_next_chain(chainfile) != 0 ) { 
      print "Error: invalid chain "$0 > "/dev/stderr"
      exit 1
    }
    if ( chain::tchr!=chr || ( chain::t0 > pos1 && chain::t0 > pos2 ) || ( chain::t1 < pos1 && chain::t1 < pos2 ) ) { continue }
    if ( ( gotpos=slide_chain(pos1, pos2) ) >= npos ) {
      break
    }
  }
  close(chainfile)
  if ( gotpos < npos ) {
    print "Warning: position "chr":"(1 in pos ? pos[1]+1 : pos[2])" not covered by the chain file" > "/dev/stderr"
    return ""
  }
  return_bed()
}

