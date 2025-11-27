@namespace "fa2chain"

function run(){
help="\
    Convert aligned FASTA into a chain file. \n\
    First sequence should come from new (query) genome, other from old (ref) genome.  \n\
    Sequences names should follow faidx output: 'chr:start-end' with '/rc' on end if reverse-complemented. \n\
    Arguments are supplied positionally (old ref, new ref, FASTA/s to convert) and/or using named options. \n\
    EXAMPLE: \n\tulitka fa2chain -1 old.fa -2 new.fa file.fa > file.chain \n\
    GNU PARALLEL EXAMPLE: \n\t parallel ulitka fa2chain -1 old.fa -2 new.fa {} ::: *.fa > file.chain"

  arg::add_argument("1", "old", 0, "old reference FASTA")
  arg::add_argument("2", "new", 0, "new reference FASTA")
  arg::add_argument("f", "fasta", 0, "FASTA file to be converted")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  if (ARGC < 4) { print help; exit 1 }
  if ( "old" in arg::args ) { 
    old=arg::args["old"] 
  } else {
    old=arg::nonargs[1]
    delete arg::nonargs[1]
  }
  if ( "new" in arg::args ) { 
    new=arg::args["new"] 
  } else {
    new=arg::nonargs[2]
    delete arg::nonargs[2]
  }
  if ( "fasta" in arg::args ) { fastas[arg::args["fasta"]]=nfastas=1 }
  for(i in arg::nonargs) {
    fastas[arg::nonargs[i]]=1
    nfastas++
  }
  ulitka::assert( old != "", "old ref, new ref, FASTA required" )
  ulitka::assert( new != "", "old ref, new ref, FASTA required" )
  ulitka::assert( nfastas > 0, "old ref, new ref, FASTA required" )

  exit main()
}

function main() {
  fai1 = old".fai"
  fai2 = new".fai"
  if ( system("test -f " fai1) != 0 ) {
    print "FAI index for "old" not found, making one now" > "/dev/stderr"
    faidx::faindex(old)
  }
  while ( getline < fai1 > 0) {
    oldlen[$1]=$2
  }
  close(fai1)
  if ( system("test -f " fai2) != 0 ) {
    print "FAI index for "new" not found, making one now" > "/dev/stderr"
    faidx::faindex(new)
  }
  while ( getline < fai2 > 0) {
    newlen[$1]=$2
  }
  close(fai2)

  for (fasta in fastas) {
    delete dash
    fnr=lastline=l=0
    header=qry=""
    while ( (r=getline < fasta > 0) || !lastline) {
      fnr++
      if ($0 ~ /^>/ || (r<=0 && lastline=1) ) {
        if (fnr==1) { 
          qry=substr($0,2) 
        } else {
          if (header != "") {
            printout()
            # Cleanup for next seq
            l/=2
          }
          header=substr($0,2) SUBSEP qry
        }
      } else {
        while ( i=index($0,"-") ) {
          l+=i
          dash[l]=1
          $0=substr($0,i+1)
        }
        l+=length($0)
      }
    }
    close(fasta)
  }
}

function printout(    start, i, dd, size, out, dtdq, hang_start, h, hh) {
  start=1
  for(i=1;i<=l/2;i++) {
    dd=2*dash[i]+1*dash[l/2+i] # 0 if both non-gap, 1 if gap in 2nd, 2 if gap in 1st, 3 if both gap
    delete dash[l/2+i]
    if ( dd==0 ) {
      start=0
      if (size > 0 && (dtdq[1] + dtdq[2] > 0) ) {
        out=out "\n" size " " int(dtdq[2]) " " int(dtdq[1])
        size=0
      }
      dtdq[1]=dtdq[2]=0
      size++
    } else if ( dd != 3 ) {
      if ( start ) {
        hang_start[dd]++
      } else {
        dtdq[dd]++
      }
    }
  }
  out=out "\n" size
  split(header, h, SUBSEP)

  faidx::parse_region(h[1])

  if (!(faidx::chr in oldlen)) {
    print "Error: chromosome "faidx::chr" not found in genome "old
    exit 1
  }

  oldstrand=faidx::strand
  faidx::start+=hang_start[2]
  faidx::end-=dtdq[2]
  hh="chain 9999 "faidx::chr" "oldlen[faidx::chr]" + "faidx::start" "faidx::end

  faidx::parse_region(h[2])

  if (!(faidx::chr in newlen)) {
    print "Error: chromosome "faidx::chr" not found in genome "new
    exit 1
  }

  if (faidx::strand != oldstrand) { faidx::strand="-" } # bcftools liftoff only takes plus-strand target features
  faidx::start+=hang_start[1]
  faidx::end-=dtdq[1]
  hh=hh" "faidx::chr" "newlen[faidx::chr]" "faidx::strand" "faidx::start" "faidx::end" "h[2]

  if (out != "\n") {
    print hh
    print substr(out,2)
    print ""
  }
}
