@namespace "faidx"

function run() {
  help="\
      Reimplementation of `samtools faidx` in AWK. \n\
      Currently only supports FASTA files and a subset of arguments. \n\
      Arguments are supplied positionally (first reference, then region/s) and/or using named options. \n\
      EXAMPLE: \n\tulitka faidx file.fa \n\tulitka faidx file.fa chr1:1-100"
  arg::add_argument("f", "fasta", 0, "reference FASTA file")
  arg::add_argument("r", "reg", 0, "chr:start-end or BED file")
  arg::add_argument("i", "reverse-complement", 1, "reverse complement sequences")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  if ( "fasta" in arg::args ) { 
    fasta=arg::args["fasta"] 
  } else {
    fasta=arg::nonargs[1]
    delete arg::nonargs[1]
  }
  if ( "reg" in arg::args ) { regs[arg::args["reg"]]=nregs=1 }
  for(i in arg::nonargs) {
    regs[arg::nonargs[i]]=1
    nregs++
  }

  ulitka::assert( fasta != "", "FASTA file required" )
  if ( nregs > 0 ) {
    for (reg in regs) {
      while ( getline < reg > 0 ) {
        isfile=1
        strand="" # can add BED column for strand here if needed
        if (arg::args["reverse-complement"] || strand=="-") {
          seq = revcomp(query(fasta, $1, $2+1, $3))
          print ">" regg "/rc"
          for (i=1;i<=length(seq);i+=60) {
            print substr(seq,i,60) 
          }
        } else {
          query(fasta, $1, $2+1, $3, "print")
        }
      }
      if ( ! isfile ) {
        parse_region(reg)
        if (arg::args["reverse-complement"] || strand=="-") {
          seq = revcomp(query(fasta, chr, start, end))
          print ">" reg "/rc"
          for (i=1;i<=length(seq);i+=60) {
            print substr(seq,i,60) 
          }
        } else {
          query(fasta, chr, start, end, "print") # any nonnull value to "doprint" triggers printing
        }
      }
    }
    exit 0
  } else {
    exit faindex(fasta)
  }
}

function faindex(fa, fai) {
  if ( fa=="" ) { 
    print "Could not detect FASTA input" > "/dev/stderr"; return 1 
  }
  if (fai=="") { fai = fa".fai" }
  while ( (r=getline < fa > 0) || !lastline ) {
    nr++
    bytes+=length($0)+1
    if ($0~/^>/ || (r<=0 && lastline=1) ) {
      if ( nr>1 ) { 
        if ( nr-lastnr > 1 ) {
          lenline=len/(nr-lastnr-1)
        } else {
          lenline=0
        }
        if ( int(lenline) != lenline ) {
          print "Error: variable line length "lenline" between lines "lastnr+1"-"nr-1 > "/dev/stderr"
          return 1
        }
        len+=lendelta
        if (lenline == 0) { lenline = len }
        print seq"\t"len"\t"lastbytes"\t"lenline"\t"lenline+1 > fai
        len=lendelta=0
      }
      seq=substr($0,2)
      lastnr=nr
      lastbytes=bytes
    } else {
      len+=lendelta
      lendelta=length($0)
    }
  }
  close(fa)
}

function query(fa, chr, start, end, doprint, fai,
                     found, skip, brkpt, tail, last, out, nr) {
  if (fai=="") { fai = fa".fai" }
  if ( system("test -f " fai) != 0 ) {
    print "FAI index not found, making one now" > "/dev/stderr"
    faindex(fa)
  }
  while ( getline < fai > 0 ) {
    if ( $1 != chr ) {
      skip += 1 + int( $2/ $4 ) + (($2 % $4) > 0)
    } else {
      if (int(start) > int($2)) { 
        print "Warning: start coordinate " start " exceeds sequence length "$2 > "/dev/stderr"
        break
      }
      if (int(end) > int($2)) { 
        print "Warning: end coordinate " end " truncated to sequence length "$2 > "/dev/stderr"
        end=$2 
      }
      found=1
      skip += 1 + int( start / $4 )
      brkpt = start % $4
      tail = end % $4
      last = skip + 1 + int((end / $4))
      if ( doprint!="" ) {
        printf("%s", ">"chr":"start+1"-"end)
        while ( getline < fa > 0 && ++nr <= last ) {
          if (nr > skip) {
            if (nr != skip+1) { printf("%.*s", brkpt, $0) }
            if (nr != last) { printf("\n%s", substr($0, brkpt+1)) 
            } else if (nr == skip+1) { tail-=brkpt; printf("\n%s", substr($0,brkpt+1,tail))
            } else { printf("\n%s", substr($0,1,tail)) }
          }
        }
        close(fa)
        if (tail) { print "" }
        out = 0
        break
      } else {
        # out=">"chr":"start+1"-"end"\n"
        while ( getline < fa > 0 && ++nr <= last ) {
          if (nr > skip) {
            if (nr != skip+1) { out = out substr($0, 1, brkpt) }
            if (nr != last) { out = out substr($0, brkpt+1) 
            } else if (nr == skip+1) { tail-=brkpt; out = out substr($0,brkpt+1,tail)
            } else { out = out substr($0,1,tail) }
          }
        }
        close(fa)
        break
      }
    }
  }
  close(fai)
  if ( !found ) {
    print "Warning: region "chr":"start+1"-"end" not found in "fa > "/dev/stderr"
  }
  return out
}

function revcomp(string,    i, n, out, s) {
  if ( _revcompinit == 0 ) {
    split("*.ACGT.*", n, "")
    for (i=1;i<7;i++) { _revcomp[n[i]]=n[9-i] }
    split("NnacgtnN", n, "")
    for (i=1;i<7;i++) { _revcomp[n[i]]=n[9-i] }
    _revcompinit=1
  }
  split(string, s, "")
  for (i in s) {
    ulitka::assert( s[i] in _revcomp, "Cannot compute reverse complement of character " s[i] )
    out = _revcomp[s[i]] out
  }
  return out
}

function parse_region(reg,    sep1, sep2) {
  # Read chr:start or chr:start-end(/rc) notation and populate chr, start, end, strand variables
  if (reg ~ /\/rc$/ ) { 
    strand="-"
    reg=substr(reg, 1, length(reg)-3)
  } else {
    strand="+"
  }
  sep1=match(reg, /:[0-9]+-?[0-9]+$/)
  sep2=index(substr(reg, sep1+1), "-")
  chr=substr(reg,1,sep1-1)
  if (sep2 > 0) {
    start=substr(reg,sep1+1,sep2-1)
    end=substr(reg,sep1+sep2+1)
  } else {
   start = end = substr(reg,sep1+1)
  }
  if (start < 1) { start=1 }
  if (start - end > 0) {
    print "Error: start coordinate exceeds end coordinate in "reg > "/dev/stderr"
    exit 1
  }
  if (chr=="" || start=="" || end=="") {
    print "Error: could not parse region notation "reg > "/dev/stderr" 
    exit 1
  }
  start--
}
