@namespace "lift"

function run(){
help="\
    Updates coordinates in a VCF file using a chain file. \n\
    A primitive liftover tool that only works well for SNPs. \n\
    Requires tabix in PATH. The VCF must be bgzipped and tabixed. \n\
    Arguments are supplied positionally (first chain, then VCF) and/or using named options. \n\
    Make sure to polish (reheader, sort) the resulting VCF before indexing! \n\
    EXAMPLE: \n\tulitka lift file.chain file.vcf.gz | bgzip > file_liftover.vcf.gz \n\
    GNU PARALLEL EXAMPLE: \n\t parallel --pipepart -a file.chain --recend '\\n\\n' ulitka lift - file.vcf.gz | bgzip > file_liftover.vcf.gz"

  arg::add_argument("v", "vcf", 0, "VCF file to be lifted over")
  arg::add_argument("c", "chain", 0, "chain file")
  arg::add_argument("f", "fasta", 0, "reference FASTA file")
  arg::parse_args(2, help)
  arg::parse_nonargs()

  if ( system("tabix --version > /dev/null") != 0 ) {
       print "Can't run tabix"; exit 1
  }
  if ( "vcf" in arg::args ) { 
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
  exit main()
}

function main() {
  for (c in chains) {
    while (getline < c > 0) {

      # Parse next chain into several variables (tchr, t0, t1, qchr, q0, q1, samestrand, idx) and arrays (start, end, delta)
      if ( chain::parse_next_chain(c) != 0 ) { continue }

      # get the relevant reference chunk
      if ( "fasta" in arg::args ) { 
        seq = faidx::query(arg::args["fasta"], chain::qchr, chain::q0+1, chain::q1)
      }

      cmd="tabix "vcf" "chain::tchr":"chain::t0+1"-"chain::t1

      # Get relevant VCF lines, replace CHR, POS, and revcomp REF and ALT if needed
      while ( cmd | getline vcfline > 0 ) {
        if (vcfline ~ /^#/) { 
          print vcfline
          continue
        }
        nf=split(vcfline, v, "\t")
        # Currently skipping indels because they get misaligned if the chain is inverted
        if ( v[5] ~ /[^,][^,]/ ) { continue }
        for (id=0;id<=idx;id++) {
          if ( v[2] > chain::end[id] ) { delete chain::delta[id] }
          else if ( v[2] >= chain::start[id] && v[2] < chain::end[id] ) {
            if ( !chain::samestrand ) { 
              newalt=""
              v[4]=revcomp(v[4])
              split(v[5],alt,",")
              for(i in alt) {
                newalt= newalt "," revcomp(alt[i]) 
              }
              v[5]=substr(newalt,2)
              v[2]-=2*(v[2]-chain::end[id]) + length(v[4]) + 1
            }
            # Compare REF against new reference and change REF/ALT/GT if needed
            if ( (ref=substr( seq, v[2]+chain::delta[id]-chain::q0, length(v[4]) )) != "" && toupper(ref) != toupper(v[4]) ) {
              whereal=0
              # print "Warning: site "v[2]+delta[id]" does not match new ref: "v[4]">"ref > "/dev/stderr"
              if (v[5]==".") {
                v[5]=v[4]
                newal=1
              } else if ( (whereal = index(v[5], ref)) > 0 ) {
                # Need more strict regex search for single-base ref match here if using indels
                v[5]=substr(v[5],1,whereal-1) v[4] substr(v[5],whereal+1)
                whereal = (whereal+1)/2
                newal=0
              } else {
                v[5]=v[5]","v[4]
                newal=gsub(/,/,",",v[5])+1
              }
              v[4]=ref
              for(i=10;i<=nf;i++) {
                gtend=index( v[i], ":" )
                if ( gtend==0 ) { gtend=length(v[i])+1 }
                for (c=1; c<gtend; c+=2) {
                  al=substr(v[i],c,1)
                  if ( al == whereal ) {
                    v[i]=substr(v[i],1,c-1) newal substr(v[i],c+1)
                  } else if ( al == 0 ) {
                    v[i]=substr(v[i],1,c-1) whereal substr(v[i],c+1)
                  }
                }
              }
            }
            printf("%s\t%s", chain::qchr, v[2]+chain::delta[id])
            for(i=3;i<=nf;i++) { printf("\t%s", v[i]) }
            print ""
            break
          }
        }
      }
      close(cmd)
    }
    close(c)
  }
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

