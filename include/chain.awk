@namespace "chain"

# Some helper functions to work with chain files
# Todo: identify and use lchains

function parse_next_chain(chainfile,
                          direction, q, t, x, y, z) {

  # Memorize deltas and breakpoints: CPU-heavy but RAM-light
  if ( $1!="chain" ) { return 1 }

  # Variables to be returned
  tchr=$3
  t0=$6
  t1=$7
  qchr=$8
  q0=$11
  q1=$12
  samestrand=( $5 == $10 ) # same as checking if $10=="+"
  idx=0 # number of chain segments

  # Local vars
  direction=( samestrand ? 1 : -1 )
  q=( samestrand ? q0 : q1 )
  t=t0

  delete delta
  delete start
  delete end
  while ( getline < chainfile > 0 ) {
    if ( $0=="" ) { break }
    delta[idx]=q-t
    start[idx]=t
    end[idx]=t+$1
    t+=$1+$2
    q+=direction*($1+$3)
    idx++
  }
  return 0

  # 1-to-1 correspondence array: simple but memory-heavy
  # delete corr
  # while ( getline < chainfile > 0 ) {
  #   if ( $0=="" ) { break }
  #   while ($1--) {
  #     corr[t]=q++
  #     t+=direction
  #   }
  #   t+=direction*$2
  #   q+=$3
  # }
}

