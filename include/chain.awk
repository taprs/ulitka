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

  # Local vars
  idx=0
  direction=( samestrand ? 1 : -1 )
  q=q0+1
  t=( samestrand ? t0+1 : t1 )

  delete delta
  delete start
  delete end
  while ( getline < chainfile > 0 ) {
    if ( $0=="" ) { break }
    idx++
    delta[idx]=q-t
    start[idx]=(samestrand ? t : t-$1+1)
    end[idx]=(samestrand ? t+$1 : t+1)
    t+=direction*($1+$2)
    q+=$1+$3
  }
  # order from smallest to largest start
  if ( !samestrand ) {
    for (i=0;i<idx/2;i++) {
      x=delta[i]
      y=start[i]
      z=end[i]
      delta[i]=delta[idx-i]
      start[i]=start[idx-i]
      end[i]=end[idx-i]
      delta[idx-i]=x
      start[idx-i]=y
      end[idx-i]=z
    }
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

