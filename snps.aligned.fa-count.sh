#!/bin/bash 


#printf "name\tN\tA\tC\t-\tn\ta\tT\tG\tc\tX\tt\tg\n"

>&2 printf "name\tN\tA\tC\t-\tn\ta\tT\tG\tc\tX\tt\tg\n"
printf $(basename $1)"\t"

awk '
function count() {
    n=length(string);
    gsub(/[A]+/,"",string); m=length(string); c["A"]+=n-m; n=m
    gsub(/[C]+/,"",string); m=length(string); c["C"]+=n-m; n=m
    gsub(/[G]+/,"",string); m=length(string); c["G"]+=n-m; n=m
    gsub(/[T]+/,"",string); m=length(string); c["T"]+=n-m; n=m
    gsub(/[N]+/,"",string); m=length(string); c["N"]+=n-m; n=m

    gsub(/[a]+/,"",string); m=length(string); c["a"]+=n-m; n=m
    gsub(/[c]+/,"",string); m=length(string); c["c"]+=n-m; n=m
    gsub(/[g]+/,"",string); m=length(string); c["g"]+=n-m; n=m
    gsub(/[t]+/,"",string); m=length(string); c["t"]+=n-m; n=m
    gsub(/[n]+/,"",string); m=length(string); c["n"]+=n-m; n=m

    gsub(/[X]+/,"",string); m=length(string); c["X"]+=n-m; n=m
    gsub(/[-]+/,"",string); m=length(string); c["-"]+=n-m; n=m
}
BEGIN{RS="\n>"; FS="\n"}
{
  #printf $1"\t"
  string=substr($0,length($1)); count()
  for(i in c) printf c[i]"\t"; print ""
  delete c; string=""
}
' ${1}/snps.aligned.fa
