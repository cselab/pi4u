#!/bin/bash
# transpose_prop.sh: this script is used to transpose the propapagated samples.
# http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash
# The output file can be displayed with the display_prop tool

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $1
