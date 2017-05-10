awk '{if($1>0) print NR" "$0}' ~/Dropbox/PsTonB/Data/coeffs331.dat | awk '{sum=sum+$2; print sum}'
awk '{if($1>0) print NR" "$0}' ~/Dropbox/PsTonB/Data/coeffs331.dat
awk '{if($1>0) print NR" "$0}' ~/Dropbox/TonB/Data/coeffs331.dat | awk '{sum=sum+$2; print sum}'
awk '{if($1>0) print NR" "$0}' ~/Dropbox/TonB/Data/coeffs331.dat
