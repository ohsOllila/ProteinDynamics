awk '{if($1>0) printf "%12.3f %12.3f %12.3f\n", NR, $1, $2}' ~/Dropbox/PsTonB/Data/coeffs331.dat | awk '{sum=sum+$2; print sum}'
awk '{if($1>0) printf "%12.3f %12.3f %12.3f\n", NR, $1, $2}' ~/Dropbox/PsTonB/Data/coeffs331.dat
awk '{if($1>0) printf "%12.3f %12.3f %12.3f\n", NR, $1, $2}' ~/Dropbox/TonB/Data/coeffs331.dat | awk '{sum=sum+$2; print sum}'
awk '{if($1>0) printf "%12.3f %12.3f %12.3f\n", NR, $1, $2}' ~/Dropbox/TonB/Data/coeffs331.dat
