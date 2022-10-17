cd /home/scilab/Desktop/shankar_fifth/4residuals/
for file in *.nc


do 
ferret << eof
use $file
let ures1 = ures
let vres1 = vres
save/file=residuals_$file/clobber ures1,vres1
sp mv residuals_$file /home/scilab/Desktop/shankar_fifth/4residuals/
eof
done
