cat /home/samuli/Dropbox/PsTonB/Data/OPfromPLATEAU.dat | awk '{if($1==331){print "0 "$2; print "100  "$2}}' > /home/samuli/Dropbox/PsTonB/Data/OPres331.dat
cat /home/samuli/Dropbox/PsTonB/Data/OPfromPLATEAUwithOPC.dat | awk '{if($1==331){print "0 "$2; print "100  "$2}}' > /home/samuli/Dropbox/PsTonB/Data/OPres331withOPC.dat
cat /home/samuli/Dropbox/PsTonB/Data/OPfromPLATEAU_T298K.dat | awk '{if($1==331){print "0 "$2; print "100  "$2}}' > /home/samuli/Dropbox/PsTonB/Data/OPres331_298K.dat
