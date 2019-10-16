#

 cat     MSD-chem-day209-230-252.readout.txt | cut -f 1-6 | awk '
 BEGIN{
  FS = "\t";
  OFS = "";
  ORS = "";
      }
  {
   if(NR % 4 == 1)
    print "name <- c(\"" $1 "\")" "\n";
   else
   {
    print "a" (NR % 4)  " <- c(";
    for(i = 1; i < NF; i++) ###pathway day209-6-elements 8th
    {
     print $i ", ";
    }
    print $NF ")" "\n";
   }
   if(NR % 4 == 0)
   {
    print "cor <- t.test(log10(a2 + 1), c(log10(a3+1),log10(a0+1))) \n";
    print "name\n";
    print "cor$p.value\n";
    print "cor$estimate\n";
   } 
  }
 '

