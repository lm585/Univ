#

echo "dna <- c(2.7513561, 2.787743772, 2.724357804, 2.830075166, 2.688241796, NA, 2.517591731, 2.624694531, 2.508664363, 2.700097705, 2.351216345, 2.880642126, 2.321805484, 2.577836341, 2.452246575, 2.696793085, 2.22297645, 2.84460142, 2.392696953, 2.463445032, 2.507180977, 2.523746467, 2.503654519, 2.621902961)"

cat  cd4-3dates-d26.txt | awk '
 BEGIN{
  FS = "\t";
  OFS = "";
  ORS = "";
      }
  {
   print "name <- c(\"" $1 "\")" "\n";
   print "path <- c(";
   for(i = 2; i < NF; i++) ###pathway day209-6-elements 8th
   {
    print $i ", ";
   }
   print $NF ")" "\n";
   print "cor <- cor.test(dna, log10(path+1)) \n";
   print "name\n";
   print "cor$p.value\n";
   print "cor$estimate\n";
  }
 '

