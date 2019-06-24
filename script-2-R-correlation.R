#

# date209 / 230 / 252
msd="msd <- c(-0.253969298, -1.277844198, -0.01461541, -0.607024862, -0.993804929, -0.447623826, 1.366856173, -0.208389096, 0.44166147, -0.310936981, 0.521554855, 2.018746443, 2.623535655, -0.459012058, -0.79218561, -0.105814457, 0.955473083, 1.138275568)"

echo "$msd"
cat pbmc-616-pathways-Malika-upAndDown-preRank-hyp-gene.dChip.txt | awk '
 BEGIN{
  FS = "\t";
  OFS = "";
  ORS = "";
      }
  NR > 1 {
   print "name <- c(\"" $1 "\")" "\n";
   print "path <- c(";
   for(i = 8; i < NF; i++) ###pathway day209-6-elements 8th
   {
    print $i ", ";
   }
   print $NF ")" "\n";
   print "cor <- cor.test(msd, path) \n";
   print "name\n";
   print "cor$p.value\n";
   print "cor$estimate\n";
  }
 '

