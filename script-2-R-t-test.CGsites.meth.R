#
#TN.w	RA.w	RARO.w	TEM.w	TN.c	RA.c	RARO.c	TEM.c
#1	0.8	0.714286	0.64	0.9	0.857143	0.714286	0.6



ls | grep "^ENSG" | while read f 
do 
 wc -l $f; 
done | awk '
 {
  if($1 > 0)
    print $2;
 }
 ' > temp-ENSG-files

cat temp-ENSG-files | head -3 | while read file
do
 echo "name <- c(\"$file\")";
 cut -f 1 $file > temp-s1
 cut -f 5 $file >> temp-s1
 cut -f 4 $file > temp-s2
 nl=`wc -l temp-s2 | awk '{print $1}'`
 cut -f 8 $file >> temp-s2
 echo "nl = $nl"
 
 cat temp-s1 | awk '
 BEGIN{
  FS = "\t";
  OFS = "";
  ORS = "";
      }
 {
  if(NR == 1)
  {
   print "a <- c(";
  }
  print $1 ", "
 }
 END{print "\n"}' | sed 's/, $/)/'

 cat temp-s2 | awk '
 BEGIN{
  FS = "\t";
  OFS = "";
  ORS = "";
      }
 {
  if(NR == 1)
  {
   print "b <- c(";
  }
  print $1 ", "
 }
 END{print "\n"}' | sed 's/, $/)/'

 echo "cor <- t.test(a, b, paired = TRUE)";
 echo "name";
 echo "nl";
 echo "mean1 <- mean(a)"
 echo "mean2 <- mean(b)"
 echo "cor\$p.value";
 echo "mean1";
 echo "mean2";
 echo "cor\$estimate";
done
