#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

struct chrSeq
{
 string name;
 vector<char> seq;
 vector<string> strSeq; //A C G N T AC  ACG
};

struct pileup
{
 string chr;
 int pos;
 char base;
 int cover;
 string readStr;
 string quality;
 void clear();
};

struct param
{
 int coverage;
 int indBaseCover;
 double indBasePerc;
 double delPerc;
 double insMajor;
};

void setChrSeq(ifstream & fasta);
void setN_ToSampleSeq();
bool line2pile(string line, pileup & pile);
string pileup2baseProfile( const pileup & pile, const param & par);
int parseInDelStr(string mapStr, int i, string & ins);
bool updateSampleSeq(const pileup& pile, string baseRep);

vector<chrSeq> glob_chr;
int glob_chr_index = 0;

int main(int argc, char * argv[])
{
 if(argc != 10)
 {
  cerr << argv[0] << "    ref-genome-fasta   output-file-name pileup-file sample-name   coverage_cutoff(2)  individual-base-cover(2)   individual-base-perc(0.3) indel-perc(0.5, 0.8, etc, [8])  insertion-majority-perc(0.5, etc. [9])\n";
  cerr << "deletion: 80% of all reads are deletions.\n";
  cerr << "insertion:  80% of all reads are deletions. Dorminant insertion >= 80% of all insertions.\n";
  return 1;
 }

 int lineNum = 0;
 ifstream input;
 ofstream output;
 param par;
 string line, baseRep;
 pileup pile;
  
 par.coverage = atoi(argv[5]);
 par.indBaseCover = atoi(argv[6]);
 par.indBasePerc = atof(argv[7]);
 par.delPerc = atof(argv[8]);
 par.insMajor = atof(argv[9]);

 input.open(argv[1], ios::in);
 if(!input)
   {
    cerr << argv[1] << " cannot be opened for reading!" << endl;
    return 1;
   }
 setChrSeq(input);
 input.close();

 for(int i = 0; i < glob_chr.size(); i++)
 {
  cout << glob_chr[i].name << '\t' << glob_chr[i].seq.size() << endl;
 }

 setN_ToSampleSeq();
 
 input.open(argv[3], ios::in);
 if(!input)
   {
    cerr << argv[3] << " pile-up-file cannot be opened for reading!" << endl;
    return 1;
   }
 getline(input, line);
 while(!input.eof())
 {
  if(lineNum % 10000000 == 0)
    cout <<  argv[3] << " pile-up-file line: " << lineNum << endl;
  lineNum++;
  bool b = line2pile(line, pile);
  if(b)
  {
   baseRep = pileup2baseProfile(pile, par);
   if(baseRep != "")
   {
     if(updateSampleSeq(pile, baseRep) == false)
     {
      cerr <<  argv[3] << " pile-up-file" << pile.chr << '\t' << pile.pos << '\t' << pile.base;
      cerr << '\t' << "not agree with reference chromosome." << endl;
      return 1;
     }
   }
  }
  getline(input, line);
 } 
 input.close();

 output.open(argv[2], ios::out);
 if(!output)
 {
  cerr << argv[2] << " cannot be written to.\n";
  return 1;
 }
 output <<   argv[4] << endl;
 for(int i = 0; i < glob_chr.size(); i++)
 {
  for(int j = 0; j < glob_chr[i].strSeq.size(); j++)
    output <<   glob_chr[i].strSeq[j] << endl;
 }

 output.close();
 
 return 0;
}

void setChrSeq(ifstream & fasta)
{
 string line, str;
 char * unitStr, * lineCstr;
 vector<string> lineStrList;
 chrSeq myChrSeq;

 getline(fasta, line);
 while(!fasta.eof())
 {
  lineStrList.clear();
  lineCstr = new char [line.length() + 10];
  strcpy(lineCstr, line.c_str());
  unitStr = strtok(lineCstr, " \t\n");
  while(unitStr != NULL)
  {
   str = unitStr;
   lineStrList.push_back(str);
   unitStr = strtok(NULL, " \t\n");
  }

  if(lineStrList.size() > 0)
  {
   if(lineStrList[0][0] == '>')
   {
    if(myChrSeq.seq.size() > 0)
    {
      glob_chr.push_back(myChrSeq);
    }
    myChrSeq.name = lineStrList[0].substr(1);
    myChrSeq.seq.clear();
   }
   else
   {
    for(int i = 0; i < lineStrList.size(); i++)
      for(int j = 0; j < lineStrList[i].length(); j++)
        myChrSeq.seq.push_back(lineStrList[i][j]); 
   }
  }
  getline(fasta, line);
 }
 if(myChrSeq.seq.size() > 0)
 {
  glob_chr.push_back(myChrSeq);
 }
 return;
}

void setN_ToSampleSeq()
{
 string str = "N";
 for(int i = 0; i < glob_chr.size(); i++)
 {
  for(int j = 0; j < glob_chr[i].seq.size(); j++)
    glob_chr[i].strSeq.push_back(str);
 }
}

bool line2pile(string line, pileup & pile)
{
 int numTab = 0;
 vector<int> tabPos;

 for(int i = 0; i < line.length(); i++)
 {
  if(line[i] == '\t')
  {
   tabPos.push_back(i);
  }
 }
 if(tabPos.size() < 5)
   return false;
 int len, coverage;
 char refBase;
 string mapStr, str;
 
 //chrosome 1st field
 len = tabPos[0] ;
 pile.chr = line.substr(0, len);
 //chr position 2nd field, between 1st and 2nd tab
 len = tabPos[1] - tabPos[0] - 1;
 str = line.substr(tabPos[0] + 1, len);
 pile.pos = atoi(str.c_str());
 //ref base is following the 2nd tab
 refBase = line[tabPos[1] + 1 ];
 pile.base = refBase;
 //coverage between 3rd tab and 4th tab
 len = tabPos[3] - tabPos[2] - 1;
 str = line.substr(tabPos[2] + 1, len);
 coverage = atoi(str.c_str());
 pile.cover = coverage;
 //length should be 4th tab and 5th tab
 len = tabPos[4] - tabPos[3] - 1;
 mapStr = line.substr(tabPos[3] + 1, len);
 pile.readStr = mapStr;
 return true;
}

string pileup2baseProfile( const pileup & pile, const param & par)
{
 char letter, refBase;
 string mapStr = "", insStr = "", acg = "ACGT"; 
 vector<string> insList;
 int numIn = 0, numDel = 0;
 int numMatches = 0, mutCover = pile.cover, insCover = pile.cover, delCover = pile.cover;
 map<char, int> baseCount;
 map<string, int> insDormMap;
 bool res = false;

 for(int i = 0; i < acg.length(); i++)
   baseCount[acg[i]] = 0;
 refBase = toupper(pile.base);
 for(int i = 0; i < pile.readStr.length(); i++)
   mapStr += toupper(pile.readStr[i]);

 for(int i = 0; i < mapStr.length(); i++)
 {
  letter = mapStr[i];
  if(letter == '^')
  {
   i++; //skip the char following ^
   delCover--;
  }
  else if (letter == '$')
  {
   insCover--;
   delCover--;
  }
  else if (letter == '+' ) //+20TTTTTTTTTTTT....TG
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
   insList.push_back(insStr);
   numIn++;
  }
  else if (letter == '-' )
  {
   i = parseInDelStr(mapStr, i, insStr); //i points to the last 'agctn'
  }
  else if (letter == '*' )
  {
   numDel++;
  }
  else if (letter == ',' || letter == '.')
  {
   numMatches++;
  }
  else if (letter == 'N' || letter == 'n')
  {
   mutCover--;
  }
  else if(acg.find(letter) != string::npos) //ACGT acgt
  {
   baseCount[letter]++;
  }
  else; //other symbol
 }
 string resStr = "";
 if(mutCover >= par.coverage)
 {
  //populate bases
  baseCount[refBase] = numMatches;
  for(int i = 0; i < acg.length(); i++)
    if(baseCount[acg[i]] >= par.indBaseCover &&
        (baseCount[acg[i]] / (double) mutCover) >= par.indBasePerc) 
    {
     resStr += acg[i];
    }
  if(resStr.length() >= 3)
  {
    resStr = "";
  }
  else if(resStr.length() == 2)
  {
   char tempCh = resStr[0];
   if(resStr[0] == refBase)
   {
    tempCh = resStr[1];
   }  
   //resStr[1] == refBase
   //resStr has no refBase
   resStr = tempCh; //garantee one char
  }
  else; //resStr.length() == 1 0 
 }

 double indelThresh = par.delPerc; //0.8;
 if(delCover  >=  par.coverage && numDel / (double) delCover >= indelThresh ) 
 {
  resStr = '-';
 }

 if(false) //(insCover >= par.coverage)
 {
  if(numIn > 0)
  {
   for(int i = 0; i < insList.size(); i++)
   {
    if(insDormMap.count(insList[i]) == 0)
    {
     insDormMap[insList[i] ] = 1;
    }
    else
      insDormMap[insList[i] ]++;
   }
   int tempInt = 0;
   string tempStr;
   map<string, int>::const_iterator iter;
   for( iter = insDormMap.begin(); iter != insDormMap.end(); iter++)
   {
     if(iter->second > tempInt)
     {
       tempInt = iter->second;
       tempStr = iter->first;
     }
   }
   if(numIn  / (double) insCover >= indelThresh 
      && tempInt / (double) numIn >= par.insMajor)
      resStr = resStr + "+" + tempStr;
  }
 }

 return resStr;
}

int parseInDelStr(string mapStr, int i, string & ins)  //+20TTTTTTTTTTTT....TG; -3AAC
{
   int begin ;
   i++; //the first digit
   string strIn(1, mapStr[i]);
   i++;
   while(isdigit(mapStr[i]))
   {
    strIn += mapStr[i];
    i++;
   }
   //exit while loop, i points to the first 'acgtn'
   begin = i;
   int numIn = atoi(strIn.c_str());
   i = i - 1 + numIn; //i points to the last 'agctn'
   int length = i + 1 - begin;
   ins = mapStr.substr(begin, length);
   return i;
}

/*
 * glob_chr_index 
 * first, pile.chr is compared with glob_chr_index, if matching, updateSampleSeq
 * If not matching, compare with ++glob_chr_index, until found or the end of all chromosomes, 
 * if still not found, give an error message, exit the program
 */
bool updateSampleSeq(const pileup& pile, string baseRep)
{
 bool res = false;
 if(pile.chr == glob_chr[glob_chr_index].name)
 {
  res = true;
 }
 else
 {
  for(++glob_chr_index; glob_chr_index < glob_chr.size(); ++glob_chr_index)
  {
   if(pile.chr == glob_chr[glob_chr_index].name)
   {
    res = true;
    break;
   }
  }
 }
 if(res == true)
 {
  if(glob_chr[glob_chr_index].strSeq.size() >= pile.pos
       && toupper(glob_chr[glob_chr_index].seq[pile.pos-1]) == toupper(pile.base)) //pos 1 based
   {
    glob_chr[glob_chr_index].strSeq[pile.pos-1] = baseRep;
    return true;
   }
 }
 return false;
}

