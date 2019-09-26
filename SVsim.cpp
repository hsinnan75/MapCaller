#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <ctype.h>
#include <cstring>
#include <cmath>
#include <ctime>
#include <map>
#include <vector>
#include <algorithm>

#define DOM			1000000
#define MutationBlock	3000
#define SNP_rate	3000 // SNPs # per 1M base
#define sIND_rate	200 // small indels # per 1M base (1~10bp)
#define lIND_rate	50 // large indels # per 1M base (11~50bp)
#define TraLoc_rate 1	// 2 translocation # per 1M base
#define Inv_rate	1 //1	 // inversion # per 1M base
#define CNV_rate	1 //1	 // CNV # per 1M base (1000bp)

using namespace std;

typedef struct
{
	int gPos;
	int mtype;
	string ori_seq;
	string mut_seq;
} SV_t;

FILE *vcf_fd, *mut_fd, *info_fd;
int iSNP, isIND, ilIND, iTraLoc, iInv, iCNV;

static const char ReverseMap[255] =
{
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
	'\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
	'\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
	'\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', /*  90 -  99 */
	'\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
	'N',  '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', /* 110 - 119 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
	'\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

int myrandom(int num)
{
	if (num == 0 || num == 1)
		return 0;
	else
		return (int)(rand() >> 1) % num;
}

bool CheckSeq(int gPos, int mLen, string& seq)
{
	bool bRet = true;

	for (int i = 0; i < mLen; gPos++, i++)
	{
		if (seq[gPos] == 'N')
		{
			bRet = false;
			break;
		}
	}
	return bRet;
}

void GetComplementarySeq(int len, string& seq, string& rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i < j; i++, j--)
	{
		rseq[i] = ReverseMap[(int)seq[j]];
		rseq[j] = ReverseMap[(int)seq[i]];
	}
	if (i == j) rseq[i] = ReverseMap[(int)seq[i]];
}

string GenRandomInsertion(int len)
{
	string str;
	for (int i = 0; i < len; i++)
	{
		switch (myrandom(4))
		{
		case 0: str.push_back('A'); break;
		case 1: str.push_back('C'); break;
		case 2: str.push_back('G'); break;
		case 3: str.push_back('T'); break;
		}
	}
	return str;
}

void OutputMutantSeq(string chr, string ref_seq, map<int64_t, SV_t>& StrVarMap)
{
	string mut_seq;
	int gPos1, gPos2, chrLen, mutLen;
	map<int64_t, SV_t>::iterator StrVarMapIter;

	gPos1 = gPos2 = 0; chrLen = (int)ref_seq.length();
	for (StrVarMapIter = StrVarMap.begin(); StrVarMapIter != StrVarMap.end(); StrVarMapIter++)
	{
		gPos2 = StrVarMapIter->second.gPos;
		if (gPos2 > gPos1)
		{
			if (StrVarMapIter->second.mtype == 0) fprintf(vcf_fd, "%s	%d	.	%s	%s	30	PASS	SVTYPE=SUBSTITUTE\n", chr.c_str(), StrVarMapIter->second.gPos + 1, (char*)StrVarMapIter->second.ori_seq.c_str(), (char*)StrVarMapIter->second.mut_seq.c_str());
			else if (StrVarMapIter->second.mtype == 1 || StrVarMapIter->second.mtype == 2) fprintf(vcf_fd, "%s	%d	.	%s	%s	30	PASS	SVTYPE=%s\n", chr.c_str(), StrVarMapIter->second.gPos + 1, (char*)StrVarMapIter->second.ori_seq.c_str(), (char*)StrVarMapIter->second.mut_seq.c_str(), StrVarMapIter->second.ori_seq.length() < StrVarMapIter->second.mut_seq.length() ? "INSERT" : "DELETE");
			else if (StrVarMapIter->second.mtype == 3) fprintf(vcf_fd, "%s	%d	.	%c	<TRANSLOCATION>	30	PASS	SVTYPE=BND\n", chr.c_str(), StrVarMapIter->second.gPos + 1, StrVarMapIter->second.ori_seq[0], (int)StrVarMapIter->second.mut_seq.length());
			else if (StrVarMapIter->second.mtype == 4) fprintf(vcf_fd, "%s	%d	.	%c	<INV>	30	PASS	size=%d;SVTYPE=INVERSION\n", chr.c_str(), StrVarMapIter->second.gPos + 1, StrVarMapIter->second.ori_seq[0], (int)StrVarMapIter->second.mut_seq.length());
			else if (StrVarMapIter->second.mtype == 5) fprintf(vcf_fd, "%s	%d	.	%dx	%dx	30	PASS	SVTYPE=CNV\n", chr.c_str(), StrVarMapIter->second.gPos + 1, 1, (int)StrVarMapIter->second.mut_seq.length() / (int)StrVarMapIter->second.ori_seq.length());

			mut_seq.append(ref_seq.substr(gPos1, gPos2 - gPos1));
			mut_seq.append(StrVarMapIter->second.mut_seq);
			gPos1 = (gPos2 + (int)StrVarMapIter->second.ori_seq.length());
		}
		else
		{
			if (StrVarMapIter->second.mtype == 0) iSNP--;
			else if (StrVarMapIter->second.mtype == 1) isIND--;
			else if (StrVarMapIter->second.mtype == 2) ilIND--;
			else if (StrVarMapIter->second.mtype == 3) iTraLoc--;
			else if (StrVarMapIter->second.mtype == 4) iInv--;
			else if (StrVarMapIter->second.mtype == 5) iCNV--;
		}
	}
	if (gPos1 < chrLen) mut_seq.append(ref_seq.substr(gPos1, chrLen - gPos1));

	mutLen = (int)mut_seq.length();
	fprintf(stderr, "\tMutatnt (%s): len = %d (ori = %d)\n", chr.c_str(), mutLen, chrLen);

	fprintf(mut_fd, ">%s_mut\n", chr.c_str());
	for (gPos1 = 0; gPos1 < mutLen; gPos1 += 70) fprintf(mut_fd, "%s\n", mut_seq.substr(gPos1, 70).c_str());
}

void GenMutantSeq(string chr, string& ref_seq)
{
	SV_t sv;
	map<int64_t, SV_t>  StrVarMap;
	int i, ref_len, gPos, mPos, mLen, mType, Dup;

	fprintf(stderr, "Generate mutant sequence\n");
	ref_len = (int)ref_seq.length(); 
	for (gPos = 0; gPos < ref_len; gPos++)
	{
		if (gPos % 10000 == 0) fprintf(stderr, "\rScan %d / %d", gPos, ref_len);
		if (ref_seq[gPos] == 'N') continue;

		if (myrandom(DOM) < SNP_rate)
		{
			sv.mtype = 0; sv.gPos = gPos; sv.ori_seq = ref_seq.substr(gPos, 1);
			switch (ref_seq[gPos])
			{
			case 'A': sv.mut_seq = "T"; break;
			case 'C': sv.mut_seq = "G"; break;
			case 'G': sv.mut_seq = "C"; break;
			case 'T': sv.mut_seq = "A"; break;
			}
			iSNP++; StrVarMap.insert(make_pair(sv.gPos, sv)); gPos+=30;
		}
		else if (myrandom(DOM) < sIND_rate)
		{
			sv.mtype = 1; sv.gPos = gPos; mLen = 1; while (mLen < 10 && myrandom(10) == 0) mLen++;
			if (myrandom(2)) // ins
			{
				sv.ori_seq = ref_seq.substr(gPos, 1);
				sv.mut_seq = sv.ori_seq + GenRandomInsertion(mLen);
			}
			else // del
			{
				sv.mut_seq = ref_seq.substr(gPos, 1);
				sv.ori_seq = ref_seq.substr(gPos, mLen + 1);
				gPos += mLen;
			}
			isIND++; StrVarMap.insert(make_pair(sv.gPos, sv)); gPos+=30;
		}
		else if (myrandom(DOM) < lIND_rate)
		{
			sv.mtype = 2; sv.gPos = gPos; mLen = 11; while (mLen < 30 && myrandom(10) < 7) mLen++;
			if (myrandom(2)) // ins
			{
				sv.ori_seq = ref_seq.substr(gPos, 1);
				sv.mut_seq = sv.ori_seq + GenRandomInsertion(mLen);
			}
			else // del
			{
				sv.mut_seq = ref_seq.substr(gPos, 1);
				sv.ori_seq = ref_seq.substr(gPos, mLen + 1);
				gPos += mLen;
			}
			ilIND++; StrVarMap.insert(make_pair(sv.gPos, sv)); gPos+=30;
		}
		else if (myrandom(DOM) < Inv_rate)
		{
			sv.mtype = 4; sv.gPos = gPos; mLen = myrandom(1000) + 1000;
			if (gPos + mLen < ref_len)
			{
				sv.ori_seq = ref_seq.substr(gPos, mLen);
				sv.mut_seq.resize(mLen); GetComplementarySeq(mLen, sv.ori_seq, sv.mut_seq);
				iInv++; StrVarMap.insert(make_pair(sv.gPos, sv));
				gPos += mLen;
			}
		}
		else if (myrandom(DOM) < TraLoc_rate && myrandom(2))
		{
			sv.mtype = 3; mLen = myrandom(1000) + 1000; mPos = gPos + myrandom(1000) + 10000;
			if ((mPos + mLen) < ref_len)
			{
				sv.gPos = gPos;
				sv.ori_seq = ref_seq.substr(gPos, mLen);
				sv.mut_seq = ref_seq.substr(mPos, mLen);
				StrVarMap.insert(make_pair(sv.gPos, sv));

				sv.gPos = mPos;
				sv.ori_seq = ref_seq.substr(mPos, mLen);
				sv.mut_seq = ref_seq.substr(gPos, mLen);
				StrVarMap.insert(make_pair(sv.gPos, sv));

				iTraLoc += 2; gPos += mLen;
				for (i = 0; i < mLen; i++, mPos++) ref_seq[mPos] = 'N';
			}
		}
		else if (myrandom(DOM) < CNV_rate)
		{
			sv.mtype = 5; sv.gPos = gPos; mLen = myrandom(1000) + 300;
			if (gPos + mLen < ref_len && CheckSeq(gPos, mLen, ref_seq))
			{
				Dup = myrandom(100) % 8 + 2;
				sv.ori_seq = ref_seq.substr(gPos, mLen); sv.mut_seq = "";
				for (; Dup > 0; Dup--) sv.mut_seq += sv.ori_seq;
				iCNV++; StrVarMap.insert(make_pair(sv.gPos, sv));
				gPos += mLen;
			}
		}
	}
	fprintf(stderr, "\rScan %d / %d\n", gPos, ref_len);
	OutputMutantSeq(chr, ref_seq, StrVarMap);
}

int main(int argc, char* argv[])
{
	fstream file;
	string fname, chr, seq, str;

	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s ref_seq\n", argv[0]);
		exit(0);
	}
	srand((unsigned int)time(NULL));
	iSNP = isIND = ilIND = iTraLoc = iInv = iCNV = 0;

	fname = ((string)argv[1]).substr(0, ((string)argv[1]).find_last_of('.')) + ".vcf";
	vcf_fd = fopen((char*)fname.c_str(), "w");
	fprintf(vcf_fd, "##maf version=1\n");

	fname = ((string)argv[1]).substr(0, ((string)argv[1]).find_last_of('.')) + ".mut";
	mut_fd = fopen((char*)fname.c_str(), "w");

	fname = ((string)argv[1]).substr(0, ((string)argv[1]).find_last_of('.')) + ".info";
	info_fd = fopen((char*)fname.c_str(), "w");

	file.open(argv[1], ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		if (str[0] == '>')
		{
			if (seq != "") GenMutantSeq(chr, seq);
			chr = str.substr(1); seq = "";
		}
		else seq.append(str);
	}
	if (seq != "") GenMutantSeq(chr, seq);

	fprintf(stderr, "SNP=%d, sIND=%d, lIND=%d, Translocation=%d, Inversion=%d, CNV=%d\n", iSNP, isIND, ilIND, iTraLoc, iInv, iCNV);
	fprintf(info_fd, "SNP=%d, sIND=%d, lIND=%d, Translocation=%d, Inversion=%d, CNV=%d\n", iSNP, isIND, ilIND, iTraLoc, iInv, iCNV);
	fclose(info_fd); fclose(vcf_fd); fclose(mut_fd);

	return 0;
}
