#include "structure.h" 

static const char ReverseMap[255] =
{
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
	'\0', '\0', '\0', '\0', '\0', '-', '\0', '\0', '\0', '\0', /*  40 -  49 */
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

void GetComplementarySeq(int len, char* seq, char* rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i<j; i++, j--)
	{
		rseq[i] = ReverseMap[(int)seq[j]];
		rseq[j] = ReverseMap[(int)seq[i]];
	}
	if (i == j) rseq[i] = ReverseMap[(int)seq[i]];

	rseq[len] = '\0';
}

void SelfComplementarySeq(int len, char* seq)
{
	int i, j;
	char aa1, aa2;

	for (j = len - 1, i = 0; i<j; i++, j--)
	{
		aa1 = seq[i]; aa2 = seq[j];
		seq[i] = ReverseMap[(int)aa2]; seq[j] = ReverseMap[(int)aa1];
	}
	if (i == j) seq[i] = ReverseMap[(int)seq[i]];
}

int CalFragPairNonIdenticalBases(int len, char* frag1, char* frag2)
{
	int i, c;

	for (c = 0, i = 0; i < len; i++) if (frag1[i] != frag2[i]) c++;

	return c;
}

void ShowFragmentPair(char* ReadSeq, FragPair_t& fp)
{
	string frag1, frag2;
	frag1.resize(fp.rLen); strncpy((char*)frag1.c_str(), ReadSeq + fp.rPos, fp.rLen);
	frag2.resize(fp.gLen); strncpy((char*)frag2.c_str(), RefSequence + fp.gPos, fp.gLen);
	printf("FragmentPair:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\n\n", frag1.c_str(), fp.rPos, fp.rPos + fp.rLen - 1, fp.rLen, frag2.c_str(), (long long)fp.gPos, (long long)(fp.gPos + fp.gLen - 1), fp.gLen);
}

void ShowSimplePairInfo(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::const_iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->rLen == 0 && iter->gLen == 0) continue;
		else
		{
			//printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)iter->gPos, (long long)(iter->gPos + iter->gLen - 1), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			if (iter->gPos < GenomeSize) printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)iter->gPos, (long long)(iter->gPos + iter->gLen - 1), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			else printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s (Rev)\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)(TwoGenomeSize - (iter->gPos + iter->gLen)), (long long)(TwoGenomeSize - 1 - iter->gPos), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			
			if (iter->bSimple)
			{
				char* str = new char[iter->gLen + 1]; str[iter->gLen] = '\0';
				strncpy(str, RefSequence + iter->gPos, iter->gLen);
				if (iter->gPos >= GenomeSize) SelfComplementarySeq(iter->gLen, str);
				printf("\t\t%s\n", str);
				delete[] str;
			}
			else //if (iter->rLen > 0 || iter->gLen > 0)
			{
				printf("\t\t%s\n\t\t%s\n", iter->aln1.c_str(), iter->aln2.c_str());
			}
		}
	}
	printf("\n"); fflush(stdout);
}

void ShowSeedLocationInfo(int64_t MyPos)
{
	int64_t gPos;

	map<int64_t, int>::const_iterator iter = ChrLocMap.lower_bound(MyPos);
	if (MyPos < GenomeSize) gPos = MyPos - ChromosomeVec[iter->second].FowardLocation;
	else gPos = iter->first - MyPos;
	printf("\t\tChr [%s, %ld]\n", ChromosomeVec[iter->second].name, gPos);
}

int64_t GetAlignmentBoundary(int64_t gPos)
{
	map<int64_t, int>::iterator iter = ChrLocMap.lower_bound(gPos);

	return iter->first;
}

bool CheckFragValidity(FragPair_t FragPair)
{
	map<int64_t, int>::iterator iter1, iter2;

	iter1 = ChrLocMap.lower_bound(FragPair.gPos);
	iter2 = ChrLocMap.lower_bound(FragPair.gPos + FragPair.gLen - 1);

	return (iter1->first == iter2->first);
}

bool CheckAlignmentValidity(vector<FragPair_t>& FragPairVec)
{
	map<int64_t, int>::iterator iter1, iter2;

	iter1 = ChrLocMap.lower_bound(FragPairVec.begin()->gPos);
	iter2 = ChrLocMap.lower_bound(FragPairVec.rbegin()->gPos + FragPairVec.rbegin()->gLen - 1);

	return (iter1->first == iter2->first);
}

Coordinate_t DetermineCoordinate(int64_t gPos)
{
	Coordinate_t coor;

	if (gPos < GenomeSize)
	{
		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = gPos + 1;
		}
		else
		{
			map<int64_t, int>::iterator iter = ChrLocMap.lower_bound(gPos);
			coor.ChromosomeIdx = iter->second;
			coor.gPos = gPos + 1 - ChromosomeVec[coor.ChromosomeIdx].FowardLocation;
		}
	}
	else
	{
		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = TwoGenomeSize - gPos;
		}
		else
		{
			map<int64_t, int>::iterator iter = ChrLocMap.lower_bound(gPos);
			coor.gPos = iter->first - gPos + 1; coor.ChromosomeIdx = iter->second;
		}
	}
	return coor;
}

int GetProfileColumnSize(MappingRecord_t& Profile)
{
	return Profile.A + Profile.C + Profile.G + Profile.T;
}

void ShowProfileColumn(int64_t gPos)
{
	int cov = GetProfileColumnSize(MappingRecordArr[gPos]);
	printf("%lld[%c]: cov=%d [A=%d C=%d G=%d T=%d] dup=%d\n", gPos, RefSequence[gPos], cov, MappingRecordArr[gPos].A, MappingRecordArr[gPos].C, MappingRecordArr[gPos].G, MappingRecordArr[gPos].T, MappingRecordArr[gPos].multi_hit);
}

void ShowVariationProfile(int64_t begin_pos, int64_t end_pos)
{
	int64_t gPos = (begin_pos + end_pos) / 2;
	Coordinate_t coor = DetermineCoordinate(gPos);
	if (end_pos >= GenomeSize) end_pos = GenomeSize - 1;
	printf("%s-%lld\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos);
	for (int64_t gPos = begin_pos; gPos <= end_pos; gPos++) ShowProfileColumn(gPos);
	printf("\n\n");
}

void ShowIndSeq(int64_t begin_pos, int64_t end_pos)
{
	for (map<int64_t, map<string, uint16_t> >::iterator iter = InsertSeqMap.begin(); iter != InsertSeqMap.end(); iter++)
	{
		if (iter->first >= begin_pos && iter->first <= end_pos)
		{
			for (map<string, uint16_t>::iterator SeqMapIter = iter->second.begin(); SeqMapIter != iter->second.end(); SeqMapIter++)
				fprintf(stdout, "INS:%lld	[%s] freq=%d\n", (long long)iter->first, (char*)SeqMapIter->first.c_str(), SeqMapIter->second);
		}
	}
	for (map<int64_t, map<string, uint16_t> >::iterator iter = DeleteSeqMap.begin(); iter != DeleteSeqMap.end(); iter++)
	{
		if (iter->first >= begin_pos && iter->first < end_pos)
		{
			for (map<string, uint16_t>::iterator SeqMapIter = iter->second.begin(); SeqMapIter != iter->second.end(); SeqMapIter++)
				fprintf(stdout, "DEL:%lld	%d	[%s]\n", (long long)iter->first, SeqMapIter->second, (char*)SeqMapIter->first.c_str());
		}
	}
}
