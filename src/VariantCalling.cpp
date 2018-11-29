#include "structure.h"

#define MaxQscore 30
#define FrequencyThr 0.33
#define minCNVlength 200
#define INV_TNL_ThrRatio 0.5
#define var_SUB 0 // substitution
#define var_INS 1 // insertion
#define var_DEL 2 // deletion
#define var_UMR 3 // unmapped region
#define var_CNV 4 // copy number variation
#define var_INV 5 // inversion
#define var_TNL 6 // translocation
#define var_NIL 255

typedef struct
{
	int64_t gPos;
	uint16_t left_score;
	uint16_t rigt_score;
} BreakPoint_t;

FILE *outFile;
vector<VarPos_t> VarPosVec;
static pthread_mutex_t Lock;
vector<BreakPoint_t> BreakPointCanVec;
extern map<int64_t, uint16_t> BreakPointMap;
 // TranslocationSiteVec;
extern uint32_t avgCov, lowCov, minCov, maxCov, avgReadLength;

extern bool CompByDiscordPos(const DiscordPair_t& p1, const DiscordPair_t& p2);

bool CompByAlnDist(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
{
	return p1.dist < p2.dist;
}

bool CompByDist(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
{
	return p1.dist < p2.dist;
}

bool CompByFirstCoor(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
{
	return p1.gPos1 < p2.gPos1;
}

bool CompByVarPos(const VarPos_t& p1, const VarPos_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.type < p2.type;
	else return p1.gPos < p2.gPos;
}

pair<char, int> GetMaxItemInProfileColumn(MappingRecord_t& Profile)
{
	pair<char, int> p = make_pair('N', 0);

	if (Profile.A > p.second) p = make_pair('A', Profile.A);
	if (Profile.C > p.second) p = make_pair('C', Profile.C);
	if (Profile.G > p.second) p = make_pair('G', Profile.G);
	if (Profile.T > p.second) p = make_pair('T', Profile.T);

	return p;
}

int CheckDiploid(int cov, MappingRecord_t& Profile)
{
	int hybrid_cov, thr;
	
	if((thr = (int)(cov * FrequencyThr)) < 5) thr = 5;
	
	hybrid_cov = 0;
	if (Profile.A > thr) hybrid_cov += Profile.A;
	if (Profile.C > thr) hybrid_cov += Profile.C;
	if (Profile.G > thr) hybrid_cov += Profile.G;
	if (Profile.T > thr) hybrid_cov += Profile.T;

	return hybrid_cov;
}

int GetIndMapFrq(map<int64_t, map<string, uint16_t> >::iterator IndMapIter, int& ind_length)
{
	int freq = 0, max_freq = 0;

	ind_length = 0;
	for (map<string, uint16_t>::iterator iter = IndMapIter->second.begin(); iter != IndMapIter->second.end(); iter++)
	{
		if (max_freq < iter->second)
		{
			max_freq = iter->second;
			ind_length = (int)iter->first.length();
		}
		freq += iter->second;
	}
	return freq;
}

float CheckRefBaseFreq(int cov, char base, MappingRecord_t& Profile)
{
	int freq;
	switch (base)
	{
	case 'A': freq = Profile.A; break;
	case 'C': freq = Profile.C; break;
	case 'G': freq = Profile.G; break;
	case 'T': freq = Profile.T; break;
	default: freq = 0;
	}
	return 1.*freq / cov;
}

int CalAvgCov(int64_t gPos)
{
	int64_t p;
	int i, cov_sum = 0;

	for (i = 1; i <= 25; i++)
	{
		if ((p = gPos + i) < GenomeSize) cov_sum += GetProfileColumnSize(MappingRecordArr[p]);
		if ((p = gPos - i) > 0) cov_sum += GetProfileColumnSize(MappingRecordArr[p]);
	}
	return cov_sum /= 50;
}

uint8_t CalQualityScore(int a, int b)
{
	uint8_t qs;
	if (a >= b) qs = MaxQscore;
	else if ((qs= -100 * log10((1.0 - (1.0*a / b)))) > MaxQscore) qs = MaxQscore;

	return qs;
}

void *IdentifyVariants(void *arg)
{
	VarPos_t VarPos;
	pair<char, int> p;
	int64_t gPos, end;
	vector<int64_t> myDupVec;
	vector<VarPos_t> MyVarPosVec;
	map<int64_t, map<string, uint16_t> >::iterator IndMapIter;
	int i, n, head_idx, num, cov, avg_cov, pre_cov, diploid_cov, ins_freq, ins_len, del_freq, del_len, dup_len, weight, tid = *((int*)arg);

	myDupVec.push_back(0);
	gPos = (tid == 0 ? 0 : (GenomeSize / iThreadNum)*tid); end = (tid == iThreadNum - 1 ? GenomeSize : (GenomeSize / iThreadNum)*(tid + 1));
	
	for (pre_cov = 0; gPos < end; gPos++)
	{
		cov = GetProfileColumnSize(MappingRecordArr[gPos]); //MappingRecordArr[gPos].count;

		if ((cov + MappingRecordArr[gPos].multi_hit) >= maxCov) myDupVec.push_back(gPos);

		if ((IndMapIter = InsertSeqMap.find(gPos)) != InsertSeqMap.end()) ins_freq = GetIndMapFrq(IndMapIter, ins_len); else ins_freq = 0;
		if ((IndMapIter = DeleteSeqMap.find(gPos)) != DeleteSeqMap.end()) del_freq = GetIndMapFrq(IndMapIter, del_len); else del_freq = 0;

		//if (ins_freq >= MinBaseDepth && (cov + ins_freq) >= minCov && (weight = ins_freq * (int)ceil(ins_len / 5.0)) >= (int)(cov*0.3))
		if (ins_freq >= MinBaseDepth)
		{
			weight = ins_freq * (int)ceil(ins_len / 5.0);
			VarPos.type = var_INS; VarPos.gPos = gPos; VarPos.DP = cov; VarPos.NS = ins_freq;
			//if (weight >= cov || (VarPos.qscore = -100 * log10((1.0 - (1.0*weight / cov)))) > MaxQscore) VarPos.qscore = MaxQscore;
			VarPos.qscore = CalQualityScore(weight, cov);
			if (VarPos.qscore >= MinVarConfScore) MyVarPosVec.push_back(VarPos);
		}
		//if (del_freq >= MinBaseDepth && (pre_cov + del_freq) >= minCov && (weight = del_freq * (int)ceil(del_len / 5.0)) >= (int)(pre_cov*0.3))
		if (del_freq >= MinBaseDepth)
		{
			weight = del_freq * (int)ceil(del_len / 5.0);
			VarPos.type = var_DEL; VarPos.gPos = gPos; VarPos.DP = cov; VarPos.NS = del_freq;
			//if (weight >= pre_cov || (VarPos.qscore = -100 * log10((1.0 - (1.0*weight / pre_cov)))) > MaxQscore) VarPos.qscore = MaxQscore;
			VarPos.qscore = CalQualityScore(weight, pre_cov);
			if (VarPos.qscore >= MinVarConfScore) MyVarPosVec.push_back(VarPos);
		}
		if ((cov >= MinBaseDepth || cov >= minCov) && cov > ins_freq && cov > del_freq)
		{
			if ((p = GetMaxItemInProfileColumn(MappingRecordArr[gPos])).first != RefSequence[gPos])
			{
				VarPos.gPos = gPos; VarPos.type = var_SUB; VarPos.DP = cov; VarPos.NS = p.second;
				//if (p.second >= (int)(cov*0.8) || ((VarPos.qscore = -100 * log10((1.0 - (1.0*p.second / cov))))) > MaxQscore) VarPos.qscore = MaxQscore;
				VarPos.qscore = CalQualityScore(p.second, cov);
				if (VarPos.qscore >= MinVarConfScore) MyVarPosVec.push_back(VarPos);
			}
			else if ((diploid_cov = CheckDiploid(cov, MappingRecordArr[gPos])) > p.second)
			{
				VarPos.gPos = gPos; VarPos.type = var_SUB; VarPos.DP = cov; VarPos.NS = diploid_cov;
				//if (diploid_cov >= (int)(cov*0.8) || (VarPos.qscore = -100 * log10((1.0 - (1.0*diploid_cov / cov)))) > MaxQscore) VarPos.qscore = MaxQscore;
				VarPos.qscore = CalQualityScore(diploid_cov, cov);
				if (VarPos.qscore >= MinVarConfScore) MyVarPosVec.push_back(VarPos);
			}
		}
		pre_cov = cov;
	}
	myDupVec.push_back(TwoGenomeSize);

	//Process myDupVec
	for (n = 0, num = (int)myDupVec.size(), head_idx = 0, i = 1; i < num; i++)
	{
		if (myDupVec[i] - myDupVec[i - 1] > FragmentSize) // break
		{
			if ((dup_len = (uint32_t)(myDupVec[i - 1] - myDupVec[head_idx] + 1)) >= minCNVlength && (i - head_idx) > (dup_len >> 1))
			{
				n++;
				VarPos.type = var_CNV;
				VarPos.gPos = myDupVec[head_idx]; // first_pos;
				VarPos.gEnd = myDupVec[i - 1]; // last_pos;
				//VarPos.qscore = (uint8_t)(-10 * log10((1.0 - (1.0*(i - head_idx) / VarPos.freq))));
				VarPos.qscore = CalQualityScore((i - head_idx), VarPos.gEnd);
				//printf("cnv: %lld - %lld (n = %d, size=%d, qscore=%d)\n", myDupVec[head_idx], myDupVec[i - 1], i - head_idx, VarPos.freq, VarPos.qscore);
				MyVarPosVec.push_back(VarPos);
			}
			head_idx = i;
		}
	}
	if (n > 0) inplace_merge(MyVarPosVec.begin(), MyVarPosVec.end() - n, MyVarPosVec.end(), CompByVarPos);

	if ((num = (int)MyVarPosVec.size()) > 0)
	{
		pthread_mutex_lock(&Lock);
		copy(MyVarPosVec.begin(), MyVarPosVec.end(), back_inserter(VarPosVec)); inplace_merge(VarPosVec.begin(), VarPosVec.end() - num, VarPosVec.end(), CompByVarPos);
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
}

void ShowMetaInfo()
{
	fprintf(outFile, "##fileformat=VCFv4.3\n");
	fprintf(outFile, "##reference=%s\n", IndexFileName);
	fprintf(outFile, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n");
	fprintf(outFile, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n");
	for (int i = 0; i < iChromsomeNum; i++) fprintf(outFile, "##Contig=<ID=%s,length=%d>\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	fprintf(outFile, "##INFO=<ID=TYPE,Type=String,Description=\"The type of allele, either SUBSTITUTE, INSERT, DELETE, DUPLICATION,or BND.\">\n");
	fprintf(outFile, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
}

void IdentifyBreakPointCandidates()
{
	BreakPoint_t bp;
	uint32_t total_freq;
	pair<int64_t, uint16_t> p;

	BreakPointMap.insert(make_pair(TwoGenomeSize, 0)); total_freq = 0; p = make_pair(0, 0);
	for (map<int64_t, uint16_t>::iterator iter = BreakPointMap.begin(); iter != BreakPointMap.end(); iter++)
	{
		if (iter->first - p.first > avgReadLength)
		{
			//printf("Break!\n");
			if (total_freq > lowCov)
			{
				bp.gPos = p.first;
				bp.left_score = bp.rigt_score = 0;
				//printf("\tAdd BP_can: %lld (score=%d)\n", (long long)bp.gPos, total_freq);
				BreakPointCanVec.push_back(bp);
			}
			p.first = iter->first;
			total_freq = p.second = iter->second;
		}
		else
		{
			total_freq += iter->second;
			if (p.second < iter->second)
			{
				p.first = iter->first;
				p.second = iter->second;
			}
		}
		//printf("Pos=%lld freq=%d\n", (long long)iter->first, iter->second);
	}
}

int CalRegionCov(int64_t begPos, int64_t endPos)
{
	int cov = 0;
	int64_t gPos;

	for (gPos = begPos; gPos <= endPos; gPos++) cov += GetProfileColumnSize(MappingRecordArr[gPos]);

	if (endPos >= begPos) return (int)(cov / (endPos - begPos + 1));
	else return 0;
}

void IdentifyTranslocations()
{
	int64_t gPos;
	VarPos_t VarPos;
	vector<int64_t> vec;
	DiscordPair_t DiscordPair;
	vector<DiscordPair_t>::iterator Iter1, Iter2;
	uint32_t i, j, n, TNLnum, num, score, LCov, RCov, Lscore, Rscore;

	//for (Iter1 = TranslocationSiteVec.begin(); Iter1 != TranslocationSiteVec.end(); Iter1++)  printf("Pos=%lld Dist=%lld\n", (long long)Iter1->gPos, (long long)Iter1->dist);
	for (num = (int)BreakPointCanVec.size(), TNLnum = i = 0; i < num; i++)
	{
		gPos = BreakPointCanVec[i].gPos; LCov = CalRegionCov(gPos - FragmentSize, gPos - (avgReadLength >> 1));
		DiscordPair.gPos = gPos - FragmentSize; Iter1 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos - (avgReadLength >> 1); Iter2 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000))/*, printf("gPos=%lld, dist=%lld\n", Iter1->gPos, Iter1->dist)*/; 
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize);
		for (n = (int)vec.size(), Lscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Lcan_%d: %lld (score=%d)\n", j + 1, vec[j], score);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Lscore) Lscore = score;
				score = 1;
			}
			else score++;
		}
		if (Lscore < minCov || Lscore < (int)(LCov*INV_TNL_ThrRatio)) Lscore = 0;

		RCov = CalRegionCov(gPos, gPos + FragmentSize);
		DiscordPair.gPos = gPos; Iter1 = upper_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos + FragmentSize; Iter2 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000))/*, printf("gPos=%lld, dist=%lld\n", Iter1->gPos, Iter1->dist)*/; 
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize);
		for (n = (int)vec.size(), Rscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Rcan_%d: %lld (score=%d)\n", j + 1, vec[j], score);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Rscore) Rscore = score;
				score = 1;
			}
			else score++;
		}
		if (Rscore < minCov || Rscore < (int)(RCov*INV_TNL_ThrRatio)) Rscore = 0;

		//printf("TNL_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			TNLnum++;
			VarPos.gPos = gPos;
			VarPos.type = var_TNL;
			VarPos.NS = Lscore > Rscore ? Lscore : Rscore;
			//VarPos.qscore = (VarPos.freq > avgCov ? MaxQscore : -100 * log10((1.0 - (1.0*VarPos.freq / avgCov))));
			VarPos.qscore = CalQualityScore(VarPos.NS, avgCov);
			VarPosVec.push_back(VarPos);
		}
	}
	if (TNLnum > 0) inplace_merge(VarPosVec.begin(), VarPosVec.end() - TNLnum, VarPosVec.end(), CompByVarPos);
}

void IdentifyInversions()
{
	int64_t gPos;
	VarPos_t VarPos;
	vector<int64_t> vec;
	DiscordPair_t DiscordPair;
	vector<DiscordPair_t>::iterator Iter1, Iter2;
	uint32_t i, j, n, LCov, RCov, INVnum, num, score, Lscore, Rscore;

	//for (Iter1 = InversionSiteVec.begin(); Iter1 != InversionSiteVec.end(); Iter1++) printf("Pos=%lld Dist=%lld\n", (long long)Iter1->gPos, (long long)Iter1->dist);
	for (num = (int)BreakPointCanVec.size(), INVnum = i = 0; i < num; i++)
	{
		gPos = BreakPointCanVec[i].gPos; LCov = CalRegionCov(gPos - FragmentSize, gPos - (avgReadLength >> 1));
		DiscordPair.gPos = gPos - FragmentSize; Iter1 = lower_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos - (avgReadLength >> 1); Iter2 = lower_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);

		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000))/*, printf("gPos=%lld, dist=%lld\n", Iter1->gPos, Iter1->dist)*/;
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize);
		for (n = (int)vec.size(), Lscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Lcan_%d: dist=%lld\n", j + 1, vec[j]);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Lscore) Lscore = score;
				score = 1;
			}
			else score++;
		}
		if (Lscore < minCov || Lscore < (int)(LCov*INV_TNL_ThrRatio)) Lscore = 0;

		RCov = CalRegionCov(gPos, gPos + FragmentSize);
		DiscordPair.gPos = gPos; Iter1 = upper_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos + FragmentSize; Iter2 = lower_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000))/*, printf("gPos=%lld, dist=%lld\n", Iter1->gPos, Iter1->dist)*/;
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize);
		for (n = (int)vec.size(), Rscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Rcan_%d: dist=%lld\n", j + 1, vec[j]);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Rscore) Rscore = score;
				score = 1;
			}
			else score++;
		}
		if (Rscore < minCov || Rscore < (int)(RCov*INV_TNL_ThrRatio)) Rscore = 0;

		//printf("INV_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			//printf("YES!\n");
			INVnum++;
			VarPos.gPos = gPos;
			VarPos.NS = Lscore > Rscore ? Lscore : Rscore;
			VarPos.type = var_INV;
			//VarPos.qscore = (score > avgCov ? MaxQscore : -100 * log10((1.0 - (1.0*VarPos.freq / avgCov))));
			VarPos.qscore = CalQualityScore(VarPos.NS, avgCov);
			VarPosVec.push_back(VarPos);
		}
	}
	if (INVnum > 0) inplace_merge(VarPosVec.begin(), VarPosVec.end() - INVnum, VarPosVec.end(), CompByVarPos);
}

void GenVariantCallingFile()
{
	int64_t gPos;
	string ALTstr;
	pair<char, int> p;
	Coordinate_t coor;
	map<int, int> VarNumMap;
	int i, len, num, thr, cov, dupN;
	map<string, uint16_t>::iterator IndSeqMapIter;

	outFile = fopen(VcfFileName, "w"); ShowMetaInfo();

	for (num = (int)VarPosVec.size(), i = 0; i < num; i++)
	{
		gPos = VarPosVec[i].gPos; coor = DetermineCoordinate(gPos); ALTstr.clear();

		if (VarPosVec[i].type == var_SUB)
		{
			VarNumMap[var_SUB]++;
			thr = GetProfileColumnSize(MappingRecordArr[gPos]) * FrequencyThr; p = GetMaxItemInProfileColumn(MappingRecordArr[gPos]);
			
			ALTstr.push_back(p.first);
			if (MappingRecordArr[gPos].A > thr && p.first != 'A') ALTstr += ",A";
			if (MappingRecordArr[gPos].C > thr && p.first != 'C') ALTstr += ",C";
			if (MappingRecordArr[gPos].G > thr && p.first != 'G') ALTstr += ",G";
			if (MappingRecordArr[gPos].T > thr && p.first != 'T') ALTstr += ",T";
			fprintf(outFile, "%s	%d	.	%c	%s	%d	%s	DP=%d;NS=%d;TYPE=SUBSTITUTE\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[VarPosVec[i].gPos], ALTstr.c_str(), VarPosVec[i].qscore, (VarPosVec[i].qscore >= 10 ? "PASS" : "q10"), VarPosVec[i].DP, VarPosVec[i].NS);
		}
		else if (VarPosVec[i].type == var_INS)
		{
			thr = (int)(VarPosVec[i].NS*FrequencyThr); ALTstr.clear();
			for (IndSeqMapIter = InsertSeqMap[gPos].begin(); IndSeqMapIter != InsertSeqMap[gPos].end(); IndSeqMapIter++)
			{
				if (IndSeqMapIter->second > thr) ALTstr += (RefSequence[gPos - 1] + IndSeqMapIter->first + ",");
			}
			if ((len = (int)ALTstr.length()) > 0)
			{
				VarNumMap[var_INS]++;
				ALTstr[len - 1] = '\0';
				fprintf(outFile, "%s	%d	.	%c	%s	%d	%s	DP=%d;NS=%d;TYPE=INSERT\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos - 1, RefSequence[gPos - 1], ALTstr.c_str(), VarPosVec[i].qscore, (VarPosVec[i].qscore >= 10 ? "PASS" : "q10"), VarPosVec[i].DP, VarPosVec[i].NS);
			}
		}
		else if (VarPosVec[i].type == var_DEL)
		{
			thr = (int)(VarPosVec[i].NS*FrequencyThr); ALTstr.clear();
			for (IndSeqMapIter = DeleteSeqMap[gPos].begin(); IndSeqMapIter != DeleteSeqMap[gPos].end(); IndSeqMapIter++)
			{
				if (IndSeqMapIter->second > thr) ALTstr += (RefSequence[gPos - 1] + IndSeqMapIter->first + ",");
			}
			if ((len = (int)ALTstr.length()) > 0)
			{
				VarNumMap[var_DEL]++;
				ALTstr[len - 1] = '\0';
				fprintf(outFile, "%s	%d	.	%s	%c	%d	%s	DP=%d;NS=%d;TYPE=DELETE\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos - 1, ALTstr.c_str(), RefSequence[gPos - 1], VarPosVec[i].qscore, (VarPosVec[i].qscore >= 10 ? "PASS" : "q10"), VarPosVec[i].DP, VarPosVec[i].NS);
			}
		}
		else if (VarPosVec[i].type == var_CNV)
		{
			for (cov = 0; gPos <= VarPosVec[i].gEnd; gPos++) cov += (GetProfileColumnSize(MappingRecordArr[gPos]) + MappingRecordArr[gPos].multi_hit);
			len = VarPosVec[i].gEnd - VarPosVec[i].gPos + 1; dupN = (int)(1.0*cov / len / avgCov + .5);
			if (dupN > 1)
			{
				VarNumMap[var_CNV]++;
				fprintf(outFile, "%s	%d	.	.	.	30	PASS	%dx;size=%d;TYPE=DUPLICATION\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, dupN, len);
			}
		}
		else if (VarPosVec[i].type == var_TNL)
		{
			VarNumMap[var_TNL]++;
			fprintf(outFile, "%s	%d	.	%c	<TRANSLOCATION>	30	PASS	DP=%d;NS=%d;TYPE=BND\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos - 1], VarPosVec[i].DP, VarPosVec[i].NS);
		}
		else if (VarPosVec[i].type == var_INV)
		{
			VarNumMap[var_INV]++;
			fprintf(outFile, "%s	%d	.	%c	<INVERSION>	30	PASS	DP=%d;NS=%d;TYPE=BND\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos - 1], VarPosVec[i].DP, VarPosVec[i].NS);
		}
		//else if (VarPosVec[i].type == var_UMR)
		//{
		//	VarNumMap[var_UMR]++;
		//	fprintf(outFile, "%s	%d	.	.	.	%d	%s	size=%d;SVTYPE=UnMapped\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, VarPosVec[i].qscore, (VarPosVec[i].qscore >= 10 ? "PASS" : "q10"), VarPosVec[i].freq);
		//}
	}
	std::fclose(outFile);
	fprintf(stderr, "\t%d(snp); %d(ins); %d(del); %d(trans); %d(inversion); %d(repeats)\n\n", VarNumMap[var_SUB], VarNumMap[var_INS], VarNumMap[var_DEL], VarNumMap[var_TNL] >> 1, VarNumMap[var_INV] >> 1, VarNumMap[var_CNV]);
}

void VariantCalling()
{
	int i, *ThrIDarr;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	ThrIDarr = new int[iThreadNum];  for (i = 0; i < iThreadNum; i++) ThrIDarr[i] = i;

	if (avgCov < 5) fprintf(stderr, "Too Low sequence coverage to identify sample variations\n");
	else
	{
		//iThreadNum = 1;
		fprintf(stderr, "Identify all variants...\n");
		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyVariants, &ThrIDarr[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
		IdentifyBreakPointCandidates();
		if (BreakPointCanVec.size() > 0 && InversionSiteVec.size() > 0) IdentifyInversions();
		if (BreakPointCanVec.size() > 0 && TranslocationSiteVec.size() > 0) IdentifyTranslocations();
		fprintf(stderr, "Write all the predicted sample variations to file [%s]...\n", VcfFileName); GenVariantCallingFile();
	}
	fprintf(stderr, "SV calling has be done in %lld seconds.\n", (long long)(time(NULL) - StartProcessTime));
	delete[] ThrIDarr; delete[] ThreadArr;
}
