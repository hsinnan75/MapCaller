#include "structure.h"

#define MaxQscore 30
#define BlockSize 100
#define INV_TNL_ThrRatio 0.5
#define var_SUB 0 // substitution
#define var_INS 1 // insertion
#define var_DEL 2 // deletion
#define var_INV 3 // inversion
#define var_TNL 4 // translocation
#define var_CNV 5 // copy number variation
#define var_UMR 6 // unmapped region
#define var_NIL 255

typedef struct
{
	int64_t gPos;
	uint16_t left_score;
	uint16_t rigt_score;
} BreakPoint_t;

int BlockNum;
FILE *outFile;
int* BlockDepthArr;
vector<Variant_t> VariantVec;
static pthread_mutex_t Lock;
vector<BreakPoint_t> BreakPointCanVec;
extern map<int64_t, uint16_t> BreakPointMap;
 // TranslocationSiteVec;
extern float FrequencyThr;
extern uint32_t avgReadLength;

extern bool CompByDiscordPos(const DiscordPair_t& p1, const DiscordPair_t& p2);

bool CompByAlnDist(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
{
	return p1.dist < p2.dist;
}

bool CompByDist(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
{
	return p1.dist < p2.dist;
}

//bool CompByFirstCoor(const CoordinatePair_t& p1, const CoordinatePair_t& p2)
//{
//	return p1.gPos1 < p2.gPos1;
//}

bool CompByVarPos(const Variant_t& p1, const Variant_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.VarType < p2.VarType;
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

//int GetIndMapFrq(map<int64_t, map<string, uint16_t> >::iterator IndMapIter, string& ind_str)
//{
//	int freq = 0, max_freq = 0;
//	for (map<string, uint16_t>::iterator iter = IndMapIter->second.begin(); iter != IndMapIter->second.end(); iter++)
//	{
//		freq += iter->second;
//		if (max_freq < iter->second)
//		{
//			ind_str = iter->first;
//			max_freq = iter->second;
//		}
//	}
//	return freq;
//}

int GetPointIndFreq(map<string, uint16_t>& IndMap)
{
	int n = 0;
	for (map<string, uint16_t>::iterator iter = IndMap.begin(); iter != IndMap.end(); iter++) n += iter->second;
	return n;
}

int GetAreaIndFrequency(int64_t gPos, map<int64_t, map<string, uint16_t> >& IndMap, string& ind_str)
{
	int64_t max_pos;
	int freq = 0, max_freq = 0;
	map<string, uint16_t>::iterator it;
	map<int64_t, map<string, uint16_t> >::iterator iter1, iter2;
	
	ind_str.clear();
	for (iter1 = IndMap.lower_bound(gPos - 5), iter2 = IndMap.upper_bound(gPos + 5); iter1 != iter2; iter1++)
	{
		if (abs(iter1->first - gPos) <= 5)
		{
			for (it = iter1->second.begin(); it != iter1->second.end(); it++)
			{
				freq += it->second;
				if (max_freq < it->second)
				{
					ind_str = it->first; 
					max_freq = it->second;
					max_pos = iter1->first;
				}
				else if (max_freq == it->second && it->first.length() > ind_str.length())
				{
					ind_str = it->first;
					max_pos = iter1->first;
				}
			}
		}
	}
	if (gPos == max_pos) return freq;
	else return 0;
}

uint8_t CalQualityScore(int a, int b)
{
	uint8_t qs;
	if (a >= b) qs = MaxQscore;
	else if ((qs= -100 * log10((1.0 - (1.0*a / b)))) > MaxQscore) qs = MaxQscore;

	return qs;
}

void *CalBlockReadDepth(void *arg)
{
	int64_t gPos, end_gPos;
	int bid, end_bid, sum, tid = *((int*)arg);

	bid = (tid == 0 ? 0 : (int)(BlockNum / iThreadNum)*tid); 
	end_bid = (tid == iThreadNum - 1 ? BlockNum : (BlockNum / iThreadNum)*(tid + 1));
	for (; bid < end_bid; bid++)
	{
		gPos = (int64_t)bid*BlockSize; if ((end_gPos = gPos + BlockSize) > GenomeSize) end_gPos = GenomeSize;
		for (sum = 0; gPos < end_gPos; gPos++) sum += GetProfileColumnSize(MappingRecordArr[gPos]);
		if (sum > 0) BlockDepthArr[bid] = sum / BlockSize;
	}
	return (void*)(1);
}

bool CheckDiploidFrequency(int cov, vector<pair<char, int> >& vec)
{
	int sum = vec[0].second + vec[1].second;
	if (sum >= (int)(cov*0.9)) return true;
	else return false;
}

//int CalPresentNT(MappingRecord_t& profile)
//{
//	int n = 0;
//	if (profile.A > 0) n++;
//	if (profile.C > 0) n++;
//	if (profile.G > 0) n++;
//	if (profile.T > 0) n++;
//
//	return n;
//}

map<int64_t, bool> LoadObservedPos()
{
	int64_t p;
	string str;
	fstream file;
	stringstream ss;
	map<int64_t, bool> m;

	file.open("a.txt", ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		ss.clear(); ss.str(str); ss >> p;
		//m.insert(make_pair(p - 1, (str[str.length()-1] == '+' ? true : false)));
		m.insert(make_pair(p, (str.find('+') != string::npos ? true : false)));
		if (m.size() == 500) break;
	}
	file.close();
	return m;
}

bool CheckObsMap(int64_t gPos, map<int64_t, bool>& obs_map)
{
	map<int64_t, bool>::iterator iter1, iter2;

	iter1 = obs_map.lower_bound(gPos - 10);
	iter2 = obs_map.lower_bound(gPos + 10);
	if (abs(iter1->first - gPos) <= 10 || abs(iter2->first - gPos) <= 10) return true;
	else return false;
}

bool CheckBreakPoints(int64_t gPos)
{
	map<int64_t, uint16_t>::iterator iter1, iter2;
	
	iter1 = BreakPointMap.lower_bound(gPos - 10);
	iter2 = BreakPointMap.lower_bound(gPos + 10);
	if (abs(iter1->first - gPos) <= 10 || abs(iter2->first - gPos) <= 10) return true;
	else return false;
}

void ShowMetaInfo()
{
	fprintf(outFile, "##fileformat=VCFv4.3\n");
	fprintf(outFile, "##reference=%s\n", IndexFileName);
	fprintf(outFile, "##source=MapCaller %s\n", VersionStr);
	fprintf(outFile, "##CommandLine=<%s>\n", CmdLine.c_str());
	fprintf(outFile, "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allel depth\">\n");
	fprintf(outFile, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n");
	fprintf(outFile, "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency\">\n");
	fprintf(outFile, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(outFile, "##FILTER=<ID=q%d,Description=\"Quality below %d\">\n", MinVarConfScore, MinVarConfScore);
	fprintf(outFile, "##INFO=<ID=TYPE,Type=String,Description=\"The type of allele, either SUBSTITUTE, INSERT, DELETE, or BND.\">\n");
	for (int i = 0; i < iChromsomeNum; i++) fprintf(outFile, "##Contig=<ID=%s,length=%d>\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
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
			if (total_freq > MinAlleleFreq)
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
	Variant_t Variant;
	vector<int64_t> vec;
	DiscordPair_t DiscordPair;
	vector<DiscordPair_t>::iterator Iter1, Iter2;
	uint32_t i, j, n, TNLnum, num, score, LCov, RCov, cov_thr, Lscore, Rscore;

	//for (Iter1 = TranslocationSiteVec.begin(); Iter1 != TranslocationSiteVec.end(); Iter1++)  printf("Pos=%lld Dist=%lld\n", (long long)Iter1->gPos, (long long)Iter1->dist);
	for (num = (int)BreakPointCanVec.size(), TNLnum = i = 0; i < num; i++)
	{
		gPos = BreakPointCanVec[i].gPos; LCov = CalRegionCov(gPos - FragmentSize, gPos - (avgReadLength >> 1));
		cov_thr = BlockDepthArr[(int)(gPos / BlockSize)] >> 1;
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
		if (Lscore < cov_thr || Lscore < (int)(LCov*INV_TNL_ThrRatio)) Lscore = 0;

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
		if (Rscore < cov_thr || Rscore < (int)(RCov*INV_TNL_ThrRatio)) Rscore = 0;

		//printf("TNL_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			TNLnum++;
			Variant.gPos = gPos;
			Variant.VarType = var_TNL;
			Variant.NS = Lscore > Rscore ? Lscore : Rscore;
			Variant.qscore = CalQualityScore(Variant.NS, cov_thr);
			VariantVec.push_back(Variant);
		}
	}
	if (TNLnum > 0) inplace_merge(VariantVec.begin(), VariantVec.end() - TNLnum, VariantVec.end(), CompByVarPos);
}

void IdentifyInversions()
{
	int64_t gPos;
	Variant_t Variant;
	vector<int64_t> vec;
	DiscordPair_t DiscordPair;
	vector<DiscordPair_t>::iterator Iter1, Iter2;
	uint32_t i, j, n, LCov, RCov, cov_thr, INVnum, num, score, Lscore, Rscore;

	//for (Iter1 = InversionSiteVec.begin(); Iter1 != InversionSiteVec.end(); Iter1++) printf("Pos=%lld Dist=%lld\n", (long long)Iter1->gPos, (long long)Iter1->dist);
	for (num = (int)BreakPointCanVec.size(), INVnum = i = 0; i < num; i++)
	{
		gPos = BreakPointCanVec[i].gPos; LCov = CalRegionCov(gPos - FragmentSize, gPos - (avgReadLength >> 1));
		cov_thr = BlockDepthArr[(int)(gPos / BlockSize)] >> 1;
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
		if (Lscore < cov_thr || Lscore < (int)(LCov*INV_TNL_ThrRatio)) Lscore = 0;

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
		if (Rscore < cov_thr || Rscore < (int)(RCov*INV_TNL_ThrRatio)) Rscore = 0;

		//printf("INV_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			INVnum++;
			Variant.gPos = gPos;
			Variant.NS = Lscore > Rscore ? Lscore : Rscore;
			Variant.VarType = var_INV;
			Variant.qscore = CalQualityScore(Variant.NS, cov_thr);
			VariantVec.push_back(Variant);
		}
	}
	if (INVnum > 0) inplace_merge(VariantVec.begin(), VariantVec.end() - INVnum, VariantVec.end(), CompByVarPos);
}

bool CheckNearbyVariant(int i, int num)
{
	bool bRet = false;
	if (i == 0)
	{
		if (VariantVec[i + 1].gPos - VariantVec[i].gPos <= 10) bRet = true;
	}
	else if (i == num - 1)
	{
		if (VariantVec[i].gPos - VariantVec[i - 1].gPos <= 10) bRet = true;
	}
	else
	{
		if ((VariantVec[i + 1].gPos - VariantVec[i].gPos) <= 10 || (VariantVec[i].gPos - VariantVec[i - 1].gPos) <= 10) bRet = true;
	}
	return bRet;
}

void GenVariantCallingFile()
{
	int64_t gPos;
	float AlleleFreq;
	char failstr[64];
	Coordinate_t coor;
	vector<int> VarNumVec(256);
	int i, n, len, num, thr, cov, dupN;
	map<string, uint16_t>::iterator IndSeqMapIter;

	outFile = fopen(VcfFileName, "w"); ShowMetaInfo();

	sprintf(failstr, "q%d", MinVarConfScore);
	for (num = (int)VariantVec.size(), i = 0; i < num; i++)
	{
		gPos = VariantVec[i].gPos; coor = DetermineCoordinate(gPos);

		if (VariantVec[i].VarType == var_SUB)
		{
			if (VariantVec[i].NS < 10 && CheckNearbyVariant(i, num)) continue;
			VarNumVec[var_SUB]++; 
			fprintf(outFile, "%s	%d	.	%c	%s	%d	%s	DP=%d;AD=%d;AF=%.3f;GT=%s;TYPE=SUBSTITUTE\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[VariantVec[i].gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, (VariantVec[i].qscore >= MinVarConfScore ? "PASS" : failstr), VariantVec[i].DP, VariantVec[i].NS, 1.0*VariantVec[i].NS / VariantVec[i].DP, (VariantVec[i].GenoType ? "0|1": "1|1"));
		}
		else if (VariantVec[i].VarType == var_INS)
		{
			if (VariantVec[i].NS < 5 && CheckNearbyVariant(i, num)) continue;
			VarNumVec[var_INS]++; AlleleFreq = 1.0*VariantVec[i].NS / VariantVec[i].DP;
			fprintf(outFile, "%s	%d	.	%c	%c%s	%d	%s	DP=%d;AD=%d;AF=%.3f;TYPE=INSERT\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos], RefSequence[gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, (VariantVec[i].qscore >= MinVarConfScore ? "PASS" : failstr), VariantVec[i].DP, VariantVec[i].NS, AlleleFreq);
		}
		else if (VariantVec[i].VarType == var_DEL)
		{
			if (VariantVec[i].NS < 5 && CheckNearbyVariant(i, num)) continue;
			VarNumVec[var_DEL]++; AlleleFreq = 1.0*VariantVec[i].NS / VariantVec[i].DP;
			fprintf(outFile, "%s	%d	.	%c%s	%c	%d	%s	DP=%d;AD=%d;AF=%.3f;TYPE=DELETE\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos], VariantVec[i].ALTstr.c_str(), RefSequence[gPos], VariantVec[i].qscore, (VariantVec[i].qscore >= MinVarConfScore ? "PASS" : failstr), VariantVec[i].DP, VariantVec[i].NS, AlleleFreq);
		}
		else if (VariantVec[i].VarType == var_TNL)
		{
			VarNumVec[var_TNL]++;
			fprintf(outFile, "%s	%d	.	%c	<TRANSLOCATION>	30	PASS	DP=%d;NS=%d;TYPE=BND\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos - 1], VariantVec[i].DP, VariantVec[i].NS);
		}
		else if (VariantVec[i].VarType == var_INV)
		{
			VarNumVec[var_INV]++;
			fprintf(outFile, "%s	%d	.	%c	<INVERSION>	30	PASS	DP=%d;NS=%d;TYPE=BND\n", ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, RefSequence[gPos - 1], VariantVec[i].DP, VariantVec[i].NS);
		}
	}
	std::fclose(outFile);
	fprintf(stderr, "\t%d(snp); %d(ins); %d(del); %d(trans); %d(inversion)\n\n", VarNumVec[var_SUB], VarNumVec[var_INS], VarNumVec[var_DEL], VarNumVec[var_TNL] >> 1, VarNumVec[var_INV] >> 1);
}

void *IdentifyVariants(void *arg)
{
	bool bShow = false;
	Variant_t Variant;
	int64_t gPos, end;
	float BaseFreq[4];
	unsigned char ref_base;
	string ins_str, del_str;
	vector<pair<char, int> > vec;
	vector<Variant_t> MyVariantVec;
	map<int64_t, map<string, uint16_t> >::iterator IndMapIter;
	int i, n, cov, cov_thr, freq_thr, ins_thr, del_thr, ins_freq, ins_len, del_freq, del_len, tid = *((int*)arg);

	if (bDebugMode)
	{
		gPos = (tid == 0 ? 0 : (ChromosomeVec[0].len / iThreadNum)*tid);
		end = (tid == iThreadNum - 1 ? ChromosomeVec[0].len : (ChromosomeVec[0].len / iThreadNum)*(tid + 1));
	}
	else
	{
		gPos = (tid == 0 ? 0 : (GenomeSize / iThreadNum)*tid);
		end = (tid == iThreadNum - 1 ? GenomeSize : (GenomeSize / iThreadNum)*(tid + 1));
	}
	map<int64_t, bool>::iterator it; map<int64_t, bool> obs_map;
	if (bDebugMode) obs_map = LoadObservedPos();
	for (; gPos < end; gPos++)
	{
		if (nst_nt4_table[RefSequence[gPos]] != 4)
		{
			cov = GetProfileColumnSize(MappingRecordArr[gPos]);
			if ((cov_thr = BlockDepthArr[(int)(gPos / BlockSize)] >> 1) < 5) cov_thr = 5;
			if ((ins_thr = (int)(cov_thr*0.25)) < MinIndFreq) ins_thr = MinIndFreq;
			if ((del_thr = (int)(cov_thr*0.35)) < MinIndFreq) del_thr = MinIndFreq;

			if (bDebugMode && (it = obs_map.find(gPos)) != obs_map.end()) bShow = true; else bShow = false;
			ins_freq = GetAreaIndFrequency(gPos, InsertSeqMap, ins_str); del_freq = GetAreaIndFrequency(gPos, DeleteSeqMap, del_str);

			if (bShow)
			{
				pthread_mutex_lock(&Lock);
				ShowVariationProfile(gPos - 5, gPos + 5); fflush(stdout);
				ShowIndSeq(gPos - 10, gPos + 10); fflush(stdout);
				printf("%cpos=%lld, cov=%d, cov_thr=%d, freq_thr=%d, ins_freq=%d / %d, del_freq=%d / %d\n\n", (it->second == true ? '*' : '!'), gPos, cov, cov_thr, (int)(cov*FrequencyThr), ins_freq, ins_thr, del_freq, del_thr);
				fflush(stdout);
				pthread_mutex_unlock(&Lock);
			}
			//if (ins_freq >= MinIndFreq)
			if (ins_freq >= ins_thr)
			{
				Variant.gPos = gPos; Variant.VarType = var_INS; Variant.DP = (cov_thr << 1); Variant.NS = ins_freq; 
				if (Variant.DP < Variant.NS) Variant.DP = Variant.NS; Variant.ALTstr = ins_str;
				if ((Variant.qscore = (int)(15.0 * ins_freq / cov_thr)) > 30) Variant.qscore = 30;
				/*if (Variant.qscore >= MinVarConfScore)*/
				MyVariantVec.push_back(Variant);
			}
			//if (del_freq >= MinIndFreq)
			if (del_freq >= del_thr)
			{
				Variant.gPos = gPos - 1; Variant.VarType = var_DEL; Variant.DP = (cov_thr << 1); Variant.NS = del_freq;
				if (Variant.DP < Variant.NS) Variant.DP = Variant.NS; Variant.ALTstr = del_str;
				if ((Variant.qscore = (int)(15.0 * del_freq / cov_thr)) > 30) Variant.qscore = 30;
				/*if (Variant.qscore >= MinVarConfScore)*/
				MyVariantVec.push_back(Variant);
			}
			//SUB
			if (cov >= cov_thr && cov > ins_freq && cov > del_freq)
			{
				vec.clear(); ref_base = nst_nt4_table[RefSequence[gPos]]; freq_thr = cov*(bSomatic ? 0.02 : FrequencyThr);
				if (freq_thr < MinAlleleFreq) freq_thr = MinAlleleFreq; if (freq_thr > cov_thr) freq_thr = cov_thr;
				if (ref_base != 0 && MappingRecordArr[gPos].A >= freq_thr) vec.push_back(make_pair('A', MappingRecordArr[gPos].A));
				if (ref_base != 1 && MappingRecordArr[gPos].C >= freq_thr) vec.push_back(make_pair('C', MappingRecordArr[gPos].C));
				if (ref_base != 2 && MappingRecordArr[gPos].G >= freq_thr) vec.push_back(make_pair('G', MappingRecordArr[gPos].G));
				if (ref_base != 3 && MappingRecordArr[gPos].T >= freq_thr) vec.push_back(make_pair('T', MappingRecordArr[gPos].T));
				if (vec.size() == 1)
				{
					Variant.gPos = gPos; Variant.VarType = var_SUB; Variant.DP = cov; Variant.NS = vec[0].second;
					if (vec[0].second >= (cov - freq_thr)) Variant.GenoType = 0;
					else Variant.GenoType = 1;

					Variant.ALTstr = vec[0].first;
					Variant.qscore = bSomatic ? (int)(1.0* Variant.NS / (cov*0.01)) : (int)(10 * (1.0* Variant.NS / (cov*FrequencyThr)));
					if (Variant.qscore > 30) Variant.qscore = 30;
					/*if (Variant.qscore >= MinVarConfScore) */
					MyVariantVec.push_back(Variant);
				}
				else if (vec.size() == 2 && CheckDiploidFrequency(cov, vec))
				{
					Variant.gPos = gPos; Variant.VarType = var_SUB; Variant.DP = cov; Variant.NS = vec[0].second + vec[1].second;
					Variant.GenoType = 1;

					Variant.ALTstr.resize(3); Variant.ALTstr[0] = vec[0].first; Variant.ALTstr[1] = ',';  Variant.ALTstr[1] = vec[1].first;
					Variant.qscore = bSomatic ? (int)(1.0* Variant.NS / (cov*0.01)) : (int)(10 * (1.0* Variant.NS / (cov*FrequencyThr)));
					if (Variant.qscore > 30) Variant.qscore = 30;
					/*if (Variant.qscore >= MinVarConfScore) */
					MyVariantVec.push_back(Variant);
				}
			}
		}
	}
	if ((n = (int)MyVariantVec.size()) > 0)
	{
		pthread_mutex_lock(&Lock);
		copy(MyVariantVec.begin(), MyVariantVec.end(), back_inserter(VariantVec)); inplace_merge(VariantVec.begin(), VariantVec.end() - n, VariantVec.end(), CompByVarPos);
		pthread_mutex_unlock(&Lock);
	}
	return (void*)(1);
}

void VariantCalling()
{
	int i, *ThrIDarr;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	ThrIDarr = new int[iThreadNum];  for (i = 0; i < iThreadNum; i++) ThrIDarr[i] = i;

	//iThreadNum = 1;
	BlockNum = (int)(GenomeSize / BlockSize); if (((int64_t)BlockNum * BlockSize) < GenomeSize) BlockNum += 1; 
	BlockDepthArr = new int[BlockNum]();
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, CalBlockReadDepth, &ThrIDarr[i]);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	fprintf(stderr, "Identify all variants...\n"); fflush(stderr);
	//iThreadNum = 1;
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyVariants, &ThrIDarr[i]);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	//IdentifyBreakPointCandidates();
	//if (BreakPointCanVec.size() > 0 && InversionSiteVec.size() > 0) IdentifyInversions();
	//if (BreakPointCanVec.size() > 0 && TranslocationSiteVec.size() > 0) IdentifyTranslocations();
	fprintf(stderr, "\tWrite all the predicted sample variations to file [%s]...\n", VcfFileName); GenVariantCallingFile();

	fprintf(stderr, "SV calling has be done in %lld seconds.\n", (long long)(time(NULL) - StartProcessTime));
	delete[] ThrIDarr; delete[] ThreadArr; delete[] BlockDepthArr;
}

