#include "structure.h"

#define MaxQscore 30
#define BlockSize 100
#define BreakPointFreqThr 3
#define INV_TNL_ThrRatio 0.5
#define Genotype_Ratio	0.50
#define var_SUB 0 // substitution
#define var_INS 1 // insertion
#define var_DEL 2 // deletion
#define var_INV 3 // inversion
#define var_TNL 4 // translocation
#define var_CNV 5 // copy number variation
#define var_UMR 6 // unmapped region
#define var_NOR 10 // normal region (for gVCF)
#define var_MON	11 // monomorphic
#define var_NIL 255

const char* GenotypeLabel[] = {"*", "0", "1", "0|0", "0|1", "1|1", "1|2"};

typedef struct
{
	int64_t gPos;
	uint16_t left_score;
	uint16_t rigt_score;
} BreakPoint_t;

FILE *outFile;
int* BlockDepthArr;
vector<int> VarNumVec(256);
int BlockNum, iTotalVarNum;
vector<Variant_t> VariantVec;
vector<BreakPoint_t> BreakPointCanVec;
extern map<int64_t, uint16_t> BreakPointMap;

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

bool CompByVarPos(const Variant_t& p1, const Variant_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.VarType < p2.VarType;
	else return p1.gPos < p2.gPos;
}

int GetPointIndFreq(map<string, uint16_t>& IndMap)
{
	int n = 0;
	for (map<string, uint16_t>::iterator iter = IndMap.begin(); iter != IndMap.end(); iter++) n += iter->second;
	return n;
}

int GetAreaIndFrequency(int64_t gPos, map<int64_t, map<string, uint16_t> >& IndMap, string& ind_str)
{
	int64_t max_pos = 0;
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
	if (sum >= (int)(cov*Genotype_Ratio)) return true;
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
	fprintf(outFile, "##fileformat=VCFv4.2\n");
	fprintf(outFile, "##reference=%s\n", (RefFileName != NULL ? RefFileName: IndexFileName));
	fprintf(outFile, "##source=MapCaller %s\n", VersionStr);
	fprintf(outFile, "##command_line=\"%s\"\n", CmdLine.c_str());
	fprintf(outFile, "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Number of reads with start coordinate at this position.\">\n");
	fprintf(outFile, "##INFO=<ID=NTFREQ,Number=4,Type=Integer,Description=\"base depth\">\n");
	fprintf(outFile, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Last position(inclusive) of the reported block\">\n");
	fprintf(outFile, "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snv, ins, del, or BP(breakpoint).\">\n");
	if (bGVCF) fprintf(outFile, "##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth in gVCF output block.\">\n");
	fprintf(outFile, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n");
	fprintf(outFile, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n");
	fprintf(outFile, "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles\">\n");
	fprintf(outFile, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(outFile, "##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=\"Count of reads in F1R2 pair orientation supporting each allele\">\n");
	fprintf(outFile, "##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=\"Count of reads in F2R1 pair orientation supporting each allele\">\n");
	fprintf(outFile, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf(outFile, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
	//fprintf(outFile, "##FILTER=<ID=LowDepth,Description=\"Read depth < %d\">\n", MinReadDepth);
	fprintf(outFile, "##FILTER=<ID=REF,Description=\"Genotyping model thinks this site is reference.\">\n");
	fprintf(outFile, "##FILTER=<ID=BreakPoint,Description=\"It is predicted as a breakpoint\">\n");
	fprintf(outFile, "##FILTER=<ID=DUP,Description=\"Duplicated regions(>=%dbp).\">\n", MinCNVsize);
	fprintf(outFile, "##FILTER=<ID=Gaps,Description=\"Region without any read alignment(>=%dbp).\">\n", MinUnmappedSize);
	fprintf(outFile, "##FILTER=<ID=q10,Description=\"Confidence score below 10\">\n");
	if (bFilter) fprintf(outFile, "##FILTER=<ID=bad_haplotype,Description=\"Variants with variable frequencies on same haplotype\">\n");
	if (bFilter) fprintf(outFile, "##FILTER=<ID=str_contraction,Description=\"Variant appears in repetitive region\">\n");
	for (int i = 0; i < iChromsomeNum; i++) fprintf(outFile, "##contig=<ID=%s,length=%d>\n", ChromosomeVec[i].name, ChromosomeVec[i].len);
	fprintf(outFile, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	unknown\n");
}

void IdentifyBreakPointCandidates()
{
	BreakPoint_t bp;
	uint32_t total_freq;
	pair<int64_t, uint16_t> p;

	BreakPointMap.insert(make_pair(TwoGenomeSize, 0)); total_freq = 0; p = make_pair(0, 0);
	for (map<int64_t, uint16_t>::iterator iter = BreakPointMap.begin(); iter != BreakPointMap.end(); iter++)
	{
		if (iter->first - p.first > avgReadLength) // break
		{
			if (total_freq >= BreakPointFreqThr)
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

	if (begPos < 0) begPos = 0; if (endPos > GenomeSize) endPos = GenomeSize - 1;
	if (endPos < begPos) return 0;
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
		gPos = BreakPointCanVec[i].gPos;
		
		LCov = CalRegionCov(gPos - FragmentSize, gPos - (avgReadLength >> 1));
		cov_thr = BlockDepthArr[(int)(gPos / BlockSize)] >> 1;
		DiscordPair.gPos = gPos - FragmentSize; Iter1 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos - (avgReadLength >> 1); Iter2 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		if (Iter1 == TranslocationSiteVec.end() || Iter2 == TranslocationSiteVec.end()) continue;
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000)); 
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize); n = (int)vec.size();
		for (Lscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Lcan_%d: %lld (score=%d)\n", j + 1, vec[j], score);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Lscore) Lscore = score;
				score = 1;
			}
			else score++;
		}
		if (Lscore < cov_thr || Lscore < (int)(LCov*INV_TNL_ThrRatio)) continue;

		RCov = CalRegionCov(gPos, gPos + FragmentSize);
		DiscordPair.gPos = gPos; Iter1 = upper_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos + FragmentSize; Iter2 = lower_bound(TranslocationSiteVec.begin(), TranslocationSiteVec.end(), DiscordPair, CompByDiscordPos);
		if (Iter1 == TranslocationSiteVec.end() || Iter2 == TranslocationSiteVec.end()) continue;
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000));
		sort(vec.begin(), vec.end()); vec.push_back(TwoGenomeSize); n = (int)vec.size();
		for (Rscore = 0, score = j = 1; j < n; j++)
		{
			//printf("Rcan_%d: %lld (score=%d)\n", j + 1, vec[j], score);
			if (vec[j] - vec[j - 1] > 1)
			{
				if (score > Rscore) Rscore = score;
				score = 1;
			}
			else score++;
		}
		if (Rscore < cov_thr || Rscore < (int)(RCov*INV_TNL_ThrRatio)) continue;

		//printf("TNL_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			TNLnum++;
			Variant.gPos = gPos;
			Variant.VarType = var_TNL;
			Variant.DP = GetProfileColumnSize(MappingRecordArr[gPos]);
			Variant.AD_alt = Lscore > Rscore ? Lscore : Rscore;
			Variant.qscore = CalQualityScore(Variant.AD_alt, cov_thr);
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
		if (Iter1 == InversionSiteVec.end() || Iter2 == InversionSiteVec.end()) continue;
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000));
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
		if (Lscore < cov_thr || Lscore < (int)(LCov*INV_TNL_ThrRatio)) continue;

		RCov = CalRegionCov(gPos, gPos + FragmentSize);
		DiscordPair.gPos = gPos; Iter1 = upper_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);
		DiscordPair.gPos = gPos + FragmentSize; Iter2 = lower_bound(InversionSiteVec.begin(), InversionSiteVec.end(), DiscordPair, CompByDiscordPos);
		if (Iter1 == InversionSiteVec.end() || Iter2 == InversionSiteVec.end()) continue;
		vec.clear(); for (; Iter1 != Iter2; Iter1++) vec.push_back((Iter1->dist / 1000));
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
		if (Rscore < cov_thr || Rscore < (int)(RCov*INV_TNL_ThrRatio)) continue;

		//printf("INV_can =%lld (Cov=%d vs %d): Lscore=%d Rscore=%d\n", (long long)gPos, LCov, RCov, Lscore, Rscore);
		if (Lscore > 0 && Rscore > 0)
		{
			INVnum++;
			Variant.gPos = gPos;
			Variant.DP = GetProfileColumnSize(MappingRecordArr[gPos]);
			Variant.AD_alt = Lscore > Rscore ? Lscore : Rscore;
			Variant.VarType = var_INV;
			Variant.qscore = CalQualityScore(Variant.AD_alt, cov_thr);
			VariantVec.push_back(Variant);
		}
	}
	if (INVnum > 0) inplace_merge(VariantVec.begin(), VariantVec.end() - INVnum, VariantVec.end(), CompByVarPos);
}

bool CheckNearbyVariant(int i, int dist)
{
	bool bRet = false;
	if (i == 0)
	{
		if (VariantVec[i + 1].gPos - VariantVec[i].gPos <= dist) bRet = true;
	}
	else if (i == iTotalVarNum - 1)
	{
		if (VariantVec[i].gPos - VariantVec[i - 1].gPos <= dist) bRet = true;
	}
	else
	{
		if ((VariantVec[i + 1].gPos - VariantVec[i].gPos) <= dist || (VariantVec[i].gPos - VariantVec[i - 1].gPos) <= dist) bRet = true;
	}
	return bRet;
}

bool CheckBadHaplotype(int i, int dist)
{
	int j, diff;
	bool bRet = false;

	for (j = i + 1; j < iTotalVarNum; j++)
	{
		if (VariantVec[j].gPos - VariantVec[i].gPos > dist) break;
		if (VariantVec[j].VarType == 0)
		{
			diff = abs(VariantVec[i].AD_alt - VariantVec[j].AD_alt);
			if (diff > 5 && (VariantVec[i].AD_alt > VariantVec[j].AD_alt ? VariantVec[i].AD_alt >> 2 : VariantVec[j].AD_alt >> 2)) bRet = true;
			break;
		}
	}
	for (j = i - 1; j >= 0; j--)
	{
		if (VariantVec[i].gPos - VariantVec[j].gPos > dist) break;
		if (VariantVec[j].VarType == 0)
		{
			diff = abs(VariantVec[i].AD_alt - VariantVec[j].AD_alt);
			if (diff > 10 && (VariantVec[i].AD_alt > VariantVec[j].AD_alt ? (int)(VariantVec[i].AD_alt * 0.33) : (int)(VariantVec[j].AD_alt * 0.33))) bRet = true;
			break;
		}
	}
	return bRet;
}

void ShowNeighboringProfile(int64_t gPos, Coordinate_t coor)
{
	int64_t i, p;

	for (i = -5; i <= 5; i++)
	{
		if (i == 0) printf("*");
		p = gPos + i;
		coor = DetermineCoordinate(p);
		printf("%s-%lld\t\t%d\t%d\t%d\t%d\tR=%d\tdepth=%d\n", ChromosomeVec[coor.ChromosomeIdx].name, (long long)coor.gPos, (int)MappingRecordArr[p].A, (int)MappingRecordArr[p].C, (int)MappingRecordArr[p].G, (int)MappingRecordArr[p].T, (int)MappingRecordArr[p].multi_hit, GetProfileColumnSize(MappingRecordArr[p]));
	}
	printf("\n\n");
}

string DetermineFileter(int VarIdx)
{
	string filter_str="";

	//if (VariantVec[VarIdx].DP < MinReadDepth) filter_str += "LowDepth;";
	if (VariantVec[VarIdx].qscore < 10) filter_str += "q10;";
	else if (VariantVec[VarIdx].VarType == var_SUB && VariantVec[VarIdx].AD_alt < 10 && CheckNearbyVariant(VarIdx, 10)) filter_str += "q10;";
	else if ((VariantVec[VarIdx].VarType == var_INS || VariantVec[VarIdx].VarType == var_DEL) && VariantVec[VarIdx].AD_alt < 5 && CheckNearbyVariant(VarIdx, 10)) filter_str += "q10;";

	if (bFilter)
	{
		if ((int)MappingRecordArr[VariantVec[VarIdx].gPos].multi_hit > (int)(GetProfileColumnSize(MappingRecordArr[VariantVec[VarIdx].gPos])*0.05)) filter_str += "str_contraction;";
		if (CheckBadHaplotype(VarIdx, 100)) filter_str += "bad_haplotype;";
	}
	if (filter_str == "") filter_str = "PASS";
	else filter_str.resize(filter_str.length() - 1);

	return filter_str;
}

void GenVariantCallingFile()
{
	int i;
	float AlleleFreq;
	Coordinate_t coor;
	string filter_str;
	int64_t gPos, gPosEnd;
	map<string, uint16_t>::iterator IndSeqMapIter;

	outFile = fopen(VcfFileName, "w"); ShowMetaInfo();

	//sort(VariantVec.begin(), VariantVec.end(), CompByVarPos);
	for (i = 0; i < iTotalVarNum; i++)
	{
		gPos = VariantVec[i].gPos; coor = DetermineCoordinate(gPos);

		if(VariantVec[i].VarType < 3) filter_str = DetermineFileter(i);
		else filter_str = ".";

		if (VariantVec[i].VarType == var_SUB)
		{
			VarNumVec[var_SUB]++; AlleleFreq = 1.0*VariantVec[i].AD_alt / VariantVec[i].DP;
			//fprintf(outFile, "%s	%d	.	%c	%s	%d	%s	DP=%d;AD=%d;RC=%d;AF=%.3f;NTFREQ=%d,%d,%d,%d;GT=%s;TYPE=snv\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[VariantVec[i].gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, filter_str.c_str(), VariantVec[i].DP, VariantVec[i].AD, (int)MappingRecordArr[gPos].readCount, 1.0*VariantVec[i].AD / VariantVec[i].DP, (int)MappingRecordArr[gPos].A, (int)MappingRecordArr[gPos].C, (int)MappingRecordArr[gPos].G, (int)MappingRecordArr[gPos].T, GenotypeLabel[VariantVec[i].GenoType]);
			fprintf(outFile, "%s	%d	.	%c	%s	%d	%s	RC=%d;NTFREQ=%d,%d,%d,%d;TYPE=snv	GT:GQ:DP:AD:AF:F1R2:F2R1	%s:%d:%d:%d,%d:%.2f:%d,%d:%d,%d\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[VariantVec[i].gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, filter_str.c_str(), (int)MappingRecordArr[gPos].readCount, (int)MappingRecordArr[gPos].A, (int)MappingRecordArr[gPos].C, (int)MappingRecordArr[gPos].G, (int)MappingRecordArr[gPos].T, GenotypeLabel[VariantVec[i].GenoType], VariantVec[i].qscore, VariantVec[i].DP, VariantVec[i].AD_ref, VariantVec[i].AD_alt, AlleleFreq, MappingRecordArr[gPos].F1, MappingRecordArr[gPos].R2, MappingRecordArr[gPos].F2, MappingRecordArr[gPos].R1);
		}
		else if (VariantVec[i].VarType == var_INS)
		{
			VarNumVec[var_INS]++; AlleleFreq = 1.0*VariantVec[i].AD_alt / VariantVec[i].DP;
			//fprintf(outFile, "%s	%d	.	%c	%c%s	%d	%s	DP=%d;AD=%d;RC=%d;AF=%.3f;GT=%s;TYPE=ins\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], RefSequence[gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, filter_str.c_str(), VariantVec[i].DP, VariantVec[i].AD, (int)MappingRecordArr[gPos].readCount, AlleleFreq, GenotypeLabel[VariantVec[i].GenoType]);
			fprintf(outFile, "%s	%d	.	%c	%c%s	%d	%s	RC=%d;TYPE=ins	GT:GQ:DP:AD:AF:F1R2:F2R1	%s:%d:%d:%d,%d:%.2f:%d,%d:%d,%d\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], RefSequence[gPos], VariantVec[i].ALTstr.c_str(), VariantVec[i].qscore, filter_str.c_str(), (int)MappingRecordArr[gPos].readCount, GenotypeLabel[VariantVec[i].GenoType], VariantVec[i].qscore, VariantVec[i].DP, VariantVec[i].AD_ref, VariantVec[i].AD_alt, AlleleFreq, MappingRecordArr[gPos].F1, MappingRecordArr[gPos].R2, MappingRecordArr[gPos].F2, MappingRecordArr[gPos].R1);
		}
		else if (VariantVec[i].VarType == var_DEL)
		{
			VarNumVec[var_DEL]++; AlleleFreq = 1.0*VariantVec[i].AD_alt / VariantVec[i].DP;
			//fprintf(outFile, "%s	%d	.	%c%s	%c	%d	%s	DP=%d;AD=%d;RC=%d;AF=%.3f;GT=%s;TYPE=del\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], VariantVec[i].ALTstr.c_str(), RefSequence[gPos], VariantVec[i].qscore, filter_str.c_str(), VariantVec[i].DP, VariantVec[i].AD, (int)MappingRecordArr[gPos].readCount, AlleleFreq, GenotypeLabel[VariantVec[i].GenoType]);
			fprintf(outFile, "%s	%d	.	%c%s	%c	%d	%s	RC=%d;TYPE=del	GT:GQ:DP:AD:AF:F1R2:F2R1	%s:%d:%d:%d,%d:%.2f:%d,%d:%d,%d\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], VariantVec[i].ALTstr.c_str(), RefSequence[gPos], VariantVec[i].qscore, filter_str.c_str(), (int)MappingRecordArr[gPos].readCount, GenotypeLabel[VariantVec[i].GenoType], VariantVec[i].qscore, VariantVec[i].DP, VariantVec[i].AD_ref, VariantVec[i].AD_alt, AlleleFreq, MappingRecordArr[gPos].F1, MappingRecordArr[gPos].R2, MappingRecordArr[gPos].F2, MappingRecordArr[gPos].R1);
		}
		else if (VariantVec[i].VarType == var_TNL)
		{
			VarNumVec[var_TNL]++;
			fprintf(outFile, "%s	%d	.	%c	<TNL>	30	BreakPoint	TYPE=BP	GT:GQ:DP:AD	*:*:*:*\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos]);
		}
		else if (VariantVec[i].VarType == var_INV)
		{
			VarNumVec[var_INV]++;
			fprintf(outFile, "%s	%d	.	%c	<INV>	30	BreakPoint	TYPE=BP	GT:GQ:DP:AD	*:*:*:*\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos]);
		}
		else if (VariantVec[i].VarType == var_CNV)
		{
			//gPosEnd = ChromosomeVec[coor.ChromosomeIdx].FowardLocation + ChromosomeVec[coor.ChromosomeIdx].len - 1;
			if (VariantVec[i].DP >= MinCNVsize) fprintf(outFile, "%s	%d	.	%c	<*>	0	DUP	END=%d	GT:GQ:DP:AD	*:*:*:*\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], (int)(coor.gPos + VariantVec[i].DP - 1));
		}
		else if (VariantVec[i].VarType == var_UMR)
		{
			if(VariantVec[i].DP >= MinUnmappedSize) fprintf(outFile, "%s	%d	.	%c	<*>	0	Gaps	END=%d	GT:GQ:DP:AD	*:*:*:*\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], (int)(coor.gPos + VariantVec[i].DP - 1));
		}
		else if (VariantVec[i].VarType == var_NOR)
		{
			gPosEnd = ChromosomeVec[coor.ChromosomeIdx].FowardLocation + ChromosomeVec[coor.ChromosomeIdx].len - 1;
			if (i + 1 < iTotalVarNum && VariantVec[i + 1].gPos < gPosEnd) gPosEnd = VariantVec[i + 1].gPos - 1;
			fprintf(outFile, "%s	%d	.	%c	<*>	0	REF	END=%d;DP=%d;MIN_DP=%d	GT:GQ:DP:AD	*:*:*:*\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], (int)DetermineCoordinate(gPosEnd).gPos, VariantVec[i].DP, VariantVec[i].AD_alt);
		}
		else if (VariantVec[i].VarType == var_MON)
		{
			fprintf(outFile, "%s	%d	.	%c	.	0	REF	DP=%d;RC=%d;NTFREQ=%d,%d,%d,%d	GT:F1R2:F2R1	%s:%d,%d:%d,%d\n", ChromosomeVec[coor.ChromosomeIdx].name, (int)coor.gPos, RefSequence[gPos], VariantVec[i].DP, (int)MappingRecordArr[gPos].readCount, (int)MappingRecordArr[gPos].A, (int)MappingRecordArr[gPos].C, (int)MappingRecordArr[gPos].G, (int)MappingRecordArr[gPos].T, GenotypeLabel[VariantVec[i].GenoType], MappingRecordArr[gPos].F1, MappingRecordArr[gPos].R2, MappingRecordArr[gPos].F2, MappingRecordArr[gPos].R1);
		}
	}
	std::fclose(outFile);
}

bool CheckNeighboringCoverage(int64_t gPos, int cov)
{
	int p, c, diff = 0;

	for (p = -5; p <= 5; p++)
	{
		if (p == 0) continue;
		c = GetProfileColumnSize(MappingRecordArr[gPos + p]);
		diff += abs(cov - c);
	}
	diff /= 10;
	if (diff >= 10 || (diff > 1 && diff >= (int)cov*0.1)) return false;
	else return true;
}

uint16_t GetRefCount(unsigned char ref_base, int64_t gPos)
{
	switch (ref_base)
	{
	case 0: return (uint16_t)MappingRecordArr[gPos].A;
	case 1: return (uint16_t)MappingRecordArr[gPos].C;
	case 2: return (uint16_t)MappingRecordArr[gPos].G;
	case 3: return (uint16_t)MappingRecordArr[gPos].T;
	default: return 0;
	}
}

uint8_t DetermineGenotype(int cov, int alt_read_count, int alt_num)
{
	uint8_t genotype = 0;
	if (iPloidy == 1)
	{
		if (alt_read_count < (int)(cov*Genotype_Ratio)) genotype = 1; //GT=0
		else genotype = 2; // GT=1
	}
	else if (iPloidy == 2)
	{
		if (alt_num == 0) genotype = 3; // GT=0/0
		else if (alt_num == 1)
		{
			if (alt_read_count < (int)(cov*Genotype_Ratio)) genotype = 4; // GT=0/1
			else genotype = 5; // GT=1/1
		}
		else if (alt_num == 2) genotype = 6; // GT=1/2
	}
	return genotype;
}

void *IdentifyVariants(void *arg)
{
	bool bNormal;
	Variant_t Variant;
	int64_t gPos, end;
	unsigned char ref_base;
	string ins_str, del_str;
	vector<pair<char, int> > vec;
	vector<Variant_t> MyVariantVec;
	map<int64_t, map<string, uint16_t> >::iterator IndMapIter;
	int n, gap, dup, cov, cov_thr, freq_thr, ins_thr, del_thr, ins_freq, del_freq, tid = *((int*)arg);

	gPos = (tid == 0 ? 0 : (GenomeSize / iThreadNum)*tid);
	end = (tid == iThreadNum - 1 ? GenomeSize : (GenomeSize / iThreadNum)*(tid + 1));

	gap = dup = 0;
	for (; gPos < end; gPos++)
	{
		cov = GetProfileColumnSize(MappingRecordArr[gPos]);
		bNormal = true; ref_base = nst_nt4_table[(unsigned short)RefSequence[gPos]];
		//if (bSomatic && (MappingRecordArr[gPos].multi_hit > (int)(cov*0.05))) continue;
		if ((cov_thr = BlockDepthArr[(int)(gPos / BlockSize)] >> 1) < MinAlleleDepth) cov_thr = MinAlleleDepth;
		if (bSomatic && cov_thr > MinAlleleDepth) cov_thr = MinAlleleDepth;

		if ((ins_thr = (int)(cov_thr*0.25)) < MinAlleleDepth) ins_thr = MinAlleleDepth;
		if ((del_thr = (int)(cov_thr*0.35)) < MinAlleleDepth) del_thr = MinAlleleDepth;
		ins_freq = GetAreaIndFrequency(gPos, InsertSeqMap, ins_str); del_freq = GetAreaIndFrequency(gPos, DeleteSeqMap, del_str);

		if (ins_freq >= ins_thr)
		{
			Variant.gPos = gPos; Variant.VarType = var_INS; Variant.DP = BlockDepthArr[(int)(gPos / BlockSize)]; Variant.AD_alt = ins_freq;
			if (Variant.DP < Variant.AD_alt) Variant.DP = Variant.AD_alt; Variant.ALTstr = ins_str;
			Variant.AD_ref = Variant.DP - Variant.AD_alt; Variant.GenoType = DetermineGenotype(Variant.DP, Variant.AD_alt, 1);
			Variant.qscore = (int)(100.0*Variant.AD_alt / cov);
			//printf("ins@%lld: dp=%d ad=%d, gt=%d\n", gPos, Variant.DP, Variant.AD_alt, Variant.GenoType);
			//if ((Variant.qscore = (int)(1.0 * MaxQscore * ins_freq / Variant.DP)) > MaxQscore) Variant.qscore = MaxQscore;
			bNormal = false; MyVariantVec.push_back(Variant);
		}
		if (del_freq >= del_thr)
		{
			Variant.gPos = gPos; Variant.VarType = var_DEL; Variant.DP = BlockDepthArr[(int)(gPos / BlockSize)]; Variant.AD_alt = del_freq;
			if (Variant.DP < Variant.AD_alt) Variant.DP = Variant.AD_alt; Variant.ALTstr = del_str; Variant.ALTstr.resize((n = (int)del_str.length())); //strncpy((char*)Variant.ALTstr.c_str(), RefSequence + gPos + 1, n);
			Variant.AD_ref = Variant.DP - Variant.AD_alt; Variant.GenoType = DetermineGenotype(Variant.DP, Variant.AD_alt, 1);
			Variant.qscore = (int)(100.0*Variant.AD_alt / cov);
			//printf("del@%lld: dp=%d ad=%d, gt=%d\n", gPos, Variant.DP, Variant.AD_alt, Variant.GenoType);
			//if ((Variant.qscore = (int)(1.0 * MaxQscore * del_freq / Variant.DP)) > MaxQscore) Variant.qscore = MaxQscore;
			bNormal = false; MyVariantVec.push_back(Variant);
		}
		//SUB
		if (cov >= cov_thr)
		{
			vec.clear(); freq_thr = (int)ceil(cov*(bSomatic ? 0.01 : FrequencyThr));

			if (freq_thr < MinAlleleDepth) freq_thr = MinAlleleDepth;

			if (ref_base != 0 && (int)MappingRecordArr[gPos].A >= freq_thr) vec.push_back(make_pair('A', (int)MappingRecordArr[gPos].A));
			if (ref_base != 1 && (int)MappingRecordArr[gPos].C >= freq_thr) vec.push_back(make_pair('C', (int)MappingRecordArr[gPos].C));
			if (ref_base != 2 && (int)MappingRecordArr[gPos].G >= freq_thr) vec.push_back(make_pair('G', (int)MappingRecordArr[gPos].G));
			if (ref_base != 3 && (int)MappingRecordArr[gPos].T >= freq_thr) vec.push_back(make_pair('T', (int)MappingRecordArr[gPos].T));
			
			Variant.AD_ref = GetRefCount(ref_base, gPos);
			if (vec.size() == 1)
			{
				Variant.gPos = gPos; Variant.VarType = var_SUB; Variant.DP = (uint16_t)cov; Variant.AD_alt = (uint16_t)vec[0].second;
				if ((Variant.GenoType = DetermineGenotype(cov, Variant.AD_alt, 1)) != 0)
				{
					Variant.ALTstr = vec[0].first; //Variant.qscore = (int)(30.0*Variant.AD_alt / cov);
					Variant.qscore = bSomatic ? (int)(35.0 * Variant.AD_alt / (cov*0.05)) : (int)(35.0*Variant.AD_alt / cov);
					bNormal = false; MyVariantVec.push_back(Variant);
				}
			}
			else if (vec.size() ==  2 && CheckDiploidFrequency(cov, vec))
			{
				Variant.gPos = gPos; Variant.VarType = var_SUB; Variant.GenoType = 1; Variant.DP = (uint16_t)cov; Variant.AD_alt = (uint16_t)(vec[0].second + vec[1].second);
				if ((Variant.GenoType = DetermineGenotype(cov, Variant.AD_alt, 2)) != 0)
				{
					Variant.ALTstr.resize(3); Variant.ALTstr[0] = vec[0].first; Variant.ALTstr[1] = ',';  Variant.ALTstr[2] = vec[1].first;
					Variant.qscore = bSomatic ? (int)(35.0 * Variant.AD_alt / (cov*0.05)) : (int)(35.0*Variant.AD_alt / cov);
					bNormal = false; MyVariantVec.push_back(Variant);
				}
			}
		}
		if (cov == 0 && MappingRecordArr[gPos].multi_hit == 0) bNormal=false, gap++;
		else if(gap > 0)
		{
			if (gap > MinUnmappedSize)
			{
				Variant.VarType = var_UMR; 	Variant.gPos = gPos - gap; Variant.DP = gap;
				MyVariantVec.push_back(Variant);
			}
			gap = 0;
		}
		if (cov == 0 && MappingRecordArr[gPos].multi_hit > 0) bNormal = false, dup++;
		else if (dup > 0)
		{
			if (dup > MinCNVsize)
			{
				Variant.VarType = var_CNV; Variant.gPos = gPos - dup; Variant.DP = dup;
				MyVariantVec.push_back(Variant);
			}
			dup = 0; 
		}
		if (bGVCF && bNormal && cov > 0)
		{
			if (MyVariantVec.size() == 0 || MyVariantVec.rbegin()->VarType != var_NOR)
			{
				Variant.qscore = 0; Variant.gPos = gPos; Variant.VarType = var_NOR; Variant.DP = Variant.AD_alt = cov; Variant.ALTstr.clear();
				MyVariantVec.push_back(Variant);
			}
			else
			{
				if (MyVariantVec.rbegin()->AD_alt > cov) MyVariantVec.rbegin()->AD_alt = (uint16_t)cov;
			}
		}
		if (bMonomorphic && bNormal && cov > 0)
		{
			Variant.qscore = 0; Variant.gPos = gPos; Variant.VarType = var_MON; Variant.DP = cov; Variant.GenoType = DetermineGenotype(cov, 0, 0); Variant.ALTstr.clear();
			Variant.AD_ref = GetRefCount(ref_base, gPos);
			MyVariantVec.push_back(Variant);
		}
		//printf("%lld: cov=%d, cnv=%d, umr=%d\n", gPos, cov, dup, gap); ShowProfileColumn(gPos);
	}
	if ((n = (int)MyVariantVec.size()) > 0)
	{
		sort(MyVariantVec.begin(), MyVariantVec.end(), CompByVarPos);
		pthread_mutex_lock(&VarLock);
		copy(MyVariantVec.begin(), MyVariantVec.end(), back_inserter(VariantVec)); inplace_merge(VariantVec.begin(), VariantVec.end() - n, VariantVec.end(), CompByVarPos);
		pthread_mutex_unlock(&VarLock);
	}
	return (void*)(1);
}

void RemoveConsecutiveGenomicVariant()
{
	vector<Variant_t>::iterator iter, next_iter;

	for (iter = VariantVec.begin(), next_iter = iter + 1; next_iter != VariantVec.end(); iter++, next_iter++)
	{
		if (iter->VarType == var_NOR && next_iter->VarType == var_NOR)
		{
			iter = VariantVec.erase(next_iter);
			next_iter = iter + 1;
		}
	}
}

void VariantCalling()
{
	FILE *log;
	int i, *ThrIDarr;
	time_t t = time(NULL);
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	ThrIDarr = new int[iThreadNum];  for (i = 0; i < iThreadNum; i++) ThrIDarr[i] = i;

	log = fopen(LogFileName, "a");

	//if (ObserveBegPos != -1) printf("Profile[%lld-%lld]\n", (long long)ObserveBegPos, (long long)ObserveEndPos), ShowVariationProfile(ObserveBegPos, ObserveEndPos);

	BlockNum = (int)(GenomeSize / BlockSize); if (((int64_t)BlockNum * BlockSize) < GenomeSize) BlockNum += 1; 
	BlockDepthArr = new int[BlockNum]();
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, CalBlockReadDepth, &ThrIDarr[i]);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

	fprintf(log, "Identify all variants (min_alt_allele_depth=%d)...\n", MinAlleleDepth); fflush(stderr);
	fprintf(stderr, "Identify all variants (min_alt_allele_depth=%d)...\n", MinAlleleDepth); fflush(stderr);

	iThreadNum = 1;
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, IdentifyVariants, &ThrIDarr[i]);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	if (bGVCF) RemoveConsecutiveGenomicVariant();

	// Identify structural variants
	IdentifyBreakPointCandidates();
	if (BreakPointCanVec.size() > 0 && InversionSiteVec.size() > 0) IdentifyInversions();
	if (BreakPointCanVec.size() > 0 && TranslocationSiteVec.size() > 0) IdentifyTranslocations();

	iTotalVarNum = (int)VariantVec.size();
	fprintf(log, "\tWrite all the predicted sample variations to file [%s]...\n", VcfFileName); 
	fprintf(stderr, "\tWrite all the predicted sample variations to file [%s]...\n", VcfFileName);
	GenVariantCallingFile();
	fprintf(log, "\t%d(snp); %d(ins); %d(del); %d(trans); %d(inversion)\n", VarNumVec[var_SUB], VarNumVec[var_INS], VarNumVec[var_DEL], VarNumVec[var_TNL] >> 1, VarNumVec[var_INV] >> 1);
	fprintf(stderr, "\t%d(snp); %d(ins); %d(del); %d(trans); %d(inversion)\n", VarNumVec[var_SUB], VarNumVec[var_INS], VarNumVec[var_DEL], VarNumVec[var_TNL] >> 1, VarNumVec[var_INV] >> 1);

	fprintf(log, "variant calling has been done in %lld seconds.\n", (long long)(time(NULL) - t));
	fprintf(stderr, "variant calling has been done in %lld seconds.\n", (long long)(time(NULL) - t));

	fclose(log);

	delete[] ThrIDarr; delete[] ThreadArr; delete[] BlockDepthArr;
}
