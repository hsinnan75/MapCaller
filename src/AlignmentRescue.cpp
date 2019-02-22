#include "structure.h"

static pthread_mutex_t Lock;

//AlnCan_t IdentifyBestAlnCan(vector<FragPair_t>& SimplePairVec)
//{
//	AlnCan_t AlnCan;
//	int i, j, score, score_thr, HeadIdx, num = (int)SimplePairVec.size();
//
//	sort(SimplePairVec.begin(), SimplePairVec.end(), CompByPosDiff);
//
//	//HeadIdx = 0; AlnCan.score = 0; score = SimplePairVec[0].rLen;
//	//for (i = 0, j = 1; j < num; i++, j++)
//	//{
//	//	if (SimplePairVec[j].rPos == SimplePairVec[i].rPos || abs(SimplePairVec[j].PosDiff - SimplePairVec[i].PosDiff) > MaxPosDiff)
//	//	{
//	//		if (score > AlnCan.score)
//	//		{
//	//			AlnCan.score = score; AlnCan.FragPairVec.clear();
//	//			copy(SimplePairVec.begin() + HeadIdx, SimplePairVec.begin() + j, back_inserter(AlnCan.FragPairVec));
//	//		}
//	//		HeadIdx = j;  score = SimplePairVec[j].rLen;
//	//	}
//	//	else score += SimplePairVec[j].rLen;
//	//}
//	//if (score > AlnCan.score)
//	//{
//	//	AlnCan.score = score; AlnCan.FragPairVec.clear();
//	//	copy(SimplePairVec.begin() + HeadIdx, SimplePairVec.begin() + j, back_inserter(AlnCan.FragPairVec));
//	//}
//	return AlnCan;
//}
AlnCan_t IdentifyBestAlnCan(vector<FragPair_t>& SimplePairVec)
{
	bool *UsedArr;
	AlnCan_t AlnCan;
	int i, j, num, score;

	AlnCan.score = 0; num = (int)SimplePairVec.size();
	for (num = (int)SimplePairVec.size(), i = 0; i < num;)
	{
		score = SimplePairVec[i].rLen;
		for (j = i + 1; j < num; j++)
		{
			if (SimplePairVec[j].PosDiff == SimplePairVec[i].PosDiff) score += SimplePairVec[j].rLen;
			else break;
		}
		if (j - i >= 1 && score > AlnCan.score)
		{
			AlnCan.score = score;
			AlnCan.FragPairVec.clear();
			copy(SimplePairVec.begin() + i, SimplePairVec.begin() + j, back_inserter(AlnCan.FragPairVec));
		}
		i = j;
	}
	//pthread_mutex_lock(&Lock);
	//printf("score=%d\n", AlnCan.score);
	//printf("Simple Pairs\n"); ShowSimplePairInfo(SimplePairVec);
	//printf("BestAlnCan\n"); ShowSimplePairInfo(AlnCan.FragPairVec);
	//printf("\n\n");
	//pthread_mutex_unlock(&Lock);

	return AlnCan;
}

int AlignmentRescue(uint32_t EstDist, ReadItem_t& read1, ReadItem_t& read2)
{
	char* RefSeg;
	AlnCan_t AlnCan;
	int64_t left_end, right_end;
	vector<KmerPair_t> KmerPairVec;
	vector<AlnCan_t>::iterator iter;
	vector<FragPair_t> SimplePairVec;
	vector<KmerItem_t> KmerVec1, KmerVec2;
	int i, j, slen, thr, score1, score2, num1, num2, iFixStrategy, nPaired = 0;

	for (score1 = 0, iter = read1.AlnCanVec.begin(); iter != read1.AlnCanVec.end(); iter++) if (iter->score > score1) score1 = iter->score;
	for (score2 = 0, iter = read2.AlnCanVec.begin(); iter != read2.AlnCanVec.end(); iter++) if (iter->score > score2) score2 = iter->score;

	if (score1 < (int)(read1.rlen>>2) && score2 < (int)(read2.rlen>>2)) return nPaired;
	else if (score1 - score2 > (int)(read2.rlen>>2)) iFixStrategy = 1;
	else if (score2 - score1 > (int)(read1.rlen>>2)) iFixStrategy = 2;
	else iFixStrategy = 3;

	//if (bDebugMode) printf("Rescue strategy = %d: score1=%d score2=%d\n", iFixStrategy, score1, score2);

	num1 = (int)read1.AlnCanVec.size(); num2 = (int)read2.AlnCanVec.size();
	if (iFixStrategy == 1 || iFixStrategy == 3) // map read2 with read1's AlnCans;
	{
		KmerVec1 = CreateKmerVecFromReadSeq(read2.rlen, read2.seq); thr = score1 >> 1;

		for (iter = read1.AlnCanVec.begin(); iter != read1.AlnCanVec.end(); iter++)
		{
			if (iter->score < thr || iter->PairedAlnCanIdx != -1) continue;
			left_end = iter->FragPairVec[0].PosDiff;
			right_end = iter->FragPairVec[0].PosDiff + EstDist + read2.rlen;
			if (right_end > TwoGenomeSize) right_end = TwoGenomeSize;
			if (ChrLocMap.lower_bound(left_end)->second != ChrLocMap.lower_bound(right_end)->second) continue;

			if ((slen = right_end - left_end) < read2.rlen) continue;
			RefSeg = RefSequence + left_end; KmerVec2 = CreateKmerVecFromReadSeq(slen, RefSeg);

			KmerPairVec = IdentifyCommonKmers(slen, KmerVec1, KmerVec2);
			SimplePairVec = GenerateSimplePairsFromCommonKmers(10, left_end, KmerPairVec);
			if (SimplePairVec.size() == 0) continue; else AlnCan = IdentifyBestAlnCan(SimplePairVec);
			if (AlnCan.score > score2)
			{
				//printf("Identify a AlnCan(read2)\n"); ShowSimplePairInfo(AlnCan.FragPairVec);
				nPaired++;
				iter->PairedAlnCanIdx = num2++;
				AlnCan.PairedAlnCanIdx = iter - read1.AlnCanVec.begin();
				read2.AlnCanVec.push_back(AlnCan);
			}
		}
	}
	if (iFixStrategy == 2 || iFixStrategy == 3) // map read1 with read2's AlnCans;
	{
		KmerVec1 = CreateKmerVecFromReadSeq(read1.rlen, read1.seq); thr = score2 >> 1;
		for (iter = read2.AlnCanVec.begin(); iter != read2.AlnCanVec.end(); iter++)
		{
			if (iter->score < thr || iter->PairedAlnCanIdx != -1) continue;
			left_end = iter->FragPairVec[0].PosDiff - EstDist;
			right_end = iter->FragPairVec[0].PosDiff + read1.rlen;
			if (right_end > TwoGenomeSize) right_end = TwoGenomeSize;
			if (ChrLocMap.lower_bound(left_end)->second != ChrLocMap.lower_bound(right_end)->second) continue;
			if ((slen = right_end - left_end) < read1.rlen) continue;
			RefSeg = RefSequence + left_end; KmerVec2 = CreateKmerVecFromReadSeq(slen, RefSeg);

			KmerPairVec = IdentifyCommonKmers(slen, KmerVec1, KmerVec2);
			SimplePairVec = GenerateSimplePairsFromCommonKmers(10, left_end, KmerPairVec);
			if (SimplePairVec.size() == 0) continue; else AlnCan = IdentifyBestAlnCan(SimplePairVec);
			if (AlnCan.score > score1)
			{
				//printf("Identify a AlnCan(read1)\n"); ShowSimplePairInfo(AlnCan.FragPairVec);
				nPaired++;
				iter->PairedAlnCanIdx = num1++;
				AlnCan.PairedAlnCanIdx = iter - read2.AlnCanVec.begin();
				read1.AlnCanVec.push_back(AlnCan);
			}
		}
	}
	//if (bDebugMode) printf("Identify %d resuced alignments\n", nPaired);

	return nPaired;
}
