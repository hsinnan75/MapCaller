#include "structure.h"
#define MinSeqIdy 0.5
#define MinAlnBlcokSize 5

float CalFragAlnSeqIdy(string& aln1, string& aln2)
{
	int i, len = (int)aln1.length(), n = 0, mis = 0;

	for (i = 0; i < len; i++)
	{
		if (aln1[i] != '-' && aln2[i] != '-')
		{
			n++;
			if (aln1[i] != aln2[i]) mis++;
		}
	}
	if (n != len && mis > 0) return 0;
	else if(n > 0) return 1.0*(n - mis) / n;
	else return 1;
}

bool CompByReadPos(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.rPos == p2.rPos) return p1.gPos < p2.gPos;
	else return p1.rPos < p2.rPos;
}

void RemoveNullFragPairs(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end();)
	{
		if (iter->rLen == 0) iter = FragPairVec.erase(iter);
		else iter++;
	}
}

bool RemoveOverlaps(vector<FragPair_t>& FragPairVec)
{
	bool bOverlap = false;
	int i, j, overlap_size, num = (int)FragPairVec.size();

	for (i = 0, j = 1; j < num; i++, j++)
	{
		if (FragPairVec[i].rPos == FragPairVec[j].rPos)
		{
			bOverlap = true;
			FragPairVec[i].rLen = FragPairVec[i].gLen = 0;
		}
		else if (FragPairVec[i].gPos >= FragPairVec[j].gPos || (FragPairVec[i].gPos + FragPairVec[i].gLen) > FragPairVec[j].gPos)
		{
			// shrink FragPairVec[i]
			bOverlap = true;
			overlap_size = FragPairVec[i].gPos + FragPairVec[i].gLen - FragPairVec[j].gPos;
			if ((FragPairVec[i].rLen -= overlap_size) < 0) FragPairVec[i].rLen = 0;
			if ((FragPairVec[i].gLen -= overlap_size) < 0) FragPairVec[i].gLen = 0;
		}
	}
	//if (bOverlap && bDebugMode)
	//{
	//	printf("After RemoveOverlaps\n");
	//	ShowSimplePairInfo(FragPairVec);
	//}
	return bOverlap;
}

void IdentifyNormalPairs(int rlen, vector<FragPair_t>& FragPairVec)
{
	FragPair_t FragPair;
	int i, j, rGaps, gGaps, num = (int)FragPairVec.size();

	FragPair.bSimple = false;

	for (i = 0, j = 1; j < num; i++, j++)
	{
		if ((rGaps = FragPairVec[j].rPos - (FragPairVec[i].rPos + FragPairVec[i].rLen)) < 0) rGaps = 0;
		if ((gGaps = FragPairVec[j].gPos - (FragPairVec[i].gPos + FragPairVec[i].gLen)) < 0) gGaps = 0;

		if (rGaps > 0 || gGaps > 0)
		{
			FragPair.rPos = FragPairVec[i].rPos + FragPairVec[i].rLen;
			FragPair.gPos = FragPairVec[i].gPos + FragPairVec[i].gLen;
			FragPair.PosDiff = FragPair.gPos - FragPair.rPos;
			FragPair.rLen = rGaps; FragPair.gLen = gGaps;
			FragPairVec.push_back(FragPair);
			//if (bDebugMode) printf("insert a normal pair: r[%d-%d] g[%lld-%lld] and r[%d-%d] g[%lld-%lld]: r[%d-%d] g[%lld-%lld]\n", FragPairVec[i].rPos, FragPairVec[i].rPos + FragPairVec[i].rLen - 1, FragPairVec[i].gPos, FragPairVec[i].gPos + FragPairVec[i].gLen - 1, FragPairVec[j].rPos, FragPairVec[j].rPos + FragPairVec[j].rLen - 1, FragPairVec[j].gPos, FragPairVec[j].gPos + FragPairVec[j].gLen - 1, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.gPos, FragPair.gPos + FragPair.gLen - 1);
		}
	}
	if ((int)FragPairVec.size() > num) inplace_merge(FragPairVec.begin(), FragPairVec.begin() + num, FragPairVec.end(), CompByReadPos);

	// Check missing blocks at both ends
	if (FragPairVec[0].rPos > 0)
	{
		FragPair.rPos = 0;
		FragPair.gPos = FragPair.PosDiff = FragPairVec[0].PosDiff;
		FragPair.rLen = FragPair.gLen = FragPairVec[0].rPos;
		FragPairVec.insert(FragPairVec.begin(), FragPair);
	}
	num = (int)FragPairVec.size();
	if (num > 0 && (FragPairVec[num - 1].rPos + FragPairVec[num - 1].rLen) < rlen)
	{
		FragPair.rPos = FragPairVec[num - 1].rPos + FragPairVec[num - 1].rLen;
		FragPair.gPos = FragPairVec[num - 1].gPos + FragPairVec[num - 1].gLen;
		FragPair.PosDiff = FragPairVec[num - 1].PosDiff;
		FragPair.rLen = FragPair.gLen = rlen - FragPair.rPos;
		FragPairVec.push_back(FragPair);
	}
}

void RemoveEmptyFragPairs(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end();)
	{
		if (iter->rLen == 0 && iter->gLen == 0) iter = FragPairVec.erase(iter);
		else iter++;
	}
}

bool CheckAlnCanCoverage(int rlen, vector<FragPair_t>& FragPairVec)
{
	bool bChecked = true;
	int TotalCovLength = 0;

	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->rLen > 0) TotalCovLength += iter->rLen;
	}
	if (TotalCovLength != rlen) bChecked = false;

	return bChecked;
}

int CalFragPairMismatches(int len, string& str1, string& str2)
{
	int i, mismatch;

	for (mismatch = i = 0; i < len; i++)
	{
		if (str1[i] != str2[i]) mismatch++;
	}
	return mismatch;
}

int CalFragPairMatches(int len, string& str1, string& str2)
{
	int i, match;

	for (match = i = 0; i < len; i++)
	{
		if (str1[i] == str2[i]) match++;
	}
	return match;
}

void ProcessNormalPair(char* seq, FragPair_t& fp)
{
	//if (fp.rLen > 0 && fp.gLen > 0)
	{
		int n;

		if (fp.rLen > 0)
		{
			fp.aln1.resize(fp.rLen);
			strncpy((char*)fp.aln1.c_str(), seq + fp.rPos, fp.rLen);
		}
		else fp.aln1.assign(fp.gLen, '-');

		if (fp.gLen > 0)
		{
			fp.aln2.resize(fp.gLen);
			strncpy((char*)fp.aln2.c_str(), RefSequence + fp.gPos, fp.gLen);
		}
		else fp.aln2.assign(fp.rLen, '-');

		//if (fp.rLen == 0 || fp.gLen == 0)
		//{
		//	printf("%s\n%s\n", fp.aln1.c_str(), fp.aln2.c_str());
		//}

		if (fp.gPos >= GenomeSize) // reverse sequence
		{
			if (fp.rLen > 0) SelfComplementarySeq(fp.rLen, (char*)fp.aln1.c_str());
			if (fp.gLen > 0) SelfComplementarySeq(fp.gLen, (char*)fp.aln2.c_str());
		}
		if (fp.rLen > 0 && fp.gLen > 0 && (fp.rLen != fp.gLen || ((n = CalFragPairMismatches(fp.rLen, fp.aln1, fp.aln2)) > 1 && n >= (int)(fp.rLen*0.2))))
		{
			if (NW_ALG) nw_alignment(fp.rLen, fp.aln1, fp.gLen, fp.aln2);
			else ksw2_alignment(fp.rLen, fp.aln1, fp.gLen, fp.aln2);
			//if (bDebugMode) ShowFragmentPair(seq, fp), printf("Alignment for normal pair:\n%s\n\%s\n\n", fp.aln1.c_str(), fp.aln2.c_str());
		}
	}
}

bool CheckLocalAlignmentQuality(FragPair_t& fp)
{
	int i, n, mis, len, AlnType, iStatus;

	AlnType = -1; len = (int)fp.aln1.length(); n = mis = iStatus = 0;
	for (i = 0; i < len; i++)
	{
		if (fp.aln1[i] == '-') // del
		{
			if (AlnType != 0)
			{
				AlnType = 0;
				iStatus++;
			}
		}
		else if (fp.aln2[i] == '-') // ins
		{
			if (AlnType != 1)
			{
				AlnType = 1; 
				iStatus++;
			}
		}
		else // type 2
		{
			n++; if (fp.aln1[i] != fp.aln2[i]) mis++;
			if (AlnType != 2)
			{
				AlnType = 2;
				iStatus++;
			}
		}
	}
	if (iStatus >= 4 || (mis >= 3 && mis >= (int)(n*0.3)))
	{
		//if (bDebugMode) printf("BadAlignment\n%s\n%s\nIdy=%.4f\n", fp.aln1.c_str(), fp.aln2.c_str(), CalFragAlnSeqIdy(fp.aln1, fp.aln2));
		return false;
	}
	else return true;
}

int EvaluateAlignmentScore(vector<FragPair_t>& FragPairVec)
{
	int len, score = 0;
	vector<FragPair_t>::iterator iter;

	for (iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->bSimple) score += iter->rLen;
		else if ((len = (int)iter->aln1.length()) > 0) score += CalFragPairMatches(len, iter->aln1, iter->aln2);
	}
	return score;
}

int FindMisMatchNumber(vector<FragPair_t>& FragPairVec)
{
	int i, len, mismatch = 0;
	vector<FragPair_t>::iterator iter;
	for (iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (!iter->bSimple)
		{
			for (len = (int)iter->aln1.length(), i = 0; i < len; i++)
			{
				if (iter->aln1[i] != iter->aln2[i] && iter->aln1[i] != '-' && iter->aln2[i] != '-') mismatch++;
			}
		}
	}
	return mismatch;
}

void RemoveHeadingGaps(bool bFirstFragPair, FragPair_t& FragPair)
{
	int j, len, Rshrink = 0, Gshrink = 0;

	for (len = (int)FragPair.aln1.length(), j = 0; j < len; j++)
	{
		if (FragPair.aln1[j] == '-') Gshrink++;
		else if (FragPair.aln2[j] == '-') Rshrink++;
		else break;
	}
	if (j > 0)
	{
		//if (bDebugMode) printf("Before:\n%s\n%s\n", FragPair.aln1.c_str(), FragPair.aln2.c_str());
		FragPair.aln1 = FragPair.aln1.substr(j);
		FragPair.aln2 = FragPair.aln2.substr(j);
		//if (bDebugMode) printf("After:\n%s\n%s\n\n",FragPair.aln1.c_str(), FragPair.aln2.c_str());
		FragPair.rLen -= Rshrink; FragPair.gLen -= Gshrink; 
		if(bFirstFragPair) FragPair.rPos += Rshrink, FragPair.gPos += Gshrink;
	}
}

void RemoveTailingGaps(bool bFirstFragPair, FragPair_t& FragPair)
{
	int j, len, Rshrink = 0, Gshrink = 0;

	for (len = (int)FragPair.aln1.length(), j = len - 1; j >= 0; j--)
	{
		if (FragPair.aln1[j] == '-') Gshrink++;
		else if (FragPair.aln2[j] == '-') Rshrink++;
		else break;
	}
	if (++j < len)
	{
		//if(bDebugMode) printf("Before:\n%s\n%s\n", FragPair.aln1.c_str(), FragPair.aln2.c_str());
		FragPair.aln1.resize(j);
		FragPair.aln2.resize(j);
		//if (bDebugMode) printf("After:\n%s\n%s\n\n", FragPair.aln1.c_str(), FragPair.aln2.c_str());
		FragPair.rLen -= Rshrink, FragPair.gLen -= Gshrink;
		if (bFirstFragPair) FragPair.rPos += Rshrink, FragPair.gPos += Gshrink;
	}
}

bool ProduceReadAlignment(ReadItem_t& read)
{
	bool bHead, bTail;
	int i, FragPairNum, TailIdx;
	vector<AlnCan_t>::iterator iter;

	for (iter = read.AlnCanVec.begin(); iter != read.AlnCanVec.end(); iter++)
	{
		if (iter->score == 0) continue;

		sort(iter->FragPairVec.begin(), iter->FragPairVec.end(), CompByReadPos);
		if (RemoveOverlaps(iter->FragPairVec)) RemoveNullFragPairs(iter->FragPairVec);
		IdentifyNormalPairs(read.rlen, iter->FragPairVec);
		if (CheckAlignmentValidity(iter->FragPairVec) == false)
		{
			iter->score = 0;
			continue;
		}
		//if (CheckAlnCanCoverage(read.rlen, iter->FragPairVec) == false)
		//{
		//	pthread_mutex_lock(&Lock);
		//	printf("read: %s\n", read.header);
		//	ShowSimplePairInfo(iter->FragPairVec);
		//	pthread_mutex_unlock(&Lock);
		//	iter->score = 0;
		//	continue;
		//}
		//if (bDebugMode)
		//{
		//	printf("Process aln_can#%d: score=%d\n", iter - read.AlnCanVec.begin() + 1, iter->score);
		//	ShowSimplePairInfo(iter->FragPairVec);
		//}
		//bQualityCheck = read.rlen < 150 ? true : false;
		bHead = bTail = true; FragPairNum = (int)iter->FragPairVec.size(), TailIdx = FragPairNum - 1;
		for (i = 0; i < FragPairNum; i++)
		{
			if (!iter->FragPairVec[i].bSimple)
			{
				ProcessNormalPair(read.seq, iter->FragPairVec[i]);
				if (i == 0)
				{
					if(iter->FragPairVec[i].gPos < GenomeSize) RemoveHeadingGaps(true, iter->FragPairVec[i]);
					else RemoveTailingGaps(true, iter->FragPairVec[i]);

					if (iter->FragPairVec[i].aln1.length() >= MinAlnBlcokSize && CheckLocalAlignmentQuality(iter->FragPairVec[i]) == false)
					{
						//printf("read:%s\n", read.header); ShowSimplePairInfo(iter->FragPairVec);
						bHead = false;
						iter->FragPairVec[i].rLen = iter->FragPairVec[i].gLen = 0;
						iter->FragPairVec[i].aln1.clear(); iter->FragPairVec[i].aln2.clear();
						iter->FragPairVec[i].rPos = iter->FragPairVec[i + 1].rPos;
						iter->FragPairVec[i].gPos = iter->FragPairVec[i + 1].gPos;
					}
				}
				else if (i == TailIdx)
				{
					if (iter->FragPairVec[i].gPos < GenomeSize) RemoveTailingGaps(false, iter->FragPairVec[i]);
					else RemoveHeadingGaps(false, iter->FragPairVec[i]);

					if (iter->FragPairVec[i].aln1.length() >= MinAlnBlcokSize && CheckLocalAlignmentQuality(iter->FragPairVec[i]) == false)
					{
						//printf("read:%s\n", read.header); ShowSimplePairInfo(iter->FragPairVec);
						bTail = false;
						iter->FragPairVec[i].rLen = iter->FragPairVec[i].gLen = 0;
						iter->FragPairVec[i].aln1.clear(); iter->FragPairVec[i].aln2.clear();
					}
				}
				else
				{
					if (iter->FragPairVec[i].rLen >= MinAlnBlcokSize && iter->FragPairVec[i].gLen >= MinAlnBlcokSize && CheckLocalAlignmentQuality(iter->FragPairVec[i]) == false)
					{
						iter->score = 0;
						break;
					}
				}
			}
		}
		//Coordinate_t coor = GetAlnCoordinate(iter->FragPairVec.begin()->gPos < GenomeSize ? true : false, iter->FragPairVec);
		//if (coor.ChromosomeIdx == 0 && coor.gPos >= ObserveBegPos && coor.gPos + read.rlen < ObserveEndPos)
		//{
		//	//Display alignments
		//	pthread_mutex_lock(&Lock);
		//	printf("read = %s, score = %d / %d\n", read.header, EvaluateAlignmentScore(iter->FragPairVec), read.rlen);
		//	ShowSimplePairInfo(iter->FragPairVec);
		//	pthread_mutex_unlock(&Lock);
		//}
		//if (bDebugMode)
		//{
		//	printf("Done mapping: head(%s), tail(%s)\n", bHead ? "Yes" : "No", bTail ? "Yes" : "No");
		//	ShowSimplePairInfo(iter->FragPairVec);
		//}
		if (iter->score == 0) continue;
		else if (!bHead && !bTail) iter->score = 0;
		else
		{
			iter->score = EvaluateAlignmentScore(iter->FragPairVec);
			if (iter->score == 0) continue;
			if (iter->score < (int)(read.rlen*0.95) && FindMisMatchNumber(iter->FragPairVec) > (int)(read.rlen*0.05)) iter->score = 0;
			else
			{
				iter->orientation = (iter->FragPairVec[0].gPos < GenomeSize ? true : false);
				if (!iter->orientation) reverse(iter->FragPairVec.begin(), iter->FragPairVec.end());

				if (iter->score > read.AlnSummary.score)
				{
					read.AlnSummary.score = iter->score;
					read.AlnSummary.BestAlnCanIdx = (int)(iter - read.AlnCanVec.begin());
				}
				else if (iter->score > read.AlnSummary.sub_score) read.AlnSummary.sub_score = iter->score;
			}
		}
	}
	for (iter = read.AlnCanVec.begin(); iter != read.AlnCanVec.end(); iter++) if (iter->score < read.AlnSummary.score) iter->score = 0;
	//if (ObserveBegPos != -1)
	//{
	//	for (iter = read.AlnCanVec.begin(); iter != read.AlnCanVec.end(); iter++)
	//	{
	//		//Display alignments
	//		if (iter->score == 0) continue;
	//		Coordinate_t coor = GetAlnCoordinate(iter->FragPairVec.begin()->gPos < GenomeSize ? true : false, iter->FragPairVec);
	//		if (coor.ChromosomeIdx == 0 && coor.gPos >= ObserveBegPos && coor.gPos + read.rlen < ObserveEndPos)
	//		{
	//			//Display alignments
	//			//pthread_mutex_lock(&Lock);
	//			printf("read: %s, score=%d (%d/%d) len=%d, PairedIdx=%d\n\n", read.header, iter->score, read.AlnSummary.score, read.AlnSummary.sub_score, read.rlen, iter->PairedAlnCanIdx);
	//			ShowSimplePairInfo(iter->FragPairVec);
	//			//pthread_mutex_unlock(&Lock);
	//		}
	//	}
	//}
	return (read.AlnSummary.score > 0);
}
