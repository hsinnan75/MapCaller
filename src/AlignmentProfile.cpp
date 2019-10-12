#include "structure.h"

#define shift 10
#define MinBreakPointSize 20

map<int64_t, uint16_t> BreakPointMap;
map<int64_t, map<string, uint16_t> > InsertSeqMap, DeleteSeqMap;

bool CheckIndStrOccu(string& IndSeq, map<int64_t, map<string, uint16_t> >::iterator iter)
{
	bool bChecked = false;
	map<string, uint16_t>::iterator SeqIter;

	for (SeqIter = iter->second.begin(); SeqIter != iter->second.end(); SeqIter++)
	{
		if (SeqIter->first == IndSeq)
		{
			SeqIter->second++;
			bChecked = true;
			break;
		}
	}
	return bChecked;
}

int CheckMismatch(vector<FragPair_t>& FragPairVec)
{
	int i, len, mis = 0;
	vector<FragPair_t>::iterator iter;
	for (iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (!iter->bSimple && iter->rLen > 0 && iter->gLen > 0)
		{
			len = (int)iter->aln1.length();
			for (i = 0; i < len; i++) if (iter->aln1[i] != iter->aln2[i] && iter->aln1[i] != '-' && iter->aln2[i] != '-') mis++;
		}
	}
	return mis;
}

void UpdateProfile(ReadItem_t* read, vector<AlnCan_t>& AlnCanVec)
{
	int64_t gPos;
	string IndSeq;
	int i, j, frag_len, ext_len, rPos, num;
	map<int64_t, map<string, uint16_t> >::iterator ind_iter, lower_iter, upper_iter;

	for (vector<AlnCan_t>::iterator iter = AlnCanVec.begin(); iter != AlnCanVec.end(); iter++)
	{
		if (iter->score == 0) continue;

		num = iter->FragPairVec.size();
		if (iter->FragPairVec.begin()->rLen == 0 && iter->FragPairVec.begin()->gLen == 0)
		{
			//if (iter->FragPairVec[1].rPos > MinBreakPointSize)
			if (iter->FragPairVec.begin()->rPos > MinBreakPointSize)
			{
				gPos = iter->FragPairVec.begin()->gPos;
				if (gPos < GenomeSize) BreakPointMap[gPos]++;
				else BreakPointMap[(TwoGenomeSize - 1 - gPos)]++;
			}
			if (iter->FragPairVec.begin()->rPos > MaxClipSize) continue;
		}
		if (iter->FragPairVec.rbegin()->rLen == 0 && iter->FragPairVec.rbegin()->gLen == 0)
		{
			if ((read->rlen - iter->FragPairVec.rbegin()->rPos) > MinBreakPointSize)
			{
				gPos = iter->FragPairVec.rbegin()->gPos;
				if (gPos < GenomeSize) BreakPointMap[gPos]++;
				else BreakPointMap[TwoGenomeSize - 1 - gPos]++;
			}
			if ((read->rlen - iter->FragPairVec.rbegin()->rPos) > MaxClipSize) continue;
		}
		//if (iter->FragPairVec.begin()->rPos > MaxClipSize || (read->rlen - iter->FragPairVec.rbegin()->rPos - iter->FragPairVec.rbegin()->rLen) > MaxClipSize) continue;
		//if (bSomatic && CheckMismatch(iter->FragPairVec) > 2) continue;
		//if (bSomatic && iter->PairedAlnCanIdx == -1) continue;

		if (iter->orientation) gPos = iter->FragPairVec.begin()->gPos;
		else gPos = TwoGenomeSize - (iter->FragPairVec.begin()->gPos + iter->FragPairVec.begin()->gLen);
		if (MappingRecordArr[gPos].readCount < iMaxDuplicate) MappingRecordArr[gPos].readCount++;
		else continue;

		if (iter->orientation) //iter->FragPairVec[0].gPos < GenomeSize
		{
			for (i = 0; i < num; i++)
			{
				rPos = iter->FragPairVec[i].rPos; gPos = iter->FragPairVec[i].gPos;
				if (iter->FragPairVec[i].bSimple)
				{
					//printf("gPos[%lld-%lld]=%d / %lld\n", gPos, gPos + iter->FragPairVec[i].gLen - 1, iter->FragPairVec[i].gLen, GenomeSize);
					for (j = 0; j < iter->FragPairVec[i].rLen; j++, rPos++, gPos++)
					{
						switch (read->seq[rPos])
						{
						case 'A': MappingRecordArr[gPos].A++; break;
						case 'C': MappingRecordArr[gPos].C++; break;
						case 'G': MappingRecordArr[gPos].G++; break;
						case 'T': MappingRecordArr[gPos].T++; break;
						}
					}
				}
				else if (iter->FragPairVec[i].gLen == 0) // ins
				{
					InsertSeqMap[gPos - 1][iter->FragPairVec[i].aln1]++;
					//printf("%lld: %s\n", gPos - 1, iter->FragPairVec[i].aln1.c_str());
					//ShowSimplePairInfo(iter->FragPairVec);
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					DeleteSeqMap[gPos-1][iter->FragPairVec[i].aln2]++;
					//printf("%lld: %s\n", (long long)(gPos - 1), iter->FragPairVec[i].aln2.c_str());
					//ShowSimplePairInfo(iter->FragPairVec);
				}
				else
				{
					for (frag_len = (int)iter->FragPairVec[i].aln1.length(), j = 0; j < frag_len;)
					{
						if (iter->FragPairVec[i].aln2[j] == '-') // ins
						{
							ext_len = 1; while (iter->FragPairVec[i].aln2[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln1.substr(j, ext_len);
							InsertSeqMap[gPos - 1][IndSeq]++;
							//printf("%lld: %s\n", gPos - 1, IndSeq.c_str());
							//ShowSimplePairInfo(iter->FragPairVec);
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							DeleteSeqMap[gPos-1][IndSeq]++;
							//printf("%lld: %s\n", (long long)(gPos - 1), IndSeq.c_str());
							//ShowSimplePairInfo(iter->FragPairVec);
							j += ext_len; gPos += ext_len;
						}
						else
						{
							switch (iter->FragPairVec[i].aln1[j])
							{
							case 'A': MappingRecordArr[gPos].A++; break;
							case 'C': MappingRecordArr[gPos].C++; break;
							case 'G': MappingRecordArr[gPos].G++; break;
							case 'T': MappingRecordArr[gPos].T++; break;
							}
							j++; rPos++; gPos++;
						}
					}
				}
			}
		}
		else
		{
			for (i = 0; i < num; i++)
			{
				if (iter->FragPairVec[i].bSimple)
				{
					rPos = iter->FragPairVec[i].rPos; gPos = TwoGenomeSize - 1 - iter->FragPairVec[i].gPos;
					for (j = 0; j < iter->FragPairVec[i].rLen; j++, rPos++, gPos--)
					{
						switch (read->seq[rPos])
						{
						case 'A': MappingRecordArr[gPos].T++; break;
						case 'C': MappingRecordArr[gPos].G++; break;
						case 'G': MappingRecordArr[gPos].C++; break;
						case 'T': MappingRecordArr[gPos].A++; break;
						}
					}
				}
				else if (iter->FragPairVec[i].gLen == 0) // ins
				{
					gPos = TwoGenomeSize - iter->FragPairVec[i].gPos;
					InsertSeqMap[gPos - 1][iter->FragPairVec[i].aln1]++;
					//printf("%lld: %s\n", gPos - 1, iter->FragPairVec[i].aln1.c_str());
					//ShowSimplePairInfo(iter->FragPairVec);
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					gPos = (TwoGenomeSize - iter->FragPairVec[i].gPos - iter->FragPairVec[i].gLen);
					DeleteSeqMap[gPos - 1][iter->FragPairVec[i].aln2]++;
					//printf("%lld: %s\n", (long long)(gPos - 1), iter->FragPairVec[i].aln2.c_str());
					//ShowSimplePairInfo(iter->FragPairVec);
				}
				else
				{
					rPos = read->rlen - (iter->FragPairVec[i].rPos + iter->FragPairVec[i].rLen);
					gPos = TwoGenomeSize - (iter->FragPairVec[i].gPos + iter->FragPairVec[i].gLen);
					for (frag_len = (int)iter->FragPairVec[i].aln1.length(), j = 0; j < frag_len;)
					{
						if (iter->FragPairVec[i].aln2[j] == '-') // ins
						{
							ext_len = 1; while (iter->FragPairVec[i].aln2[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln1.substr(j, ext_len);
							InsertSeqMap[gPos - 1][IndSeq]++;
							//printf("%lld: %s\n", (long long)(gPos - 1), IndSeq.c_str());
							//ShowSimplePairInfo(iter->FragPairVec);
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							DeleteSeqMap[gPos - 1][IndSeq]++;
							//printf("%lld: %s\n", (long long)(gPos - 1), IndSeq.c_str());
							//ShowSimplePairInfo(iter->FragPairVec);
							j += ext_len; gPos += ext_len;
						}
						else
						{
							switch (iter->FragPairVec[i].aln1[j])
							{
							case 'A': MappingRecordArr[gPos].A++; break;
							case 'C': MappingRecordArr[gPos].C++; break;
							case 'G': MappingRecordArr[gPos].G++; break;
							case 'T': MappingRecordArr[gPos].T++; break;
							}
							j++; rPos++; gPos++;
						}
					}
				}
			}
		}
	}
	//if (bShow)
	//{
	//	pthread_mutex_lock(&Lock);
	//	printf("%s\n%s\n", read->header, read->seq);
	//	ShowFragPairCluster(AlnCanVec);
	//	pthread_mutex_unlock(&Lock);
	//}
}

void UpdateMultiHitCount(ReadItem_t* read, vector<AlnCan_t>& AlnCanVec)
{
	int64_t gPos, gPosEnd;
	vector<AlnCan_t>::iterator iter;
	vector<FragPair_t>::iterator FragIter;

	for (iter = AlnCanVec.begin(); iter != AlnCanVec.end(); iter++)
	{
		if (iter->score != 0)
		{
			if (iter->orientation) gPos = iter->FragPairVec.begin()->gPos, gPosEnd = iter->FragPairVec.rbegin()->gPos + iter->FragPairVec.rbegin()->gLen;
			else gPos = TwoGenomeSize - (iter->FragPairVec.begin()->gPos + iter->FragPairVec.begin()->gLen), gPosEnd = TwoGenomeSize - iter->FragPairVec.rbegin()->gPos;
			//printf("%lld - %lld (len=%d / %d)\n", gPos, gPosEnd - 1, gPosEnd - gPos, read->rlen);
			for (; gPos < gPosEnd; gPos++) if (MappingRecordArr[gPos].multi_hit < 65535) MappingRecordArr[gPos].multi_hit++;
		}
	}
}
