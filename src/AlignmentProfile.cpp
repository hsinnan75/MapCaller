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

void UpdateProfile(ReadItem_t* read, vector<AlnCan_t>& AlnCanVec)
{
	int64_t gPos;
	string IndSeq;
	bool bChecked, bShow = false;
	int i, j, frag_len, ext_len, TailIdx, rPos, num;
	map<int64_t, map<string, uint16_t> >::iterator ind_iter, lower_iter, upper_iter;

	for (vector<AlnCan_t>::iterator iter = AlnCanVec.begin(); iter != AlnCanVec.end(); iter++)
	{
		if (iter->score == 0) continue;

		num = iter->FragPairVec.size(); TailIdx = num - 1;
		if (iter->FragPairVec[0].rLen == 0 && iter->FragPairVec[0].gLen == 0)
		{
			if (iter->FragPairVec[1].rPos > MinBreakPointSize)
			{
				gPos = iter->FragPairVec[0].gPos;
				if (gPos < GenomeSize) BreakPointMap[gPos]++;
				else BreakPointMap[(TwoGenomeSize - 1 - gPos)]++;
			}
			continue;
		}
		if (iter->FragPairVec[TailIdx].rLen == 0 && iter->FragPairVec[TailIdx].gLen == 0)
		{
			if ((read->rlen - iter->FragPairVec[TailIdx].rPos) > MinBreakPointSize)
			{
				gPos = iter->FragPairVec[TailIdx].gPos;
				if (gPos < GenomeSize) BreakPointMap[gPos]++;
				else BreakPointMap[TwoGenomeSize - 1 - gPos]++;
			}
			continue;
		}
		//if (bSomatic && iter->score < (read->rlen - 2)) continue;
		if (iter->FragPairVec[0].gPos < GenomeSize)
		{
			gPos = iter->FragPairVec[0].gPos;
			if (MappingRecordArr[gPos].iFirst == 10) continue;
			MappingRecordArr[gPos].iFirst++;

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
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					DeleteSeqMap[gPos][iter->FragPairVec[i].aln2]++;
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
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							DeleteSeqMap[gPos][IndSeq]++;
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
			gPos = TwoGenomeSize - iter->FragPairVec.rbegin()->gPos - iter->FragPairVec.rbegin()->gLen;
			if (MappingRecordArr[gPos].iFirst == 10) continue;
			MappingRecordArr[gPos].iFirst++;

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
					gPos = TwoGenomeSize - 1 - iter->FragPairVec[i].gPos;
					InsertSeqMap[gPos][iter->FragPairVec[i].aln1]++;
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					gPos = (TwoGenomeSize - iter->FragPairVec[i].gPos - iter->FragPairVec[i].gLen);
					DeleteSeqMap[gPos][iter->FragPairVec[i].aln2]++;
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
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							DeleteSeqMap[gPos][IndSeq]++;
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
