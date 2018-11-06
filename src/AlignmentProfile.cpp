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

		if (iter->FragPairVec[0].gPos < GenomeSize)
		{
			//ShowSimplePairInfo(iter->FragPairVec);
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
				else if (iter->FragPairVec[i].rLen == 0 && iter->FragPairVec[i].gLen == 0) // breakpoint
				{
					if ((i == 0 && iter->FragPairVec[i + 1].rPos > MinBreakPointSize) || (i == TailIdx && read->rlen - (iter->FragPairVec[i - 1].rPos + iter->FragPairVec[i - 1].rLen) > MinBreakPointSize))
					{
						BreakPointMap[gPos]++;
					}
				}
				else if (iter->FragPairVec[i].gLen == 0) // ins
				{
					IndSeq.resize(iter->FragPairVec[i].rLen); strncpy((char*)IndSeq.c_str(), read->seq + rPos, iter->FragPairVec[i].rLen);
					IndSeq[iter->FragPairVec[i].rLen] = '\0';
					bChecked = false; lower_iter = InsertSeqMap.lower_bound(gPos - shift); upper_iter = InsertSeqMap.upper_bound(gPos + shift);
					//if (gPos > ObserveBegPos && gPos < ObserveEndPos) printf("lower_Bound=%lld, upper_Bound=%lld\n", lower_iter->first, upper_iter->first);
					for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
					{
						if (ind_iter != InsertSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
						{
							bChecked = true;
							break;
						}
					}
					if (!bChecked) InsertSeqMap[gPos][IndSeq]++;
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					IndSeq.resize(iter->FragPairVec[i].gLen); strncpy((char*)IndSeq.c_str(), RefSequence + gPos, iter->FragPairVec[i].gLen);
					IndSeq[iter->FragPairVec[i].gLen] = '\0';
					bChecked = false; lower_iter = DeleteSeqMap.lower_bound(gPos - shift); upper_iter = DeleteSeqMap.upper_bound(gPos + shift);
					for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
					{
						if (ind_iter != DeleteSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
						{
							bChecked = true;
							break;
						}
					}
					if (!bChecked) DeleteSeqMap[gPos][IndSeq]++;
				}
				else
				{
					for (frag_len = (int)iter->FragPairVec[i].aln1.length(), j = 0; j < frag_len;)
					{
						if (iter->FragPairVec[i].aln2[j] == '-') // ins
						{
							ext_len = 1; while (iter->FragPairVec[i].aln2[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln1.substr(j, ext_len);
							bChecked = false; lower_iter = InsertSeqMap.lower_bound(gPos - shift); upper_iter = InsertSeqMap.upper_bound(gPos + shift);
							//if (gPos > ObserveBegPos && gPos < ObserveEndPos) printf("lower_Bound=%lld, upper_Bound=%lld\n", lower_iter->first, upper_iter->first);
							for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
							{
								if (ind_iter != InsertSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
								{
									bChecked = true;
									break;
								}
							}
							if (!bChecked) InsertSeqMap[gPos][IndSeq]++;
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							bChecked = false; lower_iter = DeleteSeqMap.lower_bound(gPos - shift); upper_iter = DeleteSeqMap.upper_bound(gPos + shift);
							for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
							{
								if (ind_iter != DeleteSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
								{
									bChecked = true;
									break;
								}
							}
							if (!bChecked) DeleteSeqMap[gPos][IndSeq]++;
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
				else if (iter->FragPairVec[i].rLen == 0 && iter->FragPairVec[i].gLen == 0) // breakpoint
				{
					if ((i == 0 && iter->FragPairVec[i + 1].rPos > MinBreakPointSize) || (i == TailIdx && read->rlen - (iter->FragPairVec[i - 1].rPos + iter->FragPairVec[i - 1].rLen) > MinBreakPointSize))
					{
						BreakPointMap[(TwoGenomeSize - 1 - iter->FragPairVec[i].gPos)]++;
					}
				}
				else if (iter->FragPairVec[i].gLen == 0) // ins
				{
					IndSeq.resize(iter->FragPairVec[i].rLen); strncpy((char*)IndSeq.c_str(), read->seq + iter->FragPairVec[i].rPos, iter->FragPairVec[i].rLen);
					SelfComplementarySeq(iter->FragPairVec[i].rLen, (char*)IndSeq.c_str()); IndSeq[iter->FragPairVec[i].rLen] = '\0';
					gPos = TwoGenomeSize - 1 - iter->FragPairVec[i].gPos;
					bChecked = false; lower_iter = InsertSeqMap.lower_bound(gPos - shift); upper_iter = InsertSeqMap.upper_bound(gPos + shift);
					//if (gPos > ObserveBegPos && gPos < ObserveEndPos) printf("lower_Bound=%lld, upper_Bound=%lld\n", lower_iter->first, upper_iter->first);
					for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
					{
						if (ind_iter != InsertSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
						{
							bChecked = true;
							break;
						}
					}
					if (!bChecked) InsertSeqMap[gPos][IndSeq]++;
				}
				else if (iter->FragPairVec[i].rLen == 0) // del
				{
					IndSeq.resize(iter->FragPairVec[i].gLen); strncpy((char*)IndSeq.c_str(), RefSequence + iter->FragPairVec[i].gPos, iter->FragPairVec[i].gLen);
					SelfComplementarySeq(iter->FragPairVec[i].gLen, (char*)IndSeq.c_str()); IndSeq[iter->FragPairVec[i].gLen] = '\0';
					gPos = (TwoGenomeSize - iter->FragPairVec[i].gPos - iter->FragPairVec[i].gLen);
					bChecked = false; lower_iter = DeleteSeqMap.lower_bound(gPos - shift); upper_iter = DeleteSeqMap.upper_bound(gPos + shift);
					for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
					{
						if (ind_iter != DeleteSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
						{
							bChecked = true;
							break;
						}
					}
					if (!bChecked) DeleteSeqMap[gPos][IndSeq]++;
				}
				else
				{
					rPos = read->rlen - (iter->FragPairVec[i].rPos + iter->FragPairVec[i].rLen);
					gPos = TwoGenomeSize - (iter->FragPairVec[i].gPos + iter->FragPairVec[i].gLen);
					//SelfComplementarySeq((int)iter->FragPairVec[i].aln1.length(), (char*)iter->FragPairVec[i].aln1.c_str());
					//SelfComplementarySeq((int)iter->FragPairVec[i].aln2.length(), (char*)iter->FragPairVec[i].aln2.c_str());

					for (frag_len = (int)iter->FragPairVec[i].aln1.length(), j = 0; j < frag_len;)
					{
						if (iter->FragPairVec[i].aln2[j] == '-') // ins
						{
							ext_len = 1; while (iter->FragPairVec[i].aln2[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln1.substr(j, ext_len);
							bChecked = false; lower_iter = InsertSeqMap.lower_bound(gPos - shift); upper_iter = InsertSeqMap.upper_bound(gPos + shift);
							for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
							{
								if (ind_iter != InsertSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
								{
									bChecked = true;
									break;
								}
							}
							if (!bChecked) InsertSeqMap[gPos][IndSeq]++;
							j += ext_len; rPos += ext_len;
						}
						else if (iter->FragPairVec[i].aln1[j] == '-') // del
						{
							ext_len = 1; while (iter->FragPairVec[i].aln1[j + ext_len] == '-') ext_len++;
							IndSeq = iter->FragPairVec[i].aln2.substr(j, ext_len);
							bChecked = false; lower_iter = DeleteSeqMap.lower_bound(gPos - shift); upper_iter = DeleteSeqMap.upper_bound(gPos + shift);
							for (ind_iter = lower_iter; ind_iter != upper_iter; ind_iter++)
							{
								if (ind_iter != DeleteSeqMap.end() && (ind_iter->first - gPos) < shift && CheckIndStrOccu(IndSeq, ind_iter))
								{
									bChecked = true;
									break;
								}
							}
							if (!bChecked) DeleteSeqMap[gPos][IndSeq]++;
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
		if (iter->score == 0) continue;

		for (FragIter = iter->FragPairVec.begin(); FragIter != iter->FragPairVec.end(); FragIter++)
		{
			if (FragIter->gLen == 0) continue;
			if (FragIter->gPos < GenomeSize) gPos = FragIter->gPos, gPosEnd = gPos + FragIter->gLen;
			else gPos = TwoGenomeSize - (FragIter->gPos + FragIter->gLen), gPosEnd = gPos + FragIter->gLen;
			for (; gPos < gPosEnd; gPos++) MappingRecordArr[gPos].multi_hit++;
		}
	}
}
