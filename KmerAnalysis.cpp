#include "structure.h"

bool CompByKmerID(const KmerItem_t& p1, const KmerItem_t& p2)
{
	return p1.wid < p2.wid;
}

bool CompByKmerRPos(const KmerPair_t& p1, const KmerPair_t& p2)
{
	if (p1.rPos == p2.rPos) return p1.gPos < p2.gPos;
	else return p1.rPos < p2.rPos;
}

bool CompByKmerGPos(const KmerPair_t& p1, const KmerPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}

bool CompByKmerPosDiff(const KmerPair_t& p1, const KmerPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

uint32_t CreateKmerID(const char* seq, short pos)
{
	uint32_t i, id, end_pos = pos + KmerSize;

	for (id = 0, i = pos; i < end_pos; i++) id = (id << 2) + nst_nt4_table[(int)seq[i]];

	return id;
}

string DecodeWordID(uint32_t id)
{
	string str;
	int i, index;

	for (i = 0; i < KmerSize; i++)
	{
		index = (id & 3); //id % 4;
		switch (index)
		{
		case 0: str.push_back('A'); break;
		case 1: str.push_back('C'); break;
		case 2: str.push_back('G'); break;
		default: str.push_back('T'); break;
		}
		id = (id >> 2); // id/=4
	}
	reverse(str.begin(), str.end());

	return str;
}

vector<KmerItem_t> CreateKmerVecFromReadSeq(int len, char* seq)
{
	vector<KmerItem_t> vec;
	KmerItem_t WordPosition;
	uint32_t count, head, tail;

	tail = 0; count = 0;

	while (count < KmerSize && tail < len)
	{
		if (seq[tail++] != 'N') count++;
		else count = 0;
	}
	if (count == KmerSize) // found the first kmer
	{
		WordPosition.pos = (head = tail - KmerSize); WordPosition.wid = CreateKmerID(seq, head);
		vec.push_back(WordPosition);

		for (head += 1; tail < len; head++, tail++)
		{
			if (seq[tail] != 'N')
			{
				WordPosition.pos = head;
				WordPosition.wid = ((WordPosition.wid & KmerPower) << 2) + nst_nt4_table[(int)seq[tail]];
				vec.push_back(WordPosition);
			}
			else
			{
				// find next kmer without 'N'
				count = 0; tail++;
				while (count < KmerSize && tail < len)
				{
					if (seq[tail++] != 'N') count++;
					else count = 0;
				}
				if (count == KmerSize)
				{
					WordPosition.pos = (head = tail - KmerSize); WordPosition.wid = CreateKmerID(seq, head);
					vec.push_back(WordPosition);
				}
				else break;
			}
		}
		sort(vec.begin(), vec.end(), CompByKmerID);
	}
	return vec;
}

vector<KmerPair_t> IdentifyCommonKmers(int MaxShift, vector<KmerItem_t>& vec1, vector<KmerItem_t>& vec2)
{
	int i, num;
	uint32_t wid;
	KmerPair_t KmerPair;
	vector<KmerPair_t> KmerPairVec;
	vector<KmerItem_t>::iterator KmerPtr;

	num = (int)vec1.size();
	for (i = 0; i < num; i++)
	{
		wid = vec1[i].wid; KmerPtr = lower_bound(vec2.begin(), vec2.end(), vec1[i], CompByKmerID);
		while (KmerPtr != vec2.end() && KmerPtr->wid == wid)
		{
			if ((KmerPtr->pos >= vec1[i].pos && KmerPtr->pos - vec1[i].pos < MaxShift) || ((KmerPtr->pos < vec1[i].pos && vec1[i].pos - KmerPtr->pos < MaxShift)))
			{
				KmerPair.rPos = vec1[i].pos;
				KmerPair.PosDiff = (KmerPair.gPos = KmerPtr->pos) - KmerPair.rPos;
				KmerPairVec.push_back(KmerPair);
			}
			KmerPtr++;
		}
	}
	sort(KmerPairVec.begin(), KmerPairVec.end(), CompByKmerPosDiff);

	return KmerPairVec;
}

vector<FragPair_t> GenerateSimplePairsFromCommonKmers(int thr, int64_t gPos, vector<KmerPair_t>& KmerPairVec)
{
	FragPair_t FragPair;
	int i, j, l, PosDiff, n_pos, num;
	vector<FragPair_t> SimplePairVec;

	FragPair.bSimple = true;

	num = (int)KmerPairVec.size();
	for (i = 0; i < num;)
	{
		//printf("KmerPair %d: %d PosDiff = %d\n", i+1, KmerPairVec[i].rPos, KmerPairVec[i].PosDiff);
		for (PosDiff = KmerPairVec[i].PosDiff, n_pos = KmerPairVec[i].rPos + 1, j = i + 1; j < num; j++)
		{
			if (KmerPairVec[j].rPos != n_pos || KmerPairVec[j].PosDiff != PosDiff) break;
			//printf("KmerPair: %d PosDiff = %d\n", KmerPairVec[j].rPos, KmerPairVec[j].PosDiff);
			n_pos++;
		}
		if ((l = KmerSize + (j - 1 - i)) >= thr)
		{
			//printf("Seed found! r=%d, l=%d\n\n", KmerPairVec[i].rPos, l);
			FragPair.rPos = KmerPairVec[i].rPos;
			FragPair.gPos = KmerPairVec[i].gPos + gPos;
			FragPair.PosDiff = KmerPairVec[i].PosDiff + gPos;
			FragPair.rLen = FragPair.gLen = l;
			SimplePairVec.push_back(FragPair);
		}
		i = j;
	}
	return SimplePairVec;
}
//
//vector<SeedPair_t> GenerateSimplePairsFromFragmentPair(int MaxDist, int len1, char* frag1, int len2, char* frag2)
//{
//	vector<int> rPosVec;
//	vector<KmerPair_t> KmerPairVec;
//	vector<SeedPair_t> SimplePairVec;
//	vector<KmerItem_t> KmerVec1, KmerVec2;
//
//	KmerVec1 = CreateKmerVecFromReadSeq(len1, frag1); KmerVec2 = CreateKmerVecFromReadSeq(len2, frag2);
//	KmerPairVec = IdentifyCommonKmers(MaxDist, KmerVec1, KmerVec2);
//	SimplePairVec = GenerateSimplePairsFromCommonKmers(8, KmerPairVec);
//
//	sort(SimplePairVec.begin(), SimplePairVec.end(), CompByGenomePos);
//	//if (bDebugMode) printf("Raw seeds from common kmers:\n"), ShowSeedInfo(SimplePairVec);
//
//	return SimplePairVec;
//}

