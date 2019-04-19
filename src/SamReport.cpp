#include <cmath>
#include "structure.h"

#define MAPQ_COEF 30
#define Max_MAPQ  60

void SetSingledAlignmentFlag(ReadItem_t& read)
{
	int i;

	if (read.AlnSummary.score > read.AlnSummary.sub_score || bUnique == false) // unique mapping or bUnique=false
	{
		i = read.AlnSummary.BestAlnCanIdx;
		read.AlnCanVec[i].SamFlag = (read.AlnCanVec[i].orientation ? 0 : 0x10);
	}
	else if (read.AlnSummary.score > 0)
	{
		for (i = 0; i < (int)read.AlnCanVec.size(); i++)
		{
			if (read.AlnCanVec[i].score > 0) read.AlnCanVec[i].SamFlag = (read.AlnCanVec[i].orientation ? 0 : 0x10);
		}
	}
	else read.AlnCanVec[i].SamFlag = 0x4;
}

void SetPairedAlignmentFlag(ReadItem_t& read1, ReadItem_t& read2)
{
	int i, j;

	//read1
	if (read1.AlnSummary.score > read1.AlnSummary.sub_score) // unique mapping 
	{
		i = read1.AlnSummary.BestAlnCanIdx;
		read1.AlnCanVec[i].SamFlag = 0x41;  // first read
		read1.AlnCanVec[i].SamFlag |= (read1.AlnCanVec[i].orientation ? 0x20 : 0x10);

		if ((j = read1.AlnCanVec[i].PairedAlnCanIdx) != -1 && read2.AlnCanVec[j].score > 0) read1.AlnCanVec[i].SamFlag |= 0x2;// reads are mapped in a proper pair
		else
		{
			read1.AlnCanVec[i].SamFlag |= (read1.AlnCanVec[i].orientation ? 0x10 : 0x20);
			read1.AlnCanVec[i].SamFlag |= 0x8; // next segment unmapped
		}

	}
	else if (read1.AlnSummary.score > 0)
	{
		for (i = 0; i < (int)read1.AlnCanVec.size(); i++)
		{
			if (read1.AlnCanVec[i].score > 0)
			{
				read1.AlnCanVec[i].SamFlag = 0x41; // read1 is the first read in a pair
				read1.AlnCanVec[i].SamFlag |= (read1.AlnCanVec[i].orientation ? 0x20 : 0x10);
				if ((j = read1.AlnCanVec[i].PairedAlnCanIdx) != -1 && read2.AlnCanVec[j].score > 0) read1.AlnCanVec[i].SamFlag |= 0x2;// reads are mapped in a proper pair
				else read1.AlnCanVec[i].SamFlag |= 0x8; // next segment unmapped
			}
		}
	}
	//read2
	if (read2.AlnSummary.score > read2.AlnSummary.sub_score) // unique mapping or bUnique=false
	{
		j = read2.AlnSummary.BestAlnCanIdx;
		read2.AlnCanVec[j].SamFlag = 0x81; // read2 is the second read in a pair
		read2.AlnCanVec[j].SamFlag |= (read2.AlnCanVec[j].orientation ? 0x10 : 0x20);
		if ((i = read2.AlnCanVec[j].PairedAlnCanIdx) != -1 && read1.AlnCanVec[i].score > 0) read2.AlnCanVec[j].SamFlag |= 0x2;// reads are mapped in a proper pair
		else
		{
			read2.AlnCanVec[j].SamFlag |= (read2.AlnCanVec[j].orientation ? 0x20 : 0x10);
			read2.AlnCanVec[j].SamFlag |= 0x8; // next segment unmapped
		}
	}
	else if (read2.AlnSummary.score > 0)
	{
		for (j = 0; j < (int)read2.AlnCanVec.size(); j++)
		{
			if (read2.AlnCanVec[j].score > 0)
			{
				read2.AlnCanVec[j].SamFlag = 0x81; // read2 is the second read in a pair
				read2.AlnCanVec[j].SamFlag |= (read2.AlnCanVec[j].orientation ? 0x10 : 0x20);
				if ((i = read2.AlnCanVec[j].PairedAlnCanIdx) != -1 && read1.AlnCanVec[i].score > 0) read2.AlnCanVec[j].SamFlag |= 0x2;// reads are mapped in a proper pair
				else read2.AlnCanVec[j].SamFlag |= 0x8; // next segment unmapped
			}
		}
	}
}

int EvaluateMAPQ(ReadItem_t& read)
{
	int mapq;

	if (read.AlnSummary.score == 0 || read.AlnSummary.score == read.AlnSummary.sub_score) mapq = 0;
	else
	{
		if (read.AlnSummary.sub_score == 0 || read.AlnSummary.score - read.AlnSummary.sub_score > 5) mapq = Max_MAPQ;
		else
		{
			mapq = (int)(MAPQ_COEF * (1 - (float)(read.AlnSummary.score - read.AlnSummary.sub_score) / read.AlnSummary.score)*log(read.AlnSummary.score) + 0.4999);
			if (mapq > Max_MAPQ) mapq = Max_MAPQ;
		}
	}
	return mapq;
}

string ReverseCIGAR(string CIGAR)
{
	string RevCIGAR;
	int len, pos1, pos2;

	len = (int)CIGAR.length();

	for (pos1 = 0, pos2 = 1; pos1<len; pos2++)
	{
		if (isalpha(CIGAR[pos2]))
		{
			RevCIGAR.insert(0, CIGAR.substr(pos1, pos2 - pos1 + 1));
			pos1 = pos2 + 1;
		}
	}
	return RevCIGAR;
}

Coordinate_t GetAlnCoordinate(bool orientation, vector<FragPair_t>& FragPairVec)
{
	Coordinate_t coor;
	vector<FragPair_t>::iterator iter;

	if (orientation)
	{
		for (iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
		{
			if (iter->gLen > 0)
			{
				coor = DetermineCoordinate(iter->gPos);
				break;
			}
		}
	}
	else
	{
		for (iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
		{
			if (iter->gLen > 0)
			{
				coor = DetermineCoordinate(iter->gPos + iter->gLen - 1);
				break;
			}
		}
	}
	return coor;
}

int CalCoveredBase(string CIGAR)
{
	string tmp;
	int len, pos1, pos2, count = 0;

	len = (int)CIGAR.length();
	for (pos1 = 0, pos2 = 1; pos1 < len; pos2++)
	{
		if (isalpha(CIGAR[pos2]))
		{
			if (CIGAR[pos2] == 'M' || CIGAR[pos2] == 'I' || CIGAR[pos2] == 'S')
			{
				tmp = CIGAR.substr(pos1, pos2 - pos1);
				count += atoi(tmp.c_str());
			}
			pos1 = pos2 + 1;
		}
	}
	return count;
}

string GenerateCIGARstring(int rlen, bool orientation, vector<FragPair_t>& FragPairVec)
{
	string CIGAR;
	int i, j, len, num, c;
	char state = ' ', buf[10];

	if (!FragPairVec[0].bSimple)
	{
		if (orientation)
		{
			if (FragPairVec[0].rPos != 0)
			{
				sprintf(buf, "%dS", FragPairVec[0].rPos);
				CIGAR += buf;
			}
		}
		else
		{
			if ((c = rlen - (FragPairVec[0].rPos + FragPairVec[0].rLen)) > 0)
			{
				sprintf(buf, "%dS", c);
				CIGAR += buf;
			}
		}
	}
	for (num = (int)FragPairVec.size(), c = i = 0; i < num; i++)
	{
		if (FragPairVec[i].bSimple)
		{
			if (state != 'M')
			{
				if (c > 0)
				{
					sprintf(buf, "%d%c", c, state);
					CIGAR += buf;
				}
				state = 'M'; c = 0;
			}
			c += FragPairVec[i].rLen;
		}
		else if (FragPairVec[i].aln1.length() > 0)
		{
			len = (int)FragPairVec[i].aln1.length();
			for (j = 0; j < len; j++)
			{
				if (FragPairVec[i].aln1[j] == '-')
				{
					if (state != 'D')
					{
						if (c > 0)
						{
							sprintf(buf, "%d%c", c, state);
							CIGAR += buf;
						}
						state = 'D'; c = 0;
					}
					c++;
				}
				else if (FragPairVec[i].aln2[j] == '-')
				{
					if (state != 'I')
					{
						if (c > 0)
						{
							sprintf(buf, "%d%c", c, state);
							CIGAR += buf;
						}
						state = 'I'; c = 0;
					}
					c++;
				}
				else
				{
					if (state != 'M')
					{
						if (c > 0)
						{
							sprintf(buf, "%d%c", c, state); 
							CIGAR += buf;
						}
						state = 'M'; c = 0;
					}
					c++;
				}
			}
		}
		else if (FragPairVec[i].rLen > 0) // insertion
		{
			if (state != 'I')
			{
				if (c > 0)
				{
					sprintf(buf, "%d%c", c, state);
					CIGAR += buf;
				}
				state = 'I'; c = 0;
			}
			c += FragPairVec[i].rLen;
		}
		else if (FragPairVec[i].gLen > 0) // deletion
		{
			if (state != 'D')
			{
				if (c > 0)
				{
					sprintf(buf, "%d%c", c, state);
					CIGAR += buf;
				}
				state = 'D'; c = 0;
			}
			c += FragPairVec[i].gLen;
		}
	}
	if (c > 0)
	{
		sprintf(buf, "%d%c", c, state);
		CIGAR += buf;
	}
	if ((i = num - 1) > 0 && !FragPairVec[i].bSimple)
	{
		if (orientation)
		{
			if ((c = rlen - (FragPairVec[i].rPos + FragPairVec[i].rLen)) > 0)
			{
				sprintf(buf, "%dS", c);
				CIGAR += buf;
			}
		}
		else
		{
			if (FragPairVec[i].rPos != 0)
			{
				sprintf(buf, "%dS", FragPairVec[i].rPos);
				CIGAR += buf;
			}
		}
	}
	//if ((c = CalCoveredBase(CIGAR)) != rlen)
	//{
	//	printf("coverage = %d / %d\n", c, rlen);
	//	ShowSimplePairInfo(FragPairVec);
	//	printf("Final CIGAR=%s\n", CIGAR.c_str());
	//}
	return CIGAR;
}

void GetReverseQualityStr(int len, char* qual, char* rqual)
{
	int i, j;
	for (i = 0, j = len - 1; j != 0; i++, j--) rqual[j] = qual[i];
}

void GenerateSingleSamStream(ReadItem_t& read, vector<string>& SamStreamVec)
{
	int len;
	char* buffer = NULL;

	if (read.rlen < 1000) buffer = (char*)malloc((10000));
	else buffer = (char*)malloc((read.rlen * 10));

	if (read.AlnSummary.score == 0)
	{
		len = sprintf(buffer, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read.header, read.seq, (FastQFormat ? read.qual : "*"));
		SamStreamVec.push_back(buffer);
	}
	else
	{
		int i, mapq;
		string CIGAR;
		Coordinate_t coor;
		char *seq, *rseq, *rqual;

		SetSingledAlignmentFlag(read);

		mapq = EvaluateMAPQ(read); seq = read.seq; rseq = rqual = NULL;
		for (i = read.AlnSummary.BestAlnCanIdx; i < (int)read.AlnCanVec.size(); i++)
		{
			if (read.AlnCanVec[i].score == read.AlnSummary.score)
			{
				if(!read.AlnCanVec[i].orientation && rseq == NULL)
				{
					rseq = new char[read.rlen + 1]; 
					GetComplementarySeq(read.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = new char[read.rlen + 1]; rqual[read.rlen] = '\0';
						GetReverseQualityStr(read.rlen, read.qual, rqual);
					}
				}
				CIGAR = GenerateCIGARstring(read.rlen, read.AlnCanVec[i].orientation, read.AlnCanVec[i].FragPairVec);
				coor = GetAlnCoordinate(read.AlnCanVec[i].orientation, read.AlnCanVec[i].FragPairVec);
				len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read.header, read.AlnCanVec[i].SamFlag, ChromosomeVec[coor.ChromosomeIdx].name, coor.gPos, mapq, CIGAR.c_str(), (read.AlnCanVec[i].orientation ? seq : rseq), (FastQFormat ? (read.AlnCanVec[i].orientation ? read.qual : rqual) : "*"), read.rlen - read.AlnCanVec[i].score, read.AlnSummary.score, read.AlnSummary.sub_score);
				SamStreamVec.push_back(buffer);
				if (bUnique) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			if (FastQFormat) delete[] rqual;
			rseq = rqual = NULL;
		}
	}
	free(buffer);
}

void GeneratePairedSamStream(ReadItem_t& read1, ReadItem_t& read2, vector<string>& SamStreamVec)
{
	string CIGAR, rqual;
	Coordinate_t coor1, coor2;
	int i, j, mapq, dist, len, SamFlag;
	char *seq, *rseq, *buffer = NULL;

	len = (read1.rlen > read2.rlen ? read1.rlen : read2.rlen);
	if (len < 1000) buffer = (char*)malloc((10000));
	else buffer = (char*)malloc((len * 10));

	SetPairedAlignmentFlag(read1, read2);

	if (read1.AlnSummary.score == 0)
	{
		SamFlag = 0x1; // read1 is the first read in a pair
		SamFlag |= 0x4; // segment unmapped
		SamFlag |= 0x40; // second fragment
		if (read2.AlnSummary.score == 0) SamFlag |= 0x8; // next segment unmapped
		else if (read2.AlnCanVec.size() > 0)
		{
			SamFlag |= (read2.AlnCanVec[read2.AlnSummary.BestAlnCanIdx].orientation ? 0x10 : 0x20);
			SamFlag |= (read2.AlnCanVec[read2.AlnSummary.BestAlnCanIdx].orientation ? 0x20 : 0x10);
		}
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read1.header, SamFlag, read1.seq, (FastQFormat ? read1.qual: "*"));
		SamStreamVec.push_back(buffer);
	}
	else
	{
		mapq = EvaluateMAPQ(read1); seq = read1.seq; rseq = NULL;
		for (i = read1.AlnSummary.BestAlnCanIdx; i < (int)read1.AlnCanVec.size(); i++)
		{
			//fprintf(stderr, "check %d / %d, score = %d\n", i + 1, (int)read1.AlnCanVec.size(), read1.AlnCanVec[i].score);
			if (read1.AlnCanVec[i].score == read1.AlnSummary.score)
			{
				if (!read1.AlnCanVec[i].orientation && rseq == NULL)
				{
					rseq = new char[read1.rlen + 1]; GetComplementarySeq(read1.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read1.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				CIGAR = GenerateCIGARstring(read1.rlen, read1.AlnCanVec[i].orientation, read1.AlnCanVec[i].FragPairVec);
				coor1 = GetAlnCoordinate(read1.AlnCanVec[i].orientation, read1.AlnCanVec[i].FragPairVec);
				if ((j = read1.AlnCanVec[i].PairedAlnCanIdx) != -1 && read2.AlnCanVec[j].score == read2.AlnSummary.score)
				{
					coor2 = GetAlnCoordinate(read2.AlnCanVec[j].orientation, read2.AlnCanVec[j].FragPairVec);
					dist = (int)(coor2.gPos - coor1.gPos + (read1.AlnCanVec[i].orientation ? read2.rlen : 0 - read1.rlen));
					len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read1.header, read1.AlnCanVec[i].SamFlag, ChromosomeVec[coor1.ChromosomeIdx].name, coor1.gPos, mapq, CIGAR.c_str(), coor2.gPos, dist, (read1.AlnCanVec[i].orientation ? seq : rseq), (FastQFormat ? (read1.AlnCanVec[i].orientation ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.AlnCanVec[i].score, read1.AlnSummary.score, read1.AlnSummary.sub_score);
				}
				else len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read1.header, read1.AlnCanVec[i].SamFlag, ChromosomeVec[coor1.ChromosomeIdx].name, coor1.gPos, mapq, CIGAR.c_str(), (read1.AlnCanVec[i].orientation ? seq : rseq), (FastQFormat ? (read1.AlnCanVec[i].orientation ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.AlnCanVec[i].score, read1.AlnSummary.score, read1.AlnSummary.sub_score);
				SamStreamVec.push_back(buffer);
				if (bUnique) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}
	if (read2.AlnSummary.score == 0)
	{
		SamFlag = 0x1;
		SamFlag |= 0x4; // segment unmapped
		SamFlag |= 0x80; // second fragment
		if(read1.AlnSummary.score == 0) SamFlag |= 0x8; // next segment unmapped
		else if (read1.AlnCanVec.size() > 0)
		{
			SamFlag |= (read1.AlnCanVec[read1.AlnSummary.BestAlnCanIdx].orientation ? 0x10 : 0x20);
			SamFlag |= (read1.AlnCanVec[read1.AlnSummary.BestAlnCanIdx].orientation ? 0x20 : 0x10);
		}
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0", read2.header, SamFlag, read2.seq, (FastQFormat ? read2.qual : "*"));
		SamStreamVec.push_back(buffer);
	}
	else
	{
		mapq = EvaluateMAPQ(read2); seq = read2.seq; rseq = NULL;
		for (j = read2.AlnSummary.BestAlnCanIdx; j < (int)read2.AlnCanVec.size(); j++)
		{
			if (read2.AlnCanVec[j].score == read2.AlnSummary.score)
			{
				if (!read2.AlnCanVec[j].orientation && rseq == NULL)
				{
					rseq = new char[read2.rlen + 1]; GetComplementarySeq(read2.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read2.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				CIGAR = GenerateCIGARstring(read2.rlen, read2.AlnCanVec[j].orientation, read2.AlnCanVec[j].FragPairVec);
				coor2 = GetAlnCoordinate(read2.AlnCanVec[j].orientation, read2.AlnCanVec[j].FragPairVec);
				if ((i = read2.AlnCanVec[j].PairedAlnCanIdx) != -1 && read1.AlnCanVec[i].score == read1.AlnSummary.score)
				{
					coor1 = GetAlnCoordinate(read1.AlnCanVec[i].orientation, read1.AlnCanVec[i].FragPairVec);
					dist = 0 - (int)(coor2.gPos - coor1.gPos + (read1.AlnCanVec[i].orientation ? read2.rlen : 0 - read1.rlen));
					len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read2.header, read2.AlnCanVec[j].SamFlag, ChromosomeVec[coor2.ChromosomeIdx].name, coor2.gPos, mapq, CIGAR.c_str(), coor1.gPos, dist, (read2.AlnCanVec[j].orientation ? seq : rseq), (FastQFormat ? (read2.AlnCanVec[j].orientation ? read2.qual : rqual.c_str()) : "*"), read2.rlen - read2.AlnCanVec[j].score, read2.AlnSummary.score, read2.AlnSummary.sub_score);
				}
				else len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d", read2.header, read2.AlnCanVec[j].SamFlag, ChromosomeVec[coor2.ChromosomeIdx].name, coor2.gPos, mapq, CIGAR.c_str(), (read2.AlnCanVec[j].orientation ? seq : rseq), (FastQFormat ? (read2.AlnCanVec[j].orientation ? read2.qual : rqual.c_str()) : "*"), read2.rlen - read2.AlnCanVec[j].score, read2.AlnSummary.score, read2.AlnSummary.sub_score);
				SamStreamVec.push_back(buffer);
				if (bUnique) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}
	free(buffer);
}
