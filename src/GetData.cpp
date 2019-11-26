#include "structure.h"

int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;

//int IdentifyHeaderBoundary(char* str, int len)
//{
//	int i;
//
//	for (i = 1; i < len; i++)
//	{
//		if (str[i] == ' ' || str[i] == '/' || str[i] == '\t') return i;
//	}
//	return i;
//}

int IdentifyHeaderBegPos(char* str, int len)
{
	for (int i = 1; i < len; i++)
	{
		if (str[i] != '>' && str[i] != '@') return i;
	}
	return len - 1;
}

int IdentifyHeaderEndPos(char* str, int len)
{
	if (len > 100) len = 100;
	for (int i = 1; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '/' || isprint(str[i]) == 0) return i;
	}
	return len - 1;;
}

bool CheckReadFormat(const char* filename)
{
	char buf[1];
	gzFile file = gzopen(filename, "rb");
	gzread(file, buf, 1); gzclose(file);

	if (buf[0] == '@') return true; // fastq
	else return false;
}

ReadItem_t GetNextEntry(FILE *file)
{
	int p1, p2;
	ssize_t len;
	size_t size = 0;
	ReadItem_t read;
	char *buffer = NULL;

	read.header = read.seq = read.qual = NULL; read.rlen = 0;
	if ((len = getline(&buffer, &size, file)) != -1)
	{
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, buffer + 1, len); read.header[len] = '\0';
		p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
		if (FastQFormat)
		{
			if ((read.rlen = getline(&buffer, &size, file)) != -1)
			{
				read.seq = new char[read.rlen]; strncpy(read.seq, buffer, read.rlen);
				getline(&buffer, &size, file); getline(&buffer, &size, file); 
				read.qual = new char[read.rlen]; strncpy(read.qual, buffer, read.rlen);
				read.rlen -= 1; read.seq[read.rlen] = '\0'; read.qual[read.rlen] = '\0';
			}
			else read.rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read.rlen = (int)seq.length()) > 0)
			{
				read.seq = new char[read.rlen + 1];
				strcpy(read.seq, (char*)seq.c_str());
				read.seq[read.rlen] = '\0';
			}
		}
	}
	free(buffer);

	return read;
}

int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr)
{
	char* rseq;
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;
		iCount++;

		if (bSepLibrary) ReadArr[iCount] = GetNextEntry(file2);
		else if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;

		if (bPairEnd)
		{
			rseq = new char[ReadArr[iCount].rlen + 1];
			GetComplementarySeq(ReadArr[iCount].rlen, ReadArr[iCount].seq, rseq);
			copy(rseq, rseq + ReadArr[iCount].rlen, ReadArr[iCount].seq); delete[] rseq;
			if (FastQFormat)
			{
				string rqual = ReadArr[iCount].qual; reverse(rqual.begin(), rqual.end());
				copy(rqual.c_str(), rqual.c_str() + ReadArr[iCount].rlen, ReadArr[iCount].qual);
			}
		}
		iCount++;
		if (iCount == ReadChunkSize) break;
	}
	return iCount;
}

ReadItem_t gzGetNextEntry(gzFile file)
{
	int len;
	int p1, p2;
	ReadItem_t read;
	char buffer[1024];

	read.header = read.seq = read.qual = NULL; read.rlen = 0;
	if (gzgets(file, buffer, 1024) != NULL)
	{
		//len = IdentifyHeaderBoundary(buffer, strlen(buffer)- 1) - 1;
		//read.header = new char[len + 1]; strncpy(read.header, buffer+1, len); read.header[len] = '\0';
		len = strlen(buffer);
		if (len > 0 && (buffer[0] == '@' || buffer[0] == '>'))
		{
			p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
			read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
			gzgets(file, buffer, 1024); read.rlen = strlen(buffer) - 1;
			read.seq = new char[read.rlen + 1]; read.seq[read.rlen] = '\0';
			strncpy(read.seq, buffer, read.rlen);

			if (FastQFormat)
			{
				gzgets(file, buffer, 1024); gzgets(file, buffer, 1024);
				read.qual = new char[read.rlen + 1]; read.qual[read.rlen] = '\0';
				strncpy(read.qual, buffer, read.rlen);
			}
		}
	}
	return read;
}

int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr)
{
	char* rseq;
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		iCount++;
		if (bSepLibrary) ReadArr[iCount] = gzGetNextEntry(file2);
		else if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		if (bPairEnd)
		{
			rseq = new char[ReadArr[iCount].rlen + 1];
			GetComplementarySeq(ReadArr[iCount].rlen, ReadArr[iCount].seq, rseq);
			copy(rseq, rseq + ReadArr[iCount].rlen, ReadArr[iCount].seq); delete[] rseq;
			if (FastQFormat)
			{
				string rqual = ReadArr[iCount].qual; reverse(rqual.begin(), rqual.end());
				copy(rqual.c_str(), rqual.c_str() + ReadArr[iCount].rlen, ReadArr[iCount].qual);
			}
		}
		iCount++;
		if (iCount == ReadChunkSize) break;
	}
	return iCount;
}


bool CheckBWAIndexFiles(string IndexPrefix)
{
	fstream file;
	string filename;
	bool bChecked=true;

	filename = IndexPrefix + ".ann"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();
	
	filename = IndexPrefix + ".amb"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	filename = IndexPrefix + ".pac"; file.open(filename.c_str(), ios_base::in);
	if(!file.is_open()) return false; else file.close();

	return bChecked;
}
