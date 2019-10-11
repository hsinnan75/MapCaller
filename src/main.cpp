#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "0.9.9.19";

string CmdLine;
float FrequencyThr;
uint8_t iMaxDuplicate;
time_t StartProcessTime;
MappingRecord_t* MappingRecordArr = NULL;
vector<string> ReadFileNameVec1, ReadFileNameVec2;
int64_t ObservGenomicPos, ObserveBegPos, ObserveEndPos;
char *RefSequence, *IndexFileName, *SamFileName, *VcfFileName;
int iThreadNum, iPloidy, FragmentSize, MaxMisMatches, MaxClipSize, MinAlleleDepth, MinIndFreq, MinVarConfScore, MinUnmappedSize;
bool bDebugMode, bFilter, bPairEnd, bUnique, bSAMoutput, bSAMFormat, bGVCF, bMonomorphic, bVCFoutput, bSomatic, gzCompressed, FastQFormat, NW_ALG;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "MapCaller v%s\n\n", VersionStr);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>]\n\n", program);
	fprintf(stderr, "Options: -i STR        BWT_Index_Prefix\n");
	fprintf(stderr, "         -f            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -f2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -t INT        number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -size         sequencing fragment size [%d]\n", FragmentSize);
	fprintf(stderr, "         -ad INT       minimal ALT allele count [%d]\n", MinAlleleDepth);
	fprintf(stderr, "         -dup INT      maximal PCR duplicates [%d]\n", iMaxDuplicate);
	fprintf(stderr, "         -maxmm INT    maximal mismatches in read alignment [%d]\n", MaxMisMatches);
	fprintf(stderr, "         -maxclip INT  maximal clip size at either ends [%d]\n", MaxClipSize);
	fprintf(stderr, "         -sam          SAM output filename [NULL]\n");
	fprintf(stderr, "         -bam          BAM output filename [NULL]\n");
	fprintf(stderr, "         -alg STR      gapped alignment algorithm (option: nw|ksw2)\n");
	fprintf(stderr, "         -vcf          VCF output filename [%s]\n", VcfFileName);
	fprintf(stderr, "         -gvcf         GVCF mode [false]\n");
	fprintf(stderr, "         -monomorphic  report all loci which do not have any potential alternates.\n");
	fprintf(stderr, "         -min_gap INT  the minimal gap(unmapped) size to be reported [%d].\n", MinUnmappedSize);
	fprintf(stderr, "         -ploidy INT   number of sets of chromosomes in a cell (1:monoploid, 2:diploid) [%d]\n", iPloidy);
	fprintf(stderr, "         -m            output multiple alignments\n");
	fprintf(stderr, "         -somatic      detect somatic mutations [false]\n");
	fprintf(stderr, "         -no_vcf       No VCF output [false]\n");
	fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stderr, "         -filter       apply variant filters (under test) [false]\n");
	fprintf(stderr, "         -v            version\n");
	fprintf(stderr, "\n");
}

bool CheckOutputFileName(char *FileName)
{
	struct stat s;
	bool bRet = true;

	if (stat(FileName, &s) == 0)
	{
		if (s.st_mode & S_IFDIR)
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is a directory!\n", FileName);
		}
		else if (s.st_mode & S_IFREG)
		{
		}
		else
		{
			bRet = false;
			fprintf(stderr, "Warning: %s is not a regular file!\n", FileName);
		}
	}
	return bRet;
}

bool CheckInputFiles(vector<string>& ReadFileNameVec)
{
	struct stat s;
	string filetype;
	bool bRet = true;

	for (vector<string>::iterator iter = ReadFileNameVec.begin(); iter != ReadFileNameVec.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file:[%s]\n", (char*)iter->c_str());
			break;
		}
		else
		{
			filetype = iter->substr(iter->find_last_of('.') + 1);
			for (string::iterator ii = filetype.begin(); ii != filetype.end(); ii++) *ii = tolower(*ii);
			if (filetype != "fq" && filetype != "fa" && filetype != "fastq" && filetype != "fasta" && filetype != "gz")
			{
				bRet = false;
				fprintf(stderr, "Wrong file type:[%s]\n", iter->c_str());
				break;
			}
		}
	}
	return bRet;
}

void ReadLibInput(const char* LibFileName)
{
	fstream file;
	stringstream ss;
	string str, fn1, fn2;

	file.open(LibFileName, ios_base::in);
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;
		if (str[0] == '#') continue;
		ss.clear(); ss.str(str); ss >> fn1 >> fn2;
		if (fn1 != "") ReadFileNameVec1.push_back(fn1);
		if (fn2 != "") ReadFileNameVec2.push_back(fn2);
	}
	file.close();
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	bGVCF = false;
	iPloidy = 2;
	iThreadNum = 16;
	bPairEnd = false;
	bDebugMode = false;
	bUnique = true;
	bFilter = false;
	NW_ALG = true;
	FastQFormat = true;
	bSAMoutput = false;
	bSAMFormat = true;
	bSomatic = false;
	bVCFoutput = true;
	gzCompressed = false;
	bMonomorphic = false;

	MinIndFreq = 5;
	MaxClipSize = 5;
	MaxMisMatches = 5;
	iMaxDuplicate = 5;
	FragmentSize = 500;
	MinAlleleDepth = 5;
	FrequencyThr = 0.2;
	MinVarConfScore = 10;
	MinUnmappedSize = 50;
	VcfFileName = (char*)"output.vcf";
	RefSequence = IndexFileName = SamFileName = NULL;
	ObservGenomicPos = ObserveBegPos = ObserveEndPos = -1;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "update") == 0)
	{
		system("git fetch;git merge origin/master master;make");
		exit(0);
	}
	else
	{
		for (CmdLine = argv[0], i = 1; i < argc; i++) CmdLine += " " + (string)argv[i];

		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i" && i + 1 < argc) IndexFileName = argv[++i];
			else if (parameter == "-f")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec1.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-f2")
			{
				while (++i < argc && argv[i][0] != '-') ReadFileNameVec2.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-lib" && i + 1 < argc)
			{
				ReadLibInput(argv[++i]);
			}
			else if (parameter == "-t" && i + 1 < argc)
			{
				if ((iThreadNum = atoi(argv[++i])) <= 0)
				{
					fprintf(stderr, "Warning! The thread number should be positive!\n");
					iThreadNum = 4;
				}
			}
			else if (parameter == "-dup" && i + 1 < argc)
			{
				if (atoi(argv[++i]) <= 255) iMaxDuplicate = (int8_t)atoi(argv[i]);
				else fprintf(stderr, "Warning! The PCR-duplicate range is [1-255]!\n");
			}
			else if (parameter == "-filter") bFilter = true;
			else if (parameter == "-size" && i + 1 < argc) FragmentSize = atoi(argv[++i]);
			else if (parameter == "-ad" && i + 1 < argc) MinAlleleDepth = atoi(argv[++i]);
			else if (parameter == "-ind" && i + 1 < argc) MinIndFreq = atoi(argv[++i]);
			else if (parameter == "-ploidy"  && i + 1 < argc)
			{
				if ((iPloidy = atoi(argv[++i])) > 2)
				{
					iPloidy = 2;
					fprintf(stderr, "Warning! MapCaller only supports monoploid and diploid!\n");
				}
			}
			else if ((parameter == "-sam") && i + 1 < argc)
			{
				bSAMoutput = true;
				bSAMFormat = true;
				SamFileName = argv[++i];
			}
			else if ((parameter == "-bam") && i + 1 < argc)
			{
				bSAMoutput = true;
				bSAMFormat = false;
				SamFileName = argv[++i];
			}
			else if ((parameter == "-alg") && i + 1 < argc)
			{
				str = argv[++i];
				if (str == "ksw2") NW_ALG = false;
				else NW_ALG = true; //nw
			}
			else if (parameter == "-maxmm" && i + 1 < argc) MaxMisMatches = atoi(argv[++i]);
			else if (parameter == "-maxclip" && i + 1 < argc) MaxClipSize = atoi(argv[++i]);
			else if (parameter == "-min_gap" && i + 1 < argc) MinUnmappedSize = atoi(argv[++i]);
			else if (parameter == "-vcf" && i + 1 < argc) VcfFileName = argv[++i];
			else if (parameter == "-gvcf") bGVCF = true;
			else if (parameter == "-monomorphic") bMonomorphic = true;
			else if (parameter == "-no_vcf") bVCFoutput = false;
			else if (parameter == "-somatic") bSomatic = true;
			else if (parameter == "-pair" || parameter == "-p") bPairEnd = true;
			else if (parameter == "-obs" && i + 1 < argc) ObservGenomicPos = atoi(argv[++i]);
			else if (parameter == "-obr" && i + 2 < argc)
			{
				ObserveBegPos = atoi(argv[++i]);
				ObserveEndPos = atoi(argv[++i]);
				fprintf(stderr, "obr[%lld - %lld]\n", (long long) ObserveBegPos, (long long)ObserveEndPos);
			}
			else if (parameter == "-m") bUnique = false;
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			else if (parameter == "-v" || parameter == "--version")
			{
				fprintf(stderr, "MapCaller v%s\n\n", VersionStr);
				exit(0);
			}
			else
			{
				fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
				ShowProgramUsage(argv[0]);
				exit(0);
			}
		}
		if (bGVCF && bMonomorphic) bGVCF = false;
		if (iMaxDuplicate <= 0) iMaxDuplicate = 10;
		//fprintf(stderr, "Read1:\n"); for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
		//fprintf(stderr, "Read2:\n"); for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
		if (ReadFileNameVec1.size() == 0)
		{
			fprintf(stderr, "Warning! Please specify a valid read input!\n");
			ShowProgramUsage(argv[0]);
			exit(0);
		}
		if (ReadFileNameVec2.size() > 0 && ReadFileNameVec1.size() != ReadFileNameVec2.size())
		{
			fprintf(stderr, "Warning! Paired-end reads input numbers do not match!\n");
			fprintf(stderr, "Read1:\n"); for (vector<string>::iterator iter = ReadFileNameVec1.begin(); iter != ReadFileNameVec1.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			fprintf(stderr, "Read2:\n"); for (vector<string>::iterator iter = ReadFileNameVec2.begin(); iter != ReadFileNameVec2.end(); iter++) fprintf(stderr, "\t%s\n", (char*)iter->c_str());
			exit(0);
		}
		if (CheckInputFiles(ReadFileNameVec1) == false || CheckInputFiles(ReadFileNameVec2) == false) exit(0);
		if (SamFileName != NULL && CheckOutputFileName(SamFileName) == false) exit(0);
		if (VcfFileName != NULL && CheckOutputFileName(VcfFileName) == false) exit(0);
		//if (MinAlleleFreq > MinBaseDepth) MinAlleleFreq = MinBaseDepth;
		//fprintf(stderr, "AD=%d\n", MinAlleleFreq);

		if (IndexFileName != NULL && CheckBWAIndexFiles(IndexFileName)) RefIdx = bwa_idx_load(IndexFileName);
		else
		{
			fprintf(stderr, "Warning! Please specify a valid reference index!\n");
			ShowProgramUsage(argv[0]);
			exit(0);
		}
		if (RefIdx == 0) fprintf(stderr, "\n\nError! Index files are corrupt!\n");
		else
		{
			Refbwt = RefIdx->bwt;
			RestoreReferenceInfo();

			if (GenomeSize <= 0)
			{
				fprintf(stderr, "Reference genome is empty\n");
				exit(1);
			}
			if (bVCFoutput)
			{
				fprintf(stderr, "Initialize the alignment profile...\n");
				MappingRecordArr = new MappingRecord_t[GenomeSize]();
			}
			StartProcessTime = time(NULL);
			Mapping();
			if (bVCFoutput) VariantCalling();
			bwa_idx_destroy(RefIdx);
			if (RefSequence != NULL) delete[] RefSequence;
			if (MappingRecordArr != NULL) delete[] MappingRecordArr;

			fprintf(stderr, "All done! It took %lld seconds to complete the data analysis.\n", (long long)(time(NULL)- StartProcessTime));
		}
	}
	return 0;
}
