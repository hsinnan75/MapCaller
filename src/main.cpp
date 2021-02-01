#include "structure.h"

#ifdef __cplusplus
extern "C"
{
	int bwa_idx_build(const char *fa, const char *prefix);
}
#endif

bwt_t *Refbwt;
bwaidx_t *RefIdx;
const char* VersionStr = "0.9.9.41";

string CmdLine;
int iChromsomeNum;
uint8_t iMaxDuplicate;
time_t StartProcessTime;
map<string, int> ChrIdMap;
map<int64_t, int> PosChrIdMap;
map<int64_t, bool> KnowSiteMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;
float FrequencyThr, MaxMisMatchRate;
MappingRecord_t* MappingRecordArr = NULL;
vector<string> ReadFileNameVec1, ReadFileNameVec2;
int64_t ObservGenomicPos, ObserveBegPos, ObserveEndPos;
pthread_mutex_t LibraryLock, ProfileLock, OutputLock, VarLock;
char *RefSequence, *RefFileName, *KnownSiteFileName, *IndexFileName, *SamFileName, *VcfFileName, *LogFileName, *sample_id;
int iThreadNum, MaxPosDiff, iPloidy, FragmentSize, MaxClipSize, MinReadDepth, MinAlleleDepth, MinVarConfScore, MinCNVsize, MinUnmappedSize;
bool bDebugMode, bFilter, bPairEnd, bUnique, bSAMoutput, bSAMFormat, bGVCF, bMonomorphic, bVCFoutput, bSomatic, gzCompressed, FastQFormat, NW_ALG;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "MapCaller v%s\n\n", VersionStr);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>]\n\n", program);
	fprintf(stderr, "Options: -i STR        BWT_Index_Prefix\n");
	fprintf(stderr, "         -r STR        Reference filename (format:fa)\n");
	fprintf(stderr, "         -f            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -f2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stderr, "         -t INT        number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -size         sequencing fragment size [%d]\n", FragmentSize);
	fprintf(stderr, "         -indel INT	maximal indel size [%d]\n", MaxPosDiff);
	//fprintf(stderr, "         -dp INT       minimal read depth [%d]\n", MinReadDepth);
	fprintf(stderr, "         -ad INT       minimal ALT allele count [%d]\n", MinAlleleDepth);
	fprintf(stderr, "         -dup INT      maximal PCR duplicates [%d]\n", iMaxDuplicate);
	fprintf(stderr, "         -maxmm FLOAT  maximal mismatch rate in read alignment [%.2f]\n", MaxMisMatchRate);
	fprintf(stderr, "         -maxclip INT  maximal clip size at either ends [%d]\n", MaxClipSize);
	fprintf(stderr, "         -sam          SAM output filename [NULL]\n");
	fprintf(stderr, "         -bam          BAM output filename [NULL]\n");
	fprintf(stderr, "         -alg STR      gapped alignment algorithm (option: nw|ksw2)\n");
	fprintf(stderr, "         -vcf          VCF output filename [%s]\n", VcfFileName);
	fprintf(stderr, "         -gvcf         GVCF mode [false]\n");
	fprintf(stderr, "         -log STR      log filename [%s]\n", LogFileName);
	fprintf(stderr, "         -monomorphic  report all loci which do not have any potential alternates.\n");
	fprintf(stderr, "         -min_cnv INT  the minimal cnv size to be reported [%d].\n", MinCNVsize);
	fprintf(stderr, "         -min_gap INT  the minimal gap(unmapped) size to be reported [%d].\n", MinUnmappedSize);
	fprintf(stderr, "         -ploidy INT   number of sets of chromosomes in a cell (1:monoploid, 2:diploid) [%d]\n", iPloidy);
	fprintf(stderr, "         -m            output multiple alignments\n");
	fprintf(stderr, "         -somatic      detect somatic mutations [false]\n");
	fprintf(stderr, "         -no_vcf       No VCF output [false]\n");
	fprintf(stderr, "         -p            paired-end reads are interlaced in the same file\n");
	fprintf(stderr, "         -filter       apply variant filters (under test) [false]\n");
	fprintf(stderr, "         -id STR       assign sample id\n");
	fprintf(stderr, "         -v            version\n");
	fprintf(stderr, "\n");
}

string MakeRefIdx(char* RefFileName)
{
	string str;
	srand((unsigned int)time(NULL)); str.resize(10);
	for (int i = 0; i < 10; i++) str[i] = (unsigned short)(rand() >> 3) % 26 + 97;

	return str;
}

bool CheckOutputFileName(char *FileName)
{
	struct stat s;
	bool bRet = true;
	int i, len = strlen(FileName);

	if (strcmp(FileName, "-") != 0)
	{
		for (i = 0; i < len; i++)
		{
			if (isalnum(FileName[i]) || FileName[i] == '/' || FileName[i] == '.' || FileName[i] == '_' || FileName[i] == '-');
			else
			{
				bRet = false;
				fprintf(stderr, "Warning: [%s] is not a valid filename!\n", FileName);
				break;
			}
		}
		if (stat(FileName, &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				bRet = false;
				fprintf(stderr, "Warning: %s is a directory!\n", FileName);
			}
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
	string parameter, str, random_prefix;

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

	//MinIndFreq = 5;
	MaxClipSize = 5;
	MinCNVsize = 50;
	MaxPosDiff = 30;
	MinReadDepth = 20;
	iMaxDuplicate = 5;
	FragmentSize = 500;
	MinAlleleDepth = 5;
	FrequencyThr = 0.2;
	MinVarConfScore = 10;
	MinUnmappedSize = 50;
	MaxMisMatchRate = 0.05;
	sample_id = (char*)"unknown";
	LogFileName = (char*)"job.log";
	VcfFileName = (char*)"output.vcf";
	ObservGenomicPos = ObserveBegPos = ObserveEndPos = -1;
	RefSequence = RefFileName = IndexFileName = SamFileName = KnownSiteFileName = NULL;

	if (argc == 1 || strcmp(argv[1], "-h") == 0) ShowProgramUsage(argv[0]);
	else if (strcmp(argv[1], "update") == 0)
	{
		system("git fetch;git merge origin/master master;make");
		exit(0);
	}
	else if (strcmp(argv[1], "index") == 0)
	{
		if(argc == 4) bwa_idx_build(argv[2], argv[3]);
		else
		{
			fprintf(stderr, "usage: %s index ref.fa prefix\n", argv[0]);
		}
		exit(0);
	}
	else
	{
		for (CmdLine = argv[0], i = 1; i < argc; i++) CmdLine += " " + (string)argv[i];

		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i" && i + 1 < argc) IndexFileName = argv[++i];
			else if (parameter == "-r" && i + 1 < argc) RefFileName = argv[++i];
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
				if (atoi(argv[++i]) <= 15) iMaxDuplicate = (int8_t)atoi(argv[i]);
				else fprintf(stderr, "Warning! The PCR-duplicate range is [1-15]!\n");
			}
			else if (parameter == "-filter") bFilter = true;
			else if (parameter == "-id" && i + 1 < argc) sample_id = argv[++i];
			else if (parameter == "-label" && i + 1 < argc) sample_id = argv[++i];
			else if (parameter == "-size" && i + 1 < argc) FragmentSize = atoi(argv[++i]);
			else if (parameter == "-indel" && i + 1 < argc)
			{
				if ((MaxPosDiff = atoi(argv[++i])) > 100)
				{
					MaxPosDiff = 100;
					fprintf(stderr, "Warning! The maximal indel size is 100!\n");
				}
			}
			else if (parameter == "-min_cnv" && i + 1 < argc) MinCNVsize = atoi(argv[++i]);
			else if (parameter == "-min_gap" && i + 1 < argc) MinUnmappedSize = atoi(argv[++i]);
			//else if (parameter == "-dp" && i + 1 < argc) MinReadDepth = atoi(argv[++i]);
			else if (parameter == "-ad" && i + 1 < argc) MinAlleleDepth = atoi(argv[++i]);
			//else if (parameter == "-ind" && i + 1 < argc) MinIndFreq = atoi(argv[++i]);
			else if (parameter == "-ploidy" && i + 1 < argc)
			{
				if ((iPloidy = atoi(argv[++i])) > 2)
				{
					iPloidy = 2;
					fprintf(stderr, "Warning! MapCaller only supports monoploid and diploid!\n");
				}
			}
			//else if (parameter == "-known_sites" && i + 1 < argc)
			//{
			//	KnownSiteFileName = argv[++i];
			//}
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
			else if (parameter == "-log" && i + 1 < argc) LogFileName = argv[++i];
			else if ((parameter == "-alg") && i + 1 < argc)
			{
				str = argv[++i];
				if (str == "ksw2") NW_ALG = false;
				else NW_ALG = true; //nw
			}
			else if (parameter == "-maxmm" && i + 1 < argc) MaxMisMatchRate = atof(argv[++i]);
			else if (parameter == "-maxclip" && i + 1 < argc) MaxClipSize = atoi(argv[++i]);
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
		if (iMaxDuplicate <= 0 || iMaxDuplicate > 15) iMaxDuplicate = 15;

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

		if (strcmp(LogFileName, "job.log")!= 0 && CheckOutputFileName(LogFileName) == false) exit(0);
		if (SamFileName != NULL && CheckOutputFileName(SamFileName) == false) exit(0);
		if (VcfFileName != NULL && CheckOutputFileName(VcfFileName) == false) exit(0);

		if (RefFileName != NULL)
		{
			random_prefix = MakeRefIdx(RefFileName);
			IndexFileName = (char*)random_prefix.c_str();
			bwa_idx_build(RefFileName, IndexFileName);
		}
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
			//if (KnownSiteFileName != NULL) LoadKnownSites(KnownSiteFileName);

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
			pthread_mutex_init(&VarLock, NULL); pthread_mutex_init(&OutputLock, NULL); pthread_mutex_init(&LibraryLock, NULL); pthread_mutex_init(&ProfileLock, NULL);

			StartProcessTime = time(NULL);
			FILE *log = fopen(LogFileName, "a"); fprintf(log, "%s\n[CMD]", string().assign(80, '*').c_str()); for (i = 0; i < argc; i++) fprintf(log, " %s", argv[i]); fprintf(log, "\n\n"); fclose(log);

			Mapping();
			if (bVCFoutput) VariantCalling();

			bwa_idx_destroy(RefIdx);
			if (RefSequence != NULL) delete[] RefSequence;
			if (MappingRecordArr != NULL) delete[] MappingRecordArr;
			if (RefFileName != NULL)
			{
				random_prefix = "rm -f " + random_prefix + "*";
				system(random_prefix.c_str());
			}
			log = fopen(LogFileName, "a"); 
			fprintf(log, "All done! It took %lld seconds to complete the data analysis.\n\n\n", (long long)(time(NULL) - StartProcessTime)); fclose(log);
			fprintf(stderr, "All done! It took %lld seconds to complete the data analysis.\n", (long long)(time(NULL)- StartProcessTime));
		}
	}
	return 0;
}
