#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm>
#include <ctime>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <inttypes.h>

#define KmerSize 8
#define KmerPower 0x3FFF

#define MaxPosDiff 30
#define MinSeedLength 16
#define ReadChunkSize 4000

using namespace std;

typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	bwtint_t x[3];
} bwtintv_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
} bwtSearchResult_t;

typedef struct
{
	uint32_t wid; // word id
	uint32_t pos; // occurrence position
} KmerItem_t;

typedef struct
{
	int PosDiff;
	uint32_t rPos;
	uint32_t gPos;
} KmerPair_t;

typedef struct
{
	int len; // chromosome length
	char* name; // chromosome name
	int64_t FowardLocation;
	int64_t ReverseLocation;
} Chromosome_t;

typedef struct
{
	int64_t gPos;
	int ChromosomeIdx;
} Coordinate_t;

typedef struct
{
	bool bSimple;
	int rPos; // read position
	int64_t gPos; // genome position
	int rLen; // read block size
	int gLen; // genome block size
	int64_t PosDiff; // gPos-rPos
	string aln1; // read fragment alignment
	string aln2; // genomic fragment alignment
} FragPair_t;

typedef struct
{
	int score;
	int SamFlag;
	string CIGAR;
	bool orientation;
	int PairedAlnCanIdx;
	vector<FragPair_t> FragPairVec;
} AlnCan_t;

typedef struct
{
	int BestAlnCanIdx;
	int score;
	int sub_score;
} AlnSummary_t;

typedef struct
{
	int rlen;
	char* seq;
	char* qual;
	char* header;
	uint8_t* EncodeSeq;
	AlnSummary_t AlnSummary;
	vector<AlnCan_t> AlnCanVec;
} ReadItem_t;

typedef struct
{
	uint16_t A;
	uint16_t C;
	uint16_t G;
	uint16_t T;
	uint32_t multi_hit;
} MappingRecord_t;

typedef struct
{
	int idx1;
	int idx2;
	int p_score;
} PairedReads_t;

typedef struct
{
	int64_t dist;
	int64_t gPos1;
	int64_t gPos2;
} CoordinatePair_t;

typedef struct
{
	int64_t gPos;
	int64_t dist;
} DiscordPair_t;

typedef struct
{
	uint8_t type;
	int64_t gPos;
	uint8_t qscore;
	uint32_t gEnd;
	uint16_t NS; // Number of samples with data
	uint16_t DP; // read depth
} VarPos_t;

// Global variables
extern bwt_t *Refbwt;
extern bwaidx_t *RefIdx;
extern uint32_t avgDist;
extern const char* VersionStr;
extern vector<string> ReadVec;
extern time_t StartProcessTime;
extern map<int64_t, int> ChrLocMap;
extern unsigned char nst_nt4_table[256];
extern MappingRecord_t* MappingRecordArr;
extern vector<Chromosome_t> ChromosomeVec;
extern vector<CoordinatePair_t> DistantPairVec;
extern vector<string> ReadFileNameVec1, ReadFileNameVec2;
extern int64_t GenomeSize, TwoGenomeSize, ObservGenomicPos;
extern char *RefSequence, *IndexFileName, *SamFileName, *VcfFileName;
extern bool bDebugMode, bPairEnd, bUnique, gzCompressed, FastQFormat, bSAMoutput, bVCFoutput;
extern int iThreadNum, iChromsomeNum, WholeChromosomeNum, ChromosomeNumMinusOne, FragmentSize, MinBaseDepth, MinVarConfScore, ObserveBegPos, ObserveEndPos;

extern vector<DiscordPair_t> InversionSiteVec, TranslocationSiteVec;
extern map<int64_t, map<string, uint16_t> > InsertSeqMap, DeleteSeqMap;

// GetData.cpp
extern bool CheckReadFormat(const char* filename);
extern bool CheckBWAIndexFiles(string IndexPrefix);
extern bool CheckReadFile(char* filename, bool& bReadFormat);
extern int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr);
extern int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr);

// VariantCalling.cpp
extern void VariantCalling();

// ReadMapping.cpp
extern void Mapping();
extern void ShowFragPairCluster(vector<AlnCan_t>& AlnCanVec);
extern bool CompByPosDiff(const FragPair_t& p1, const FragPair_t& p2);
extern vector<AlnCan_t> SimplePairClustering(int rlen, vector<FragPair_t>& SimplePairVec);

// ReadAlignment.cpp
extern bool ProduceReadAlignment(ReadItem_t& read);

// AlignmentRescue.cpp
extern int AlignmentRescue(uint32_t EstDist, ReadItem_t& read1, ReadItem_t& read2);

// AlignmentProfile.cpp
extern void UpdateProfile(ReadItem_t* read, vector<AlnCan_t>& AlnCanVec);
extern void UpdateMultiHitCount(ReadItem_t* read, vector<AlnCan_t>& AlnCanVec);

// SamReport.cpp
extern void GenerateSingleSamStream(ReadItem_t& read, vector<string>& SamStreamVec);
extern void GeneratePairedSamStream(ReadItem_t& read1, ReadItem_t& read2, vector<string>& SamStreamVec);

// tools.cpp
extern void ShowProfileColumn(int64_t gPos);
extern void ShowSeedLocationInfo(int64_t MyPos);
extern int64_t GetAlignmentBoundary(int64_t gPos);
extern bool CheckFragValidity(FragPair_t FragPair);
extern void SelfComplementarySeq(int len, char* rseq);
extern Coordinate_t DetermineCoordinate(int64_t gPos);
extern int GetProfileColumnSize(MappingRecord_t& Profile);
extern void ShowFragmentPair(char* ReadSeq, FragPair_t& fp);
extern void ShowSimplePairInfo(vector<FragPair_t>& FragPairVec);
extern void GetComplementarySeq(int len, char* seq, char* rseq);
extern bool CheckAlignmentValidity(vector<FragPair_t>& FragPairVec);
extern void ShowVariationProfile(int64_t begin_pos, int64_t end_pos);

// bwt_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);

// bwt_search.cpp
extern void BWT_Check(uint8_t* seq, int start, int stop);
extern bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop);
extern vector<FragPair_t> SimplePairRescue(int64_t LowerBound, int64_t UpperBound, ReadItem_t& read);

// KmerAnalysis.cpp
extern vector<KmerItem_t> CreateKmerVecFromReadSeq(int len, char* seq);
extern vector<KmerPair_t> IdentifyCommonKmers(int MaxShift, vector<KmerItem_t>& vec1, vector<KmerItem_t>& vec2);
extern vector<FragPair_t> GenerateSimplePairsFromCommonKmers(int thr, int64_t gPos, vector<KmerPair_t>& KmerPairVec);
//extern vector<SeedPair_t> GenerateSimplePairsFromFragmentPair(int MaxDist, int len1, char* frag1, int len2, char* frag2);

// nw_alignment.cpp
extern void nw_alignment(int m, string& s1, int n, string& s2);
