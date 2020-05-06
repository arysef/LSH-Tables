#include "LSH.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string>
#include <map>
#include <algorithm>
#include "MurMurHash3.h"
#include "config.h"
#include <chrono>
using namespace std;

KSEQ_INIT(int, read)
int len = 32;

struct comparison_info
{
    /* data */
    int dist_from_query;
    string dna_sequence;
};

static int min3(int x, int y, int z);
static int leven_distance(std::string str1, std::string str2);
//uhh before I do this, I need to figure out how to store a sequence of strings
static vector<string> brute_topk(int k, std::string query, vector<string> dna_strings);
static uint32_t getSequenceMinHash(string sequence, uint32_t seed, uint32_t subseq_len);
static unsigned int *getDNAMinhashes(string dna_sequence, int numHashes, unsigned int * dHashes);
static void insert_into_lsh(string dna_sequence, uint32_t string_idx, int numHashes, LSH *lsh_table);
static vector<uint32_t> getVectorMinhashes(string dna_sequence, uint32_t seed, int numHashes);

int main() {
    int test_edit;

    FILE* fp;
    kseq_t *seq;
    int n = 0, slen = 0, qlen = 0;
    fp = fopen("SRR000001.fastq", "r");
    seq = kseq_init(fileno(fp));
    int count = 0;
    std::string seq0;
    uint32_t idx = 0;
    map<uint32_t, string> idx_to_string_map;
   
    //Intializing LSH table
    LSH *lsh2 = new LSH();
    cout << "Done allocating memory for LSH table" << "\n";
    //vector to store all the dna sequences to be processed 
    vector<std::string> dna_arr;
    unsigned int seq_hashes[4];
    unsigned int seq_top10[10];
    // Reads FASTQ file sequence by sequence
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    while (kseq_read(seq) >= 0) {
        
        string str_seq = string(seq->seq.s);
        //mapping an index to the read string and updating the index -- needed for Leven distance
        idx_to_string_map.insert({idx, str_seq});
        idx += 1;

        // heartbeat for insertion
        if (count % 10000 == 0){
            uint32_t min = getSequenceMinHash(str_seq, 0, 10);
            uint32_t full_hash;
            int len = 32;
            MurmurHash3_x86_32(&str_seq, len, 0, &full_hash);
            cout << count << " " << full_hash << " " << min << "\n";
        }
        insert_into_lsh(str_seq, count, 4, lsh2);
        count += 1;
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // below are the test values that we used 
    unsigned int idxs[10] =  {308597, 319039, 464183, 394249, 439143, 442117, 18615, 153262, 260764, 82983};
    cout << "Experiment begins" << "\n";
    // Retrieving top 10 for each test value 
    for(int i = 0; i < 10; i++) {
        uint32_t idx = idxs[i];
        cout << "test";
        string seq = idx_to_string_map[idx];
        seq_hashes[0] = getSequenceMinHash(seq, 0, 6);
        seq_hashes[1] = getSequenceMinHash(seq, 1, 6);
        seq_hashes[2] = getSequenceMinHash(seq, 2, 6);
        seq_hashes[3] = getSequenceMinHash(seq, 3, 6);
        lsh2->top_k(1, 10, seq_hashes, seq_top10);
        
        // Uncomment below to print top k for each sequence in test values
        //cout << "Top 10 for Index: " << idx << "\n";
        for (int j = 0; j < 10; j++) {
            //cout << seq_top10[j] << "\n";
        }
    }

   std::cout << "Time to insert = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    kseq_destroy(seq);
    fclose(fp);
    return 0;

}

/*
 * Requires: 
 * dna_sequence is a dna sequence in string form
 * string_idx is the corresponding integer index of the string
 * lsh_table is the lsh data structure where we'll be inserting the string
 * 
 * Effect: 
 * Gets the hashes of the string and inserts the index of the string with the corresponding hashes into the lsh table   
 **/
static void insert_into_lsh(string dna_sequence, uint32_t string_idx, int numHashes, LSH *lsh_table) {
   unsigned int sequence_hashes[numHashes];
    //get the hashes of the string
    getDNAMinhashes(dna_sequence, numHashes, sequence_hashes);
    // insert the index with the hashes into the lsh table 
    lsh_table->insert((unsigned int)string_idx, sequence_hashes);
}

/*
 * Requires: 
 * dna sequence to be hashed
 * number of times to hash the string
 * integer which is the length of ngram used for the MinHashes
 * 
 * Effect: 
 * returns vector with num_hashes MinHashes for the sequence
 * this vector will have length of len(num_hashes)   
 **/
static unsigned int *getDNAMinhashes(string dna_sequence, int numHashes, unsigned int * dHashes) {
    for (int i = 0; i < numHashes; i++)
        dHashes[i] = getSequenceMinHash(dna_sequence, i, NGRAM_LEN);
    return  dHashes;
}

/*
 * Requires: 
 * vector with sequences in it
 * integer which represents number of hashes per sequence ( the total amount needed to be hashed into the L tables)
 * integer which is the length of ngram used for the MinHashes
 * 
 * Effect: 
 * returns vector with num_hashes MinHashes for each sequence that was input in sequences 
 * this vector will have length of len(sequences)*num_hashes   
 **/
static vector<uint32_t> getMinHashes(vector<string> sequences, int num_hashes) {
    vector<uint32_t> seqMinHashes;
    for (string seq: sequences)  {
        for (int i = 0; i < num_hashes; i++)
            seqMinHashes.push_back(getSequenceMinHash(seq, i, NGRAM_LEN));
    }
    return seqMinHashes;
}

/*
 * Requires: 
 * sequence is string dna sequence to take MinHash of
 * seed is seed to be used for hashing
 * subseq_len is the length of the n-gram used to take MinHash
 * 
 * Effect: 
 * Returns 32bit unsigned int MinHash of sequence
 **/
static uint32_t getSequenceMinHash(string sequence, uint32_t seed, uint32_t subseq_len) {
    uint32_t min = UINT_MAX;
    uint32_t subseq_hash;
    // Iterates through all ngrams in sequence and returns lowest hash value for given seed
    for(int i = 0; i < sequence.length()-subseq_len; i++) {
        string temp = sequence.substr(i, i+subseq_len);
        uint32_t cur;
        MurmurHash3_x86_32(&temp, len, seed, &cur);
        if (cur < min)
            min = cur;
    }
    return min % (1 << RANGE_POW);
} 



// simple comparison operator that will put the smaller things first
bool compareByDist(const comparison_info  &a, const comparison_info &b)
{
    return a.dist_from_query < b.dist_from_query;
}

/*
* Requires:
*   k is the number of strings to return
*   query is the string in question in which we want to find the k strings that are closest to it.
*   dna_strings is a vector of strings that we will compare to the query
*
* Effect:
*   Finds the top k similar dna strings to the query. 
*   if multiple dna sequences have the same distance from the query, we only take one of them. Need to ammend this tomorrow.
*   Returns of vector of those topk dna sequences
*/
static vector<string> brute_topk(int k, std::string query, vector<string> dna_strings) 
{

    vector<comparison_info> dna_info_vec;
    int cnt = 0; // will be used to tell me if we've inserted the k things in the vector already
    //iterating over the dna sequences
    for (auto p_dna_seq = dna_strings.begin(); p_dna_seq != dna_strings.end(); ++p_dna_seq) {
        int dist_val = leven_distance(query, *p_dna_seq);
        if (query.compare(*p_dna_seq) != 0) { //if the strings aren't equal
        //new struct
        auto dna_struct = new comparison_info; //returns a pointer to a comparison_info struct
        //distance of query and dna_seq
        dna_struct->dist_from_query = dist_val;
        dna_struct->dna_sequence = *p_dna_seq;
        dna_info_vec.push_back(*dna_struct);
       }
    }
    // now sorting the information in place
    std::sort(dna_info_vec.begin(), dna_info_vec.end(), compareByDist);

    vector<string>topk_vec;
    int cnt2 = 0;
    //putting only the top k dna sequences 
    for(auto it = dna_info_vec.begin(); cnt2 < k; ++it) {
        topk_vec.push_back(it->dna_sequence);
        cnt2 += 1;
    }

    return (topk_vec);
}

//Find the minimum value of 3 integers
int 
min3(int x, int y, int z) 
{ 
    return std::min(std::min(x, y), z); 
} 

/*
* Requires:
*   str1 and str2 are both genomic sequences represented as strings
*   m is the index of the character we're looking at in str1
*   n is the index of the character we're looking at in str2
*
* Effect:
*   Computes the Levenshtein distance between the two strings(the minimum number of changes to turn str1 into str2)
*   returns the levenshtein distance as an integer
*/

int leven_distance(std::string str1, std::string str2)
{
    //dynamic programming with storing
    int m;
    int n;
    m = (int)strlen(str1.c_str());
    n = (int)strlen(str2.c_str());
    //Initializing my dp matrix
   
    int dp[m + 1][n + 1]; //an extra box for when looking at the empty string prefix of either str1 or str2 or both

    //filling matrix
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++){
            //if i =0, first string is empty, we need to insert all the characters of string 2 to string 1
            if (i == 0) {
               // printf("i is 0 \n");
                dp[i][j] = j;
                continue;
            }
            //if j =0, first string is empty, we need to insert all the characters of string 2 to string 1
            if (j == 0) {
                dp[i][j] = i;
                continue;
            }
            //if the previous character of str1 and str2 are equal, 
            if (str1[i-1] == str2[j-1]) {//remember that the dp[i][j] represents the distance between str1 up to index i-1 and str2 up to index j-1
                dp[i][j] = dp[i-1][j-1];
                continue;
            }
            // if the previous characters are different, consider all the possibilities
            dp[i][j] = 1 + min3(dp[i-1][j-1], //replace
                                    dp[i-1][j], //remove
                                    dp[i][j-1]); //insert
        }
    }
    return (dp[m][n]);
}
