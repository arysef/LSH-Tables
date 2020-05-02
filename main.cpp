#include "LSH.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string>
#include <map>
#include <algorithm>
using namespace std;
KSEQ_INIT(int, read)

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
int main() {

    LSH *lsh = new LSH();


    unsigned int i[] = {1, 1, 1, 2, 2, 2, 3, 3, 3};

    unsigned int h[] = {0, 1, 2, 3, 0, 1, 2, 3, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2,
                        1, 0, 1, 3, 2, 2, 1, 2, 3, 3, 3, 3, 2, 1, 1, 1, 2, 3};

    lsh->insert(9, i, h);

    //lsh->view();

    unsigned int q[] = {0, 1, 2, 3, 1, 2, 0, 1};

    unsigned int r[10];

    lsh->top_k(2, 5, q, r);

    
    // for (int i; i < 10; i++) {
    //     std::cout << r[i] << "\n";
    // }
    
    int test_edit;

    //test_edit = leven_distance("oluwapelumi","back");
    //printf("should say 10: %d \n", test_edit); 
    FILE* fp;
    kseq_t *seq;
    int n = 0, slen = 0, qlen = 0;
    fp = fopen("SRR000001.fastq", "r");
    seq = kseq_init(fileno(fp));
    int count = 0;
    std::string seq0;
    std::string seq1;
    
    //I can use a vector to dynamically allocate space for an array
    // will store all the dna strings 
    vector<std::string> dna_arr;
    // Reads FASTQ file sequence by sequence
    while (kseq_read(seq) >= 0) {
        std::string str_seq = std::string(seq->seq.s);

        if (count < 20) {
        dna_arr.push_back(str_seq);
        }

        //if statements to use to check that the distance metric works for the actual dna data
        if (count == 0) {
            seq0 = str_seq;
        } 
        if (count == 1) {
            seq1 = str_seq;
            //doing the data comparison
            //printf("seq0: %s\n", seq0.c_str());
            //printf("seq1: %s\n", seq1.c_str());
            //std::cout << "dna_seq distance: " <<leven_distance(seq0, seq1) << "\n";
        }
        //Trying to get basic minHash of sequence 
        int min = INT_MAX;
        int sub_len = 5; // n-gram, so if sublen is 3 we would be MinHashing based on trigram
        for(int i = 0; i < str_seq.length()-sub_len; i=i+sub_len) {
            unsigned int cur = std::hash<std::string>{}(str_seq.substr(i, i+sub_len));
            if (cur < min)
                min = cur;
        }
        //This is just sanity check: we would expect MinHash of sequence to in most cases be smaller than hash of sequence itself
        //std::cout << std::hash<std::string>{}(str_seq) << " " << min << "\n";
        count += 1;
    }
    //printf("%d\t%d\t%d\n", n, slen, qlen);
    brute_topk(5, seq1, dna_arr);
    kseq_destroy(seq);
    fclose(fp);
    return 0;
    
}

//comparison operator that will put the smaller things first
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
    //okay I'm going to change this all up and create a bunch of structs with the dna string information along with the distance information
    //I'll put the distance and dna info in a struct, and put them all in a vector. Then I'll sort the vector by the distance from the query and take the first k elements 


    vector<comparison_info> dna_info_vec;
    int cnt = 0; // will be used to tell me if we've inserted the k things in the vector already
    //iterating over the dna sequences
    for (auto p_dna_seq = dna_strings.begin(); p_dna_seq != dna_strings.end(); ++p_dna_seq) {
        int dist_val = leven_distance(query, *p_dna_seq);
        cout << "sequence " << cnt << *p_dna_seq << "and length: " << dist_val << "\n";
        if (query.compare(*p_dna_seq) != 0) { //if the strings aren't equal
        //new struct
        auto dna_struct = new comparison_info; //returns a pointer to a comparison_info struct
        //distance of query and dna_seq
        //int dist_val = leven_distance(query, *p_dna_seq);
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
        std::cout << cnt2 << ": " << it->dna_sequence << " dist: " << it->dist_from_query << "\n";
        cnt2 += 1;
    }

    return (topk_vec);
}

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
    //std::cout << "me: " << m << "\n";
    //Initializing my dp matrix
   
   // printf("m,n: %d,%d\n", m, n);
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

//Thinking about all the comparisons we'll have to do, we may need to store the information in an array 