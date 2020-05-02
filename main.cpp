#include "LSH.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string>
KSEQ_INIT(int, read)


static int min3(int x, int y, int z);
static int leven_distance(std::string str1, std::string str2);

int main() {

    LSH *lsh = new LSH();


    unsigned int i[] = {1, 1, 1, 2, 2, 2, 3, 3, 3};

    unsigned int h[] = {0, 1, 2, 3, 0, 1, 2, 3, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2,
                        1, 0, 1, 3, 2, 2, 1, 2, 3, 3, 3, 3, 2, 1, 1, 1, 2, 3};

    lsh->insert(9, i, h);

    lsh->view();

    unsigned int q[] = {0, 1, 2, 3, 1, 2, 0, 1};

    unsigned int r[10];

    lsh->top_k(2, 5, q, r);

    for (int i; i < 10; i++) {
        std::cout << r[i] << "\n";
    }
    
    int test_edit;

    test_edit = leven_distance("oluwapelumi","back");
    printf("should say 10: %d \n", test_edit); 
    FILE* fp;
    kseq_t *seq;
    int n = 0, slen = 0, qlen = 0;
    fp = fopen("SRR000001.fastq", "r");
    seq = kseq_init(fileno(fp));
    // Reads FASTQ file sequence by sequence
    while (kseq_read(seq) >= 0){
        std::string str_seq = std::string(seq->seq.s);

        //Trying to get basic minHash of sequence 
        int min = INT_MAX;
        int sub_len = 10; // n-gram, so if sublen is 3 we would be MinHashing based on trigram
        for(int i = 0; i < str_seq.length()-sub_len; i=i+sub_len) {
            unsigned int cur = std::hash<std::string>{}(str_seq.substr(i, i+sub_len));
            if (cur < min)
                min = cur;
        }
        //This is just sanity check: we would expect MinHash of sequence to in most cases be smaller than hash of sequence itself
        std::cout << std::hash<std::string>{}(str_seq) << " " << min << "\n";
    }
    printf("%d\t%d\t%d\n", n, slen, qlen);
    kseq_destroy(seq);
    fclose(fp);
    return 0;
    
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
    std::cout << "me: " << m << "\n";
    //Initializing my dp matrix
   
    printf("m,n: %d,%d\n", m, n);
    int dp[m + 1][n + 1]; //an extra box for when looking at the empty string prefix of either str1 or str2 or both

    //filling matrix
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++){
            //if i =0, first string is empty, we need to insert all the characters of string 2 to string 1
            if (i == 0) {
                printf("i is 0 \n");
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