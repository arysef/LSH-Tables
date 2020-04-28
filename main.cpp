#include "LSH.h"
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(int, read)


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
        int sub_len = 10;
        for(int i = 0; i < str_seq.length()-sub_len; i=i+sub_len) {
            //std::cout << str_seq.substr(i, 1+sub_len);
            unsigned int cur = std::hash<std::string>{}(str_seq.substr(i, i+sub_len));
            if (cur < min)
                min = cur;
        }
        //std::cout << "\n Above is put together \n" << str_seq << "\n";
        std::cout << std::hash<std::string>{}(str_seq) << " " << min << "\n";
    }
    printf("%d\t%d\t%d\n", n, slen, qlen);
    kseq_destroy(seq);
    fclose(fp);
    return 0;
    
}