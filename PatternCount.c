/* 
 * File:   main.c
 * Author: kevincando
 *
 * Created on April 20, 2015, 10:21 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
 * 
 */
/*
int main(int argc, char** argv) {

    char str1[] = "CGTTTCATCATCCGCTGTCCGGTCCAGTCCGGTTTAGGTCCGGTAGTCCGGTGTTGTGTCCGGTAACGAGAGACCTCGTCCGGTCACCTCTGTCCGGTAAAGTCCGGTAAGGTCCGGTCTGTCCGGTGGGTTTGACCGGGCGTCCGGTCCTATCGTCCGGTAGGTCCGGTCGTCCGGTTGTCCGGTCCGGGGTCCGGTGCTGTCCGGTGTCCGGTAGTCCGGTGTCCGGTCGCCACCTGTGCGTCCGGTCGTCCGGTGGAAGCGTCCGGTCCTCCCGGTCCGGTTGTCCGGTAGCGTCCGGTAAGTCCGGTAGTCCGGTTTGTCCGGTCGTCCGGTAAGGTCCGGTTCATCGTCCGGTGGTCCGGTTGTCCGGTTGTCCGGTCTTTTGAAGTCCGGTTGTCCGGTGTCCGGTTGGTCCGGTCTGGTCCGGTAAAGTCCGGTTTTGTCCGGTGGTCCGGTGGTCCGGTAGTCCGGTGTCCGGTTAAGTCCGGTGTCCGGTATTGTCCGGTGTCCGGTTCAGGTCCGGTCCCCGTCCGGTAAGTCCGGTTGGTTGGTCCGGTCACATAGGATACAGTCCGGTGTAGTCCGGTCTAAGGTCCGGTGAAATTCAAGGTCCGGTTCGTGTCCGGTGTCCGGTGGATTGTCCGGTTGTCCGGTAGTCCGGTACGTCCGGTCCGTCCGGTTGTCCGGTGGGTCCGGTGATCAAGTCCGGTGCGTCCGGTTGTAGTCCGGTCGTCCGGTGAGTTAGGATGTCCGGTTGTCCGGTGTCCGGTGTCCGGTCCGTCCGGTCGTCCGGTGTCCGGTTAAACCGGTCCGGTAAGTCCGGTGTCCGGTAGGTCCGGTTCGTCCGGTTCATATCGTCCGGTTAGTCCGGTAAGGGTCCGGTGTCCGGTCATTCTTCCAGGGTGTCCGGTTGTCCGGTCTTCTGTCCGGTAGTGTCCGGT";
    char str2[] = "GTCCGGTGT";
    
    printf("%d",PatternCount(str1,str2));
    
    return (EXIT_SUCCESS);
}
*/
int PatternCount(char text[], char pattern[]){
    unsigned int count = 0;
    unsigned int i,j,z,flag = 0;
    
    for(i = 0 ; i < strlen(text) ; i++){
        
        if( text[i] == pattern[0] ){
            j = 1;
            z = i + 1;
            if( z + strlen(pattern) < strlen(text) )
                do{
                    if(text[z] == pattern[j])
                        flag = 1; 
                    else{ 
                        flag = 0;
                        break;
                    }
                    j++;
                    z++;
                }while( j < strlen(pattern));
        }
        if(flag){
            count++;
            flag = 0;
        }
    }
    
    return count;
}

