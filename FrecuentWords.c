/* 
 * File:   FrecuentWords.c
 * Author: kevincando
 *
 * Created on April 20, 2015, 11:53 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "list.h"
#include "queue.h"

/*
 * 
 */
char **frecuentWords(char text[], int k);
char *ReverseComplementProblem(const char DNAString[]);
int *PatternMatching(const char Pattern[], const char Genome[]);
unsigned short int SymbolToNumber(const char symbol);
unsigned long int PatternToNumber(char Pattern[]);
char *NumberToSymbol(int index);
char *NumberToPattern(int index,int k);
char *strBuilder(const char *c1,const char *c2);
long int *ComputingFrequencies(const char DNAStr[] ,int k);
char **ClumpFindingProblem(const char Genome[], unsigned int k, unsigned int L, unsigned int t);
char **BetterClumpFindingProblem(const char Genome[], unsigned int k, unsigned int L, unsigned int t);
int *Skew(const char Genome[]);
int *MinimumSkew(const char Genome[]);
int *MaximumSkew(const char Genome[]);
int HammingDistanceProblem(const char String[],const char String2[]);
int *ApproximatePatternMatching(const char Pattern[],const char Genome[],int d);
int ApproximatePatternCount(const char Text[], const char Pattern[],int d);
char **FrequentWordswithMismatchesProblem(const char Text[], int k, int d);
//static char **Neighbors(char Pattern[],int d);
char *Suffix(char *Pattern);
char *FirstSymbol(const char Pattern[]);
List *Neighbors(char Pattern[],int d);
List *FrequentWordsWithMismatches(const char Text[], int k,int d);
List *FrequentWordsWithMismatchesAndReverseComplements(const char Text[], int k,int d);
List *MotifEnumeration(List *DNA, int k ,int d);
Generic StringReadFromFile(FILE *pf);
void listRemoveDuplicates(List *L, cmpfn fn);
char *strCreateSelect(const char Text[], int posIni, int length);
int listElemExist(List *L, Generic cont, cmpfn fn);
List *listCat(List *L1,List *L2);
int DistanceBetweenPatternAndStrings(const char Pattern[], List *DNA);
char *MedianString(List *DNA, int k);
float **readMatrixFromFile(const char file[],int nrow,int ncol);
char *ProfileMostProbableKmer(const char Text[], int k, float **Profile);
List *GreedyMotifSearch(List *DNA, int k, int t);
float **createProfileMatrix(List *motifs);
int Count(const char Text[], char symbol);
int Score(List *motifs);
char *findConsensus(List *motifs);

//Segundo modulo

Queue *Composition(const char Text[], int k);
int strAlfCmp(const char str1[], const char str2[]);

int main(int argc, char** argv) {

    char str[] = "CAGTGA";
    char dna[] = "CATTCCAGTACTTCGATGATGGCGTGAAGA";
    char pattern[] = "GGCATGACTC";
    char genome[] = "CGTCCGTCGTCGTCTTACTTACGTCGTATACATACCTTACGTATACATACATACCTTACGTATACATACCGTCACCCCGTCCGTCACCCACCCACCCGTCTTACTTACGTCACCCACCCGTCGTCCGTCTTACGTCGTCCGTCACCCACCATACCGTCGTCGTCTTACGTATACCTTAATACCCGTCACCCTTACCGTCTTACTTACACCCACCCGTATACATACATACCCGTCGTCCGTCGTCACCCACCCCGTCACCCGTCACCCTTACCGTCTTACGTCCGTCGTCACCCCGTCCGTCACCCACCCTTACCGTCTTACGTCACCCTTACCGTCACCCCGTCCGTATACCTTACTTAATACCCGTATACCTTACTTAATACCGTCGTCGT";
    char String[10000] = "CAATCCAAC", 
         String2[10000] = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG";
    int k = 5, *matching, L = 500, t = 5 , d = 2;
    int index = 5437,cont = 0;
    int k2 = 0;
    char **frec, *revDna ;
    char *pattern2;
    long int num = pow(4,k), *frecArray;
    float  probMat[4][5]={{0.2,0.2,0.3,0.2,0.3},
                          {0.4,0.3,0.1,0.5,0.1},
                          {0.3,0.3,0.5,0.2,0.4},
                          {0.1,0.2,0.1,0.1,0.2}};
    float **prob2;
    List *L1,*L2;
    FILE *pf = fopen("datos.txt","w");
    FILE *pf2 = fopen("datosing.txt","r");
    FILE *pf3 = fopen("Answers.txt","r");
    //printf("%d\n",PatternCount(str,pattern));
    
    //frec = FrequentWordswithMismatchesProblem(genome,k,d);
    /**/
    
    //L1 = listReadFile("datosing.txt",StringReadFromFile);
    //prob2 = readMatrixFromFile("datosing.txt",4,k);
    //for(int i = 0 ; i < 4; i++){
    //    for(int j = 0; j < k ; j++)
      //      printf("%f ",prob2[i][j]);
    //    printf("\n");
   // }
    //fprintf(pf,"%s", ProfileMostProbableKmer(str,k,prob2));
    //fprintf(pf,"%s", MedianString(L1,k));
    
    //L2 = GreedyMotifSearch(L1,k,t);
    fscanf(pf2,"%d",&k);
    fscanf(pf2,"%s",String);
    
    L2 = Composition(String, k);
    for(NodeList *it = listGetHeader(L2); it!=NULL; it = nodeListGetNext(it)){
        fprintf(pf,"%s\n",nodeListGetCont(it));
        cont++;
    }
    printf("%d\n",cont);
    fclose(pf);
    fclose(pf2);
    printf("%d",fileEqual("datos.txt","Answers.txt",StringReadFromFile,(cmpfn)strcmp));
    //frec = Neighbors(pattern,d);
    
    //revDna = ReverseComplementProblem(dna);
    //fflush(stdin);
    //printf("%s  \n",revDna);
    
    
    //matching = PatternMatching(pattern,genome);
    //for(int i = 0; i< strlen(genome); i++)
    //    printf("%d ", matching[i]);
    /**
    num = PatternToNumber(pattern);
    printf("%ld\n",num);
    //printf("%s",pattern);
    
    pattern2 = NumberToPattern(index,k2);
    printf("%s",pattern2);
    */
    
    //frecArray = ComputingFrequencies(dna,k);
    //for(int i = 0; i < pow(4,k); i++)
     //   printf("%ld ",frecArray[i]);
    
    //frec = BetterClumpFindingProblem(genome,k,L,t);
    /**
    for(int i = 0; i < num; i++)
        if(frec[i] != NULL)
            k2++;
    
    printf("%d",k2);
    
    */
    
    //matching = MaximumSkew(dna);
    //for(int i = 0; i <= strlen(dna); i++)
      //  printf("%d ",matching[i]);
    
    //printf("\n");
    //printf("%d",HammingDistanceProblem(String,String2));
    //printf("\n");
    //matching = ApproximatePatternMatching(pattern,genome,d);
    //for(int i = 0; i <= strlen(genome); i++)
    //printf("%d ",ApproximatePatternCount("CGTGACAGTGTATGGGCATCTTT","TGT",1));
    
    return (EXIT_SUCCESS);
}

//Primer modulo

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

unsigned int maxValue(int Array[] , unsigned int size){
    int value = Array[0];
    
    for(int i = 0; i < size; i++)
        if(value < Array[i])
            value = Array[i];
    
    return value;
}

char **frecuentWords(char text[], int k){
    char    //(*FrecuentPatterns)[k] = malloc(sizeof(char)*k*(strlen(text) - k)),
            **FrecuentPatterns,
            *pattern = malloc(sizeof(char)*k),
            **FrecuentWords = malloc(sizeof(char*)*(strlen(text)-k));
    
    unsigned int    i, maxCount; 
    int             *Count = (int *)malloc(sizeof(int)*strlen(text));
    
    FrecuentPatterns = malloc(sizeof(char*)*strlen(text) - k);
    
    for(i = 0; i < strlen(text) - k; i++){
        FrecuentPatterns[i] = malloc(sizeof(char)*k);
        for(int j = 0; j < k; j++)
            FrecuentPatterns[i][j]='\0';
    }
    for(i = 0; i <= (strlen(text) - k); i++){
        for(int j = i, z=0; z < k; j++,z++)
            pattern[z] = text[j];
        Count[i] = PatternCount(text,pattern);    
    }
    maxCount = maxValue(Count,strlen(text));
    
    for(i = 0; i <= (strlen(text) - k); i++){
        if(Count[i] == maxCount){
            for(int j = i, z=0; z < k; j++,z++)
                pattern[z] = text[j];
            strcpy(FrecuentPatterns[i],pattern);
        }
    }
    
    int z = 0;
    char *temp = malloc(sizeof(char)*k);
    for(i = 0; i < k; i++)
        temp[i] = '\0';
    
    for(i = 0; i < strlen(text)-k; i++){
        FrecuentWords[i] = NULL;
        for(int j = i + 1; j < strlen(text)-k; j++)
            if(strcmp(FrecuentPatterns[i],FrecuentPatterns[j]) == 0){
                strcpy(FrecuentPatterns[i],temp);
                break;
            }
    }
    return FrecuentPatterns;
}

char *ReverseComplementProblem(const char DNAString[]){
    int len, i;
    char *DNAReverse;
    
    len = strlen(DNAString);
    
    DNAReverse = malloc(sizeof(char)*(len+1));
    
    if( DNAReverse == NULL)
        return NULL;
    
    for(i = 0; i < len; i++){
        DNAReverse[i] = '\0';
        if(DNAString[len -i -1] == 'A')
            DNAReverse[i] = 'T' ;
        if(DNAString[len -i -1] == 'T')
            DNAReverse[i] = 'A' ;
        if(DNAString[len -i -1] == 'C')
            DNAReverse[i] = 'G' ;
        if(DNAString[len -i -1] == 'G')
            DNAReverse[i] = 'C' ;
    }
    
    return DNAReverse;
}

int *PatternMatching(const char Pattern[], const char Genome[]){
    unsigned int lenPattern,lenGenome,flag = 0;
    int *values;
    
    lenPattern = strlen(Pattern);
    lenGenome = strlen(Genome);
    
    values = malloc(sizeof(int)*lenGenome);
    if(values == NULL)
        return NULL;
    
    int k=0;
    for(int i = 0; i < lenGenome; i++){
        values[i] = 0;
        flag = 0;
        if(Genome[i] == Pattern[0]){
            
            for(int j = 1,z = i+1; j < lenPattern & z < lenGenome; j++,z++)
                if(Genome[z] == Pattern[j])
                    flag = 1;
                else{
                    flag = 0;
                    break;
                }
            if(flag){   
                values[k] = i;
                k++;
            }
        }    
    }
    
    return values;
}
/*
  FrequentPatterns ← an empty set
        for i ←0 to 4k − 1
            Clump(i) ← 0
        for i ← 0 to |Genome| − L
            Text ← the string of length L starting at position i in Genome 
            FrequencyArray ← ComputingFrequencies(Text, k)
            for j ← 0 to 4k − 1
                if FrequencyArray(j) ≥ t
                    Clump(j) ← 1
        for i ← 0 to 4k − 1
            if Clump(i) = 1
                Pattern ← NumberToPattern(i, k)
                add Pattern to the set FrequentPatterns
        return FrequentPatterns
 */
char **ClumpFindingProblem(const char Genome[], unsigned int k, unsigned int L, unsigned int t){
    long int *FrecuencyArray;
    int num = pow(4,k),i = 0;
    int *Clump = malloc(sizeof(int)*num);
    char *Pattern, *text = malloc(sizeof(char)*L);
    char **FrecuentPatterns = malloc(sizeof(char*)*num);
    
    for(i = 0; i < num; i++)
        Clump[i] = 0;
    
    for(i = 0; i <= (strlen(Genome) - L); i++){
        for(int j = i, z = 0; z < L; j++,z++)
            text[z] = Genome[j];
        FrecuencyArray = ComputingFrequencies(text,k);
        for(int h = 0; h < num; h++)
            if(FrecuencyArray[h] >= t)
                Clump[h] = 1;
        free(FrecuencyArray);
    }
    //free(FrecuencyArray);
    free(text);
    
    for(i = 0; i < num; i++)
        if(Clump[i] == 1){
            Pattern = NumberToPattern(i,k);
            FrecuentPatterns[i] = Pattern;
        }
    
    return FrecuentPatterns;
}

char **BetterClumpFindingProblem(const char Genome[], unsigned int k, unsigned int L, unsigned int t){
    long int *FrecuencyArray,x,cont;
    int num = pow(4,k),i = 0;
    int *Clump = malloc(sizeof(int)*num);
    char *Pattern, *text = malloc(sizeof(char)*(L)),*firstPattern = malloc(sizeof(char)*(k)), *lastPattern = malloc(sizeof(char)*(k));
    char **FrecuentPatterns = malloc(sizeof(char*)*num);
    
    for(i = 0; i < num; i++)
        Clump[i] = 0;
    
    for(i = 0; i<L; i++)
        text[i] = Genome[i];
    //strncpy(text,Genome,L);
    FrecuencyArray = ComputingFrequencies(text,k);
    
    for(int h = 0; h < num; h++)
        if(FrecuencyArray[h] >= t)
            Clump[h] = 1;
    cont = strlen(Genome) - L;
    for(i = 1; i <= (strlen(Genome) - L); i++){
        for(int j = i-1, z = 0; z < k; j++,z++)
            firstPattern[z] = Genome[j];
        x = PatternToNumber(firstPattern);
        FrecuencyArray[x] = FrecuencyArray[x] - 1;
        for(int j = i+L-k ,z = 0; z < k; j++,z++)
            lastPattern[z] = Genome[j];
        x = PatternToNumber(lastPattern);
        FrecuencyArray[x] = FrecuencyArray[x] + 1;
        if(FrecuencyArray[x] >= t)
            Clump[x] = 1;
        
    }
    free(firstPattern);
    free(lastPattern);
    free(FrecuencyArray);
    free(text);
    
    for(i = 0; i < num; i++)
        if(Clump[i] == 1){
            Pattern = NumberToPattern(i,k);
            FrecuentPatterns[i] = Pattern;
        }
    
    free(Clump);
    return FrecuentPatterns;
}

unsigned short int SymbolToNumber(const char symbol){
    switch (symbol){
        case 'A':{
            return 0;
        }break;
        case 'C':{
            return 1;
        }break;
        case 'G':{
            return 2;
        }break;
        case 'T':{
            return 3;
        }break;
        default:{
        }
    }
    
}

unsigned long int PatternToNumber(char Pattern[]){
    int len = strlen(Pattern);
    char symbol = '\0';
    static unsigned long int num = 0;
    
    if(len == 0)
        return 0;
    
    symbol = Pattern[len-1];
    Pattern[len-1] = '\0';
    num = 4*PatternToNumber(Pattern)+SymbolToNumber(symbol);
    
    return num;
}

char *NumberToSymbol(int index){
    switch (index){
        case 0:{
            return "A";
        }break;
        case 1:{
            return "C";
        }break;
        case 2:{
            return "G";
        }break;
        case 3:{
            return "T";
        }break;
        
    }
}

char *NumberToPattern(int index,int k){
    int prefixIndex = 0, remainder = 0;
    static char tem;
    char *PrefixPattern,*symbol;
    
    if(k == 1)
        return NumberToSymbol(index);
    
    prefixIndex = index/4;
    remainder = index%4;
    
    PrefixPattern = NumberToPattern(prefixIndex,k-1);
    symbol = NumberToSymbol(remainder);
    
    return strBuilder(PrefixPattern,symbol);
}

char *strBuilder(const char *c1,const char *c2){
    char *temp = malloc(sizeof(c1)+sizeof(c2));
    
    strcat(strcat(temp,c1),c2);
    
    return temp;
}
long int *ComputingFrequencies(const char DNAStr[] ,int k){
    int num = pow(4,k),i = 0, j = 0, z = 0;
    long int x = 0;
    long int *FrecuencyArray = malloc(sizeof(long int)*num);
    char *pattern = malloc(sizeof(char)*(k+1));
    
    for(i = 0; i < num; i++){
        FrecuencyArray[i] = 0;
    }
    
    for(i = 0; i <= (strlen(DNAStr) - k); i++){
        for(j = i, z=0; z < k; j++,z++)
            pattern[z] = DNAStr[j];
        z = PatternToNumber(pattern);
        FrecuencyArray[z] = FrecuencyArray[z] + 1;
    }
    
    free(pattern);
    return FrecuencyArray;
}

int *Skew(const char Genome[]){
    int *skew = malloc(sizeof(int)*(strlen(Genome)+2));
    int length = strlen(Genome);
    
    skew[0] = 0;
    for(int i = 0; i <= length; i++){
        skew[i+1] = 0;
        
        if(Genome[i] == 'C'){
            skew[i+1] = skew[i] - 1;
            continue;
        }
        if(Genome[i] == 'G'){
            skew[i+1] = skew[i] + 1;
            continue;
        }
        
        skew[i+1] = skew[i];
    }
    
    return skew;
}

int *MinimumSkew(const char Genome[]){
    int length = strlen(Genome), *skew, value = 0, z = 0;
    int *minskew = malloc(sizeof(int)*length);
    
    skew = Skew(Genome);
    value = skew[0];
    for(int i = 0; i < length; i++){
        if(value >= skew[i])
            value = skew[i];
    }
    
    for(int i = 0; i <= length; i++){
        minskew[i] = 0;
        if(value == skew[i]){
            minskew[z] = i;
            z++;
        }
    }
    return minskew;
}
int *MaximumSkew(const char Genome[]){
    int length = strlen(Genome), *skew, value = 0, z = 0;
    int *minskew = malloc(sizeof(int)*length);
    
    skew = Skew(Genome);
    value = skew[0];
    for(int i = 0; i < length; i++){
        if(value < skew[i])
            value = skew[i];
    }
    
    for(int i = 0; i <= length; i++){
        minskew[i] = 0;
        if(value == skew[i]){
            minskew[z] = i;
            z++;
        }
    }
    return minskew;
}
int HammingDistanceProblem(const char String[],const char String2[]){
    int distance = 0;
    
    for(int i = 0; i < strlen(String); i++)
        if(String[i] != String2[i])
            distance++;
    
    return distance;
}

int *ApproximatePatternMatching(const char Pattern[],const char Genome[],int d){
    unsigned int lenPattern,lenGenome,flag = 0;
    int *values;
    
    char *PartGenome;
    
    lenPattern = strlen(Pattern);
    lenGenome = strlen(Genome);
    
    values = malloc(sizeof(int)*lenGenome);
    PartGenome = malloc(sizeof(char)*lenPattern);
    
    if(values == NULL)
        return NULL;
    
    int k=0;
    for(int i = 0; i < lenGenome; i++){
        values[i] = 0;
        flag = 0;
        PartGenome[0] = Genome[i];
        for(int j = 1,z = i+1; j < lenPattern & z < lenGenome; j++,z++)
            PartGenome[j] = Genome[z];
                //flag = 1;
        
        if(HammingDistanceProblem(PartGenome,Pattern) <= d){
            values[k] = i;
            k++;
        }
          
    }
    
    return values;
}

int ApproximatePatternCount(const char Text[], const char Pattern[],int d){
    int count = 0;
    char *pattern = malloc(sizeof(char)*strlen(Pattern));
    
    for(int i = 0; i <= (strlen(Text)-strlen(Pattern)); i++){
        for(int j = 0,z = i; j < strlen(Pattern); j++,z++)
            pattern[j] = Text[z];
        if(HammingDistanceProblem(pattern,Pattern) <= d)
            count = count + 1;
    }
    
    return count;
}

char **FrequentWordswithMismatchesProblem(const char Text[], int k, int d){
    char    //(*FrecuentPatterns)[k] = malloc(sizeof(char)*k*(strlen(text) - k)),
            **FrecuentPatterns,
            *pattern = malloc(sizeof(char)*k),
            **FrecuentWords = malloc(sizeof(char*)*(strlen(Text)-k));
    
    unsigned int    i, maxCount; 
    int             *Count = (int *)malloc(sizeof(int)*strlen(Text));
    
    FrecuentPatterns = malloc(sizeof(char*)*(strlen(Text) - k));
    
    for(i = 0; i < (strlen(Text) - k); i++){
        FrecuentPatterns[i] = NULL;
    }
    for(i = 0; i <= (strlen(Text) - k); i++){
        for(int j = i, z=0; z < k; j++,z++)
            pattern[z] = Text[j];
        Count[i] = ApproximatePatternCount(Text,pattern,d);    
    }
    maxCount = maxValue(Count,strlen(Text));
    
    for(i = 0; i <= (strlen(Text) - k); i++){
        if(Count[i] == maxCount){
            pattern = malloc(sizeof(char)*(k+1));
            for(int j = i, z=0; z <= k; j++,z++)
                pattern[z] = Text[j];
            FrecuentPatterns[i] = pattern;
        }
    }
    
    return FrecuentPatterns;
}
/**
static char **Neighbors(char Pattern[],int d){
    //char **Neighborhood2;
    static char **SuffixNeighbors = NULL;
    char *symbol = NULL, *NucleotideSet[] = {"A","C","G","T"}; 
    int length = pow(strlen(Pattern),strlen(Pattern));
    char **Neighborhood = malloc(sizeof(char*)*(length));
    
    long int i = 0;
    
    if(d == 0){
        Neighborhood[0] = Pattern;
        return Neighborhood;
    }
    for(int i = 0; i < (length); i++)
        Neighborhood[i] = NULL;
    
    if (length == 1){
        for(int i = 0; i < 4; i++)
            Neighborhood[i] = NucleotideSet[i];
        return Neighborhood;
    }
    
    SuffixNeighbors = Neighbors(Suffix(Pattern),d);
    
    //Neighborhood2 = malloc(sizeof(SuffixNeighbors)*4);
    
    i = 0;
    
    int j,z = 0;    
    for(i = 0; i < pow(strlen(Pattern)-1,strlen(Pattern)-1); i++){
        if(SuffixNeighbors[i] != NULL){
            char *text = SuffixNeighbors[i];
            char *temp = Suffix(Pattern);
            if(HammingDistanceProblem(temp,text) < d){
                for(j = 0; j < 4; j++){
                    char *temp2;
                    temp2 = strBuilder(NucleotideSet[j],text);
                    Neighborhood[z] = temp2;
                    z++;
                }
            }
            else{
                Neighborhood[z] = strBuilder(FirstSymbol(Pattern),text);
                z++;
            }
        }
    }
    
    return Neighborhood;
}
*/
char *Suffix(char *Pattern){
    int length = strlen(Pattern);
    char *suffix = malloc(sizeof(char)*length);
    
    for(int i = 1,z = 0; i < length; i++,z++){
        suffix[z+1] = '\0';
        suffix[z] = Pattern[i];
    }
    Pattern = suffix;
    
    return suffix;
}

char *FirstSymbol(const char Pattern[]){
    char *symbol = malloc(sizeof(char));
    
    *symbol = Pattern[0];
    
    return symbol;
}

List *Neighbors(char Pattern[],int d){
    //char **Neighborhood2;
    List *SuffixNeighbors = listNew();
    List *Neighborhood = listNew();
    //static char **SuffixNeighbors = NULL;
    char *symbol = NULL, *NucleotideSet[] = {"A","C","G","T"}; 
    int length = strlen(Pattern);
    //char **Neighborhood = malloc(sizeof(char*)*(length));
    
    long int i = 0;
    
    if(d == 0){
        listAddNode(Neighborhood,nodeListNew(strdup(Pattern)));
        return Neighborhood;
    }
    if (length == 1){
        for(int i = 0; i < 4; i++)
            listAddNode(Neighborhood,nodeListNew(NucleotideSet[i]));
        return Neighborhood;
    }
    
    SuffixNeighbors = Neighbors(Suffix(Pattern),d);
    
    i = 0;
    
    int j,z = 0;    
    
    for(NodeList *iterator = listGetHeader(SuffixNeighbors); 
        iterator != NULL; iterator = nodeListGetNext(iterator)){
        char *text = nodeListGetCont(iterator);
        char *temp = Suffix(Pattern);
        if(HammingDistanceProblem(temp,text) < d){
            for(j = 0; j < 4; j++){
                char *temp2;
                temp2 = strBuilder(NucleotideSet[j],text);
                listAddNode(Neighborhood,nodeListNew(temp2));
            }
        }
        else{
            listAddNode(Neighborhood,nodeListNew(strBuilder(FirstSymbol(Pattern),text)));
        }
    }
    
    return Neighborhood;
}

List *FrequentWordsWithMismatches(const char Text[], int k,int d){
    List *FrequentPatterns = listNew();
    List *Neighborhood;
    int *Close = malloc(sizeof(int)*pow(4,k));
    int *FrecuencyArray = malloc(sizeof(int)*pow(4,k));
    int index, maxCount;
    char *temp = malloc(sizeof(char)*(k + 1));
    
    for(int i = 0; i < (pow(4,k)); i++){
        Close[i] = 0;
        FrecuencyArray[i] = 0;
    }
    
    for(int i = 0; i <= (strlen(Text) - k); i++){
        for(int j = i,z = 0; (z < k) && (j <= (strlen(Text) - k)) ; j++,z++)
            temp[z] = Text[j];
        Neighborhood = Neighbors(temp,d);
        for(NodeList *it = listGetHeader(Neighborhood);
            it != NULL; it = nodeListGetNext(it)){
            
            index = PatternToNumber(nodeListGetCont(it));
            Close[index] = 1;
        }
    }
    
    for(int i = 0; i < pow(4,k);i++)
        if(Close[i] == 1){
            char *Pattern;
            Pattern = NumberToPattern(i,k);
            FrecuencyArray[i] = ApproximatePatternCount(Text,Pattern,d);
        }
    
    maxCount = maxValue(FrecuencyArray,pow(4,k));
    
    for(int i = 0; i < pow(4,k); i++)
        if(FrecuencyArray[i] == maxCount){
            char *Pattern;
            Pattern = NumberToPattern(i,k);
            listAddNode(FrequentPatterns,nodeListNew(Pattern));
        }
            
    return FrequentPatterns;
}

List *FrequentWordsWithMismatchesAndReverseComplements(const char Text[], int k,int d){
    List *FrequentPatterns = listNew();
    List *Neighborhood;
    int *Close = malloc(sizeof(int)*pow(4,k));
    int *FrecuencyArray = malloc(sizeof(int)*pow(4,k));
    int index, maxCount;
    char *temp = malloc(sizeof(char)*(k + 1));
    
    for(int i = 0; i < (pow(4,k)); i++){
        Close[i] = 0;
        FrecuencyArray[i] = 0;
    }
    
    for(int i = 0; i <= (strlen(Text) - k); i++){
        for(int j = i,z = 0; (z < k) && (j <= (strlen(Text) - k)) ; j++,z++)
            temp[z] = Text[j];
        Neighborhood = Neighbors(temp,d);
        for(NodeList *it = listGetHeader(Neighborhood);
            it != NULL; it = nodeListGetNext(it)){
            
            index = PatternToNumber(nodeListGetCont(it));
            Close[index] = 1;
        }
    }
    
    for(int i = 0; i < pow(4,k);i++)
        if(Close[i] == 1){
            char *Pattern;
            Pattern = NumberToPattern(i,k);
            FrecuencyArray[i] = ApproximatePatternCount(Text,Pattern,d) 
                                + ApproximatePatternCount(Text,ReverseComplementProblem(Pattern),d);
        }
    
    maxCount = maxValue(FrecuencyArray,pow(4,k));
    
    for(int i = 0; i < pow(4,k); i++)
        if(FrecuencyArray[i] == maxCount){
            char *Pattern;
            Pattern = NumberToPattern(i,k);
            listAddNode(FrequentPatterns,nodeListNew(Pattern));
        }
            
    return FrequentPatterns;
}

List *MotifEnumeration(List *DNA, int k ,int d){
    List *Patterns = listNew();
    int numWords = listGetSize(DNA);
    List **patternsList = malloc(sizeof(List*)*numWords);
    char *patternPrim;
    int j;
    
    int z = 0;
    for(NodeList *it = listGetHeader(DNA); it != NULL; it = nodeListGetNext(it), z++){
        char *pattern = nodeListGetCont(it);
        List *kmers = listNew();
        for(int i = 0; i <= (strlen(pattern) - k); i++){
            patternPrim = strCreateSelect(pattern,i,k);
            List *kmers2 = Neighbors(patternPrim,d);
            j = listGetSize(kmers2);
            kmers = listCat(kmers,kmers2);
            int o = listGetSize(kmers);
        }
        patternsList[z] = kmers;
        j = listGetSize(kmers);
    }

    for(NodeList *it = listGetHeader(patternsList[0]); it != NULL; it = nodeListGetNext(it)){
        int flag = 0;
        for(int i = 1; i < numWords; i++){
            char *temp = nodeListGetCont(it);
            if(listElemExist(patternsList[i], temp, (cmpfn)strcmp) == 1)
                flag = 1;
            else{
                flag = 0;
                break;
            }
        }
        if(flag)
            listAddNode(Patterns,nodeListNew(strdup(nodeListGetCont(it))));
    }
    
    listRemoveDuplicates(Patterns,(cmpfn)strcmp);
    
    return Patterns;
}

Generic StringReadFromFile(FILE *pf){
    char *lee=malloc(sizeof(char)*1000000);
    int longi;

    if(pf!=NULL) fscanf(pf,"%[^\n]",lee);
    
    longi = strlen(lee);
    lee = reallocf(lee,sizeof(char)*(longi+1));
    fgetc(pf);
    return lee;
}

void listRemoveDuplicates(List *L, cmpfn fn){
    for(NodeList *it = listGetHeader(L); it != NULL; it = nodeListGetNext(it)){
        NodeList *it2 = nodeListGetNext(it);
        if(it2 == NULL)
            it2 = listGetHeader(L);
        for(; ; it2 = nodeListGetNext(it2)){
            char *temp = nodeListGetCont(it);
            char *temp2 = nodeListGetCont(it2);
            
            if(nodeListGetNext(it2) == NULL && it != listGetLast(L))
               it2 = listGetHeader(L);
            
            if(it == it2){
                break;
            }
            
            if(fn(nodeListGetCont(it2),nodeListGetCont(it)) == 0){
                NodeList *it3 = listGetPrevious(L,it2);
                listRemoveNode(L,it2);
                it2 = it3;
            }
        }
    }
}

char *strCreateSelect(const char Text[], int posIni, int length){
    char *temp = malloc(sizeof(char)*length);
    
    for(int j = posIni, z = 0; z < length; j++,z++)
        temp[z] = Text[j];
    
    return temp;
}

int listElemExist(List *L, Generic cont, cmpfn fn){
    for(NodeList *it = listGetHeader(L); it != NULL; it = nodeListGetNext(it)){
        if(fn(cont,nodeListGetCont(it)) == 0)
            return 1;
    }
    return 0;
}

List *listCat(List *L1,List *L2){
    if(listIsEmpty(L1))
        L1 = L2;
    else{
        nodeListSetNext(listGetLast(L1),listGetHeader(L2));
        L1->last = L2->last;
    }
    return L1;
}

int DistanceBetweenPatternAndStrings(const char Pattern[], List *DNA){
    int k = strlen(Pattern);
    int distance = 0;
    int i = listGetSize(DNA);
    
    for(NodeList *it = listGetHeader(DNA); it != NULL; it = nodeListGetNext(it)){
        char *text = nodeListGetCont(it);
        int HammingDistance = INT32_MAX;
        for(int i = 0; i <= (strlen(text)-k); i++){
            char *patternPrim = strCreateSelect(text,i,k);
            if(HammingDistance > HammingDistanceProblem(Pattern,patternPrim))
                HammingDistance = HammingDistanceProblem(Pattern,patternPrim);
        }
        distance = distance + HammingDistance;
    }
    
    return distance;
}

char *MedianString(List *DNA, int k){
    int distance = INT32_MAX;
    char *Median;
    
    for(int i = 0; i <= (pow(4,k)-1); i++){
        char *Pattern = NumberToPattern(i,k);
        if( distance > DistanceBetweenPatternAndStrings(Pattern,DNA)){
            distance = DistanceBetweenPatternAndStrings(Pattern,DNA);
            Median = Pattern;
        }
    }
    
    return Median;
}

float **readMatrixFromFile(const char file[],int nrow,int ncol){
    float **Matrix = malloc(sizeof(float*)*nrow);
    int i = 0, j = 0;
    FILE *pf = fopen(file, "r");
    
    for(int i = 0 ; i < nrow; i++)
        Matrix[i] = malloc(sizeof(float)*ncol);
    
    while(!feof(pf)){
        char *temp = malloc(sizeof(char*)*ncol);
        char *pch;
        fscanf(pf,"%[^\n]",temp);
        pch = strtok(temp," ");
        j = 0;
        while(pch != NULL){
            float num = atof(pch);
            Matrix[i][j] = num;
            pch = strtok(NULL," ");
            j++;
        }    
        fgetc(pf);
        i++;
    }
    for(int i = 0 ; i < 4; i++){
        for(int j = 0; j < ncol ; j++)
            printf("%f ",Matrix[i][j]);
        printf("\n");
    }
    return Matrix;
}

char *ProfileMostProbableKmer(const char Text[], int k, float **Profile){
    int len = strlen(Text);
    char *MostProbKmer = malloc(sizeof(char)*(k+1));
    float prob = 0;
    
    for(int i = 0; i <= (len - k); i++){
        char *temp = strCreateSelect(Text,i,k);
        float tempfloat = 1;
        for(int j = 0; j < k ; j++){
            switch(temp[j]){
                case 'A':{
                    tempfloat = tempfloat*Profile[0][j];
                }break;
                case 'C':{
                    tempfloat = tempfloat*Profile[1][j];
                }break;
                case 'G':{
                    tempfloat = tempfloat*Profile[2][j];
                }break;
                case 'T':{
                    tempfloat = tempfloat*Profile[3][j];
                }break; 
            }
        }
        if(prob < tempfloat){
            prob = tempfloat;
            strcpy(MostProbKmer,temp);
        }
    }
    
    return MostProbKmer;
}

List *GreedyMotifSearch(List *DNA, int k, int t){
    List *bestMotifs = listNew(), *otherStrands = listNew(), *motifs = listNew();
    char *baseStrand;
    
    for(NodeList *it = listGetHeader(DNA); it != NULL; it = nodeListGetNext(it)){
        char *temp = strCreateSelect(nodeListGetCont(it),0,k);
        listAddNode(bestMotifs,nodeListNew(temp));
    }
    
    baseStrand = nodeListGetCont(listGetHeader(DNA));
    for(NodeList *it = nodeListGetNext(listGetHeader(DNA)); it != NULL; it = nodeListGetNext(it)){
        listAddNode(otherStrands,nodeListNew(nodeListGetCont(it)));
    }
    
    for(int i = 0; i < (strlen(baseStrand) - k); i++){
        char *temp = strCreateSelect(baseStrand,i,k);
        //List *motifs
        listAddNode(motifs,nodeListNew(temp));
        for(NodeList *it = listGetHeader(otherStrands); it != NULL; it = nodeListGetNext(it)){
            char *strand = nodeListGetCont(it);
            float **profileMatrix = createProfileMatrix(motifs);
            char *nextMotif = ProfileMostProbableKmer(strand,k,profileMatrix);
            listAddNode(motifs,nodeListNew(nextMotif));
        }
        if(Score(motifs) < Score(bestMotifs))
            bestMotifs = motifs;    
    }
    
    return bestMotifs;
}
float **createProfileMatrix(List *motifs){
    float **matrix = malloc(sizeof(float*)*4);
    int j = 0, ncols = listGetSize(motifs);
    
    for(int i = 0; i < 4; i++)
        matrix[i] = malloc(sizeof(float)*ncols);
    
    for(NodeList *it = listGetHeader(motifs); it != NULL; it = nodeListGetNext(it),j++){
        char *temp = nodeListGetCont(it);
        int a = Count(temp,'A');
        int c = Count(temp,'C');
        int g = Count(temp,'G');
        int t = Count(temp,'T');
        
        if(a == 0) matrix[0][j] = 0.0;
        else matrix[0][j] = 1/a;
        if(c == 0) matrix[1][j] = 0.0;
        else matrix[1][j] = 1/c;
        if(g == 0) matrix[2][j] = 0.0;
        else matrix[2][j] = 1/g;
        if(t == 0) matrix[3][j] = 0.0;
        else matrix[3][j] = 1/t;
        
    }
    return matrix;
}
int Count(const char Text[], char symbol){
    int cont = 0;
    
    for(int i = 0; i < strlen(Text); i++){
        if(Text[i] == symbol)
            cont++;
    }
    
    return cont;
}
int Score(List *motifs){
    //char *consensus = findConsensus(motifs);
    int score = 0;
    /**
    for(NodeList *it = listGetHeader(motifs); it != NULL; it = nodeListGetNext(it)){
        //score += HammingDistanceProblem(consensus,nodeListGetCont(it));
        char *temp = nodeListGetCont(it);
        int len = strlen(temp);
        int flag1 = 0,flag2 = 0;
        int a = Count(temp,'A');
        int c = Count(temp,'C');
        int g = Count(temp,'G');
        int t = Count(temp,'T');
        
        if(a > c)
            flag1 = a;
        else
            flag1 = c;
        if(g > t)
            flag2 = g;
        else
            flag2 = t;
        
        if(flag1 >= flag2)
            score = score + len - flag1;
        else 
            score = score + len - flag2; 
        
    }
    */
    char *consensus = findConsensus(motifs);
    for(NodeList *it = listGetHeader(motifs); it != NULL; it = nodeListGetNext(it)){
        char *motif = nodeListGetCont(it);
        score += HammingDistanceProblem(consensus,motif);
    }
    return score;
}
char *findConsensus(List *motifs){
    char *consensus,**matrix = malloc(sizeof(char*));
    int i = 0, len = 0;
    
    for(NodeList *it = listGetHeader(motifs); it!=NULL ; it = nodeListGetNext(it),i++){
        char *temp = nodeListGetCont(it);
        len = strlen(temp);
        matrix[i] = temp;
    }
    
    consensus = malloc(sizeof(char)*len);
    
    for(int j = 0; j < len; j++){
        char symbol;
        int cont,cont1,cont2,cont3;
        
        for(int z = 0; z < listGetSize(motifs); z++){
            switch(matrix[j][z]){
                case 'A':{
                    cont++;
                }break;
                case 'C':{
                    cont1++;
                }break;
                case 'G':{
                    cont2++;
                }break;
                case 'T':{
                    cont3++;
                }break;
            }
        }
        
        if((cont >= cont2)&&(cont >= cont1)&&(cont >= cont3))
            symbol = 'A';
        else if((cont1 >= cont)&&(cont1 >= cont2)&&(cont1 >= cont3))
            symbol = 'C';
        else if((cont2 >= cont)&&(cont2 >= cont1)&&(cont2 >= cont3))
            symbol = 'G';
        else if((cont3 >= cont)&&(cont3 >= cont1)&&(cont3 >= cont2))
            symbol = 'T';
        
        consensus[j] = symbol;
    }
    return consensus;
}

/* Segundo Curso de BioAlgoritmos*/

Queue *Composition(const char Text[], int k){
    Queue *composition = listNew(), *temp = listNew();
    
    for(int i = 0; i < (strlen(Text)-k); i++){
        queueEnqueue(composition,nodeListNew(strCreateSelect(Text,i,k)));
    }
        
    while(!queueIsEmpty(composition)){
        NodeList *cont = queueDequeuePriority(composition,(cmpfn)strAlfCmp);
        char *tempString = nodeListGetCont(cont);
        queueEnqueue(temp,cont);
    }
    
    composition = temp;
    
    return composition;
}

int strAlfCmp(const char str1[], const char str2[]){
    int lenStr1 = strlen(str1), lenStr2 = strlen(str2);
    int weigthStr1 = 0, weigthStr2 = 0, cont = 0;
    
    if(lenStr1 == lenStr2){
        for(int i = 0; i < lenStr1; i++){
            if(str1[i] < str2[i])
                return -1;
            else if(str1[i] > str2[i])
                return 1;
        }
    }
    
    if(lenStr1 < lenStr2){
        return 1;
    }
    
    if(lenStr1 > lenStr2){
        return -1;
    }
}

int fileEqual(char file1[], char file2[], readfn rfn, cmpfn cfn){
    List *L1 = listReadFile(file1,rfn);
    List *L2 = listReadFile(file2,rfn);
    
    int i = 0;
    
    for(NodeList *it = listGetHeader(L1); it != NULL ; it = nodeListGetNext(it)){
        for(NodeList *it2 = listGetHeader(L2); it2!= NULL; it2 = nodeListGetNext(it2)){
            char *temp = nodeListGetCont(it);
            char *temp2 = nodeListGetCont(it2);
            
            if(cfn(nodeListGetCont(it),nodeListGetCont(it2)) != 0){
                return 0;
            }
            i++;
        }
    }
    
    return 1;
}
