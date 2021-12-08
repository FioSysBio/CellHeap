#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include <wait.h>
#include <string.h>
#define MAX 30

struct stack {
    char sra[25];
    struct stack * next;
} stack;

struct stack *  push (struct stack * s, char *sra, int size){
    struct stack *newS;
    newS = (struct stack *) malloc (sizeof(struct stack));
    if(strcmp(s->sra, "-") == 0){
        memcpy(s->sra, sra, size);
        return s;
    }else{
        memcpy(newS->sra, sra, size);
        newS->next = s;
        return newS;
    }
};

struct stack * pop (struct stack *s){
    char *tmp = s->sra;
    s=s->next;
    return s;
};

struct stack * initStack (){
    struct stack *s;
    char c[] = "-";
    s = (struct stack *) malloc (sizeof(struct stack));
    memcpy(s->sra, c, sizeof(c));
    s->next=NULL;
    return s;
};

int isValid(char *fileName, char* sra) {
    char * tmp;
    char reading[500];
    char is_consistent[]=".sra' is consistent";
    int value = 0;
    int i=0;
    char c;

    FILE *file; 
    file = fopen(fileName, "r"); 
    if (file==NULL) perror ("Error opening output file for validation\n");

    tmp = (char*) malloc (sizeof(char)*(strlen(sra)+strlen(is_consistent)));
    strcpy(tmp, sra);
    strcat(tmp, is_consistent);
	printf("PROCURANDO POR %s\n", tmp);
    while(fgets(reading, 500, file)!=NULL){
        if((strstr(reading, tmp))!=NULL){
            value = 1;
            break;
        }
    }

	fclose(file); 
 
    return value;

}

int main (int argc, char *argv[] )  {

    FILE *inFile;
    FILE *fileCom;
    FILE *outputFile;
    FILE *logFile;
    FILE *myLog;
    char output_name[] = "call_thread.out";
    char c;
    char tmp[MAX];
    char *outputCommand;
    int pMl = 0;
    int pMc = 0;
    int i=0;
    int col=0;
    struct stack *s;
    pid_t pid=1;
    pid_t wpid;
    int status = 0;
    char prefetch[] = "/scratch/app/sratoolkit/2.10.5/bin/prefetch --max-size 400G ";
    char validate[] = "/scratch/app/sratoolkit/2.10.5/bin/vdb-validate ";
    char prepath[] = " -O ";
    char path[500];
    char fastqdump[] = "/scratch/app/sratoolkit/2.10.5/bin/fastq-dump --split-files --gzip  ";
    char mv[]= "mv ";
    char r1In[] = "_1.fastq.gz ";
    char r2In[] = "_2.fastq.gz ";
    char r1Out[] = "_S1_L001_R1_001.fastq.gz";
    char r2Out[] = "_S1_L001_R2_001.fastq.gz";


    if( argc == 4 ) {
        printf("The input file is %s\nThe output directory is: %s\nThe log file is: %s\n", argv[1], argv[2], argv[3]);
    }
    else if( argc > 3 ) {
        printf("Too many arguments supplied.\n");
    }
    else {
        printf("One argument expected.\n");
    }

    inFile = fopen(argv[1], "r"); 
    if (inFile==NULL) perror ("Error opening input file\n");

    logFile = fopen(argv[3], "r"); 
    if (logFile==NULL) perror ("Error log file\n");

    //outputFile = fopen(output_name, "w+");
    //if (inFile==NULL) perror ("Error opening output file\n");


    myLog = fopen("my_log.txt", "w+"); 
    if (myLog==NULL) perror ("Error opening error file\n");


    strcpy(path, argv[2]);

    if(path[strlen(argv[2])-1]!='/'){
        strcat(path, "/");
    }
    printf("%s\n\n", path);
    s = initStack();

    do{
        c = fgetc(inFile);
        if(!isspace(c) && c!= EOF){
            tmp[i]=c;
            i++;
        }
        if(isspace(c)){
            s=push(s, tmp, i);
            col++;
            //printf("size: %d, value :%s\n", i, s->sra);
            i=0;
        }
    }while(c!=EOF);

    fclose(inFile);
    //printf("Entries: %d\n", col);


    i=col;
    while(i>0 && pid!=0){
        char *str;
        char to_file[]=" >> ";
        char sufixo[]= ".val";
        char remove[] = "rm ";
        int sys=0;
        int max_tries=10;
        str = s->sra;
        i--;
        s=pop(s);

        pid=fork();
        if(pid==0){
            int valido=0;
            char complete[500] = "";

            //faz o prefetch
            strcat(complete, prefetch);
            strcat(complete, str);
            strcat(complete, prepath);
            strcat(complete, path);
            printf("command prefetch %d: %s\n", i, complete);
            fputs(complete, myLog);
            sys=system(complete);
            //try again if error occurs
            while(sys!=0 && max_tries>0){
                system("sleep 3");
                sys=system(complete);
                printf("ERRO NO PREFETCH de sra %s\n", str);
                fputs("ERRO NO PREFETCH de sra %s\n", myLog);
                max_tries--;
            }
            max_tries=10;

            //verifica se o arquivo esta valido, mandando resultado para $sra.val
            memcpy(complete, "", sizeof(""));
            strcat(complete, validate);
            strcat(complete, path);
            strcat(complete, str);
            strcat(complete, ".sra");
            //strcat(complete, to_file);
            //strcat(complete, path);
            //strcat(complete, str);
            //strcat(complete, sufixo);
            ////system(complete); 
            fputs(complete, myLog);
            printf("command validate %d: %s\n", i, complete);
            sys=system(complete);
            //try again if error occurs
            while(sys!=0 && max_tries>0){
                sys=system(complete);
                max_tries--;
            }
            max_tries=10;

            // está válido?
            // memcpy(complete, "", sizeof(""));
            // strcat(complete, str);
            // strcat(complete, sufixo);
            // printf("temporary file %d: %s\n",i,complete);
            valido = isValid(argv[3], str);

            // // apagar arquivo temporário
            // memcpy(complete, "", sizeof(""));
            // strcat(complete, remove);
            // strcat(complete, str);
            // strcat(complete, sufixo);
            // printf("remove %d: %s\n",i,complete);
            // ////system(complete);

            
            if(valido==1){
                printf("%s file is valid\n", str);
                memcpy(complete, "", sizeof(""));
                strcat(complete, fastqdump);
                strcat(complete, path);
                strcat(complete, str);
                strcat(complete, ".sra ");
                strcat(complete, "-O ");
                strcat(complete, path);
                printf("command %d: %s\n", i, complete);
                fputs(complete, myLog);
                sys=system(complete);
                //try again if error occurs
                while(sys!=0 && max_tries>0){
                    sys=system(complete);
                    max_tries--;
                }
                max_tries=10;

                memcpy(complete, "", sizeof(""));
                strcat(complete, mv);
                strcat(complete, path);
                strcat(complete, str);
                strcat(complete, r1In);
                strcat(complete, path);
                strcat(complete, str);
                strcat(complete, r1Out);
                printf("command %d: %s\n", i, complete);
                fputs(complete, myLog);
                sys=system(complete);
                //try again if error occurs
                while(sys!=0 && max_tries>0){
                    sys=system(complete);
                    max_tries--;
                }
                max_tries=10;

                memcpy(complete, "", sizeof(""));
                strcat(complete, mv);
                strcat(complete, path);
                strcat(complete, str);
                strcat(complete, r2In);
                strcat(complete, path);
                strcat(complete, str);
                strcat(complete, r2Out);
                printf("command %d: %s\n", i, complete);
                fputs(complete, myLog);
                sys=system(complete);
                //try again if error occurs
                while(sys!=0 && max_tries>0){
                    sys=system(complete);
                    max_tries--;
                }
                max_tries=10;
                }else{
                    printf("not valid (unconsistent)\n");   
                    fputs("not valid (unconsistent)\n", myLog);
                }

        }else{
            if(pid == -1){
                char error[100];
		        system("echo falha na thread ");
                memcpy(error, "", sizeof(""));
                strcat(error, str);
                strcat(error, "\n");
                fputs(error, myLog);
            }
        }
    }
    
    while((wpid = wait(&status))>0);
    
    if(pid==0){
        exit(0);
    }else{        
        return 0;
    }
}
