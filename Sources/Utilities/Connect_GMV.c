/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Program: Connect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#define MAXF  9999
#define MAXBC   10

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

int main(int argc, char *argv[])
 {
  int  i, j, k, N, s, Nsc, insc, isc;
  char base_name[80], name[80], format[80], 
       token1[80], token2[80], token3[80];
  file *fp[MAXF+1];
  int  NC[MAXF+1];
  int  NbC[MAXF+1][MAXBC+1];  
  int  visited;

  if(argc!=4)
   {
    printf("=====================================================\n");
    printf(" _________                                     __    \n");      
    printf(" \\_   ___ \\  ____   ____   ____   ____   _____/  |_        \n");
    printf(" /    \\  \\/ /  _ \\ /    \\ /    \\_/ __ \\_/ ___\\   __\\ \n");
    printf(" \\     \\___(  <_> )   |  \\   |  \\  ___/\\  \\___|  |     \n");
    printf("  \\______  /\\____/|___|  /___|  /\\___  >\\___  >__|       \n");
    printf("         \\/            \\/     \\/     \\/     \\/          \n");
    printf("                                                     \n"); 
    printf("Usage:\n");                          
    printf("  %s <base_name> <N> <format>\n", argv[0]);
    printf("  where:\n");
    printf("  <format> is either GMV, or DAT.\n");
    printf("=====================================================\n");
    return -1; 
   }

  strcpy(base_name,argv[1]);
  printf("%s\n", base_name);
  
  N=atoi(argv[2]);
  printf("%d\n", N);

  strcpy(format,argv[3]);
  printf("%s\n", base_name);
  

  for(i=0;i<strlen(base_name)-3;i++)
   {
    if( base_name[  i] == '-' &&   
        base_name[i+1] == '0' && base_name[i+1] == '0' &&
        base_name[i+2] == '0' && base_name[i+3] == '0' )
      {printf("Number starts at %d\n", i+1);
       s = i+1;}
   }

/*%%%%%%%%%%%%%%%%%%%%%%%
%  Open the basic file  %
%%%%%%%%%%%%%%%%%%%%%%%*/
    fp[0]=fopen(base_name,"w");

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Open all the remaining files  %
%  and read NC for each of them  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  for(i=1;i<=N;i++)
   {
    strcpy(name,base_name);
    if     (i <    10) sprintf(&name[s+3],"%d", i);
    else if(i <   100) sprintf(&name[s+2],"%2d",i);
    else if(i <  1000) sprintf(&name[s+1],"%3d",i);
    else if(i < 10000) sprintf(&name[s+0],"%4d",i);
    strcat(&name[s+3],&base_name[s+4]);

    printf("Opening: %s\n", name); 

    fp[i]=fopen(name,"r");
    fseek(fp[i],-9L,SEEK_END);         /* go 8 bytes before the end */  
    fscanf(fp[i],"%d", &NC[i]);  /* read that last digit */
    printf("%8d\n",NC[i]);       /* check what are you reading */
    for(j=1;j<=MAXBC;j++)
     {
      fscanf(fp[i],"%d", &NbC[i][j]);
      printf("%8d\n",NbC[i][j]);
     }
    fseek(fp[i],0L,0);           /* rewind */  
   }

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G. M. V. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  if(strcmp(format,"GMV")==0) 
   {

/*-------------+ 
|  velocities  |
+-------------*/
/* read "velocity" and "0" from all the files */
  for(i=1;i<=N;i++) 
   {
    fscanf(fp[i],"%s",token1);
    printf("==>%s\n", token1);
    fscanf(fp[i],"%s",token2); /* zero */
   }

/* write "velocity" and "0" to basic file */
  fprintf(fp[0]," %s %s\n", token1, token2);

/* read & write all U */
  for(i=1;i<=N;i++)
    for(j=1;j<=NC[i];j++)  
     {fscanf (fp[i],"%s",  token2);
      fprintf(fp[0],"%s\n",token2);}

/* read & write all V */
  for(i=1;i<=N;i++)
    for(j=1;j<=NC[i];j++)  
     {fscanf (fp[i],"%s",  token2);
      fprintf(fp[0],"%s\n",token2);}

/* read & write all W */
  for(i=1;i<=N;i++)
    for(j=1;j<=NC[i];j++)  
     {fscanf (fp[i],"%s",  token2);
      fprintf(fp[0],"%s\n",token2);}

/*------------+ 
|  variables  |
+------------*/
/* read "variable" from all the files */
  for(i=1;i<=N;i++)
   {
    fscanf(fp[i],"%s",token1);
    printf("==>%s\n", token1);
   }

/* write "variable" basic file */
  fprintf(fp[0]," %s\n", token1);

/* infinite loop through all the variables */
  for(;;)
   { 
/* read "var_name" and "0" from all the files */
    for(i=1;i<=N;i++)
     {
      fscanf(fp[i],"%s",token1);
/*      printf("%s\n", token1); */
      if(strcmp(token1,"endvars")==0) goto exit;
      fscanf(fp[i],"%s",token2); /* zero */ 
     }

/* write "var_name" and "0" to the basic file */
    fprintf(fp[0]," %s %s\n", token1, token2);

/* read & write all values */
    for(i=1;i<=N;i++)
      for(j=1;j<=NC[i];j++)  
       {fscanf( fp[i],"%s",  token2);
        fprintf(fp[0],"%s\n",token2);}
   } /* for (;;) */

exit:

  fprintf(fp[0]," endvars\n endgmv\n");
   } /* GMV */

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FLUENT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  else if(strcmp(format,"DAT")==0) 
   {
    int  numb, count, cel;
    char ch;  
    char dummy[20];

/*---------+ 
|  header  |
+---------*/
    fprintf(fp[0],"(0 \"============================\")\n");
    fprintf(fp[0],"(0 \"Created by T-Rex - Processor\")\n");
    fprintf(fp[0],"(0 \"============================\")\n");

    NC[0]=0;
    for(i=1;i<=N;i++)
      NC[0] += NC[i];
   
    for(j=1;j<=MAXBC;j++)
      NbC[0][j]=0;

    for(i=1;i<=N;i++)
      for(j=1;j<=MAXBC;j++)
        NbC[0][j] += NbC[i][j]; 

    for(j=1;j<=MAXBC;j++)
      printf("%d \n", NbC[0][j]); 

/*-------------+ 
|  velocities  |
+-------------*/
    fprintf(fp[0],"(0 \"velocity\")\n");
    fprintf(fp[0],"(300 (2 1 3 0 0 %d %d)(\n", 1, NC[0]);
    for(i=1;i<=N;i++)
     {
/* skip three lines */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* (0 "============================") */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* (0 "Created by T-Rex - Processor") */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* (0 "============================") */
/* number of scalars */
      fscanf(fp[i], "%s %d\n", dummy, &Nsc);   /* Nsc */                       
      Nsc++;                                   /* add pressure */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* to the end of line     */
/* skip two lines */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* (0 "velocity")         */
      do{ch=fgetc(fp[i]);} while(ch != '\n');  /* (300 (2 1 3 0 0 ...    */
      for(j=1;j<=NC[i];j++)
       {
        fscanf(fp[i],"%s %s %s",    token1, token2, token3);
        fprintf(fp[0],"%s %s %s\n", token1, token2, token3);
       }
/* skip two lines */
      do{ch=fgetc(fp[i]); } while(ch != '\n'); /* finish last line */
      do{ch=fgetc(fp[i]); } while(ch != '\n'); /* )) */
     }
    fprintf(fp[0], "))\n");

/*----------+
|  scalars  |
+----------*/
    printf("Nsc = %d\n", Nsc);
    for(isc=0; isc<Nsc; isc++)         /* Nsc!!!! replaced by 6 IZMJENA */
     {
      insc = 699 + isc;             /* user defined scalar */
      if(isc==0) insc = 1;          /* pressure is 1 */
      do{ch=fgetc(fp[1]); fputc(ch,fp[0]);} while(ch != '\n'); /* scalar name */
      fprintf(fp[0],"(300 (%3d 1 1 0 0 %d %d)(\n", insc, 1, NC[0]);
      for(i=1;i<=N;i++)
       {
/* skip two lines */
        if(i>1) do{ch=fgetc(fp[i]);} while(ch != '\n'); /* (0 "scalar name") */
                do{ch=fgetc(fp[i]);} while(ch != '\n'); /* (300 (1 101 ...   */
        for(j=1;j<=NC[i];j++)
         {
          fscanf(fp[i],"%s",    token1);
          fprintf(fp[0],"%s\n", token1);
         }
/* skip two lines */
        do{ch=fgetc(fp[i]);} while(ch != '\n'); /* finish last line */
        do{ch=fgetc(fp[i]);} while(ch != '\n'); /* )) */
       }
      fprintf(fp[0], "))\n");

      for(j=1;j<=MAXBC;j++)
       {
        visited = 0;      
        if(NbC[0][j] > 0)  
         {
          for(i=1;i<=N;i++)
            if(NbC[i][j] > 0)
             {
/* skip two lines */
              if(visited == 0)
               {do{ch=fgetc(fp[i]); fputc(ch, fp[0]);} while(ch != '\n');
                fprintf(fp[0],"(300 (%3d %3d 1 0 0 %d %d)(\n", insc, 100+j, 1, NbC[0][j]);
                visited++;}
              else
               {do{ch=fgetc(fp[i]);} while(ch != '\n');} /* (0 "scalar on the boundary") */
              do{ch=fgetc(fp[i]);} while(ch != '\n');    /* (300 (1 101 ... */
              for(k=1;k<=NbC[i][j];k++)
               {
                fscanf(fp[i],"%s",    token1);
                fprintf(fp[0],"%s\n", token1);
               }              
/* skip two lines */
              do{ch=fgetc(fp[i]);} while(ch != '\n');
              do{ch=fgetc(fp[i]);} while(ch != '\n'); /* )) */
             }
            fprintf(fp[0], "))\n");
         }
       }
     }

   } /* FLUENT */

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Close all the remaining files  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  for(i=1;i<=N;i++)
    fclose(fp[i]);

/*%%%%%%%%%%%%%%%%%%%%%%%%
%  Close the basic file  %
%%%%%%%%%%%%%%%%%%%%%%%%*/
  fclose(fp[0]);
  
  return 0;

 } /* main */
