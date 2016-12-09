//
//  BatchCorrectionMPI.c
//  
//
//  Created by Hemant Sharma on 6/28/13.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <stdint.h>

///////////
int
main(int argc, char *argv[]){
    
    clock_t start0, end;
    clock_t inter2, inter3, inter, start, start1;
    char *ParamFN;
    char *FolderName;
    int i=0;
    int LowNr;
    int DarkNr=1;
    int AllType=0;
    int KeepFiles=0;
    int GENumber=2;
    char *str;
    char *DarkFileStem="dark";
    char *OutputFilePath="";
    char DarkFileName[1000];
    char *FileStem="";
    char OutFileName[1000];
    char SumFileName[1000];
    FILE *fp;
    FILE *fb;
    FILE *ft;
    FILE *fc;
    FILE *fz;
    FILE *fr;
    char *FileList;
    char *BadFile="";
    char FileName[1000];
    char *BadFileName="";
    char *ext=".ge2";
    char *extcor=".cor";
    int PosHere;
    int DoNdel=0;
    int FileNumber;
    int SizeFile=2*2048*2048;
    int ElementNr=0;
    int Numels=2048*2048;
    char *aveext=".ave";
    char AveFileName[1000];
    char *sumext=".sum";
    double diftotal;
    double diffhere;
    double diffhere2;
    double diffhere3;
    int size;
    int nFrames;
    unsigned short int *Contents;
    double *Data;
    double *SumMat;
    unsigned short int *SumMat2;
    unsigned short int *AveMat2;
    double *AveMat;
    unsigned short int *SubData;
    double *Dark;
    int FrameNr=0;
    int NoExt=0;
    char *sep="/";
    char *DirName;
    int LowNrs[1000];
    int HiNrs[1000];
    int NThreads=1;
    
#ifdef _WIN32
    sep="\\";
#elif _WIN64
    sep="\\";
#endif

    printf("\n\n\n\t\t\t    Batch Correction\n\t          Written by Hemant Sharma, 1ID,APS\n       In case of questions or errors, email hsharma@anl.gov\n\n\n");
    
    getcwd(FolderName, sizeof(FolderName));
    ParamFN = argv[1];
    for (i=1; i< argc; i++) {
        str = "--ndel";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            DoNdel = atoi(argv[i+1]);
            if (DoNdel == 1) {
                KeepFiles = 1;
            }
        }
        str = "--drk";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            DarkFileStem = argv[i+1];
        }
        str = "--darknr";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            DarkNr = atoi(argv[i+1]);
        }
        str = "--ext";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            ext = argv[i+1];
        }
   //     str = "--noext";
     //   LowNr = strncmp(argv[i], str, strlen(str));
       // if (LowNr == 0) {
         //   NoExt = 1;
       // }
        str = "--inpath";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            FileStem = argv[i+1];
        }
        str = "--indir";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            FolderName = argv[i+1];
        }
        str = "--outpath";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            OutputFilePath = argv[i+1];
        }
        str = "--genum";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            GENumber = atoi(argv[i+1]);
        }
        str = "--bad";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            BadFileName = argv[i+1];
        }
        str = "--help";
        LowNr = strncmp(argv[i], str, strlen(str));
        if (LowNr == 0) {
            printf("Usage: ./BatchCorrection --lo --hi --ndel --drk --darkr --inpath --outpath --genum --fn --bad\n");
            return(0);
        }
    }
    if (GENumber == 1) {
        ext = ".ge1";
    } else if (GENumber == 3) {
        ext = ".ge3";
    } else if (GENumber == 4) {
        ext = ".ge4";
    }
    if (ext == "empty"){
        ext = "";
    }
    
    start0 = clock();

    //Read bad pixel file
    fr = fopen(BadFileName,"rb");
    if (fr==NULL){
        printf("Cannot open file: %s.\n",BadFileName);
    }
    fseek(fr, 0, SEEK_END);
    int sizeBad=ftell(fr);
    rewind(fr);
    fseek(fr, 8192, SEEK_SET);
    unsigned short int *BadFileContents;
    BadFileContents = malloc(Numels * sizeof(*BadFileContents));
    fread(BadFileContents, SizeFile, 1, fr);
    fclose(fr);
    
    //Read dark current values
    sprintf(DarkFileName, "%s%s%s_%05d%s", FolderName, sep, DarkFileStem, DarkNr, ext);
    fz = fopen(DarkFileName,"rb");
    if (fz==NULL){
        printf("Cannot open file: %s.\n",DarkFileName);
    }
    fseek(fz, 0, SEEK_END);
    int sizeDark = ftell(fz);
    rewind(fz);
    int nFramesDark=((sizeDark-8192)/SizeFile);
    printf("Number of frames in dark file %s = %d\n",DarkFileName,nFramesDark);
    fseek(fz, 8192, SEEK_SET);
    unsigned short int *DarkFileContents;
    DarkFileContents = malloc(SizeFile*nFramesDark * sizeof(*DarkFileContents));
    fread(DarkFileContents, SizeFile, nFramesDark, fz);
    fclose(fz);
    Dark = malloc(Numels * sizeof(*Dark));
    for (FrameNr = 0; FrameNr<nFramesDark; FrameNr++){
        for (ElementNr=0; ElementNr<(Numels); ElementNr++){
            Dark[ElementNr] += (((double)(DarkFileContents[(Numels*FrameNr)+ElementNr]))/((double)nFramesDark));
        }
    }
    free(DarkFileContents);
    
    //Read file
        sprintf(FileName, "%s%s%s", FolderName, sep, FileStem);
        fp = fopen(FileName,"rb");
        if (fp==NULL){
            printf("Cannot open file: %s.\n",FileName);
            return 0;
        }
        fseek(fp, 0, SEEK_END);
        size=ftell(fp);
        rewind(fp);
        fseek(fp, 8192, SEEK_SET); // Skip header.
        nFrames=(size/SizeFile);
        printf("Read file: %s\nNumber of frames in file = %i\n",FileName,nFrames);
        SumMat   = malloc(Numels * sizeof(*SumMat));
        AveMat   = malloc(Numels * sizeof(*AveMat));
        sprintf(SumFileName, "%s%s%s%s", OutputFilePath, sep, FileStem, sumext);
        sprintf(AveFileName, "%s%s%s%s", OutputFilePath, sep, FileStem, aveext);

        // Subtract background, check for bad pixels, sum, save files.
        if (KeepFiles==1){ // Will save individual files.
            start = clock();
            for (FrameNr=0; FrameNr<nFrames; FrameNr++){
                SubData  = malloc(Numels * sizeof(*SubData));
                Contents = malloc(Numels * sizeof(*Contents));
                Data     = malloc(Numels * sizeof(*Data));
                fread(Contents, SizeFile, 1, fp);
                // Convert Contents to double
                for (i=0; i<(Numels);i++){
                    Data[i] = (double)Contents[i];
                }
                for (ElementNr=0; ElementNr<(Numels); ElementNr++) {
                    Data[ElementNr] -= Dark[ElementNr];
                    if (BadFileContents[ElementNr]==2 && ElementNr > 2048 && ElementNr < (2047*2048)){
                        Data[ElementNr] = (Data[ElementNr-1] + Data[ElementNr+1] + Data[ElementNr + 2048] + Data[ElementNr - 2048])/(4.0);
                    }
                    if (Data[ElementNr] < 0 || (BadFileContents[ElementNr])%2 == 1) {
                        Data[ElementNr] = 0;
                    }
                    SumMat[ElementNr] += Data[ElementNr];
                    AveMat[ElementNr] += ((Data[ElementNr])/nFrames);
                    SubData[ElementNr] = Data[ElementNr];
                }
                // Write the .cor files.
                sprintf(OutFileName, "%s%s%s.frame%d%s", OutputFilePath, sep, FileStem, FrameNr+1, extcor);
                printf("Now writing corrected file %s, number %d of %d.\n",OutFileName,FrameNr+1,nFrames);
                fb = fopen(OutFileName,"wb");
                fwrite(SubData,SizeFile,1,fb);
                fclose(fb);
                free(SubData);
                free(Contents);
                free(Data);
            }
        } else {
            start = clock();
            for (FrameNr=0; FrameNr<nFrames; FrameNr++){
                Contents = malloc(Numels * sizeof(*Contents));
                Data     = malloc(Numels * sizeof(*Data));
                fread(Contents, SizeFile, 1, fp);
                // Convert Contents to double
                for (i=0; i<(Numels);i++){
                    Data[i] = (double)Contents[i];
                }
                for (ElementNr=0; ElementNr<(Numels); ElementNr++) {
                    Data[ElementNr] -= Dark[ElementNr];
                    if (BadFileContents[ElementNr]==2 && ElementNr > 2048 && ElementNr < (2047*2048)){
                        Data[ElementNr] = (Data[ElementNr-1] + Data[ElementNr+1] + Data[ElementNr + 2048] + Data[ElementNr - 2048])/(4.0);
                    }
                    if (Data[ElementNr] < 0 || (BadFileContents[ElementNr])%2 == 1) {
                        Data[ElementNr] = 0;
                    }
                    SumMat[ElementNr] += Data[ElementNr];
                    AveMat[ElementNr] += ((Data[ElementNr])/nFrames);
                }
                free(Contents);
                free(Data);
            }
        }
        fclose(fp);
        SumMat2 = malloc((Numels)*sizeof(*SumMat2));
        for (ElementNr=0; ElementNr<(Numels); ElementNr++) {
            SumMat2[ElementNr] = SumMat[ElementNr];
        }
        AveMat2 = malloc((Numels)*sizeof(*AveMat2));
        for (ElementNr=0; ElementNr<(Numels); ElementNr++) {
            AveMat2[ElementNr] = AveMat[ElementNr];
        }
        //Write sum and ave files.
        ft = fopen(SumFileName,"wb+");
        if (ft==NULL) {
            printf("Could not make an output file.");
        }
        printf("Writing .sum file %s\n",SumFileName);
        fwrite(SumMat2,SizeFile,1,ft);
        fclose(ft);
        fc = fopen(AveFileName,"wb");
        if (fc==NULL) {
            printf("Could not make an output file.");
        }
        printf("Writing .ave file %s\n",AveFileName);
        fwrite(AveMat2,SizeFile,1,fc);
        fclose(fc);
        free(SumMat);
        free(AveMat);
        free(SumMat2);
        free(AveMat2);
        inter = clock();
        diffhere = ((double)(inter-start))/CLOCKS_PER_SEC;
        printf("Time for this file [s]: %f\n",diffhere);
    
    printf("Done\n");
    end=clock();
    diftotal = ((double)(end-start0))/CLOCKS_PER_SEC;
    printf("\nTotal time elapsed [s] [min]: %f %f\n", diftotal, diftotal/60.0);
    return(0);
}
