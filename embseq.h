//
//  EMBSEQ.hpp
//  Graphs LDA
//
//  Created by Scott Cramer on 8/16/16.
//  Copyright © 2016 Scott Cramer. All rights reserved.
//

#ifndef EMBSEQ_h
#define EMBSEQ_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include <QGraphicsView>
#include <QVector>
#include <QTextStream>

#include "emb.h"

#define jASCII 106
#define LEFTP 40
#define RIGHTP 41
#define starASCII 42


int outputProgress(int cur, int tot, int opt);

int CompareTwoStrings(char* str1, char* str2, int nstart, int nend);
void *FillGap(void *ta);

//#define MAXEMBS 10000
#define MAXSQRS 15
#define MAXSTRLEN 1000
#define MAXLEFTDIV 100

class EMBSEQ {
public:
    int len; //number of nodes/terms
    EMB *embs; //list of nodes
    int *eord; //order index of nodes in the < ordering, so eord[0] = index of least node
    int *eordi; //order lookup eordi[n] = what order number the index n node is
    int invcomp; //whether the order lookup has been computed
    int MAXEMBS; //max number of nodes
    int *forceq1; //the queue of nodes to add, first component
    int *forceq2; //the queue of nodes to add, second component
    int *fscore; //scores of nodes is queue -- says how far removed they are from what we really want
    int nforce; //number of nodes in queue
    int nthreads; // not used

    int numiterations;
    int ILnum;


    EMBSEQ(int max);
    ~EMBSEQ();

    EMB* GetEmbs();
    int* GetEord();
    int* GetEordi();
    int GetLength();

    int InsertEmb(int ind, int num);
    int CreateILSeq(int num);
    int CreateILSeqModified(int num, int sqrnum);
    void AddToRngSeq(int ind, int start, int finish);
    void AddToNotRngSeq(int ind, int start, int finish);
    void AddToNotRngSeqOrd(int ind, int start, int finish);
    int WriteDiagram(const char* file, int opt);
    int WriteCreationAni(const char* file);
    int InsertPullback(int ind, int rind);
    int FindMaxInRng(int ind1, int ind2);
    int FindMinNotInRng( int ind1, int ind2);
    void AddSqrSeqRng (int ind);
    int FindOrderInd(int ind);
    int FindOrderIndComp(int ind);
    void ComputeScores();
    void InsertAllPullbacks();
    void InsertAllPullbacksAbove(int start, int end);
    void ComputeOrderInverse();
    void SortRng(int ind);

    int PushForward(int ind, int tind);
    int FindImage (int ind, int tind);
    int FindImage (int ind, int tind, int start, int end);
    void PushForwardBelow();
    void PushForwardBelow(int start, int end);
    void PushForwardBelowPerc(int start, int end, int perc);
    void PushForwardSeq(int above, int push);

    int ComputePushForward(int ind, int tind, int iters);
    void ComputeAllPushForwards();
    void ComputeAllPushForwardsCheck();
    void ComputeAllPushForwardsCheck(int rank);

    void SmartAlgo(int opt);
    void RecheckEmbs();
    int GetLeastCommonSqr(int ind, int tind);
    int GetSqrSeq(int ind, int* seq);
    int GetEmbRank(int ind);

    int GetSuccPreimage(int ind);
    int FindLeastSentAbove(int ind, int rind);
    int FindLeastImgAbove(int ind, int tind);
    void SetRngChangedBelow(int ind);

    int WriteSeq (const char* file);
    int LoadSeq (const char* file);

    int WillEverCompute(int ind, int tind);
    int WillEverCompute(int ind, int tind, int rank);
    void RemoveOuterParenth(char* str);
    void ProcessString(char* str);
    int CompFromNotation (char* str);
    int ForceCompFromNotation(char* str, int* hullv, int* hullnum);
    void ParseString(char* str, char* p1, char* p2);

    void ComputeImages(int ind);

    int CompareFromNotation (char* str1, char* str2);
    int FindXCoord( int ind);
    int GetPowerEqual(int ind, int pwr);

    void GetNodeString(char* str, int ind, int depth, int opt);
    void GetNodeString(char* str, int ind);
    int FindLeftDivision(int ind1, int ind2, int* res);
    void PrintLeftDivision(int ind1, int ind2);
    int FindLeftDivision(int ind1, int ind2, int* res, int depth);
    void PrepNodeStrings();
    void SetGenSqrStrs();
    void SetRightSqrStrs();
    void ReplaceEmb(int newemb, int oldemb);
    void ReplaceAllInst(int newi, int oldi);

    int ForceCompPushForward(int ind, int tind, int score);
    int ForcePushForward( int ind, int tind, int score);
    void AddToForceQ(int ind, int tind, int score);
    int ForceCompute(int ind, int tind, int *hullv, int *hullnum);
    void DoPullbacks();
    int WeedForceComp();
    void FillGaps(int opt, int threadnum);
    int CheckForPushForward(int ind, int tind);
    void PerformPushForwardChecks(int ind, int tind);
    void RunGapFill(int ind, int tind);

    void MakeHull(int* hullv, int vlen, int* imgs);
    void RestrictEmbs(int* keep, int* imgs);
    int BasicRestrict(int hull1, int* hullv, int* hullnum, int num);
    void SortForceQ();
    void SortForceQ(int start, int end);
    void MergeSortForceQ(int s1, int e1, int e2);

    void AddSqrRtsToHull(int* hullv, int* hullnum, int ind);

    void SearchForPattern();

};


struct THREADARG{
    EMBSEQ *seq;
    int i;
    int i2;

    THREADARG( EMBSEQ* s, int ind, int tind)
    {
        seq = s;
        i = ind;
        i2 = tind;
    }
};

#endif /* EMBSEQ_h */
