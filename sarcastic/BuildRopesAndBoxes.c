//
//  BuildRopesAndBoxes.cpp
//  Sadilac
//
//  Created by Darren on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//

#include "BuildRopesAndBoxes.h"

void BuildRopesAndBoxes(KdData * Node, int *RS, AABB aabb, KdData * KdTree){
    
    int Sl, Sr; //SplitPosLeft, SPlitPosRight;
    float V;
    int RSLeft[6], RSRight[6]; // Ropes for left and right side
    AABB aabbLeft, aabbRight;
    KdData * Nr, * Nl;
    
    if( KDT_ISLEAF(Node) ){
        for(int i=0; i<3; i++){
            Node->leaf.aabb.AA.cell[i] = aabb.AA.cell[i];
            Node->leaf.aabb.BB.cell[i] = aabb.BB.cell[i];
            Node->leaf.Ropes[2*i] = RS[2*i];
            Node->leaf.Ropes[2*i+1] = RS[2*i+1];
        }
    } else {
        for(int i=0; i<3; i++){
            Node->branch.aabb.AA.cell[i] = aabb.AA.cell[i];
            Node->branch.aabb.BB.cell[i] = aabb.BB.cell[i];
        }
        
        Sl = KDT_DIMENSION(Node) * 2;
        Sr = KDT_DIMENSION(Node) * 2 + 1;
        V  = Node->branch.splitPosition;
        for(int i=0; i<6; i++)RSLeft[i]=RSRight[i]=RS[i];
        
        RSLeft[Sr] = KDT_OFFSET(Node)+1; // Right child for Node
        aabbLeft = aabb;
        aabbLeft.BB.cell[KDT_DIMENSION(Node)] = V;
        Nl = &(KdTree[KDT_OFFSET(Node)]);
        
        BuildRopesAndBoxes(Nl, RSLeft, aabbLeft, KdTree);
        
        RSRight[Sl] = KDT_OFFSET(Node); // Left child for Node
        aabbRight = aabb;
        aabbRight.AA.cell[KDT_DIMENSION(Node)] = V;
        Nr = &(KdTree[KDT_OFFSET(Node)+1]);
        
        BuildRopesAndBoxes(Nr, RSRight, aabbRight, KdTree);
        
    }
    
    return ;
}
