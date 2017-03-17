//
//  buildRopesAndBoxes.cpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "buildRopesAndBoxes.hpp"
#include "buildTree.hpp"

void BuildRopesAndBoxes(kdTree::KdData * Node, int *RS, AABB aabb, kdTree::KdData * KdTree){
    
    int Sl, Sr; //SplitPosLeft, SPlitPosRight;
    float V;
    int RSLeft[6], RSRight[6]; // Ropes for left and right side
    AABB aabbLeft, aabbRight;
    kdTree::KdData * Nr, * Nl;
    
    if( KDT_ISLEAF(Node) ){
        for(int i=0; i<3; i++){
            Node->leaf.aabb.AA.cell[i] = aabb.AA.cell[i];
            Node->leaf.aabb.BB.cell[i] = aabb.BB.cell[i];
            Node->leaf.Ropes[2*i] = RS[2*i];
            Node->leaf.Ropes[2*i+1] = RS[2*i+1];
        }
    } else {
        for(int i=0; i<3; i++){
            Node->brch.aabb.AA.cell[i] = aabb.AA.cell[i];
            Node->brch.aabb.BB.cell[i] = aabb.BB.cell[i];
        }
        
        Sl = KDT_DIMENSION(Node) * 2;
        Sr = KDT_DIMENSION(Node) * 2 + 1;
        V  = Node->brch.splitPosition;
        for(int i=0; i<6; i++)RSLeft[i]=RSRight[i]=RS[i];
        
        RSLeft[Sr] = KDT_RGHTCHILD(Node); // Right child for Node
        aabbLeft = aabb;
        aabbLeft.BB.cell[KDT_DIMENSION(Node)] = V;
        Nl = &(KdTree[KDT_LEFTCHILD(Node)]);
        
        BuildRopesAndBoxes(Nl, RSLeft, aabbLeft, KdTree);
        
        RSRight[Sl] = KDT_LEFTCHILD(Node); // Left child for Node
        aabbRight = aabb;
        aabbRight.AA.cell[KDT_DIMENSION(Node)] = V;
        Nr = &(KdTree[KDT_RGHTCHILD(Node)]);
        
        BuildRopesAndBoxes(Nr, RSRight, aabbRight, KdTree);
        
    }
    
    return ;
}
