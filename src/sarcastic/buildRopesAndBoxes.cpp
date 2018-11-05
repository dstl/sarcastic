/***************************************************************************
 * 
 *           Module :  buildRopesAndBoxes.cpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 15/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Simple routine to calculate the AABB around a node in a KdTree
 *      And then attach 'ropes' from leaves at the bottom of the tree to 
 *      adjoining leaves. This means that a search algorithm doesnt have to 
 *      climb back up the tree to go down a new branch
 *
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 ***************************************************************************/

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
