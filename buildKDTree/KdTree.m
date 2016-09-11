/****************************************************************************
 *  KdTree.m
 *  Sadilac
 *
 *  Created by Darren on 23/05/2012.
 *  Copyright (c) 2012 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   23/05/2012
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
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
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 ***************************************************************************/

#import "KdTree.h"

@implementation KdTreeEvent
@synthesize ptype=_ptype;
@synthesize dim=_dim;
@synthesize tri=_tri;
@synthesize pos=_pos;

- (void) setType: (unsigned int) tp Dim: (unsigned int) d Tri: (Triangle *) t Pos: (float) p {
    _ptype = tp; _dim=d; _tri=t;_pos=p;
}
- (id) copyWithZone:(NSZone *) zone {
    KdTreeEvent * e = [[KdTreeEvent allocWithZone:zone] init];
    [e setType:_ptype Dim:_dim Tri:_tri Pos:_pos];
    return e ;
}
- (NSString *) description {
    return [NSString stringWithFormat:@"dim:%d Pos:%1.5f Type:%d tri:%@",
                    _dim, _pos, _ptype, _tri ];
}
- (BOOL) isEqual:(id)object {
    return ([self compareEvents:object] == NSOrderedSame );
}
- (NSUInteger) hash {
    NSUInteger n = [[NSNumber numberWithUnsignedInt:_dim] hash];
    NSUInteger m = [[NSNumber numberWithUnsignedInt: _pos] hash];
    NSUInteger o = [[NSNumber numberWithUnsignedInt:_ptype] hash];
    NSUInteger p = [[NSNumber numberWithUnsignedInt:(unsigned int)_tri] hash];
    return (n ^ m ^ o ^ p) ;
}
- (NSComparisonResult) compareEvents: (KdTreeEvent *) e {
    // first dim, then pos, then type, then tri
    if (_dim < [e dim] ) {
        return NSOrderedAscending;
    } else if (_dim > [e dim] ) {
        return NSOrderedDescending;
    }else{
        if (_pos < [e pos]) {
            return NSOrderedAscending;
        }else if (_pos > [e pos] ) {
            return NSOrderedDescending;
        }else{
            if (_ptype < [e ptype]) {
                return NSOrderedAscending;
            } else if ( _ptype > [e ptype] ) {
                return NSOrderedDescending;
            } else {
                if ( _tri < [e tri] ){
                    return NSOrderedAscending;
                } else if ( _tri > [e tri] ) {
                    return NSOrderedDescending;
                } else {
                    return NSOrderedSame;
                }
            }
        }
    }
}
@end

@implementation KdTreeNode
@synthesize isLeaf=_isLeaf;
@synthesize splitPosition=_splitPosition;
@synthesize leftNode=_leftNode;
@synthesize rightNode=_rightNode;
@synthesize axis=_axis;
@synthesize depth=_depth;
@synthesize triangles=_triangles;
@synthesize Events=_Events;
@synthesize cheapSide=_cheapSide;
- (id) init {
    self = [super init];
    if (self) {
        aabb = [[AABB alloc] init];
        _isLeaf = NO;
        _splitPosition = 0.0;
        _leftNode = nil;
        _rightNode = nil;
        _axis = 0;
        _depth = 0;
        _triangles = nil;
        _Events = [NSMutableArray array];
        NSSortDescriptor * dimDescriptor 
        = [NSSortDescriptor sortDescriptorWithKey:@"dim" 
                                        ascending:YES 
                                         selector:@selector(compare:)];
        NSSortDescriptor * posDescriptor 
        = [NSSortDescriptor sortDescriptorWithKey:@"pos"
                                        ascending:YES 
                                         selector:@selector(compare:)];
        NSSortDescriptor * typeDescriptor 
        = [NSSortDescriptor sortDescriptorWithKey:@"ptype" 
                                        ascending:YES 
                                         selector:@selector(compare:)];
        NSSortDescriptor * triDescriptor 
        = [NSSortDescriptor sortDescriptorWithKey:@"tri" 
                                        ascending:YES 
                                         selector:@selector(compareTriangle:)];
        sortDescriptors = [NSArray arrayWithObjects:dimDescriptor, posDescriptor, typeDescriptor, triDescriptor, nil];
    }
    return self;
}
- (id) initWithTriangleArray: (NSArray *) theTris andAABB: (AABB *) v{
    self = [super init];
    if ( self ){
        aabb           = [v copy] ;
        _isLeaf        = NO ;
        _splitPosition = 0.0 ;
        _leftNode      = nil ;
        _rightNode     = nil ;
        _axis          = 0 ;
        _depth         = 0 ;
        _triangles     = [NSArray arrayWithArray:theTris] ;
        for( int t=0; t<[_triangles count]; t++){
            [[_triangles objectAtIndex:t] setTriId:t];
        }
        _Events = [NSArray array];
        NSSortDescriptor * dimDescriptor 
            = [NSSortDescriptor sortDescriptorWithKey:@"dim" 
                                            ascending:YES 
                                             selector:@selector(compare:)];
        NSSortDescriptor * posDescriptor 
            = [NSSortDescriptor sortDescriptorWithKey:@"pos"
                                            ascending:YES 
                                             selector:@selector(compare:)];
        NSSortDescriptor * typeDescriptor 
            = [NSSortDescriptor sortDescriptorWithKey:@"ptype" 
                                            ascending:YES 
                                             selector:@selector(compare:)];
        NSSortDescriptor * triDescriptor 
            = [NSSortDescriptor sortDescriptorWithKey:@"tri" 
                                            ascending:YES 
                                            selector:@selector(compareTriangle:)];
         sortDescriptors = [NSArray arrayWithObjects:dimDescriptor, posDescriptor, typeDescriptor, triDescriptor, nil];
    }
    return self;
}
- (id) KdTreeNodeWithTriangleArray: (NSArray *) triangles andAABB: (AABB *) v{
    return [[KdTreeNode alloc] initWithTriangleArray:triangles andAABB:v];
}
- (void) setAabb:(AABB *)v { aabb = [v copy]; }
- (AABB *) Aabb { return aabb; }

// This is the main control function that recursively subdivides the KdTree
//
- (void) subDivide {
//    printf("\r  Tree depth: %d. Triangles : %lu",[self depth], (unsigned long)[[self triangles] count]);
    [self findBestSplitPlane];  // Sets Axis and plane pos for best split
    if( [self isLeaf] ){        // oh and it also works out if it should terminate
        return;
    }else{
        [self splitTreeAndEvents];
        [_leftNode subDivide];
        [_rightNode subDivide];
    }
    return ;
}

- (void) splitTreeAndEvents {
    KdTreeEvent * event;
    Triangle * tri;

    NSMapTable * leftRightBoth = [[NSMapTable alloc] initWithKeyOptions:NSPointerFunctionsObjectPointerPersonality
                                                           valueOptions:NSPointerFunctionsObjectPointerPersonality
                                                               capacity:[_Events count]];
    
    // classify all triangles to be left of, right of or overlapping
    // the event plane in a single sweep over all events
    
    for ( event in _Events){
        tri = [event tri];
        [leftRightBoth setObject:[NSNumber numberWithInt:BOTH] forKey:tri];
    }

    // Triangles are now classified into LEFT,RIGHT or BOTH

    for (event in _Events){
        if ( ([event ptype] == EVENT_END)
            && ([event dim] == _axis)
            && ([event pos] <= _splitPosition)){
            tri = [event tri];
            [leftRightBoth setObject:[NSNumber numberWithInt:LEFT] forKey:tri];
        }
        else if ( ([event ptype] == EVENT_START)
                 && ([event dim] == _axis)
                 && ([event pos] >= _splitPosition)){
            tri = [event tri];
            [leftRightBoth setObject:[NSNumber numberWithInt:RIGHT] forKey:tri];
        }
        else if ( ([event ptype] == EVENT_PLANAR)
                 && ([event dim] == _axis) ){
            if ( ([event pos] < _splitPosition) 
                || ( ([event pos] == _splitPosition) && (_cheapSide == LEFT))) {
                tri = [event tri];
                [leftRightBoth setObject:[NSNumber numberWithInt:LEFT] forKey:tri];
            }
            if ( ([event pos] > _splitPosition)
                || ( ([event pos] == _splitPosition) && (_cheapSide == RIGHT))){
                tri = [event tri];
                [leftRightBoth setObject:[NSNumber numberWithInt:RIGHT] forKey:tri];
            }
        }
    }
    
    // Triangles are now classified into LEFT,RIGHT or BOTH

    NSMutableSet * Elo = [NSMutableSet setWithCapacity:[_Events count]];
    NSMutableSet * Ero = [NSMutableSet setWithCapacity:[_Events count]];
    NSMutableSet * Ebl = [NSMutableSet setWithCapacity:[_Events count]];
    NSMutableSet * Ebr = [NSMutableSet setWithCapacity:[_Events count]];
    
    for (event in _Events){
        tri = [event tri];

        if ( [[leftRightBoth objectForKey:tri] intValue] == LEFT) {
            [Elo addObject:event];
        }else if ( [[leftRightBoth objectForKey:tri] intValue] == RIGHT){
            [Ero addObject:event];
        }else{
            int k = _axis;
            float pos = _splitPosition;
            NSArray * tArr = [NSArray arrayWithObject:tri];
            AABB * Vl = [aabb CopyAndSplitLeftDim:k Pos:pos];
            NSArray * eArrl = calculateEventsForTris(tArr, Vl);
            [Ebl addObjectsFromArray:eArrl];
            AABB * Vr = [aabb CopyAndSplitRightDim:k Pos:pos];
            NSArray * eArrr = calculateEventsForTris(tArr, Vr);
            [Ebr addObjectsFromArray:eArrr];
        }
    }
    
    [Elo unionSet:Ebl];
    [Ero unionSet:Ebr];
    
    NSArray * El = [Elo sortedArrayUsingDescriptors:sortDescriptors];
    NSArray * Er = [Ero sortedArrayUsingDescriptors:sortDescriptors];

    _leftNode = [self splitLeftWithEvents:El];
    _rightNode = [self splitRightWithEvents:Er];

    return ;
}
- (KdTreeNode *) splitLeftWithEvents: (NSArray *) events {
   
    KdTreeEvent * e;
    NSMutableSet * leftTris = [[NSMutableSet alloc] initWithCapacity:[events count]];
    AABB * v = [aabb CopyAndSplitLeftDim:_axis Pos:_splitPosition];
 
    for ( e in events) [leftTris addObject:[e tri]];

    NSArray * newTris   = [[NSArray alloc] initWithArray:[leftTris allObjects] copyItems:NO];
    NSArray * newEvents = [[NSArray alloc] initWithArray:events copyItems:YES];
    KdTreeNode * node   = [[KdTreeNode alloc] initWithTriangleArray:newTris andAABB:v];
    [node setDepth:_depth+1];
    [node setEvents: newEvents];
    return  node ;
}
- (KdTreeNode *) splitRightWithEvents:(NSArray *) events {
    KdTreeEvent * e;
    NSMutableSet * rightTris = [[NSMutableSet alloc] initWithCapacity:[events count]];
    AABB * v = [aabb CopyAndSplitRightDim:_axis Pos:_splitPosition];
    
    for (e in events) [rightTris addObject:[e tri]];
    
    NSArray * newTris  = [[NSArray alloc] initWithArray:[rightTris allObjects] copyItems:NO];
    NSArray * newEvents = [[NSArray alloc] initWithArray:events copyItems:YES];
    
    KdTreeNode * node = [[KdTreeNode alloc] initWithTriangleArray:newTris andAABB:v];
    [node setDepth:_depth+1];
    [node setEvents:newEvents];
    return node ;
}

// Find best split plane using a plane sweeping algorithm
//
- (void) findBestSplitPlane {
    // Early return if this is a leaf node
    //
    NSUInteger E = [_Events count];
    if ( ([_triangles count] == 0) || (E == 0) || (_depth == MAXDEPTH)) { _isLeaf = YES; return; }
    unsigned long Nl, Np, Nr;
    float p, mincost=9e99, bestPlanePos=0;
    int bestAxis=0, pneg, pplane, ppos, nEvents[3], eventsForDim=0, cheapside=0, mincheapside=0;
    
    // Count the number of events in event list that exist for each dimension k
    //
    KdTreeEvent * e;
    for (int k=0; k<3; k++){
        for ( e in _Events)if ( [e dim] == k)eventsForDim++;
        nEvents[k] = eventsForDim;
    }
    int i = 0;
    for (int k=0; k<3; k++){
        float cost;
        Nl = Np = 0;
        Nr = [_triangles count];
        
        // sweep over all plane candidates for this dim
        //
        while (i<nEvents[k]){
            p = [[_Events objectAtIndex:i] pos];
            pneg = pplane = ppos = 0;

            // if next event is an end event increment the pneg counter and move on
            //
            while ( (i < E ) 
                   && ([[_Events objectAtIndex:i] dim] == k)
                   && ([[_Events objectAtIndex:i] pos] == p)
                   && ([[_Events objectAtIndex:i] ptype] == EVENT_END) ) {
                pneg++;
                i++;
            }
            
            // if next event is a plane event increment the pplane counter and move on
            //
            while ( (i < E )
                   && ([[_Events objectAtIndex:i] dim] == k)
                   && ([[_Events objectAtIndex:i] pos] == p)
                   && ([[_Events objectAtIndex:i] ptype] == EVENT_PLANAR) ){
                pplane++;
                i++;
            }
            
            // if next event is a start event increment the ppos counter and move on
            //
            while ( (i < E )
                   && ([[_Events objectAtIndex:i] dim] == k)
                   && ([[_Events objectAtIndex:i] pos] == p)
                   && ([[_Events objectAtIndex:i] ptype] == EVENT_START) ){
                ppos++;
                i++;
            }
           
            // Now initialise event counts for events in the plane Np and on its right Nr
            //
            Np = pplane;
            Nr -= pplane;
            Nr -= pneg;            
           
            // Calc SAH cost here SAH(V, p, Nl, Nr, Np)
            //
            cost = [self SAHDim:k Pos:p Nleft:Nl Nright:Nr Nplane:Np] ;
            if(cost <0) {
                cost *= -1.;
                cheapside = LEFT;
            }else{
                cheapside = RIGHT;
            }
            if(cost <mincost){
                mincost = cost;
                bestPlanePos = p;
                bestAxis = k;
                if( (p == [[aabb AA] cell:k]) && (Np != 0) ){
                    mincheapside = LEFT;
                }else if( (p == [[aabb BB] cell:k]) && (Np != 0) ){
                    mincheapside = RIGHT;
                }else{
                    mincheapside = cheapside;
                }
            }
            
            // finally update the number of events on the left of the plane
            // After SAH calcs. Nr and Np get updated before the next SAH call
            Nl += ppos;
            Nl += pplane;
        }
    }
    if (mincost > INTERSECTIONCOST * [_triangles count]){
        [self setIsLeaf:YES];           // terminate this branch
    }else{
        _splitPosition = bestPlanePos;  // else set split pos and subdivide
        _axis = bestAxis;
        _cheapSide = mincheapside;

    }
    return;
}
- (float) SAHDim: (NSUInteger) dim Pos: (float) p Nleft: (NSUInteger) Nl Nright: (NSUInteger) Nr Nplane: (NSUInteger) Np {
    AABB * Vl = [aabb CopyAndSplitLeftDim:dim Pos:p];
    AABB * Vr = [aabb CopyAndSplitRightDim:dim Pos:p];
    BOOL maxSpace = YES;
    if ((Np == 0) && ([Vl isPlanarInDimension:dim] || [Vr isPlanarInDimension:dim])){
        maxSpace = NO;
    }
    float ProbL = [Vl SurfaceArea] / [aabb SurfaceArea] ;
    float ProbR = [Vr SurfaceArea] / [aabb SurfaceArea] ;
    float CostL = cost(ProbL, ProbR, Nl+Np, Nr,maxSpace);
    float CostR = cost(ProbL, ProbR, Nl, Np+Nr,maxSpace);
    if (CostL < CostR)
        return -1*CostL; // use minus flag to specify that left side is cheaper
    else
        return CostR;
}
- (NSString *) description {
    NSMutableString * s = [NSMutableString string];
    [s appendFormat:@"[%d]",_depth];
    if(_isLeaf)[s appendString:@"@ "];
    else [s appendString:@"I "];
    switch (_axis) {
        case 0: 
            [s appendString:@"| "];
            break;
        case 1:
            [s appendString:@"- "];
            break;
        case 2:
            [s appendString:@"+ "];
            break;
    }
    [s appendFormat:@"%1.5f ",_splitPosition];
    [s appendFormat:@"%ldT ",[_triangles count]];
    [s appendFormat:@"%@", aabb];
    return s;
}
- (id) copyWithZone:(NSZone *)zone{
    KdTreeNode * kdt = [[KdTreeNode alloc] initWithTriangleArray:_triangles andAABB:aabb];
    [kdt setIsLeaf:_isLeaf];
    [kdt setSplitPosition:_splitPosition];
    [kdt setLeftNode:_leftNode];
    [kdt setRightNode:_rightNode];
    [kdt setAxis:_axis];
    [kdt setDepth:_depth];
    [kdt setAabb:aabb];
    return kdt;
}
- (void) PackNode: (KdTreeNode *) node inArray: (NSMutableArray *) array {
    if ( [node isLeaf] ){
        return;
    }else{
        [array addObject:[node leftNode]];
        [array addObject:[node rightNode]];
        [node PackNode:[node leftNode] inArray:array];
        [node PackNode:[node rightNode] inArray:array];
    }
    return;
}
@end


// Implementation of KdTree
//
@implementation KdTree;
@synthesize packArray=_packArray;
//@synthesize textures=_textures;
- (id) initWithTriangles: ( NSArray * ) triangles {
    self = [super init];
    if( self ){

        AABB * v      = [[AABB alloc] initWithTriangleArray:triangles];
        _root         = [[KdTreeNode alloc] initWithTriangleArray:triangles andAABB:v];
        _packArray    = [NSMutableArray array];
        _triangleData = (triDataStruct *)malloc(sizeof(triDataStruct)*[triangles count]);
        if(_triangleData == NULL){
            NSLog(@"Error: ran out of memory for Triangles");
            exit(1);
        }
        
//        NSMutableArray * materialsFound = [[NSMutableArray alloc] init] ;
//        Triangle *tri ;
//        for (tri in triangles){
//            int triid = [tri matId];
//            int found = 0 ;
//            for(int i=0; i< [materialsFound count]; i++){
//                if( [[materialsFound objectAtIndex:i] materialID] == triid){
//                    found =1 ;
//                }
//            }
//            if (found == 0) {
//                radarMaterial * material = [radarMaterial materialWithID:[tri matId] andName:[tri materialName]];
//                [materialsFound addObject:material];
//            }
//        }
//        _textures = [NSArray arrayWithArray:materialsFound];
    }
    return self;
}

- (void) Build {
    NSArray * tris = [_root triangles];
    AABB * V = [_root Aabb]; 
    NSMutableArray * events;
    
    printf("Calculating events.....\n");
    events = calculateEventsForTris(tris, V);
    printf("  %lu events for this scene\n", [events count]);
    printf("  Sorting Events...");
    [events sortUsingSelector:@selector(compareEvents:)];
    printf("...Done!");
    [_root setEvents: [NSMutableArray arrayWithArray:events]];
    printf("  Starting tree subdivision...");
    [_root subDivide];
    for(int i=0; i<[tris count]; i++){
        triDataStruct data;
        Triangle * tri   = [tris objectAtIndex:i] ;
        TriAccel * acc   = [[TriAccel alloc] initWithTriangle:tri] ;
        data.d           = [acc d];
        data.nd_u        = [acc nd_u];
        data.nd_v        = [acc nd_v];
        data.k           = [acc k];
        data.kbu         = [acc kbu];
        data.kbv         = [acc kbv];
        data.kbd         = [acc kbd];
        data.kcu         = [acc kcu];
        data.kcv         = [acc kcv];
        data.kcd         = [acc kcd];
        data.texture     = [acc texture];
        _triangleData[i] = data ;
    }
    printf("\n");
    
    return;
}
- (void) setRoot:(KdTreeNode *)n {
    _root = [n copy];
    return ;
}
- (NSString *) description {
    
       return [NSString stringWithFormat:@"Tree: %@", _packArray];

}
- (void) PacktoFile: (NSString *) filename {
    double d, nd_u, nd_v, k , kbu , kbv;
    double kbd, kcu, kcv, kcd;
    double AAx,AAy,AAz,BBx,BBy,BBz;
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
    NSInteger ntri,nleaves,nnodes;
    KdTreeNode * node;
    unsigned int flagDimAndOffset;
    float splitPosition;
    int tex;

    FILE *fp;
    fp = fopen([filename UTF8String],"wb");
    if (fp == NULL){
        NSLog(@"Error: failed to open file %@",filename);
        exit(0);
    }
    [_packArray addObject:_root];
    [_root PackNode:_root inArray:_packArray];
    
    /* File structure
     double x, double y, double z - AA for AABB
     double x, double y, double z - BB for AABB
     uint number_of_tris
     for each triangle:
        struct accelerated_tri_structure
     end for
     uint number of leaves
     for each leaf
        uint num_of_tris_in_this_leaf
        for each tri
            uint index of tri in above tri_array
        end for
     end for
     uint number_of_nodes
     for each node
        node structure
     end for
     */
    
    printf("Writing KdTree to file %s...\n",[filename UTF8String]);
    AAx = [[[_root Aabb] AA] cell:0];
    AAy = [[[_root Aabb] AA] cell:1];
    AAz = [[[_root Aabb] AA] cell:2];
    BBx = [[[_root Aabb] BB] cell:0];
    BBy = [[[_root Aabb] BB] cell:1];
    BBz = [[[_root Aabb] BB] cell:2];
    
    fwrite(&AAx,sizeof(double),1,fp);
    fwrite(&AAy,sizeof(double),1,fp);
    fwrite(&AAz,sizeof(double),1,fp);
    fwrite(&BBx,sizeof(double),1,fp);
    fwrite(&BBy,sizeof(double),1,fp);
    fwrite(&BBz,sizeof(double),1,fp);
    
    NSLog(@"Bounding box is %@",[_root Aabb]);
    
    ntri = [[_root triangles] count];
    fwrite(&ntri,sizeof(int),1,fp);
    nleaves = 0;
    Triangle *t;
    
    printf("Packing %ld triangles\n",ntri);
    for (t in [_root triangles]){
               
        TriAccel * tri = [[TriAccel alloc] initWithTriangle:t]; 
        d    = [tri d];
        nd_u = [tri nd_u];
        nd_v = [tri nd_v];
        k    = [tri k];
        kbu  = [tri kbu];
        kbv  = [tri kbv];
        kbd  = [tri kbd];
        kcu  = [tri kcu];
        kcv  = [tri kcv];
        kcd  = [tri kcd];
        tex  = [tri texture];
                
        fwrite(&d,sizeof(double),1,fp);
        fwrite(&nd_u,sizeof(double),1,fp);
        fwrite(&nd_v,sizeof(double),1,fp);
        fwrite(&k,sizeof(double),1,fp);
        fwrite(&kbu,sizeof(double),1,fp);
        fwrite(&kbv,sizeof(double),1,fp);
        fwrite(&kbd,sizeof(double),1,fp);
        fwrite(&kcu,sizeof(double),1,fp);
        fwrite(&kcv,sizeof(double),1,fp);
        fwrite(&kcd,sizeof(double),1,fp);
        fwrite(&tex, sizeof(int), 1, fp);
    }
    
   
    // Calculate the number of leaves and then for each leaf write the index of each triangle in the leaf
    
    for (node in _packArray)if( [node isLeaf] )nleaves++;
    fwrite(&nleaves,sizeof(int),1,fp);
    printf("Packing %ld Leaves\n",nleaves);
    NSMapTable * KVCtriangles = [NSMapTable mapTableWithKeyOptions:NSPointerFunctionsObjectPointerPersonality
                                                      valueOptions:NSPointerFunctionsObjectPointerPersonality];
    unsigned int index=0;
    for (t in [_root triangles])[KVCtriangles setObject:[NSNumber numberWithUnsignedInt:index++] forKey:t];
    unsigned int indexOfTri;
    for (node in _packArray){
        if ([node isLeaf]) {
            ntri = [[node triangles] count];
            fwrite(&ntri,sizeof(int),1,fp);
            for ( t in [node triangles] ){
                indexOfTri = [(NSNumber *)[ KVCtriangles objectForKey:t] unsignedIntValue];
                fwrite(&indexOfTri,sizeof(int),1,fp);
            }
        }
    }

    // Now write the kdTree node info
    
    nnodes = [_packArray count];
    fwrite(&nnodes,sizeof(int),1,fp);
    
    printf("Packing %ld Nodes\n",nnodes);
    NSMapTable * KVCPackArray = [NSMapTable mapTableWithKeyOptions:NSPointerFunctionsObjectPointerPersonality
                                                      valueOptions:NSPointerFunctionsObjectPointerPersonality];
    index=0;
    for (node in _packArray)[KVCPackArray setObject:[NSNumber numberWithUnsignedInt:index++] forKey:node];
    
    int leafOffset = 0;
    unsigned int indexOfNode;
    for (node in _packArray){
        splitPosition = [node splitPosition];
        if ([node isLeaf]){
            int splitDim = [node axis];
            flagDimAndOffset = (unsigned int)1<<31;
            flagDimAndOffset = flagDimAndOffset | ((unsigned int)(leafOffset << 2));
            flagDimAndOffset = flagDimAndOffset | ((unsigned int)(splitDim));
            leafOffset++; 
        }else{
            KdTreeNode * firstChildNode = [node leftNode];
            indexOfNode = [(NSNumber *)[ KVCPackArray objectForKey:firstChildNode] unsignedIntValue];
            int splitDim = [node axis];
            flagDimAndOffset = (unsigned int)0<<31;
            flagDimAndOffset = flagDimAndOffset | ((unsigned int)(indexOfNode << 2 ));
            flagDimAndOffset = flagDimAndOffset | ((unsigned int)(splitDim));
        }

        fwrite(&flagDimAndOffset, sizeof(unsigned int),1,fp);
        fwrite(&splitPosition, sizeof(float),1,fp);
    }
    
    printf("Writing triangle coords\n");
    double Nx,Ny,Nz,area;
    double *gtolmtx, *ltogmtx;
    int mat ;
    for (t in [_root triangles]){
        Ax = [[t vertexAt:0] X] ;
        Ay = [[t vertexAt:0] Y] ;
        Az = [[t vertexAt:0] Z] ;
        Bx = [[t vertexAt:1] X] ;
        By = [[t vertexAt:1] Y] ;
        Bz = [[t vertexAt:1] Z] ;
        Cx = [[t vertexAt:2] X] ;
        Cy = [[t vertexAt:2] Y] ;
        Cz = [[t vertexAt:2] Z] ;
        Nx = [[t normal] X];
        Ny = [[t normal] Y];
        Nz = [[t normal] Z];
        area = [t area] ;
        gtolmtx = [t globalToLocalMatrix];
        ltogmtx = [t localToGlobalMatrix];
        mat = [t matId];
        
        fwrite(&Ax,sizeof(double),1,fp);
        fwrite(&Ay,sizeof(double),1,fp);
        fwrite(&Az,sizeof(double),1,fp);
        fwrite(&Bx,sizeof(double),1,fp);
        fwrite(&By,sizeof(double),1,fp);
        fwrite(&Bz,sizeof(double),1,fp);
        fwrite(&Cx,sizeof(double),1,fp);
        fwrite(&Cy,sizeof(double),1,fp);
        fwrite(&Cz,sizeof(double),1,fp);
        fwrite(&Nx,sizeof(double),1,fp);
        fwrite(&Ny,sizeof(double),1,fp);
        fwrite(&Nz,sizeof(double),1,fp);
        fwrite(&area,sizeof(double),1,fp);
        for (int i=0; i<9; i++)fwrite(&gtolmtx[i],sizeof(double),1,fp);
        for (int i=0; i<9; i++)fwrite(&ltogmtx[i],sizeof(double),1,fp);
        fwrite(&mat, sizeof(int), 1, fp);
    }

    fclose(fp);

}
- (void) printBoxes{
    KdTreeNode * node;
    AABB * V;
    float minx,miny,minz,maxx,maxy,maxz;
    for ( node in _packArray ){
        printf("[%d], %lu triangles\n",[node depth], [[node triangles] count]);
        V = [node Aabb];
        minx = [[V AA] cell:0];
        miny = [[V AA] cell:1];
        minz = [[V AA] cell:2];
        maxx = [[V BB] cell:0];
        maxy = [[V BB] cell:1];
        maxz = [[V BB] cell:2];
        printf(" %f,%f,%f\n",minx,miny,minz);
        printf(" %f,%f,%f\n",maxx,miny,minz);
        printf(" %f,%f,%f\n",maxx,maxy,minz);
        printf(" %f,%f,%f\n",minx,maxy,minz);
        printf(" %f,%f,%f\n",minx,miny,minz);
        printf(" %f,%f,%f\n",minx,miny,maxz);
        printf(" %f,%f,%f\n",maxx,miny,maxz);
        printf(" %f,%f,%f\n",maxx,maxy,maxz);
        printf(" %f,%f,%f\n",minx,maxy,maxz);
        printf(" %f,%f,%f\n",minx,miny,maxz);

    }
}
- (void) printSummary {
    int Nl, Nne,ntri,nt;
    double  Nat, Et, El, Ei, Ct, SAs, x, Kt, Ki;
    KdTreeNode *node; 
    
    Kt = TRAVERSALCOST ;
    Ki = INTERSECTIONCOST ;
    
    if (_packArray == nil) {
        printf("KdTree is not yet packed\n");
    }
    
    SAs = [[_root Aabb] SurfaceArea];
    
    Nl = Nne = ntri = Et = El = Ei = 0;
    for ( node in _packArray){
        Et += ([[node Aabb] SurfaceArea]) / SAs ;
        if ([node isLeaf]){ 
            x = ([[node Aabb] SurfaceArea]) / SAs ;
            El += x;
            Nl++;
            if ( (nt = (int)[[node triangles] count]) != 0) {
                Ei += ( x * nt );
                Nne++;
                ntri += nt;
            }
        }
    
    }
    Nat = (double)ntri / (double)Nne;
    Ct  = (Et * Kt) + (Ei * Ki) ;
    printf("\n            KDTree Summary\n");
    printf("==========================================\n");
    printf("  Triangles                       : %4ld\n",[[_root triangles] count]);
    printf("  Number of Nodes                 : %4ld\n",[_packArray count]);
    printf("  Number of leaves (NL)           : %4d \n",Nl);
    printf("  Number of non-empty leaves      : %4d\n",Nne);
    printf("  Ave triangles / non-empty leaf  : %4.1f\n",Nat);
    printf("  Expected Traversals             : %4.1f\n",Et);
    printf("  Expected leaves visited         : %4.1f\n",El);
    printf("  Expected triangle intersections : %4.1f\n",Ei);
    printf("------------------------------------------\n");
    printf("  Total tree traversal cost       : %4.1f\n\n",Ct);
    
    return ;
    
}

@end

void eventCheck(NSArray * events){
    KdTreeEvent * e;
    unsigned int fastmod[] = {0,1,2,0,1};
    for (e in events){
        if( [e ptype] == EVENT_PLANAR){
            unsigned int k = [e dim];
            Triangle * tri = [e tri];
            for(int v=0; v<3; v++){
                if ([[tri vertexAt:v] cell:k]
                    != [[tri vertexAt:fastmod[v+1]] cell:k] ){
                    NSLog(@"******** Event Error in events!!");
                    NSLog(@"Events: %@",events);
                    NSLog(@"Triangle: %@",tri);
                    NSLog(@"Dimension: %d",k);
                    NSLog(@"Vertex: %d",v);
                    NSLog(@"************************************\n\n");
                    exit(0);
                }
            }
        }
    }
}
float cost(float ProbLeft, float ProbRight, NSUInteger NLeft, NSUInteger NRight,BOOL maxSpace){
    float Kt = TRAVERSALCOST;       // Cost for a tree traversal step
    float Ki = INTERSECTIONCOST;    // Cost for an intersection test
    float lambda; // Parameter used to weight cost in order to bias towards empty V
    
    if ( (NLeft == 0 || NRight == 0) && maxSpace){
        lambda=0.8;
    }
    else lambda = 1.;
    
    return lambda*(Kt + Ki*(ProbLeft*NLeft + ProbRight*NRight)) ;
}
NSMutableArray * calculateEventsForTris( NSArray * tris, AABB * V) {
    // Calculate an event tree for triangles in this node
    NSUInteger T = [tris count];
    Triangle * tri;
    AABB * B;
    KdTreeEvent * e = [[KdTreeEvent alloc] init];
    NSMutableArray * events = [NSMutableArray array];
    unsigned int fastmod[] = {0,1,2,0,1};
    
    for (int k=0; k<3; k++){
        for (int t=0; t<T; t++){
            
            tri = [tris objectAtIndex:t];
            
            if ( [V isPlanarInDimension:k] ){
                if ( [tri isPlanarInDimension:k] 
                    && ((![tri isPlanarInDimension:fastmod[k+1]]) && (![tri isPlanarInDimension:fastmod[k+2]])) ){
                    [e setType:EVENT_PLANAR 
                           Dim:k 
                           Tri:tri 
                           Pos:[[tri vertexAt:0] cell:k] ];
                    [events addObject:[e copy]];
                }
                
            }else{

                B = [V copy];
                [B clipToTriangle:tri];

                if ( ( [B isPlanarInDimension:k] 
                      && ((![B isPlanarInDimension:fastmod[k+1]])  && (![B isPlanarInDimension:fastmod[k+2]])) )
                    && [tri isPlanarInDimension:k] ){
                    [e setType:EVENT_PLANAR 
                           Dim:k 
                           Tri:tri 
                           Pos: [[B AA] cell:k]];
                    [events addObject:[e copy]];
                }else{
                    [e setType:EVENT_START 
                           Dim:k 
                           Tri:tri
                           Pos:[[B AA] cell:k]];
                    [events addObject:[e copy]];
                    [e setType:EVENT_END 
                           Dim:k
                           Tri:tri 
                           Pos:[[B BB] cell:k]];
                    [events addObject:[e copy]];
                }
            }
        }
    }
    return events;
}

@implementation TriAccel
@synthesize d=_d;
@synthesize nd_u=_nd_u;
@synthesize nd_v=_nd_v;
@synthesize k=_k;
@synthesize kbu=_kbu;
@synthesize kbv=_kbv;
@synthesize kbd=_kbd;
@synthesize kcu=_kcu;
@synthesize kcv=_kcv;
@synthesize kcd=_kcd;
@synthesize texture=_texture;
- (id) initWithTriangle: (Triangle *) triangle {
    
    self = [super init];
    if (self){
        
        MVector * Aa = [triangle vertexAt:0];
        MVector * Bb = [triangle vertexAt:1];
        MVector * Cc = [triangle vertexAt:2];
        MVector * N  = [triangle normal];

        // calculate the projection dimension
        float nx = [N X];
        float ny = [N Y];
        float nz = [N Z];
        float anx = ABS(nx);
        float any = ABS(ny);
        float anz = ABS(nz);
        if( anx > any )
        if (anx > anz) _k = 0; /* X */ else _k=2; /* Z */
        else
        if ( any > anz) _k=1; /* Y */ else _k=2; /* Z */
        int u = quickmodulo[_k+1];
        int v = quickmodulo[_k+2];
        
        _nd_u = [N cell:u] / [N cell:_k];
        _nd_v = [N cell:v] / [N cell:_k];
        _d = [Aa dot:[N divide:[N cell:_k]]];
        
        MVector * b = [Cc subtract:Aa];
        MVector * c = [Bb subtract:Aa];
        
        float denom = 1./(([b cell:u]*[c cell:v]) - ([b cell:v]*[c cell:u]));
        
        _kbu     = -[b cell:v] * denom;
        _kbv     =  [b cell:u] * denom;
        _kbd     =  (([b cell:v] * [Aa cell:u]) - ([b cell:u] * [Aa cell:v])) * denom;
        _kcu     =  [c cell:v] * denom;
        _kcv     = -[c cell:u] * denom;
        _kcd     =  (([c cell:u] * [Aa cell:v]) - ([c cell:v] * [Aa cell:u])) * denom;
        _texture = [triangle matId];
        
    }
    return self;
}

@end

