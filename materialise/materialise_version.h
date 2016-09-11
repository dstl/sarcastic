//
//  materialise_version.h
//  sarcastic
//
//  Created by Darren on 29/04/2016.
//  Copyright (c) 2016 Dstl. All rights reserved.
//

#ifndef sarcastic_materialise_version_h
#define sarcastic_materialise_version_h

#ifdef MATERIALISE_REVISION
#undef MATERIALISE_REVISION
#endif
#ifdef MATERIALISE_SHORT_VERSION
#undef MATERIALISE_SHORT_VERSION
#endif
#ifdef MATERIALISE_FULL_VERSION
#undef MATERIALISE_FULL_VERSION
#endif
#ifdef MATERIALISE_VERSION_DATE
#undef MATERIALISE_VERSION_DATE
#endif
#define MATERIALISE_REVISION "125"
#define MATERIALISE_SHORT_VERSION "2.1"
#define MATERIALISE_FULL_VERSION "2.1-3-8d7fc27-dirty"
#define MATERIALISE_VERSION_DATE "2016-03-07 16:41:41 +0000"

#endif
