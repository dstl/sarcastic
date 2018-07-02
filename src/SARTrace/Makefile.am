bin_PROGRAMS = SARTrace
SARTrace_SOURCES = ../fastKdTree/AABB.cpp ../fastKdTree/clipToAABB.cpp ../fastKdTree/rayTrace.cpp ../fastKdTree/kdTreeNode.cpp ../fastKdTree/splitCandidate.cpp ../Common/tinyxml/tinystr.cpp ../Common/tinyxml/tinyxml.cpp ../Common/tinyxml/tinyxmlerror.cpp ../Common/tinyxml/tinyxmlparser.cpp ../fastKdTree/buildTree.cpp ../colladaToPlyFile/colladainterface.cpp ../Common/TriangleMesh.cpp ../sarcastic/buildRays.cpp ../sarcastic/accelerateTriangles.cpp ../sarcastic/buildRopesAndBoxes.cpp ../bircs/buildKernel.cpp ../sarcastic/reflect.cpp getSARTraceUserInput.cpp SARTrace.cpp SARTracebanner.cpp traceThreadCore.cpp ../Common/readMaterialFile.cpp ../Common/matrixMultiplication.c ../Common/OpenCLUtils.c ../sarcastic/ecef2SceneCoords.c ../Common/pmath_rng.c ../Common/boxMullerRandom.c ../Common/ranf.c

SARTrace_CPPFLAGS = -I.. -I../Common -I../fastKdTree -I../colladaToPlyFile -I../sarcastic -I../bircs -DOCLKERNELSPATH=$(abs_top_srcdir)/src/kernels
SARTrace_CXXFLAGS = -std=gnu++11

AM_CXXFLAGS = @CL_CFLAGS@
SARTrace_LDADD = ../sarclib/libsarclib.a @CL_LIBS@

