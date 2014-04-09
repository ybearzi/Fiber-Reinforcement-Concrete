CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT( concrete )

SET(CMAKE_AUTOMOC ON)

#Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

IF(CMAKE_COMPILER_IS_GNUCXX)
	ADD_DEFINITIONS(-fpermissive -std=gnu++0x)
ENDIF()

#Find Qt
#FIND_PACKAGE(Qt4 REQUIRED)
#INCLUDE(${QT_USE_FILE})

# Add Qt files
#ADD_DEFINITIONS(${QT_DEFINITIONS})

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

#SET( My_sources main.cxx )
#SET( My_hearders )
#SET( My_forms )
#SET( My_RCC )

#QT4_ADD_RESOURCES(My_RCC_srcs ${My_RCC})
#QT4_WRAP_UI( My_forms_headers ${My_forms} )
#QT4_WRAP_CPP( My_headers_moc ${My_headers} )

SET( OtherHeaders itkBinaryThinningImageFilter3D.h itkSkeletonLabeledGraphImageFilter.h itkLineShapeImageFilter.h itkBiggestComponentImageFilter.h itkSaliencyMeasureBaseImageFilter.h itkFiberSaliencyImageFilter.h Graph/itkGraphSource.h itkSkeletonToGraphFilter.h Graph/itkDefaultGraphTraits.h itkLineGraphTraits.h Graph/itkGraph.h itkMarkovChainLineGraphFilter.h itkMarkovChainGraphFilter.h Graph/itkGraphDuplicator.h itkGraphToSkeletonImageFilter.h Graph/itkGraphToImageFilter.h itkSplitAndMergeGraphFilter.h Graph/itkInPlaceGraphFilter.h itkSkeletonGraph.h itkFiberExtractionGraphFilter.h ) # All files not using Qt
SET( OtherSources main.cxx itkBinaryThinningImageFilter3D.txx itkSkeletonLabeledGraphImageFilter.hxx itkLineShapeImageFilter.hxx itkBiggestComponentImageFilter.hxx itkSaliencyMeasureBaseImageFilter.hxx itkFiberSaliencyImageFilter.hxx Graph/itkGraphSource.txx Graph/itkGraph.txx itkMarkovChainLineGraphFilter.hxx itkMarkovChainGraphFilter.hxx Graph/itkGraphDuplicator.txx itkGraphToSkeletonImageFilter.hxx Graph/itkGraphToImageFilter.txx itkSplitAndMergeGraphFilter.hxx Graph/itkInPlaceGraphFilter.txx itkSkeletonGraph.hxx itkFiberExtractionGraphFilter.hxx )


#ADD_EXECUTABLE(skull_stripping ${My_sources} ${My_headers_moc} ${My_forms_headers} ${My_RCC_srcs} ${OtherHeaders} ${OtherSources} )
ADD_EXECUTABLE(concrete ${OtherHeaders} ${OtherSources} )

#TARGET_LINK_LIBRARIES(concrete ${QT_LIBRARIES} ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(concrete ${ITK_LIBRARIES})
