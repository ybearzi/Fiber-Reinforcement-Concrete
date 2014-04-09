#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <stack>


	/*
	 * Image
	 */

#include "itkImage.h"


const unsigned int Dim = 3;
typedef float RealType;

typedef unsigned short 					PixelType;
typedef itk::Image< PixelType, Dim > 	ImageType;
typedef itk::Image< RealType, Dim > 	RealImageType;
typedef itk::ImageRegion< Dim >			ImageRegionType;

const PixelType MaxPixelValue = itk::NumericTraits< PixelType >::max();
const PixelType MinPixelValue = itk::NumericTraits< PixelType >::min();

#include "itkSkeletonGraph.h"
#include "itkLineGraphTraits.h"

typedef itk::LineGraphTraits< Dim >				GraphTraitsType;
typedef itk::SkeletonGraph< GraphTraitsType >	GraphType;

	/*
	 * Inputs and Outputs
	 */

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileReader< RealImageType > RealReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

	/*
 	 * Iterators
 	 */


	/*
	 * Filters
	 */

#include "itkBinaryThinningImageFilter3D.h"
#include "itkLineShapeImageFilter.h"
#include "itkBiggestComponentImageFilter.h"
#include "itkSaliencyMeasureBaseImageFilter.h"
#include "itkFiberSaliencyImageFilter.h"
//#include "itkSkeletonLabeledGraphImageFilter.h"
#include "itkSkeletonToGraphFilter.h"
#include "itkMarkovChainLineGraphFilter.h"
#include "itkGraphToSkeletonImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSplitAndMergeGraphFilter.h"
#include "itkFiberExtractionGraphFilter.h"

typedef itk::LineShapeImageFilter< ImageType > 							LineShapeFilterType;
typedef itk::BiggestComponentImageFilter< ImageType >					BiggestComponentFilterType;
typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType >		BinaryThinningFilterType;
typedef itk::SaliencyMeasureBaseImageFilter< ImageType, RealImageType >	SaliencyFilterType;
typedef itk::FiberSaliencyImageFilter< ImageType, RealImageType >		FiberSaliencyFilterType;
//typedef itk::SkeletonLabeledGraphImageFilter< ImageType, ImageType >	GraphFilterType;

typedef itk::ImageFileWriter< RealImageType >							RealWriterType;

typedef itk::SkeletonToGraphFilter< ImageType, GraphType, RealImageType, ImageType >	SkeletonToGraphFilterType;
typedef itk::MarkovChainLineGraphFilter< GraphType >					MarkovChainGraphFilterType;
typedef itk::GraphToSkeletonImageFilter< GraphType, RealImageType >			GraphToSkeletonFilterType;
typedef itk::SplitAndMergeGraphFilter< GraphType >						SplitAndMergeFilterType;
typedef itk::FiberExtractionGraphFilter< GraphType >					FiberExtractionFilterType;

typedef itk::MultiplyImageFilter< ImageType, RealImageType, RealImageType > MultiplyFilterType;

//void histoMap( GraphFilterType::RealListMapType phi, GraphFilterType::RealListMapType thet, GraphFilterType::LabelImagePixelType max );
void histoLabelMap( BiggestComponentFilterType::LabelMapType* labelMap, BiggestComponentFilterType::LabelMapType* labelMapContour, BiggestComponentFilterType::LabelMapType* labelMapSkeleton );
void histoGraphMap( GraphType::Pointer graph );


#include <vnl/vnl_rank.h>

int main( int argc, char* argv[] )
{
	/*
	 * Checking input parameters
	 */

	if( argc < 2 )
	{
		std::cout << "Error : Not enough input arguments" << std::endl;
		return EXIT_FAILURE;
	}

	itk::Matrix< double, 2, 2 > M;
	M(0,0) = 1;
	M(0,1) = 0.99;
	M(1,0) = M(0,1);
	M(1,1) = 1;
	std::cout << " M " << M << std::endl;
	std::cout << vnl_svd< double >( M.GetVnlMatrix() ).pinverse() << std::endl;
	std::cout << vnl_svd< double >( M.GetVnlMatrix() ).pinverse( 0 ) << std::endl;
	std::cout << vnl_svd< double >( M.GetVnlMatrix() ).pinverse( 1 ) << std::endl;
	std::cout << vnl_svd< double >( M.GetVnlMatrix() ).pinverse( 2 ) << std::endl;
	std::cout << vnl_svd< double >( M.GetVnlMatrix() ).pinverse( vnl_rank< double >( M.GetVnlMatrix() ) ) << std::endl;
	std::cout << M(0,1) * M(0,1) * vnl_svd< double >( M.GetVnlMatrix() ).pinverse( 1 ) + ( 1 - M(0,1)*M(0,1) ) * vnl_svd< double >( M.GetVnlMatrix() ).pinverse( 2 )<< std::endl;

	std::list< unsigned int > a, b, afbf, abbf, afbb, abbb;

	a.push_back( 0 );
	a.push_back( 1 );
	a.push_back( 2 );
	a.push_back( 3 );

	b.push_back( 5 );
	b.push_back( 6 );
	b.push_back( 7 );
	b.push_back( 8 );

	std::copy( a.begin(), a.end(), std::back_inserter( afbf ) );
	std::copy( a.begin(), a.end(), std::back_inserter( abbf ) );
	std::copy( a.begin(), a.end(), std::back_inserter( afbb ) );
	std::copy( a.begin(), a.end(), std::back_inserter( abbb ) );

	std::copy( b.begin(), b.end(), std::front_inserter( afbf ) );
	std::copy( b.rbegin(), b.rend(), std::front_inserter( afbb ) );
	std::copy( b.begin(), b.end(), std::back_inserter( abbf ) );
	std::copy( b.rbegin(), b.rend(), std::back_inserter( abbb ) );

	std::cout << " afbf ";
	for( std::list< unsigned int >::iterator it = afbf.begin(); it != afbf.end(); ++it )
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << " afbb ";
	for( std::list< unsigned int >::iterator it = afbb.begin(); it != afbb.end(); ++it )
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << " abbb ";
	for( std::list< unsigned int >::iterator it = abbb.begin(); it != abbb.end(); ++it )
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;
	std::cout << " abbf ";
	for( std::list< unsigned int >::iterator it = abbf.begin(); it != abbf.end(); ++it )
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;



	std::string str = std::string( argv[1] );

	WriterType::Pointer writer = WriterType::New();
	RealWriterType::Pointer realWriter = RealWriterType::New();

	/*
	 *	Reading input image
	 */

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( str + ".nrrd" );
	reader->Update();

//	RealReaderType::Pointer reader2 = RealReaderType::New();
//	reader2->SetFileName( str + "_vesselness.nrrd" );
//	reader2->Update();

	std::cout << " ***** LineShapeFilter " << std::endl;

	LineShapeFilterType::Pointer lineShapeFilter = LineShapeFilterType::New();

	lineShapeFilter->SetInput( reader->GetOutput() );
	lineShapeFilter->SetSigma( 2.0f );
	lineShapeFilter->SetExtractBrightLine( false );
	lineShapeFilter->EigenValuesExtractionOn();
	lineShapeFilter->LabelImageOn();
	lineShapeFilter->Update();

	writer->SetInput( lineShapeFilter->GetBinaryOutput() );
	writer->SetFileName( str + "_binary.nrrd" );
	writer->Update();

	writer->SetInput( lineShapeFilter->GetLabelOutput() );
	writer->SetFileName( str + "_label.nrrd" );
	writer->Update();

	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 0 ) );
	realWriter->SetFileName( str + "_eigen0.nrrd" );
	realWriter->Update();

	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 1 ) );
	realWriter->SetFileName( str + "_eigen1.nrrd" );
	realWriter->Update();

	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 2 ) );
	realWriter->SetFileName( str + "_eigen2.nrrd" );
	realWriter->Update();
	
	realWriter->SetInput( lineShapeFilter->GetVesselnessOutput() );
	realWriter->SetFileName( str + "_home_vesselness.nrrd" );
	realWriter->Update();
/*
	std::cout << " ***** Biggest Component " << std::endl;

	BiggestComponentFilterType::Pointer biggestComponentFilter = BiggestComponentFilterType::New();
	biggestComponentFilter->SetInput( reader->GetOutput() );
	biggestComponentFilter->SetMaskInput( lineShapeFilter->GetBinaryOutput() );
	biggestComponentFilter->Update();

	realWriter->SetInput( biggestComponentFilter->GetMetricOutput() );
	realWriter->SetFileName( str + "_SR_V.nrrd" );
	realWriter->Update();

	writer->SetInput( biggestComponentFilter->GetLabelOutput() );
	writer->SetFileName( str + "_connected_label.nrrd" );
	writer->Update();

	writer->SetInput( biggestComponentFilter->GetLabelContourOutput() );
	writer->SetFileName( str + "_connected_contour_label.nrrd" );
	writer->Update();

	writer->SetInput( biggestComponentFilter->GetLabelSkeletonOutput() );
	writer->SetFileName( str + "_connected_skeleton_label.nrrd" );
	writer->Update();


	histoLabelMap( biggestComponentFilter->GetLabelMapOutput(), biggestComponentFilter->GetLabelMapContourOutput(), biggestComponentFilter->GetLabelMapSkeletonOutput() );
*/
	std::cout << " ***** Skeletonization " << std::endl;

	BinaryThinningFilterType::Pointer thinningFilter = BinaryThinningFilterType::New();
	thinningFilter->SetInput( lineShapeFilter->GetOutput() );
	thinningFilter->Update();

	writer->SetInput( thinningFilter->GetOutput() );
	writer->SetFileName( str + "_skeleton.nrrd" );
	writer->Update();

	std::cout << " ***** Vesselness Skeleton " << std::endl;

	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
	multiplyFilter->SetInput1( thinningFilter->GetOutput() );
	multiplyFilter->SetInput2( lineShapeFilter->GetVesselnessOutput() );
	multiplyFilter->Update();

	realWriter->SetInput( multiplyFilter->GetOutput() );
	realWriter->SetFileName( str + "_skeleton_vesselness.nrrd" );
	realWriter->Update();

	std::cout << " ***** Skeleton To Graph Filter " << std::endl;

	SkeletonToGraphFilterType::Pointer skeletonToGraphFilter = SkeletonToGraphFilterType::New();
	skeletonToGraphFilter->SetInput( thinningFilter->GetOutput() );
	skeletonToGraphFilter->SetVesselnessInput( lineShapeFilter->GetVesselnessOutput() );
	//skeletonToGraphFilter->SetInput( reader->GetOutput() );
	//skeletonToGraphFilter->SetVesselnessInput( reader2->GetOutput() );
	skeletonToGraphFilter->Update();

	itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::Pointer iwriter = itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::New();
	iwriter->SetInput( skeletonToGraphFilter->GetIdentifierOutput() );
	iwriter->SetFileName( str + "_identifier.nrrd" );
	iwriter->Update();

	std::cout << " ***** Fiber extraction Graph Filter  " << std::endl;

	FiberExtractionFilterType::Pointer fiberExtractionFilter = FiberExtractionFilterType::New();
	fiberExtractionFilter->SetInput( skeletonToGraphFilter->GetOutput() );
	fiberExtractionFilter->Update();

/*	MarkovChainGraphFilterType::Pointer markovFilter = MarkovChainGraphFilterType::New();
	markovFilter->SetInput( skeletonToGraphFilter->GetOutput() );
	markovFilter->SetAlpha( 0.01 );
	markovFilter->SetBeta( 1 );
	markovFilter->SetGamma( 1 );
	markovFilter->Update();
*/
	GraphToSkeletonFilterType::Pointer graphToSkeletonFilter = GraphToSkeletonFilterType::New();
	graphToSkeletonFilter->SetInput( fiberExtractionFilter->GetOutput() );
	graphToSkeletonFilter->SetRegion( reader->GetOutput()->GetLargestPossibleRegion() );
	graphToSkeletonFilter->Update();

	std::cout << " Writer " << std::endl;

	realWriter->SetInput( graphToSkeletonFilter->GetOutput() );
	realWriter->SetFileName( str + "_fiber_extraction.nrrd" );
	realWriter->Update();
//	histoGraphMap( markovFilter->GetOutput() );



//	std::cout << " adress " << markovFilter << std::endl;
/*	std::cout << " size " << markovFilter->GetOutput()->GetNodeContainer()->Size() << std::endl;

	std::cout << " ***** Split and Merge GraphFilter " << std::endl;

	SplitAndMergeFilterType::Pointer splitAndMergeFilter = SplitAndMergeFilterType::New();
	splitAndMergeFilter->SetInput( markovFilter->GetOutput() );
	splitAndMergeFilter->Update();

	std::cout << " ***** Graph To Skeleton Image Filter " << std::endl;

	GraphToSkeletonFilterType::Pointer graphToSkeletonFilter2 = GraphToSkeletonFilterType::New();
	graphToSkeletonFilter2->SetInput( splitAndMergeFilter->GetOutput() );
	graphToSkeletonFilter2->SetRegion( reader->GetOutput()->GetLargestPossibleRegion() );
	graphToSkeletonFilter2->Update();



	
	std::cout << " Writer " << std::endl;

	realWriter->SetInput( graphToSkeletonFilter2->GetOutput() );
	realWriter->SetFileName( str + "_merge.nrrd" );
	realWriter->Update();
	histoGraphMap( markovFilter->GetOutput() );
*/
	std::cout << " done " << std::endl;

/*	std::cout << " ***** Saliency " << std::endl;

	FiberSaliencyFilterType::Pointer saliencyFilter = FiberSaliencyFilterType::New();
	saliencyFilter->SetInput( reader->GetOutput() );
	saliencyFilter->SetNumberOfIteration( 150 );
	saliencyFilter->SetNumberOfThreads( 4 );
	saliencyFilter->Update();

	std::cout << " ***** Making graph " << std::endl;
*/	
/*	GraphFilterType::Pointer graphFilter = GraphFilterType::New();
	graphFilter->SetInput( thinningFilter->GetOutput() );
	graphFilter->Update();

	GraphFilterType::MapType map = graphFilter->GetCostMap();
	GraphFilterType::LabelImagePixelType max = graphFilter->GetMaxLabel();

	histoMap( graphFilter->GetPhiListMap(), graphFilter->GetThetaListMap(), max );

 	writer->SetInput( graphFilter->GetLabelOutput() );
	writer->SetFileName( str + "_skeleton_label.nrrd" );
	writer->Update();
*/
	return EXIT_SUCCESS;
}

void histoGraphMap( GraphType::Pointer graph )
{
	std::ofstream file( "graph.m", std::ios::out | std::ios::trunc );

	file << "P = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << it.Value().Point << "; ";
	}

	file << "];" << std::endl << "V = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << it.Value().Vector << "; ";
	}

	file << "];" << std::endl << "D = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << std::abs( it.Value().AverageDistance - it.Value().TheoricalAverageDistance ) << "; ";
	}

	file << "];" << std::endl << "W = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << it.Value().Weight << "; ";
	}

	file << "];" << std::endl << "Vess = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << it.Value().MeanVesselness << "; ";
	}

	file << "];" << std::endl << "F = [";

	for( GraphType::NodeContainerType::Iterator it = graph->GetNodeContainer()->Begin();
			it != graph->GetNodeContainer()->End();
			++it )
	{
		file << " " << it.Value().Force << "; ";
	}



	file << "];" << std::endl;



	file.close();

}

void histoLabelMap( BiggestComponentFilterType::LabelMapType* labelMap, BiggestComponentFilterType::LabelMapType* labelMapContour, BiggestComponentFilterType::LabelMapType* labelMapSkeleton )
{
	BiggestComponentFilterType::LabelMapType::SizeValueType max = labelMap->GetNumberOfLabelObjects();

	std::ofstream file( "metricmap.m", std::ios::out | std::ios::trunc );

	file << "roundness = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMap->GetLabelObject( i )->GetRoundness() << " ";
	}

	file << " ]; " << std::endl;

	file << "elongation = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMap->GetLabelObject( i )->GetElongation() << " ";
	}

	file << " ]; " << std::endl;

	file << "V = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMap->GetLabelObject( i )->GetNumberOfPixels() << " ";
	}

	file << " ]; " << std::endl;

	file << "S = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMapContour->GetLabelObject( i )->GetNumberOfPixels() << " ";
	}

	file << " ]; " << std::endl;

	file << "R = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMap->GetLabelObject( i )->GetEquivalentSphericalRadius() << " ";
	}

	file << " ]; " << std::endl;
	
	file << "L = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		file << labelMapSkeleton->GetLabelObject( i )->GetNumberOfPixels() << " ";
	}

	file << " ]; " << std::endl;

	file << "X;";
	unsigned int j = 0;

	for( unsigned int i = 1; i <= max; ++i )
	{
		if( labelMap->HasLabel( i ) && labelMapContour->HasLabel( i ) && labelMapSkeleton->HasLabel( i ) )
		{
			++j;
			file << "X(" << i << ") = " << j << ";" << std::endl;
		}
	}





/*
	file << "Rmax = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		BiggestComponentFilterType::VectorType vector = labelMap->GetLabelObject( i )->GetEquivalentEllipsoidDiameter();
		double max = vector[0];
		for( unsigned int j = 1; j < Dim; ++j )
		{
			max = std::max( max, vector[j] );
		}
		file << max << " ";
	}

	file << " ]; " << std::endl;

	file << "Rmin = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		BiggestComponentFilterType::VectorType vector = labelMap->GetLabelObject( i )->GetEquivalentEllipsoidDiameter();
		double min = vector[0];

		std::cout << vector[0] << "  ";
		for( unsigned int j = 1; j < Dim; ++j )
		{
			std::cout << vector[j] << "  ";
			min = std::min( min, vector[j] );
		}
		std::cout << std::endl;
		file << min << " ";
	}

	file << " ]; " << std::endl;

	file << "feret = [";

	for( unsigned int i = 1; i <= max; ++i )
	{
		file << labelMapContour->GetLabelObject( i )->GetFeretDiameter() << " ";
	}

	file << " ]; " << std::endl;


*/


	file.close();

}

/*void histoMap( GraphFilterType::RealListMapType phi, GraphFilterType::RealListMapType theta, GraphFilterType::LabelImagePixelType max )
{
	std::cout << " coucou " << std::endl;
	GraphFilterType::CostType cost;
	GraphFilterType::RealListMapType::iterator it;

	std::ofstream file( "map.m", std::ios::out | std::ios::trunc );

	file << "clear all; close all;" << std::endl;
	file << "phi = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( phi.count( label ) )
		{
			it = phi.find( label );
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << *list << " ";
			}
		}
	}

	file << "];" << std::endl;

	file << "theta = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( theta.count( label ) )
		{
			it = theta.find( label );
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << *list << " ";
			}
		}
	}

	file << "];" << std::endl;

	file << "labelphi = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( phi.count( label ) )
		{
			it = phi.find( label );
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << label << " ";
			}
		}
	}

	file << "];" << std::endl;

	file << "labeltheta = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( theta.count( label ) )
		{
			it = theta.find( label );
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << label << " ";
			}
		}
	}

	file << "];" << std::endl;



	file << "positionphi = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( phi.count( label ) )
		{
			it = phi.find( label );
			unsigned int position = 0;
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << position << " ";
				++position;
			}
		}
	}

	file << "];" << std::endl;

	file << "positiontheta = [";

	for( GraphFilterType::LabelImagePixelType label = 1; label < max; ++label )
	{
		if( theta.count( label ) )
		{
			it = theta.find( label );
			unsigned int position = 0;
			for( GraphFilterType::ListRealType::iterator list = it->second.begin(); list != it->second.end(); ++list )
			{
				file << position << " ";
				++position;
			}
		}
	}

	file << "];" << std::endl;


	file << "sortphi = sort(phi);" << std::endl;
	file << "sorttheta = sort(theta);" << std::endl;


	file.close();
}*/
