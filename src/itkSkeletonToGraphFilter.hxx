#ifndef __itkSkeletonToGraphFilter_hxx
#define __itkSkeletonToGraphFilter_hxx

#include "itkSkeletonToGraphFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"

#include <iostream>

namespace itk
{

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SkeletonToGraphFilter()
{
//	std::cout << " constructor " << std::endl;

		std::cout << sizeof( NodeType ) << std::endl;
	this->ProcessObject::SetNumberOfRequiredInputs( 2 );

	m_NeighborhoodCount = 1;

	for( unsigned int i = 0; i < ImageDimension; ++i )
	{
		m_NeighborhoodCount *= 3;
	}

	m_Label = 1;

	GraphPointer output = dynamic_cast< GraphType* > ( this->MakeOutput( 0  ).GetPointer() );
	LabelImagePointer outputLabel = dynamic_cast< LabelImageType* > ( this->MakeOutput( 1  ).GetPointer() );
	IdentifierImagePointer outputIdentifier = dynamic_cast< IdentifierImageType* > ( this->MakeOutput( 2  ).GetPointer() );

	this->ProcessObject::SetNumberOfRequiredOutputs( 3 );
	this->ProcessObject::SetNthOutput( 1, outputLabel.GetPointer() );
	this->ProcessObject::SetNthOutput( 2, outputIdentifier.GetPointer() );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::~SkeletonToGraphFilter() 
{}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
DataObject::Pointer
SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::MakeOutput( unsigned int idx )
{
	switch( idx )
	{
		case 0: {
			GraphPointer outputGraph = GraphType::New();
			return dynamic_cast< DataObject * > ( outputGraph.GetPointer() );
			break; }
		case 1: {
			LabelImagePointer outputLabel = LabelImageType::New();
			return dynamic_cast< DataObject * > ( outputLabel.GetPointer() );
			break; }
		case 2: {
			IdentifierImagePointer outputIdentifier = IdentifierImageType::New();
			return dynamic_cast< DataObject * > ( outputIdentifier.GetPointer() );
			break; }
		default:
			return NULL;
			break;
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetInput( ImageType * image )
{
	this->SetNthInput( 0, const_cast< ImageType * > ( image ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetInput( unsigned int idx, const ImageType * input )
{
	this->ProcessObject::SetNthInput( idx, const_cast< ImageType * > ( input ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
const typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetInput( unsigned int idx )
{
	return dynamic_cast< const ImageType * > ( this->ProcessObject::GetInput( idx ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GraphType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetOutput( void )
{
	return dynamic_cast< GraphType * > ( this->ProcessObject::GetOutput( 0 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::LabelImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetLabelOutput( void )
{
	return dynamic_cast< LabelImageType * > ( this->ProcessObject::GetOutput( 1 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
typename SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::IdentifierImageType * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetIdentifierOutput( void )
{
	return dynamic_cast< IdentifierImageType * > ( this->ProcessObject::GetOutput( 2 ) );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GenerateData()
{
	std::cout << " ---- GenerateData " << std::endl;

	ImageConstPointer input = static_cast< const ImageType * > ( this->GetInput( 0 ) );

	typename ImageIteratorType::RadiusType radius;
	radius.Fill( 1 );

	this->GetLabelOutput()->SetRegions( input->GetLargestPossibleRegion() );
	this->GetLabelOutput()->Allocate();
	this->GetLabelOutput()->FillBuffer( 0 );
	this->GetLabelOutput()->Update();

	this->GetIdentifierOutput()->SetRegions( input->GetLargestPossibleRegion() );
	this->GetIdentifierOutput()->Allocate();
	this->GetIdentifierOutput()->FillBuffer( 0 );
	this->GetIdentifierOutput()->Update();

	this->GetOutput()->GetEdgeContainer()->Reserve( 1000000 );

	ImageIteratorType it( radius, input, input->GetLargestPossibleRegion() );
	LabelImageIteratorType lit( radius, this->GetLabelOutput(), input->GetLargestPossibleRegion() );
	IdentifierImageIteratorType iit( radius, this->GetIdentifierOutput(), input->GetLargestPossibleRegion() );
	VesselnessImageIteratorType vit( radius, this->GetVesselnessInput(), input->GetLargestPossibleRegion() );

	m_Region = input->GetLargestPossibleRegion();

	GraphPointer output = this->GetOutput();

	LabelType label = 1;

	unsigned int neighborhoodCount = std::pow( 3, ImageDimension );
	unsigned int center = neighborhoodCount / 2;

	it.GoToBegin();
	lit.GoToBegin();
	iit.GoToBegin();
	vit.GoToBegin();

	while( this->ReachBranchSkeleton( it, lit, iit, vit ) )
	{
		std::cout << " Skeleton found " << it.GetIndex() << std::endl;

		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		this->TransformSkeletonToGraph( tit, tlit, tiit, tvit );

		std::cout << " Reach " << std::endl;

	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::TransformSkeletonToGraph( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit )
{
	std::cout << std::endl;
	std::cout << " ---- TransformSkeletonToGraph " << std::endl;

	std::cout << " node1 " << std::endl;
	std::cout << " Create New Node " << std::endl;
	NodePointerType node = this->GetOutput()->CreateNewNode();
	std::cout << " Node created " << node->Identifier << std::endl;
	node->LineIndexes = IndexContainerType();

	std::cout << " max " << this->GetOutput()->GetNodeContainer()->Size() << std::endl;

	node->Weight = 1;
	node->BackLabel = 0;
	node->FrontLabel = 0;

	NodePointerType nodes[2];

	this->ConnectSelf( node );

	std::cout << " ok " << std::endl;

	ImageIteratorType tit = it;
	LabelImageIteratorType tlit = lit;
	IdentifierImageIteratorType tiit = iit;
	VesselnessImageIteratorType tvit = vit;

	std::cout << " ok " << std::endl;

	bool borderFound = false;
	this->MakeNode( tit, tlit, tiit, tvit, node, &borderFound );
	std::cout << " Debug " << node->Identifier << " " <<  node->LineIndexes.size() << std::endl;
	if( !borderFound )
	{
		borderFound = true;
		std::cout << " *** make edge reverse " << std::endl;
		this->MakeNode( it, lit, iit, tvit, node, &borderFound );
		std::cout << " Debug " << node->Identifier << " " <<  node->LineIndexes.size() << std::endl;
	}

	std::cout << " Debug " << node->Identifier;
	
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::MakeNode( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, bool * borderFound )
{ 
	std::cout << " *** MakeNode " << std::endl;
	ArrayType map; //m_NeighborhoodCount );
	unsigned int count = this->GetNeighborhoodMap( it, map );

	std::cout << " count " << count << std::endl;
	std::cout << " borderFound " << *borderFound << std::endl;
	std::cout << " Identifier " << node->Identifier << std::endl;

	std::cout << " index " << it.GetIndex() << " pixel " << it.GetCenterPixel() << std::endl;
	bool end = false;

	if( count != 2 )
	{

		// Intersection
		if( count > 2 )
		{
			std::cout << " intersection " << std::endl;

			ImageIteratorType tit = it;
			LabelImageIteratorType tlit = lit;
			IdentifierImageIteratorType tiit = iit;
			VesselnessImageIteratorType tvit = vit;
			
			// Node not existing -> creating new node
			if( !lit.GetCenterPixel() )
			{
				std::cout << " intersection not existing -> creating new nodes and connect them " << std::endl;	


				std::cout << " ---- FillIntersection " << std::endl;
				this->FillIntersection( tit, tlit, tiit, tvit, m_Label );
				std::cout << " end " << std::endl;
			}
			else
			{
				std::cout << " node existing " << std::endl;
				std::cout << " ---- ConnectNodes " << std::endl;
				std::cout<< " Identifier " << node->Identifier << " ConnectNodes " << std::endl;
				this->ConnectNodes( tit, tlit, tiit, tvit, node, m_Label, *borderFound );
				++m_Label;
				std::cout << " end " << std::endl;
			}
		}

		if( *borderFound )
		{
//			iit.SetCenterPixel( node->Identifier );
//			this->UpdateNodeAttributes( node );	
//			return
			end = true;
		}
//		std::cout << " salut " << *borderFound << " " << nodes[*borderFound] << std::endl;
//		std::cout << nodes[*borderFound]->Weight << std::endl;
//		std::cout << nodes[*borderFound]->IncomingEdges.size() << std::endl;
//		nodes[*borderFound]->LineIndexes = IndexContainerType();
//		std::cout << " coucou " << std::endl;
//		std::cout << nodes[*borderFound]->LineIndexes.size() << std::endl;
		if( !( node->LineIndexes.size() ) )
		{
			std::cout << " oui " << std::endl;
			*borderFound = !*borderFound;
		}
		else
		{
			end = true;
//			iit.SetCenterPixel( node->Identifier );
//			return;
		}
	}

	std::cout << " line " << std::endl;

	unsigned int i = 0;

	if( node->LineIndexes.size() )
	{	
		if( node->LineIndexes.front() == it.GetIndex() || node->LineIndexes.back() == it.GetIndex() )
		{
			std::cout << " if " << std::endl;
			i = 1;
		}
		else
		{
			std::cout << " else " << std::endl;
			IndexType backIndex = node->LineIndexes.back();
			IndexType frontIndex = node->LineIndexes.front();
			std::cout << " back " << backIndex << " front " << frontIndex << " index " << it.GetIndex( map[0] ) << std::endl;
			i = ( backIndex == it.GetIndex( map[0] ) || frontIndex == it.GetIndex( map[0] ) );
			if( i && count > 1 && ( backIndex == it.GetIndex( map[1] ) || frontIndex == it.GetIndex( map[1] ) ) )
			{
				std::cout << " end ";
				end = true;
				*borderFound = true;
			}
		}
	}

	if( count > 1 )
	{
		std::cout << " i 0 " << iit.GetIndex( map[0] ) << " i 1 " << iit.GetIndex( map[1] );
	}

	std::cout << " i " << i << std::endl;
	std::cout << " index " << iit.GetIndex() << std::endl;
	std::cout << " iit " << iit.GetCenterPixel() << std::endl;

	if( iit.GetCenterPixel() != node->Identifier + 1 )
	{
		std::cout << " push size " << node->LineIndexes.size() << std::endl;
		if( *borderFound )
		{
			node->LineIndexes.push_back( it.GetIndex() );
			std::cout << " first " << std::endl;
			std::cout << " vit " << vit.GetIndex() << std::endl;
			std::cout << " vit " << vit.GetCenterPixel() << std::endl;
			std::cout << " push size " << node->VesselnessValues.size() << std::endl;
			node->VesselnessValues.push_back( vit.GetCenterPixel() );
			std::cout << " second " << std::endl;
		}
		else
		{
			node->LineIndexes.push_front( it.GetIndex() );
			node->VesselnessValues.push_front( vit.GetCenterPixel() );
		}
		std::cout << " end push back " << std::endl;

		iit.SetCenterPixel( node->Identifier + 1 );
	}
	std::cout << " i " << i << " iit " << iit.GetCenterPixel() << std::endl;

	if( end )
	{
		return;
	}

	it += it.GetOffset( map[i] );
	lit += lit.GetOffset( map[i] );
	iit += iit.GetOffset( map[i] );
	vit += vit.GetOffset( map[i] );

	this->MakeNode( it, lit, iit, vit, node, borderFound );
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ConnectSelf( NodePointerType node )
{
	std::cout << " Create New edge " << std::endl;
	EdgePointerType edge = this->GetOutput()->CreateNewEdge( node, node, 1 );
	std::cout << " New Edge created " << edge->Identifier << std::endl;
	edge->ReverseEdgeIdentifier = edge->Identifier;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ConnectNodes( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, LabelType label, bool borderFound )
{
	std::cout << " Index " << it.GetIndex() << " Identifier " << iit.GetCenterPixel() << std::endl;

	ArrayType map;

	unsigned int count = this->GetNeighborhoodMap( it, map );

	if( count <= 2 )
	{
		if( iit.GetCenterPixel() )
		{
			NodePointerType tmp = this->GetOutput()->GetNodePointer( iit.GetCenterPixel() - 1 );
			if( tmp->LineIndexes.size() < 2 || iit.GetIndex() == *( ++( tmp->LineIndexes.begin() ) ) )
			{
				tmp->FrontLabel = m_Label;
			}
			else
			{
				tmp->BackLabel = m_Label;
			}
		}
	}


	if( lit.GetCenterPixel() == label )
	{
		return;
	}

	std::cout << " count " << count << std::endl;
	std::cout << " it " << it.GetCenterPixel() << std::endl;

	if( count <= 2 )
	{
		if( iit.GetCenterPixel() )
		{
			bool nodeAlreadyProcessed = false;
			NodePointerType nodetmp = this->GetOutput()->GetNodePointer( iit.GetCenterPixel()-1 );
//			std::cout << " iit " << iit.GetCenterPixel() << " id " << node->Identifier << std::endl;
//			if( borderFound )
//			{
				for( unsigned int i = 0; i < nodetmp->IncomingEdges.size(); ++i )
				{
					std::cout << " incoming " << this->GetOutput()->GetSourceNode( nodetmp->IncomingEdges[i] ).Identifier << std::endl;
//						std::cout << i << std::endl;
					if( this->GetOutput()->GetSourceNodePointer( nodetmp->IncomingEdges[i] )->Identifier == node->Identifier
					|| this->GetOutput()->GetSourceNodePointer( nodetmp->IncomingEdges[i] )->Identifier == node->Identifier )
					{
//						std::cout << " id " <<  this->GetOutput()->GetSourceNodePointer( node->IncomingEdges[i] )->Identifier << " id " <<  nodes[borderFound]->Identifier << std::endl;
//						std::cout << " nodeAlreadyProcessed " << std::endl;
						nodeAlreadyProcessed = true;
					}
				}
				for( unsigned int i = 0; i < nodetmp->OutgoingEdges.size(); ++i )
				{
					std::cout << " outgoing " << std::endl;
//						std::cout << i << std::endl;
					if( this->GetOutput()->GetSourceNodePointer( nodetmp->OutgoingEdges[i] )->Identifier == node->Identifier
					|| this->GetOutput()->GetSourceNodePointer( nodetmp->OutgoingEdges[i] )->Identifier == node->Identifier )
					{
//						std::cout << " id " <<  this->GetOutput()->GetSourceNodePointer( node->IncomingEdges[i] )->Identifier << " id " <<  nodes[borderFound]->Identifier << std::endl;
//						std::cout << " nodeAlreadyProcessed " << std::endl;
						nodeAlreadyProcessed = true;
					}
				}
//	
//			}
			

			
			if( !nodeAlreadyProcessed )
			{
				std::cout << " Node not already processed " << std::endl;
				std::cout << " edge1 " << std::endl;
				std::cout << " Create New Edge " << std::endl;
				EdgePointerType edge2 = this->GetOutput()->CreateNewEdge( nodetmp, node, 1 );
				std::cout << " Edge created " << edge2->Identifier << std::endl;
				std::cout << " Create New Edge " << std::endl;
				EdgePointerType edge1 = this->GetOutput()->CreateNewEdge( node, nodetmp, 1 );
				std::cout << " Edge created " << edge1->Identifier << std::endl;
				std::cout << " coucou " << std::endl;
				std::cout << edge1->Identifier << std::endl;
				std::cout << edge2->Identifier << std::endl;
				edge1->ReverseEdgeIdentifier = edge2->Identifier;
				std::cout << " reverse " << std::endl;
				edge2->ReverseEdgeIdentifier = edge1->Identifier;

//				EdgePointerType edges[2];
//				edges[0] = edge1;
//				edges[1] = edge2;

//				this->UpdateEdgeAttributes( edges );

			}
		}
		std::cout << " return " << std::endl;
		return;
	}

	bool keepGoing = false;
	if( lit.GetCenterPixel() != label )
	{
		keepGoing = true;
		lit.SetCenterPixel( label );
	}


	for( unsigned int i = 0; i < count; ++i )
	{
		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		tit += tit.GetOffset(map[i]);
		tlit += tit.GetOffset(map[i]);
		tiit += tit.GetOffset(map[i]);
		tvit += tit.GetOffset(map[i]);


		if( keepGoing )
		{
			this->ConnectNodes( tit, tlit, tiit, tvit, node, label, borderFound );
		}
	}

}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::FillIntersection( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, LabelType & label )
{
	std::cout << "inter " << std::endl;
	std::cout << lit.GetCenterPixel() << std::endl;
	std::cout << it.GetCenterPixel() << std::endl;
	if( lit.GetCenterPixel() )
	{
		return;
	}

	lit.SetCenterPixel( label );

	ArrayType map;

	unsigned int count = this->GetNeighborhoodMap( it, map );
	std::cout << " count " << count << std::endl;
	for( unsigned int i = 0; i < count; ++i )
	{
		std::cout << " map " << map[i] << " i " << i << std::endl;
	}

	if( count <= 2 )
	{
		return;
	}

	for( unsigned int i = 0; i < count; ++i )
	{
		std::cout << " i " << i << " count " << count << std::endl;
		ImageIteratorType tit = it;
		LabelImageIteratorType tlit = lit;
		IdentifierImageIteratorType tiit = iit;
		VesselnessImageIteratorType tvit = vit;

		std::cout << "map " << map[i] << std::endl;

		tit += tit.GetOffset(map[i]);
		tlit += tit.GetOffset(map[i]);
		tiit += tit.GetOffset(map[i]);
		tvit += tit.GetOffset(map[i]);


		std::cout << " FillIntersection " << std::endl;
		this->FillIntersection( tit, tlit, tiit, tvit, label );
	}
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
unsigned int SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetNeighborhoodMap( const ImageIteratorType & it, ArrayType & map ) const
{
	unsigned int count = 0;
	for( unsigned int i = 0; i < m_NeighborhoodCount; ++i )
	{
		if( m_Region.IsInside( it.GetIndex( i ) ) && it.GetPixel( i ) && it.GetIndex( i ) != it.GetIndex() )
		{
			map.push_back( i );
			++count;
		}
	}
	return count;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
unsigned int SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetNeighborhoodCount( const ImageIteratorType & it ) const
{
	unsigned int count = 0;
	for( unsigned int i = 0; i < m_NeighborhoodCount; ++i )
	{
		if( m_Region.IsInside( it.GetIndex( i ) ) && it.GetPixel( i ) )
		{
			++count;
		}
	}
	return count == 2 || count == 3;
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
bool SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::ReachBranchSkeleton( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit ) const
{
	std::cout << " ---- ReachBranchSkeleton " << std::endl;

	while( !it.IsAtEnd() && ( !this->GetNeighborhoodCount( it ) || ( iit.GetCenterPixel() || !it.GetCenterPixel() ) ) )
	{
		++it;
		++iit;
		++lit;
		++vit;
	}
	std::cout << " count " << this->GetNeighborhoodCount( it ) << std::endl;
	std::cout << " index " << it.GetIndex() << std::endl;
	std::cout << " lit " << lit.GetCenterPixel() << std::endl;
	std::cout << " !it " << !it.GetCenterPixel() << std::endl;
	return !it.IsAtEnd();
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
void SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::SetVesselnessInput( const TVesselnessImage* vesselness )
{
	this->SetNthInput( 1, const_cast< TVesselnessImage* > ( vesselness ) ); 
}

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage >
const TVesselnessImage * SkeletonToGraphFilter< TInputImage, TOutputGraph, TVesselnessImage, TLabelImage >
::GetVesselnessInput() const
{
	return static_cast< const TVesselnessImage * > ( this->ProcessObject::GetInput( 1 ) );
}

} // namespace itk

#endif
