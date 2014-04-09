#ifndef __itkFiberExtractionGraphFilter_hxx
#define __itkFiberExtractionGraphFilter_hxx

#include "itkFiberExtractionGraphFilter.h"

#include "itkNumericTraits.h"

namespace itk
{

template< class TGraph >
FiberExtractionGraphFilter< TGraph >
::FiberExtractionGraphFilter()
{
	m_Maximum = NumericTraits< NodeWeightType >::max();
	m_MergeThreshold = 5;
	m_SplitThreshold = 4;
	m_MaxIntersectionSet = false;
	m_MaxIntersection = 0;
}

template< class TGraph >
void FiberExtractionGraphFilter< TGraph >
::GenerateOutputInformation()
{
	Superclass::GenerateOutputInformation();
	m_Output = this->GetInput()->Copy();
}

template< class TGraph >
void FiberExtractionGraphFilter< TGraph >
::GenerateData()
{
	std::cout << " Initialize " << std::endl;
	this->Initialize();
	std::cout << " Split and Merge " << std::endl;
	this->SplitAndMerge();
	std::cout << " Graft Output " << std::endl;
	this->GraftOutput( m_Output );
}

template< class TGraph >
void FiberExtractionGraphFilter< TGraph >
::Initialize()
{
	m_Output->InitializeGraphAttributes();
	if( !m_MaxIntersectionSet )
	{
		m_MaxIntersection = m_Output->MaxIntersectionLabel();
	}
}

template< class TGraph >
void FiberExtractionGraphFilter< TGraph >
::SplitAndMerge()
{
	NodePointerType node;

	unsigned int nodesNumber = m_Output->GetNodeContainer()->Size();
	unsigned int i = 0;
	for( typename NodeContainerType::Iterator it = m_Output->GetNodeContainer()->Begin(); i < nodesNumber; ++it, ++i)
	{
		m_Output->TryToSplit( &( it.Value() ), m_SplitThreshold, m_MaxIntersection );
	}


	do
	{
		node = m_Output->GetMaximumNodeWeight( m_Maximum );
		while( m_Output->MergeSuccess( node, m_MergeThreshold ) );
		if( node != NULL )
		{
			m_Maximum = node->Weight;
		}
	} while( node != NULL );

//	m_Output->FillGaps( 0.5 );

	m_Output->UpdateNodesWeightWithoutInteraction();
}

template< class TGraph >
void FiberExtractionGraphFilter< TGraph >
::SetMaxIntersection( NodeIdentifierType inter )
{
	m_MaxIntersectionSet = true;
	m_MaxIntersection = inter;
}

} // namespace itk

#endif
