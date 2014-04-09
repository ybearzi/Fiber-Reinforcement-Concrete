#ifndef __itkFiberExtractionGraphFilter_h
#define __itkFiberExtractionGraphFilter_h

#include "Graph/itkGraphToGraphFilter.h"

namespace itk
{
template< class TGraph >
class ITK_EXPORT FiberExtractionGraphFilter : public GraphToGraphFilter< TGraph, TGraph >
{
public:
	/** Standard class typedefs. */
	typedef FiberExtractionGraphFilter				Self;
	typedef GraphToGraphFilter< TGraph, TGraph >	Superclass;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( GraphToGraphFilter, FiberExractionGraphFilter );

	/** Some convenient typedefs. */
	typedef TGraph													GraphType;
	typedef typename GraphType::Pointer								GraphPointer;
	typedef typename GraphType::GraphTraitsType						GraphTraitsType;
	typedef typename GraphType::NodeContainerType					NodeContainerType;
	typedef typename GraphType::EdgeContainerType					EdgeContainerType;
	typedef typename GraphTraitsType::NodeType						NodeType;
	typedef typename GraphTraitsType::NodePointerType				NodePointerType;
	typedef typename GraphTraitsType::EdgeType						EdgeType;
	typedef typename GraphTraitsType::EdgePointerType				EdgePointerType;
	typedef typename GraphTraitsType::VectorType					VectorType;
	typedef typename GraphTraitsType::IndexContainerType			IndexContainerType;
	typedef typename GraphTraitsType::NodeIdentifierType			NodeIdentifierType;
	typedef typename GraphTraitsType::EdgeIdentifierType			EdgeIdentifierType;
	typedef typename IndexContainerType::value_type					IndexType;
	typedef typename GraphTraitsType::ScalarContainerType			ScalarContainerType;
	typedef typename GraphTraitsType::NodeWeightType				NodeWeightType;
	typedef typename GraphTraitsType::EdgeIdentifierContainerType	EdgeIdentifierContainerType;

	static const unsigned int Dimension = IndexType::Dimension;

	itkGetConstMacro( MaxIntersection, NodeIdentifierType );
	itkGetConstMacro( MergeThreshold, double );
	itkGetConstMacro( SplitThreshold, double );

	itkSetMacro( MergeThreshold, double );
	itkSetMacro( SplitThreshold, double );
	void SetMaxIntersection( NodeIdentifierType inter );

protected:

	FiberExtractionGraphFilter();
	~FiberExtractionGraphFilter(){};

	virtual void GenerateData();

	virtual void GenerateOutputInformation();
	
	virtual void Initialize();
	virtual void SplitAndMerge();

	NodeIdentifierType m_MaxIntersection;
	double m_MergeThreshold, m_SplitThreshold;
	NodeWeightType m_Maximum;
	bool m_MaxIntersectionSet;

private:
	FiberExtractionGraphFilter( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely node implemented

	
	GraphPointer m_Output;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberExtractionGraphFilter.hxx"
#endif

#endif
