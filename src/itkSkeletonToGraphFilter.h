#ifndef __itkSkeletonToGraphFilter_h
#define __itkSkeletonToGraphFilter_h

#include "Graph/itkGraphSource.h"
#include "itkObjectFactory.h"

#include "itkFixedArray.h"

namespace itk
{

template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage = TInputImage >
class ITK_EXPORT SkeletonToGraphFilter : public GraphSource< TOutputGraph >
{
public:
	/** Standard class typedefs. */
	typedef SkeletonToGraphFilter			Self;
	typedef GraphSource< TOutputGraph >		Superclass;
	typedef SmartPointer< Self >			Pointer;
	typedef SmartPointer< const Self >		ConstPointer;

	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( SkeletonToGraphFilter, GraphSource );

	/**  Some Image related typedefs. */
	typedef TInputImage							ImageType;
	typedef typename ImageType::Pointer			ImagePointer;
	typedef typename ImageType::ConstPointer	ImageConstPointer;
	typedef typename ImageType::RegionType		RegionType;
	typedef typename ImageType::PixelType		PixelType;
	typedef typename ImageType::IndexType		IndexType;
	typedef typename ImageType::OffsetType		OffsetType;

	/** Some Graph related typedefs. */
	typedef TOutputGraph								GraphType;
	typedef typename GraphType::GraphTraitsType			GraphTraitsType;
	typedef typename GraphType::Pointer					GraphPointer;
	typedef typename GraphType::NodeType				NodeType;
	typedef typename GraphType::NodePointerType			NodePointerType;
	typedef typename GraphType::NodeIdentifierType		NodeIdentifierType;
	typedef typename GraphTraitsType::NodeWeightType	NodeWeightType;
	typedef typename GraphType::EdgeType				EdgeType;
	typedef typename GraphType::EdgePointerType			EdgePointerType;
	typedef typename GraphType::EdgeIdentifierType		EdgeIdentifierType;
	typedef typename GraphTraitsType::EdgeWeightType	EdgeWeightType;
	typedef typename GraphTraitsType::VectorType		VectorType;
	typedef typename GraphTraitsType::PointType			PointType;
	typedef typename GraphTraitsType::IndexContainerType	IndexContainerType;
	typedef EdgeIdentifierType							IdentifierType;

	/** Some Label Image related typedefs. */
	typedef TLabelImage									LabelImageType;
	typedef typename LabelImageType::Pointer			LabelImagePointer;
	typedef typename LabelImageType::ConstPointer		LabelImageConstPointer;
	typedef typename LabelImageType::PixelType			LabelType;

	static const unsigned int ImageDimension = ImageType::ImageDimension;

	/** Some Identifier related typedefs. */
	typedef Image< NodeIdentifierType, ImageDimension >	IdentifierImageType;
	typedef typename IdentifierImageType::Pointer		IdentifierImagePointer;
	typedef typename IdentifierImageType::ConstPointer	IdentifierImageConstPointer;

	/** Some Vesselness related typedefs. */
	typedef TVesselnessImage							VesselnessImageType;
	typedef typename VesselnessImageType::Pointer		VesselnessImagePointer;
	typedef typename VesselnessImageType::ConstPointer	VesselnessImageConstPointer;

	/** Iterator typedefs */
	typedef ConstNeighborhoodIterator< ImageType >				ImageIteratorType;
	typedef NeighborhoodIterator< LabelImageType >				LabelImageIteratorType;
	typedef NeighborhoodIterator< IdentifierImageType >			IdentifierImageIteratorType;
	typedef ConstNeighborhoodIterator< VesselnessImageType >	VesselnessImageIteratorType;

	typedef std::vector< unsigned int > 				ArrayType;

	/** Set the input image of this process object. */
	void SetInput( unsigned int idx, const ImageType *input );
	void SetInput( ImageType * );

	/** Get the input image of this process object. */
	const ImageType * GetInput( unsigned int idx );
	
	/** Get the output graph of this process object. */
	GraphType * GetOutput( void );

	/** Get the output label image of this process object. */
	LabelImageType * GetLabelOutput( void );

	/** Get the output identifier image of this process object. */
	IdentifierImageType * GetIdentifierOutput( void );

	void SetVesselnessInput( const VesselnessImageType * vesselness );
	const VesselnessImageType * GetVesselnessInput() const;

	/** Prepare the output */
//	virtual void GenerateOutputInformation( void );

protected:
	SkeletonToGraphFilter();
	~SkeletonToGraphFilter();
//	void PrintSelf( std::ostream& os, Indent indent ) const;

	virtual void GenerateData();

	void TransformSkeletonToGraph( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit );

	void MakeNode( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, bool * borderFound );

	void FillIntersection( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, LabelType & label );

	void ConnectNodes( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, LabelType label, bool borderFound );

	void ConnectSelf( NodePointerType node );

	unsigned int GetNeighborhoodCount( const ImageIteratorType & it ) const;

	unsigned int GetNeighborhoodMap( const ImageIteratorType & it, ArrayType & map ) const;

//	void UpdateNodeAttributes( NodePointerType * nodes );
//	void UpdateEdgeAttributes( EdgePointerType * edges );

	bool ReachBranchSkeleton( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit ) const;

	DataObject::Pointer MakeOutput( unsigned int idx );




private:
	SkeletonToGraphFilter( const SkeletonToGraphFilter & ); // purposely not implemented
	void operator=( const SkeletonToGraphFilter & );		// purposely not implemented
	
	unsigned int m_NeighborhoodCount;
	LabelType m_Label;
	RegionType m_Region;
};

} // namesapce itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSkeletonToGraphFilter.hxx"
#endif

#endif
