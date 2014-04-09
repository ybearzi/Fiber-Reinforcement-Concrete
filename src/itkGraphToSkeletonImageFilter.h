#ifndef __itkGraphToSkeletonImageFilter_h
#define __itkGraphToSkeletonImageFilter_h

#include "Graph/itkGraphToImageFilter.h"

namespace itk
{
template< class TInputGraph, class TOutputImage >
class ITK_EXPORT GraphToSkeletonImageFilter : public GraphToImageFilter< TInputGraph, TOutputImage >
{
public:	
	/** Standard class typedefs. */
	typedef GraphToSkeletonImageFilter						Self;
	typedef GraphToImageFilter< TInputGraph, TOutputImage >	Superclass;
	typedef SmartPointer< Self >							Pointer;
	typedef SmartPointer< const Self >						ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( GraphToSkeletonImageFilter, GraphToImageFilter );

	/** Some convenient typedefs. */
	typedef typename Superclass::OutputImageRegionType	OutputImageRegionType;
	typedef TInputGraph									InputGraphType;
	typedef typename InputGraphType::Pointer			InputGraphPointer;
	typedef typename InputGraphType::ConstPointer		InputGraphConstPointer;
	typedef typename InputGraphType::NodeContainerType	NodeContainerType;
	typedef typename InputGraphType::GraphTraitsType	GraphTraitsType;
	typedef typename GraphTraitsType::NodeType			NodeType;
	typedef TOutputImage								OutputImageType;
	typedef typename OutputImageType::Pointer			OutputImagePointer;
	typedef typename OutputImageType::SizeType			SizeType;
	typedef typename OutputImageType::ValueType			ValueType;
	typedef typename OutputImageType::PixelType			PixelType;

	/** ImageDimension constants */
	itkStaticConstMacro( OutputImageDimension, unsigned int, TOutputImage::ImageDimension );

	itkGetConstMacro( Region, OutputImageRegionType );
	itkSetMacro( Region, OutputImageRegionType );

protected:
	GraphToSkeletonImageFilter();
	~GraphToSkeletonImageFilter(){};
	
	virtual void GenerateData();
	virtual void GenerateOutputInformation();
	virtual DataObject::Pointer MakeOutput( unsigned int idx );

private:
	GraphToSkeletonImageFilter( const Self & ); //purposely not implemented
	void operator=( const Self & ); //purposely not implemented

	OutputImageRegionType m_Region;


};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGraphToSkeletonImageFilter.hxx"
#endif

#endif
