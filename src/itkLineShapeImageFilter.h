#ifndef __itkLineShapeImageFilter_h
#define __itkLineShapeImageFilter_h

// C++ includes

#include <iostream>

#include "itkNumericTraits.h"
#include "itkCrossHelper.h"

// Image includes

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkVector.h"

// Filter includes

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"

// Eigen values processing

#include "itkSymmetricEigenAnalysis.h"

// Functor

#include "itkUnaryFunctorImageFilter.h"

namespace itk {

template< class TInputImage, class TOutputImage = TInputImage, class TLabelImage = TOutputImage, class TVectorImage = Image< Vector< typename NumericTraits< typename TInputImage::PixelType >::RealType, TInputImage::ImageDimension >, TInputImage::ImageDimension > >
class ITK_EXPORT LineShapeImageFilter :
	public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
	/** Standard class typedefs. */
	typedef LineShapeImageFilter 								Self;
	typedef ImageToImageFilter< TInputImage, TOutputImage >		Superclass;
	typedef SmartPointer< Self >								Pointer;
	typedef SmartPointer< const Self >							ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( LineShapeImageFilter, ImageToImageFilter );


	/** Convenient typedef. */
	typedef typename Superclass::InputImageType					InputImageType;
	typedef typename Superclass::InputImageConstPointer			InputImageConstPointer;
	typedef typename Superclass::InputImagePixelType			InputImagePixelType;
	typedef typename Superclass::InputImagePointer				InputImagePointer;
	typedef typename Superclass::InputImageRegionType			InputImageRegionType;
	typedef typename Superclass::OutputImageType				OutputImageType;
	typedef typename Superclass::OutputImagePixelType			OutputImagePixelType;
	typedef typename Superclass::OutputImagePointer				OutputImagePointer;
	typedef typename Superclass::OutputImageRegionType			OutputImageRegionType;
	typedef TLabelImage											LabelImageType;
	typedef typename LabelImageType::PixelType 					LabelType;
	typedef typename LabelImageType::Pointer					LabelImagePointer;
	typedef TVectorImage										VectorImageType;
	typedef typename VectorImageType::PixelType					VectorType;

	typedef float		RealType;

	static const unsigned int ImageDimension = InputImageType::ImageDimension;

	typedef Image< RealType, ImageDimension >					RealImageType;


	OutputImageType* GetBinaryOutput();
	LabelImageType* GetLabelOutput();
	VectorImageType* GetLineDirectionOutput();
	RealImageType* GetEigenValuesOutput( unsigned int i );
	RealImageType* GetVesselnessOutput();

	/** Get const Macros */
	itkGetConstMacro( BrightLine, LabelType );
	itkGetConstMacro( DarkLine, LabelType );
	itkGetConstMacro( LabelCount, LabelType );
	itkGetConstMacro( ExtractBrightLine, bool );
	itkGetConstMacro( EigenValuesExtraction, bool );
	itkGetConstMacro( LabelImage, bool );
	itkGetConstMacro( Region, InputImageRegionType );
	itkGetConstMacro( Sigma, RealType );
	itkGetConstMacro( OutsideValue, OutputImagePixelType );
	itkGetConstMacro( InsideValue, OutputImagePixelType );
	itkGetConstMacro( DimensionsProcessed, unsigned int );

	/** Set Macros */
	itkSetMacro( ExtractBrightLine, bool );
	itkSetMacro( EigenValuesExtraction, bool );
	itkSetMacro( LabelImage, bool );
	itkSetMacro( Sigma, RealType );
	itkSetMacro( OutsideValue, OutputImagePixelType );
	itkSetMacro( InsideValue, OutputImagePixelType );
	itkSetMacro( DimensionsProcessed, unsigned int );

	void EigenValuesExtractionOn();
	void EigenValuesExtractionOff();
	void LabelImageOn();
	void LabelImageOff();

protected:

	/** Filter typedef. */
	typedef HessianRecursiveGaussianImageFilter< InputImageType >			HessianGaussianFilterType;
	typedef BinaryThresholdImageFilter< LabelImageType, OutputImageType >	BinaryThresholdFilterType;
	typedef BinaryThresholdImageFilter< RealImageType, RealImageType >		RealBinaryThresholdFilterType;
	typedef Hessian3DToVesselnessMeasureImageFilter< RealType >			VesselnessFilterType;

	/** Inside pixel typedefs. */
	typedef typename HessianGaussianFilterType::OutputImageType				TensorImageType;
	typedef typename HessianGaussianFilterType::OutputPixelType				TensorType;
	typedef typename TensorType::EigenValuesArrayType						EigenValuesArrayType;
	typedef typename TensorType::EigenVectorsMatrixType						EigenVectorsMatrixType;
	typedef Image< EigenValuesArrayType, ImageDimension >					EigenValuesArrayImageType;


	/** Matrix analysis tool typedef. */
	typedef SymmetricEigenAnalysis< TensorType, EigenValuesArrayType >			EigenAnalysisType;

	/** Functor class. */
	class AbsoluteSortedEigenValuesArrayExtractor
	{
		public:
			AbsoluteSortedEigenValuesArrayExtractor() {};
			~AbsoluteSortedEigenValuesArrayExtractor() {};
			bool operator != ( const AbsoluteSortedEigenValuesArrayExtractor & ) const { return false; }
			bool operator == ( const AbsoluteSortedEigenValuesArrayExtractor & other ) const { return !( *this != other ); }
			inline EigenValuesArrayType operator()( const TensorType & A ) const;

	};

	class SignedBasedEigenValuesClassifier
	{
		public:
			SignedBasedEigenValuesClassifier() { m_DimensionsProcessed = ImageDimension - 1; }
			~SignedBasedEigenValuesClassifier() {};
			bool operator != ( const SignedBasedEigenValuesClassifier & ) const { return false; }
			bool operator == ( const SignedBasedEigenValuesClassifier & other ) const { return *this != other; }
			void SetDimensionsProcessed( unsigned int & dimensionsProcessed ) { m_DimensionsProcessed = dimensionsProcessed; }
			unsigned int GetDimensionsProcessed() const { return m_DimensionsProcessed; }
			inline LabelType operator() ( const EigenValuesArrayType & A ) const;
		protected:
			unsigned int m_DimensionsProcessed;
	};

	class EigenValuesExtractor
	{
		public:
			EigenValuesExtractor() { m_EigenIndex = 0; }
			~EigenValuesExtractor() {};
			bool operator != ( const EigenValuesExtractor & ) const { return false; }
			bool operator == ( const EigenValuesExtractor & other ) const { return *this != other; }
			void SetEigenIndex( unsigned int & eigenIndex ) { m_EigenIndex = eigenIndex; }
			unsigned int GetEigenIndex() const { return m_EigenIndex; }
			inline RealType operator() ( const EigenValuesArrayType & A ) const;
		protected:
			unsigned int m_EigenIndex;
	};
	
	class LineDirectionVectorExtractor
	{
		public:
			LineDirectionVectorExtractor() {}
			~LineDirectionVectorExtractor() {};
			bool operator != ( const LineDirectionVectorExtractor & ) const { return false; }
			bool operator == ( const LineDirectionVectorExtractor & other ) const { return *this != other; }
			inline VectorType operator() ( const TensorType & A ) const;
	};

	class VesselnessMeasure
	{
		public:
			VesselnessMeasure() {}
			~VesselnessMeasure() {};
			bool operator != ( const VesselnessMeasure & ) const { return false; }
			bool operator == ( const VesselnessMeasure & other ) const { return *this != other; }
			inline RealType operator() ( const EigenValuesArrayType & A ) const;
	};

	class NotOperator
	{
		public:
			NotOperator() {};
			~NotOperator() {};
			bool operator != ( const NotOperator & ) const { return false; }
			bool operator == ( const NotOperator & other ) const { return *this != other; }
			inline PixelType operator() ( const PixelType & A ) const
			{
				return ~A;
			}
	};

	/** Functor typedef. */
	typedef UnaryFunctorImageFilter< TensorImageType, EigenValuesArrayImageType, AbsoluteSortedEigenValuesArrayExtractor >	EigenValuesArrayFilterType;
	typedef UnaryFunctorImageFilter< EigenValuesArrayImageType, LabelImageType, SignedBasedEigenValuesClassifier >			LabelFilterType;
	typedef UnaryFunctorImageFilter< EigenValuesArrayImageType, RealImageType, EigenValuesExtractor >						EigenValuesFilterType;
	typedef UnaryFunctorImageFilter< TensorImageType, VectorImageType, LineDirectionVectorExtractor >						LineDirectionFilterType;
	//typedef UnaryFunctorImageFilter< EigenValuesArrayImageType, RealImageType, VesselnessMeasure >							VesselnessFilterType;
	typedef UnaryFunctorImageFilter< ImageType, ImageType, NotOperator >													NotFilterType;


	LineShapeImageFilter();

	/** GenerateData function. */
	void GenerateData();

	void ComputeData();

	DataObject::Pointer MakeOutput( unsigned int idx );

	virtual ~LineShapeImageFilter() {}
//	void PrintSelf( std::ostream & s, Indent indent ) const;

	/* Sets region. */
	void SetRegion( const InputImageRegionType & region );

private:
	LineShapeImageFilter( const Self & );	// purposely not implemented
	void operator=( const Self & );			// purposely not implemented

	InputImageRegionType m_Region;
	bool m_RegionSetByUser;

	/** Label definition. */
	LabelType m_BrightLine;
	LabelType m_DarkLine;
	LabelType m_LabelCount;

	bool m_ExtractBrightLine;

	bool m_EigenValuesExtraction;
	bool m_LabelImage;

	OutputImagePixelType m_OutsideValue;
	OutputImagePixelType m_InsideValue;

	unsigned int m_DimensionsProcessed;

	RealType m_Sigma;

	class Compare
	{
		public:
			Compare() {};
			~Compare() {};
			bool operator()( const RealType x, const RealType y ) const
			{ return std::abs( x ) < std::abs( y ); }
	};

};
}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLineShapeImageFilter.hxx"
#endif

#endif
