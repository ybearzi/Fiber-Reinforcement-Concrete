#ifndef __itkLineGraphTraits_h
#define __itkLineGraphTraits_h

#include "Graph/itkDefaultGraphTraits.h"
#include <vector>
#include <list>

namespace itk
{

template< unsigned int VDimension, typename TNodeWeight = float, typename TEdgeWeight = float >
class LineGraphTraits : public DefaultGraphTraits< TNodeWeight, TEdgeWeight >
{
public:
	typedef LineGraphTraits							Self;
	typedef DefaultGraphTraits< TNodeWeight, TEdgeWeight >	Superclass;
	
	typedef Index< VDimension > 					IndexType;
	typedef TNodeWeight								NodeWeightType;
	typedef TEdgeWeight								EdgeWeightType;
	typedef typename Superclass::NodeIdentifierType	NodeIdentifierType;
	typedef typename Superclass::EdgeIdentifierType EdgeIdentifierType;
	typedef typename Superclass::EdgeIdentifierContainerType	EdgeIdentifierContainerType;

	static const unsigned int Dimension = VDimension;


	typedef std::list< IndexType >					IndexContainerType;
	typedef std::list< double >						ScalarContainerType;
	typedef std::vector< EdgeIdentifierType >		EdgeExtremityIdentifierContainerType;
	typedef Vector< float, VDimension >				VectorType;
	typedef VectorType								PointType;
	typedef typename IndexType::OffsetType			OffsetType;


	struct NodeType
	{
		NodeIdentifierType Identifier;
		NodeWeightType Weight;
		NodeWeightType Weight0;
		VectorType Vector;
		PointType Point;
		double AverageDistance;
		double TheoricalAverageDistance;
		double Force;
		double MeanVesselness;
		EdgeIdentifierContainerType IncomingEdges;
		EdgeIdentifierContainerType OutgoingEdges;
		EdgeExtremityIdentifierContainerType IncomingEdgesFront;
		EdgeExtremityIdentifierContainerType IncomingEdgesBack;
		EdgeExtremityIdentifierContainerType OutgoingEdgesFront;
		EdgeExtremityIdentifierContainerType OutgoingEdgesBack;
		IndexContainerType LineIndexes;
		ScalarContainerType VesselnessValues;
		NodeIdentifierType BackLabel;
		NodeIdentifierType FrontLabel;
	};

	typedef NodeType* NodePointerType;

	struct EdgeType
	{
		EdgeIdentifierType Identifier;
		NodeIdentifierType SourceIdentifier;
		NodeIdentifierType TargetIdentifier;
		EdgeIdentifierType ReverseEdgeIdentifier;
		EdgeWeightType Weight;
		float f;
//		float lambda;
		VectorType Vector;
		PointType Point;
		double TheoricalAverageDistance;
		double AverageDistance;
		//float LinkWeight;
		IndexContainerType LineIndexes;
	};

	typedef EdgeType* EdgePointerType;

};

}// namesapce itk

#endif
