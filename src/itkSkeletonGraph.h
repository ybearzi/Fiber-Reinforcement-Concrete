#ifndef __itkSkeletonGraph_h
#define __itkSkeletonGraph_h

#include "Graph/itkGraph.h"

namespace itk
{
template< class TGraphTraits >
class SkeletonGraph : public Graph< TGraphTraits >
{
public:
	/** Standard class typedefs. */
	typedef SkeletonGraph				Self;
	typedef Graph< TGraphTraits >		Superclass;
	typedef SmartPointer< Self >		Pointer;
	typedef SmartPointer< const Self >	ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( Graph, SkeletonGraph );

	/** Some convenient typedefs. */
	typedef TGraphTraits											GraphTraitsType;
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
	typedef typename Superclass::NodeContainerType					NodeContainerType;
	typedef typename Superclass::EdgeContainerType					EdgeContainerType;
	typedef typename GraphTraitsType::NodeWeightType				NodeWeightType;
	typedef typename GraphTraitsType::EdgeWeightType				EdgeWeightType;
	typedef typename GraphTraitsType::EdgeIdentifierContainerType	EdgeIdentifierContainerType;

	static const unsigned int Dimension = IndexType::Dimension;

	virtual inline Pointer Copy() const;
	virtual inline void InitializeGraphAttributes();
	virtual inline void InitializeOutgoingEdgesInteraction( const NodeType & node );
	virtual inline void InitializeNodeInteraction( NodePointerType node );
	virtual inline void InitializeEdge( EdgePointerType edge, const NodeType & source, const NodeType & target );
	virtual inline void InitializeNode( NodePointerType node );
	virtual inline void UpdateNode( NodePointerType node );
	virtual inline void UpdateNodesWeightWithoutInteraction();
	virtual inline NodeWeightType GetFiberLikelihood( const NodeType & node ) const;
	virtual inline void ConnectNodeExtremities();
	virtual inline void SetOutgoingEdgesWeight( const NodeType & node );
	virtual inline NodeWeightType Weight( const NodeType & node ) const;
	virtual inline double f( const NodeType & i, const NodeType & j ) const;
	virtual inline double MinimumDistance( const NodeType & i, const NodeType & j ) const;
	virtual inline double Force( const NodeType & node ) const;
	virtual inline double Force( const NodeType & i, const NodeType & j, const EdgeType & ij ) const;
	virtual inline double GeometricLikelihood( const NodeType & node ) const;
	virtual inline double ConnectionLikelyhood( const NodeType & i, const NodeType & j, const EdgeType & edge ) const;
	virtual inline double VirtualConnectionLikelyhood( const NodeType & i, const NodeType & j ) const;


	template< class TIterator >
	inline VectorType P( TIterator begin, TIterator end ) const;
	template< class TIterator >
	inline VectorType V( TIterator begin, TIterator end, const VectorType & p ) const;
	template< class TIterator >
	inline double AverageDistance( TIterator begin, TIterator end, const VectorType & p, const VectorType & v ) const;
	template< class TIterator >
	inline double Mean( TIterator begin, TIterator end ) const;
	template< class TIterator >
	inline double AverageDistanceVariance( const VectorType & point, const VectorType & vector, double AverageDistance, TIterator begin, TIterator end ) const;

	template< class TContainer, class TIterator >
	inline double MeanCurvatureDerivative( TIterator begin, TIterator end ) const;
	template< class TContainer, class TIterator >
	inline TContainer Diff( TIterator begin, TIterator end ) const;
	template< class TIterator >
	inline double DimensionalMean( TIterator begin, TIterator end ) const;


	inline void FillGaps( double threshold );
	inline NodePointerType GetMaximumNodeWeight( const NodeWeightType & threshold );
	inline IndexContainerType MergeIndexes( const NodeType & i, const NodeType & j, bool ifront, bool jfront ) const;
	inline double Distance( const IndexType & index, const VectorType & p, const VectorType & v ) const;
	inline double TheoricalAverageDistance( const VectorType & v, double norm ) const;
	inline bool IsConnectedToFront( const NodeType & i, const NodeType & j ) const;
	inline bool IsConnectedToBack( const NodeType & i, const NodeType & j ) const;
	inline bool IsLikelyConnectedToFront( const NodeType & i, const NodeType & j ) const;
	inline bool IsLikelyConnectedToBack( const NodeType & i, const NodeType & j ) const;
	inline bool IsOnlyConnectedToFront( const NodeType & i, const NodeType & j ) const;
	inline bool IsOnlyConnectedToBack( const NodeType & i, const NodeType & j ) const;
	inline bool IsSelfConnected( const NodeType & i ) const;
	inline bool AreConnected( const NodeType & i, const NodeType & j ) const;
	inline bool OutgoingEdgesSide( const NodePointerType & i, const NodePointerType & j ) const;
	inline void SynchronizeIncomingOnOutgoingEdges( NodePointerType node );
	inline void DeleteLinks( NodePointerType & i, bool front );
	inline double MaxDistance( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, unsigned int & pos ) const;
	inline bool CanBeMerged( const NodeType & source, const NodeType & target, double threshold ) const;
	inline bool CanBeSplit( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, double threshold, unsigned int & pos ) const;
	inline void Connect( NodePointerType & i, NodePointerType & j, bool iFront, bool jFront );
	inline bool MergeSuccess( NodePointerType node, const double & threshold );
	inline bool Merge( NodePointerType & i, NodePointerType & j, bool ijOutgoingSide, bool jiOutgoingSide, bool reloopingPossible = true );
	inline NodeIdentifierType MaxIntersectionLabel() const;
	NodeIdentifierType TryToSplit( NodePointerType node, double threshold, NodeIdentifierType & intersectionLabel );
	inline NodePointerType Split( NodePointerType node, unsigned int pos, NodeIdentifierType inter );



protected:

	SkeletonGraph(){};
	~SkeletonGraph(){};


private:
	SkeletonGraph( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemeted
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSkeletonGraph.hxx"
#endif

#endif
