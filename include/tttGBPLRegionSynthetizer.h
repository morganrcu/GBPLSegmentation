/*
 * tttGBPLRegions.h
 *
 *  Created on: Aug 18, 2016
 *      Author: morgan
 */

#ifndef MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLREGIONSYNTHETIZER_H_
#define MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLREGIONSYNTHETIZER_H_
#include <tttGBPL.h>
#include <itkImage.h>
#include <itkLevelSet.h>
#include <queue>
namespace ttt {
template<class TGBPL, class TClusterImage> class GBPLRegionSynthetizer: public itk::ImageSource<TClusterImage> {

public:
	typedef GBPLRegionSynthetizer Self;
	typedef itk::ImageSource<TClusterImage> Superclass;
	typedef itk::SmartPointer<GBPLRegionSynthetizer> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(Self,Superclass)

	typedef TGBPL GBPLType;
	typedef itk::Image<float,2> LevelSetImageType;
	typedef LevelSetImageType::PointType PointType;
	/** Typedef support of level set method types. */
	typedef itk::LevelSetTypeDefault<LevelSetImageType> LevelSetType;
	typedef typename LevelSetType::LevelSetPointer LevelSetPointer;
	typedef typename LevelSetType::PixelType PixelType;
	typedef typename LevelSetType::NodeType NodeType;
	typedef typename NodeType::IndexType NodeIndexType;
	typedef typename LevelSetType::NodeContainer NodeContainer;
	typedef typename LevelSetType::NodeContainerPointer NodeContainerPointer;

	typedef typename LevelSetImageType::SizeType OutputSizeType;
	typedef typename LevelSetImageType::RegionType OutputRegionType;
	typedef typename LevelSetImageType::SpacingType OutputSpacingType;
	typedef typename LevelSetImageType::DirectionType OutputDirectionType;
	typedef typename LevelSetImageType::PointType OutputPointType;
	typedef typename LevelSetImageType::OffsetType OffsetType;

	typedef itk::Image<unsigned char, 2> LabelImageType;
	typedef typename LabelImageType::Pointer LabelImagePointer;


	typedef TClusterImage ClusterImageType;
	/** Enum of Fast Marching algorithm point types. FarPoints represent far
	 * away points; TrialPoints represent points within a narrowband of the
	 * propagating front; and AlivePoints represent points which have already
	 * been processed. */
	enum LabelType {
		FarPoint = 0, AlivePoint, TrialPoint, InitialTrialPoint, OutsidePoint
	};

	class AxisNodeType: public NodeType {
	public:
		int GetAxis() const {
			return m_Axis;
		}
		void SetAxis(int axis) {
			m_Axis = axis;
		}
		const AxisNodeType & operator=(const NodeType & node) {
			this->NodeType::operator=(node);
			return *this;
		}

	private:
		int m_Axis;
	};

	void SetGBPL(const GBPLType & gbpl){
		m_GBPL=gbpl;
		this->Modified();
	}

	void SetOutputRegion(const OutputRegionType & region){
		m_Region=region;
		this->Modified();
	}
	void SetOutputSpacing(const OutputSpacingType & spacing){
		m_Spacing=spacing;
	}
	void SetOutputOrigin(const OutputPointType & origin){
		m_Origin=origin;
	}
protected:

	virtual void GenerateOutputInformation() ITK_OVERRIDE{
		this->GetOutput()->SetLargestPossibleRegion(m_Region);
		this->GetOutput()->SetSpacing(m_Spacing);
		this->GetOutput()->SetOrigin(m_Origin);
	}
	virtual void Initialize(ClusterImageType *output) {
		// allocate memory for the output buffer

		OutputRegionType requestedRegion = output->GetRequestedRegion();
		output->SetBufferedRegion(requestedRegion);
		output->Allocate();

		// cache some buffered region information

		m_StartIndex = requestedRegion.GetIndex();
		m_LastIndex = m_StartIndex + requestedRegion.GetSize();
		typename LevelSetImageType::OffsetType offset;
		offset.Fill(1);
		m_LastIndex -= offset;

		m_LabelImage = LabelImageType::New();
		// allocate memory for the PointTypeImage
		m_LabelImage->CopyInformation(output);
		m_LabelImage->SetBufferedRegion(output->GetBufferedRegion());
		m_LabelImage->Allocate();
		// set all points type to FarPoint
		m_LabelImage->FillBuffer(FarPoint);

		m_LevelSetImage =LevelSetImageType::New();
		m_LevelSetImage->CopyInformation(output);
		m_LevelSetImage->SetBufferedRegion(output->GetBufferedRegion());
		m_LevelSetImage->Allocate();
		// set all output value to infinity
		m_LevelSetImage->FillBuffer(m_LargeValue);


		// process input alive points
		AxisNodeType node;
		NodeIndexType idx;

		// make sure the heap is empty
		while (!m_TrialHeap.empty()) {
			m_TrialHeap.pop();
		}

		unsigned int n=0;

		for(auto it = m_GBPL.BeginAtoms();it!=m_GBPL.EndAtoms();++it){
			PointType center;
			center[0]=it->first[0];
			center[1]=it->first[1];
			m_LabelImage->TransformPhysicalPointToIndex(center,idx);
			m_LabelImage->SetPixel(idx, InitialTrialPoint);
			m_LevelSetImage->SetPixel(idx,-it->second.m_W);
			this->GetOutput()->SetPixel(idx, n);
			node.SetIndex(idx);
			node.SetValue(-it->second.m_W);
			m_TrialHeap.push(node);
			n++;
		}

	}

	void UpdateNeighbors(const NodeIndexType & index, ClusterImageType *output) {

		unsigned char label;

		int cluster = output->GetPixel(index);

		PointType point;
		output->TransformIndexToPhysicalPoint(index,point);
		OffsetType offset;
		for(int i=-1;i<=1;i++){
			offset[0]=i;
			for(int j=-1;j<=1;j++){
				offset[1]=j;
				if(i==j) continue;

				NodeIndexType neighIndex= index + offset;
				if(output->GetBufferedRegion().IsInside(neighIndex)){
					PointType neighPoint;
					output->TransformIndexToPhysicalPoint(neighIndex,neighPoint);
					typename GBPLType::VectorType diff = (neighPoint-point).GetVnlVector();

					auto currentAtom = m_GBPL.BeginAtoms();
					std::advance(currentAtom,cluster);
					auto ip=inner_product(diff,diff.pre_multiply(currentAtom->second.m_M));
					auto distance=std::sqrt(ip);


					auto currentDistance=m_LevelSetImage->GetPixel(neighIndex);
					auto proposedDistance=m_LevelSetImage->GetPixel(index) + distance;
					if(proposedDistance< currentDistance){
						AxisNodeType node;
						node.SetIndex(neighIndex);
						node.SetValue(proposedDistance);
						m_LevelSetImage->SetPixel(neighIndex,proposedDistance);
						output->SetPixel(neighIndex,cluster);
						m_LabelImage->SetPixel(neighIndex,TrialPoint);
						m_TrialHeap.push(node);
					}
				}
			}
		}

	}

	virtual void GenerateData() {
		if (m_NormalizationFactor < vnl_math::eps) {
			itk::ExceptionObject err(__FILE__, __LINE__);
			err.SetLocation(ITK_LOCATION);
			err.SetDescription("Normalization Factor is null or negative");
			throw err;
		}

		typename ClusterImageType::Pointer output = this->GetOutput();

		this->Initialize(output);

		// process points on the heap
		AxisNodeType node;
		double oldProgress = 0;

		this->UpdateProgress(0.0);   // Send first progress event

		// CACHE
		while (!m_TrialHeap.empty()) {
			// get the node with the smallest value
			node = m_TrialHeap.top();
			m_TrialHeap.pop();


			// does this node contain the current value ?
			m_CurrentValue = static_cast<double>(m_LevelSetImage->GetPixel(node.GetIndex()));
			auto cluster = this->GetOutput()->GetPixel(node.GetIndex());
			//std::cout << node.GetIndex() << "\t" << m_CurrentValue <<  "\t" << cluster << std::endl;

			if (node.GetValue() == m_CurrentValue) {
				// is this node already alive ?
				if (m_LabelImage->GetPixel(node.GetIndex()) != AlivePoint) {

					// set this node as alive
					m_LabelImage->SetPixel(node.GetIndex(), AlivePoint);

					// update its neighbors
					this->UpdateNeighbors(node.GetIndex(),  output);

					// Send events every certain number of points.
					const double newProgress = m_CurrentValue / m_StoppingValue;
					if (newProgress - oldProgress > 0.01)  // update every 1%
							{
						this->UpdateProgress(newProgress);
						oldProgress = newProgress;
						if (this->GetAbortGenerateData()) {
							this->InvokeEvent(itk::AbortEvent());
							this->ResetPipeline();
							itk::ProcessAborted e(__FILE__, __LINE__);
							e.SetDescription("Process aborted.");
							e.SetLocation(ITK_LOCATION);
							throw e;
						}
					}
				}
			}
		}


	}


#if 0
	typename itk::DataObject::Pointer MakeOutput(unsigned int idx){
	  itk::DataObject::Pointer output;

	  switch ( idx )
	    {
	    case 0:
	      output = ( LevelSetImageType::New() ).GetPointer();
	      break;

	    default:
	      std::cerr << "No output " << idx << std::endl;
	      output = NULL;
	      break;
	    }
	  return output.GetPointer();
	}
#endif
	GBPLRegionSynthetizer() {
		m_StoppingValue=10;
		m_CurrentValue=itk::NumericTraits<PixelType>::min();
		m_NormalizationFactor=1;
		m_LargeValue=static_cast< PixelType >( itk::NumericTraits< PixelType >::max() / 2.0 );

		m_Spacing.Fill(1.0);
		m_Origin.Fill(0.0);
		this->SetNumberOfRequiredOutputs(1);

	}
	virtual ~GBPLRegionSynthetizer() {

	}

private:
	typedef std::vector<AxisNodeType> HeapContainer;
	typedef std::greater<AxisNodeType> NodeComparer;
	typedef std::priority_queue<AxisNodeType, HeapContainer, NodeComparer> HeapType;

	HeapType m_TrialHeap;

	typename LabelImageType::Pointer m_LabelImage;

	typename LevelSetImageType::Pointer m_LevelSetImage;


	double m_NormalizationFactor;

	double m_StoppingValue;
	double m_CurrentValue;
	NodeIndexType m_StartIndex;
	NodeIndexType m_LastIndex;

	OutputRegionType m_Region;
	OutputSpacingType m_Spacing;
	OutputPointType m_Origin;
	double m_LargeValue;

	AxisNodeType m_NodesUsed[2];

	GBPLType m_GBPL;


};
}




#endif /* MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLREGIONSYNTHETIZER_H_ */
