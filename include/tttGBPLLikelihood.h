/*
 * tttGBPLLikelihood.h
 *
 *  Created on: Aug 19, 2016
 *      Author: morgan
 */

#ifndef MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLLIKELIHOOD_H_
#define MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLLIKELIHOOD_H_
#include <itkStatisticsImageFilter.h>
#include <tttGBPLRegionSynthetizer.h>
#include <itkNeighborhoodIterator.h>
#include <itkImageFileWriter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkSquareImageFilter.h>
namespace ttt{
template<class TGBPL,class TObservationImage> class GBPLLikelihood{
public:
	typedef TGBPL GBPLType;
	typedef TObservationImage ObservationImageType;

	GBPLLikelihood(){
		m_Sigma2=1.0;
	}

	double GetLogLikelihood(const TGBPL & gbpl){
		typedef ttt::GBPLRegionSynthetizer<GBPLType,ObservationImageType> RegionSynthetizerType;
		typename RegionSynthetizerType::Pointer regionSynthetizer = RegionSynthetizerType::New();
		regionSynthetizer->SetGBPL(gbpl);
		regionSynthetizer->SetOutputRegion(m_ObservedImage->GetLargestPossibleRegion());
		regionSynthetizer->SetOutputSpacing(m_ObservedImage->GetSpacing());
		regionSynthetizer->SetOutputOrigin(m_ObservedImage->GetOrigin());
		regionSynthetizer->Update();


		typename ObservationImageType::SizeType radius;
		radius.Fill(1.0);
		itk::NeighborhoodIterator<ObservationImageType> neighIterator(radius, regionSynthetizer->GetOutput(),m_ObservedImage->GetLargestPossibleRegion());

		typename ObservationImageType::Pointer skeleton = ObservationImageType::New();

		skeleton->SetRegions(m_ObservedImage->GetLargestPossibleRegion());
		skeleton->CopyInformation(m_ObservedImage);
		skeleton->Allocate();

		while (!neighIterator.IsAtEnd()) {

			bool different = false;

			auto center = neighIterator.GetIndex();
			auto label = *neighIterator.GetCenterValue();
			// Set the current pixel to white

			for (unsigned int i = 0; i < neighIterator.Size(); i++) {
				typename ObservationImageType::IndexType index = neighIterator.GetIndex(i);
				auto value = neighIterator.GetPixel(i);
				if (label != value) {
					different = true;
				}
			}
			if (different) {
				skeleton->SetPixel(center, float(1.0));
			} else {
				skeleton->SetPixel(center, float(0.0f));
			}
			++neighIterator;
		}



#if 0

		typedef itk::ImageFileWriter<ObservationImageType> WriterType;
		typename WriterType::Pointer skeletonWriter = WriterType::New();
		skeletonWriter->SetFileName("./skeleton.mha");
		skeletonWriter->SetInput(skeleton);
		skeletonWriter->Update();
#endif

		typedef itk::DiscreteGaussianImageFilter<ObservationImageType,ObservationImageType> BlurType;
		typename BlurType::Pointer blurFilter = BlurType::New();
		blurFilter->SetInput(skeleton);
		blurFilter->SetVariance(6.0);
		blurFilter->SetUseImageSpacing(false);

#if 0
		typedef itk::ImageFileWriter<ObservationImageType> WriterType;

		typename WriterType::Pointer syntethizedWriter = WriterType::New();
		syntethizedWriter->SetFileName("./synthetized.mha");
		syntethizedWriter->SetInput(blurFilter->GetOutput());
		syntethizedWriter->Update();
#endif

		typename NormalizerType::Pointer normalizer = NormalizerType::New();
		normalizer->SetInput(blurFilter->GetOutput());


		typename SubtractorType::Pointer subtractor2 = SubtractorType::New();

		subtractor2->SetInput1(normalizer->GetOutput());
		subtractor2->SetInput2(m_ObservedImage);


		typename SquareType::Pointer square = SquareType::New();
		square->SetInput(subtractor2->GetOutput());
#if 1
		typedef itk::ImageFileWriter<ObservationImageType> WriterType;

		typename WriterType::Pointer squaredWriter = WriterType::New();
		squaredWriter->SetFileName("./squared.mha");
		squaredWriter->SetInput(square->GetOutput());
		squaredWriter->Update();
#endif
		typename StatisticsType::Pointer statistics2 = StatisticsType::New();
		statistics2->SetInput(square->GetOutput());
		statistics2->Update();
		return -statistics2->GetSum()/(2.0*m_Sigma2);
	}
	inline void SetVariance(double variance){
		m_Sigma2=variance;
	}
	void SetObservedImage(const TObservationImage * observation){
		typename NormalizerType::Pointer normalizer = NormalizerType::New();
		normalizer->SetInput(observation);
		normalizer->Update();
		m_ObservedImage=normalizer->GetOutput();
		m_ObservedImage->DisconnectPipeline();

	}
protected:

private:


	typedef itk::NormalizeImageFilter<ObservationImageType,ObservationImageType> NormalizerType;
	typedef itk::SubtractImageFilter<ObservationImageType,ObservationImageType,ObservationImageType> SubtractorType;
	typedef itk::SquareImageFilter<ObservationImageType,ObservationImageType> SquareType;
	typedef itk::StatisticsImageFilter<ObservationImageType> StatisticsType;

	typename ObservationImageType::Pointer m_ObservedImage;

	double m_Sigma2;

};
}


#endif /* MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLLIKELIHOOD_H_ */
