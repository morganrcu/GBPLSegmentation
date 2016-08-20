#include <itkImage.h>
#include <itkImageFileReader.h>
#include <tttRJMCMCGBPLInference.h>


int main(int argc,char ** argv){


	typedef ttt::GBPL<double,2> GBPLType;
	typedef itk::Image<float,2> ObservedImageType;
	typedef itk::ImageFileReader<ObservedImageType> ObservationReaderType;
	typedef ttt::RJMCMCGBPLInference<GBPLType,ObservedImageType> InferenceType;

	ObservationReaderType::Pointer observationReader = ObservationReaderType::New();
	observationReader->SetFileName(argv[1]);
	observationReader->Update();

	InferenceType inference;

	inference.GetGBPLDistribution().SetLambda(400.0);
	inference.SetPLimits(0,0.0,(observationReader->GetOutput()->GetLargestPossibleRegion().GetSize(0)-1)*observationReader->GetOutput()->GetSpacing()[0]);
	inference.SetPLimits(1,0.0,(observationReader->GetOutput()->GetLargestPossibleRegion().GetSize(1)-1)*observationReader->GetOutput()->GetSpacing()[1]);
	inference.GetGBPLDistribution().SetPLimits(0,0.0,(observationReader->GetOutput()->GetLargestPossibleRegion().GetSize(0)-1)*observationReader->GetOutput()->GetSpacing()[0]);
	inference.GetGBPLDistribution().SetPLimits(1,0.0,(observationReader->GetOutput()->GetLargestPossibleRegion().GetSize(1)-1)*observationReader->GetOutput()->GetSpacing()[1]);
	inference.GetLikelihoodModel().SetVariance(0.1);
	inference.GetLikelihoodModel().SetObservedImage(observationReader->GetOutput());

	inference.Infer();
}
