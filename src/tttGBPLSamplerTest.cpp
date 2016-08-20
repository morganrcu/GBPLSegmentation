
#include <tttGBPL.h>
#include <iostream>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkLaplacianImageFilter.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <tttGBPLDistribution.h>
#include <tttGBPLRegionSynthetizer.h>
int main(int argc,char ** argv){

	typedef ttt::GBPL<double,2> GBPLType;
	GBPLType gbpl;
	typedef ttt::GBPLDistribution<GBPLType> GBPLDistributionType;

	GBPLDistributionType distribution;
	distribution.SetLambda(250.0);
	distribution.SetPLimits(0,0.0,568.0);
	distribution.SetPLimits(1,0.0,568.0);
	gbpl=distribution.Sample();


	typedef itk::Image<float,2> ImageType;
	ImageType::RegionType region;
	region.SetIndex(0,0);
	region.SetIndex(1,0);
	region.SetSize(0,568);
	region.SetSize(1,568);

	typedef ttt::GBPLRegionSynthetizer<GBPLType,ImageType> GBPLRegionExtractorType;
	GBPLRegionExtractorType::Pointer gbplRegionExtractor = GBPLRegionExtractorType::New();
	gbplRegionExtractor->SetGBPL(gbpl);
	gbplRegionExtractor->SetOutputRegion(region);
	gbplRegionExtractor->Update();
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("./gbpl.mha");
	writer->SetInput(gbplRegionExtractor->GetOutput());
	writer->Update();


	ImageType::SizeType radius;
	radius.Fill(1.0);
	itk::NeighborhoodIterator<ImageType> neighIterator(radius, gbplRegionExtractor->GetOutput(),region);


	ImageType::Pointer skeleton = ImageType::New();

	skeleton->SetRegions(region);
	skeleton->Allocate();

	while(!neighIterator.IsAtEnd()){

		bool different=false;

		auto center = neighIterator.GetIndex();
		auto label = *neighIterator.GetCenterValue();
	    // Set the current pixel to white

	    for(unsigned int i = 0; i < neighIterator.Size(); i++)	{
	    	ImageType::IndexType index = neighIterator.GetIndex(i);
	    	auto value =neighIterator.GetPixel(i);
	    	if(label != value){
	    		different=true;
	    	}
	    }
	    if(different){
	    	skeleton->SetPixel(center,float(1.0));
	    }else{
	    	skeleton->SetPixel(center,float(0.0f));
	    }
	    ++neighIterator;
	}



#if 0
	typedef ImageType::IndexType IndexType;
	typedef ImageType::PointType PointType;


	ImageType::Pointer image = ImageType::New();

	image->SetRegions(region);
	image->Allocate();

	typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

	IteratorType iterator(image,region);

	while(!iterator.IsAtEnd()){
		IndexType index=iterator.GetIndex();
		PointType point;
		image->TransformIndexToPhysicalPoint(index,point);
		double minDistance = std::numeric_limits<double>::max();
		unsigned minIndex;

		for(int n=0;n<gbpl.m_N;n++){
			auto diff = (point.GetVnlVector()-gbpl.m_P[n]);
			auto distance=std::sqrt(inner_product(diff,diff.pre_multiply(gbpl.m_M[n].as_matrix()))) - gbpl.m_W[n];

			if(distance<minDistance){
				minIndex=n;
				minDistance=distance;
			}
		}
		iterator.Set(minIndex);
		++iterator;
	}







#endif
}
