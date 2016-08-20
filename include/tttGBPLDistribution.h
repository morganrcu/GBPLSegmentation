/*
 * tttGBPLSampler.h
 *
 *  Created on: Aug 17, 2016
 *      Author: morgan
 */

#ifndef MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLDISTRIBUTION_H_
#define MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLDISTRIBUTION_H_

#include <random>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_sym_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_qr.h>
#include <tttGBPL.h>
#include <iostream>
#include <cmath>
namespace ttt{

template<class TGBPL> class GBPLDistribution {

public:
	typedef TGBPL GBPLType;
	typedef typename TGBPL::RealType RealType;
	static const unsigned int Dimension = TGBPL::Dimension;
	typedef typename TGBPL::PointType PointType;
	typedef typename TGBPL::SymmetricMatrixType SymmetricMatrixType;
	typedef vnl_matrix<RealType> MatrixType;
	typedef vnl_diag_matrix<RealType> DiagonalMatrixType;

	double PoissonLogLikelihood(const RealType x, const RealType lambda){
		return x*std::log(lambda) -std::lgamma(x+1) -lambda;
	}
	double GammaLogLikelihood(const RealType x, const RealType alpha, const RealType beta){
		return (alpha-1)*std::log(x) - x/beta  - std::lgamma(alpha) - alpha*log(beta);
	}
	double LogP(const GBPLType & sample){
		double total =0;

		total+= PoissonLogLikelihood(sample.GetNumberOfAtoms(), m_NDistribution.mean());

		total+=std::accumulate(sample.BeginAtoms(),sample.EndAtoms(),0.0,[&](double result,const typename GBPLType::AtomType & atom){
			total+=GammaLogLikelihood(atom.second.m_W,m_WDistribution.alpha(),m_WDistribution.beta());

			vnl_symmetric_eigensystem<double> eigensystem(atom.second.m_M);
			for(int d=0;d<Dimension;d++){
				total+= GammaLogLikelihood(eigensystem.D(d),m_HDistribution.alpha(),m_HDistribution.beta());
			}
			return result;
		});

		return total;
	}
	MatrixType SampleM(){
		vnl_matrix<RealType> zeroMatrix(Dimension,Dimension);
		zeroMatrix.fill(0.0);
		vnl_symmetric_eigensystem<RealType> eigensystem(zeroMatrix);
		eigensystem.V=SampleOrthogonal();
		eigensystem.D=SampleCoefficients();
		MatrixType matrix = eigensystem.recompose();
		return (matrix+matrix.transpose())/2.0;

	}
	RealType SampleW(){
		return m_WDistribution(m_Generator);
	}
	typename TGBPL::AtomType SampleAtom(){


		//Sample centroid
		typename TGBPL::PointType point(Dimension);
		for(int d=0;d<Dimension;d++){

			point[d]=m_PDistribution[d](m_Generator);
		}


		typename TGBPL::GBPLMarkType mark;
		mark.m_M=SampleM();
		mark.m_W=SampleW();

		return std::make_pair(point,mark);

	}
	TGBPL Sample(){

		TGBPL result;


		unsigned N = m_NDistribution(m_Generator);

		for(int n=0;n<N;n++){
			auto newAtom = SampleAtom();
			result.AddAtom(newAtom);
		}
		return result;
	}
#if 0
	RealType SampleTruncatedNormal(const RealType & mean, const RealType & sigma ){
		bool valid =false;
		std::normal_distribution<RealType> normal_distribution(mean,sigma);
		RealType sample=0;
		while(!valid){
			sample = normal_distribution(m_Generator);
			valid = sample>0;
		}
		return sample;
	}
#endif

	DiagonalMatrixType SampleCoefficients(){
		DiagonalMatrixType result(Dimension);
		for(int d=0;d<Dimension;d++){
			result(d)=m_HDistribution(m_Generator);
		}
		return result;
	}
	MatrixType SampleOrthogonal(){
		MatrixType matrix(Dimension,Dimension);
		std::normal_distribution<RealType> normal_distribution;
		for(int d0=0;d0<Dimension;d0++){
			for(int d1=0;d1<Dimension;d1++){
				matrix(d0,d1)=normal_distribution(m_Generator);
			}
		}
		return vnl_qr<RealType>(matrix).Q();
	}
	void SetLambda(const RealType  & lambda){
		m_NDistribution = std::poisson_distribution<unsigned>(lambda);
	}
	RealType GetLambda(){
		return m_NDistribution.mean();

	}
	void SetPLimits(const unsigned int d, const RealType & min,const RealType & max){
		m_PDistribution[d]= std::uniform_real_distribution<RealType>(min,max);
	}

	GBPLDistribution(){
		m_Generator= std::default_random_engine(0);
	}
protected:

private:
	std::poisson_distribution<unsigned> m_NDistribution;
	std::uniform_real_distribution<RealType> m_PDistribution[Dimension];
	std::gamma_distribution<RealType> m_WDistribution;
	std::gamma_distribution<RealType> m_HDistribution;
	std::default_random_engine m_Generator;

};

}


#endif /* MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPLDISTRIBUTION_H_ */
