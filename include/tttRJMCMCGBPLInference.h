/*
 * tttRJMCMCGBPLInference.h
 *
 *  Created on: Aug 19, 2016
 *      Author: morgan
 */

#ifndef MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTRJMCMCGBPLINFERENCE_H_
#define MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTRJMCMCGBPLINFERENCE_H_

#include <tttGBPLDistribution.h>
#include <tttGBPLLikelihood.h>
#include <random>

namespace ttt{

template<class TGBPL, class TObservationImage> class RJMCMCGBPLInference{

public:
	typedef TGBPL GBPLType;
	typedef TObservationImage ObservationImageType;

	typedef ttt::GBPLDistribution<GBPLType> GBPLDistributionType;
	typedef ttt::GBPLLikelihood<GBPLType,ObservationImageType> GBPLLikelihoodType;
	double m_Pmin[GBPLType::Dimension];
	double m_Pmax[GBPLType::Dimension];
	void SetPLimits(unsigned int d, double pmin,double pmax){
		m_Pmin[d]=pmin;
		m_Pmax[d]=pmax;
	}
	void Infer(){

		TGBPL init = m_GBPLDistribution.Sample();
		double logprior =m_GBPLDistribution.LogP(init);
		double loglikelihood=m_GBPLLikelihood.GetLogLikelihood(init);
		double logposterior = logprior + loglikelihood;

		m_MarkovChain.push_back(init);
		m_LogPosteriors.push_back(logposterior);

		GBPLType current  =init;
		double currentLogPosterior = logposterior;
		m_NumberOfAcceptedJumps=0;
		//while(!HasConverged()){
		while(true){

			GBPLType next = Jump(current);
			double nextLogPosterior = m_GBPLDistribution.LogP(next) + m_GBPLLikelihood.GetLogLikelihood(next);

			double acceptanceRatio = std::min(1.0,std::exp(nextLogPosterior-currentLogPosterior));

			std::uniform_real_distribution<double> uniform;
			double acceptance =uniform(m_Generator);

			if(acceptance<acceptanceRatio){
				std::cout << "Accept " << nextLogPosterior << std::endl;
				//m_MarkovChain.push_back(next);
				//m_LogPosteriors.push_back(nextLogPosterior);
				m_NumberOfAcceptedJumps++;
				current=next;
				currentLogPosterior=nextLogPosterior;
			}else{
				std::cout << "Reject " << nextLogPosterior << std::endl;

			}
		}
	}


	bool HasConverged(){

		return m_NumberOfAcceptedJumps>=m_MaxJumps;
	}
	RJMCMCGBPLInference() {
		m_MaxJumps=100000000;
		m_NumberOfAcceptedJumps=0;
	}
protected:

	std::default_random_engine m_Generator;

	GBPLType Jump(const GBPLType & current){
		GBPLType next=current;

		//1. generate random number
		std::uniform_int_distribution<unsigned> uniform(0,2);

		unsigned move = uniform(m_Generator);

		switch(move){
		case 0: // Birth
		{
			std::cout << "\t Add atom" << std::endl;
			auto newAtom = m_GBPLDistribution.SampleAtom();
			next.AddAtom(newAtom);
			break;

		}

		case 1: //Death
		{
			std::cout << "\t Remove atom" << std::endl;
			std::uniform_int_distribution<unsigned> uniform(0,current.GetNumberOfAtoms()-1);
			unsigned selected = uniform(m_Generator);
			auto toRemove = next.BeginAtoms();
			std::advance(toRemove,selected);
			next.RemoveAtom(toRemove);

			break;
		}
		case 2: //Move
		{
			std::cout << "\t Move" << std::endl;
			std::uniform_int_distribution<unsigned> uniform(0,current.GetNumberOfAtoms()-1);

			unsigned selected = uniform(m_Generator);
			auto toModify = next.BeginAtoms();
			std::advance(toModify,selected);
			std::normal_distribution<double> gaussian(0,0.5);

			std::uniform_int_distribution<unsigned> uniform2(0,1);

			unsigned mode = uniform2(m_Generator);
			switch(mode){
			case 0:
			{
				std::cout << "\t\tCentroid" << std::endl;
				toModify->first[0]+=gaussian(m_Generator);
				if(toModify->first[0]<m_Pmin[0]){
					toModify->first[0]=m_Pmin[0];
				}
				if(toModify->first[0]>m_Pmax[0]){
					toModify->first[0]=m_Pmax[0];
				}

				toModify->first[1]+=gaussian(m_Generator);
				if(toModify->first[1]<0){
					toModify->first[1]=0;
				}
				if(toModify->first[1]>m_Pmax[1]){
					toModify->first[1]=m_Pmax[1];
				}
				break;
			}
			case 1:
			{
				std::cout << "\t\tGrowth" << std::endl;
				toModify->second.m_W=m_GBPLDistribution.SampleW();
				toModify->second.m_M=m_GBPLDistribution.SampleM();
				break;
			}
			}


			#if 0
			vnl_symmetric_eigensystem<double> eigs(toModify->second.m_M);
			eigs.D[0]+=std::max(gaussian(m_Generator),0.0);
			eigs.D[1]+=std::max(gaussian(m_Generator),0.0);
			auto matrix =eigs.recompose();
			typedef typename GBPLType::SymmetricMatrixType SymmetricMatrixType;
			toModify->second.m_M=0.5*(matrix + matrix.transpose());
			#endif

			break;
		}

		}
		return next;
	}
public:
	GBPLLikelihoodType & GetLikelihoodModel(){
		return m_GBPLLikelihood;
	}
	GBPLDistributionType & GetGBPLDistribution(){
		return m_GBPLDistribution;
	}

private:


	typename ObservationImageType::Pointer m_Observation;

	GBPLDistributionType m_GBPLDistribution;
	GBPLLikelihoodType m_GBPLLikelihood;

	std::vector<GBPLType> m_MarkovChain;
	std::vector<double> m_LogPosteriors;

	unsigned int m_MaxJumps;
	unsigned int m_NumberOfAcceptedJumps;
};


};



#endif /* MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTRJMCMCGBPLINFERENCE_H_ */
