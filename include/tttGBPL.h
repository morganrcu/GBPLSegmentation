/*
 * tttGBPL.h
 *
 *  Created on: Aug 17, 2016
 *      Author: morgan
 */

#ifndef MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPL_H_
#define MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPL_H_

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_sym_matrix.h>
#include <map>
#include <algorithm>
namespace ttt{

template<class TReal, int dim> struct VectorComparator : public std::binary_function<vnl_vector_fixed<TReal,dim>, vnl_vector_fixed<TReal,dim>, bool>{
	bool operator()(const vnl_vector_fixed<TReal,dim> & a,const vnl_vector_fixed<TReal,dim> & b) const{
		return std::lexicographical_compare(a.begin(),a.end(),b.begin(),b.end());
	}
};

template<class TReal,unsigned int dim> class GBPL{
public:
	const static unsigned int Dimension = dim;
	typedef TReal RealType;

	typedef vnl_vector<TReal> PointType;
	typedef vnl_vector<TReal> VectorType;
	typedef vnl_matrix<TReal> SymmetricMatrixType;
	typedef struct GBPLMark{
		SymmetricMatrixType m_M;
		TReal m_W;
	}GBPLMarkType;

	//typedef std::map<PointType,GBPLMarkType,VectorComparator<TReal,dim>> AtomContainerType;

	//typedef typename AtomContainerType::value_type AtomType;
	typedef std::pair<PointType,GBPLMarkType> AtomType;
	typedef std::vector<AtomType> AtomContainerType;
	typedef typename AtomContainerType::iterator AtomIterator;
	typedef typename AtomContainerType::const_iterator AtomConstIterator;


	void AddAtom(const AtomType & atom){

		m_Atoms.push_back(atom);
	}

	AtomIterator BeginAtoms(){
		return m_Atoms.begin();
	}
	AtomIterator EndAtoms(){
		return m_Atoms.end();
	}
	AtomConstIterator BeginAtoms() const{
		return m_Atoms.begin();
	}
	AtomConstIterator EndAtoms() const{
		return m_Atoms.end();
	}
	void RemoveAtom(const AtomIterator & atom){
		m_Atoms.erase(atom);
	}

	unsigned int GetNumberOfAtoms() const{
		return m_Atoms.size();
	}


private:
	AtomContainerType m_Atoms;
};

}


#endif /* MODULES_FUNCTIONS_IMAGEFUNCTIONS_GBPLSEGMENTATION_INCLUDE_TTTGBPL_H_ */
