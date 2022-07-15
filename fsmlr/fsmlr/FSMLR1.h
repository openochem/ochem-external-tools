// FSMLR1.h: interface for the FSMLR class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FSMLR1_H__13D8C287_E073_408F_92BC_D07F7E44A71D__INCLUDED_)
#define AFX_FSMLR1_H__13D8C287_E073_408F_92BC_D07F7E44A71D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "common.h"

#define cDeleted 0
#define cTraining 1
#define cValidation 2
#define cPrediction 3

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

class FSMLR  
{
public:
	double F5;
	UINT8 * cSet;
	double ** X;
	double *Y;
	int nDescr;
	int nComp;
	void ReduceModel(int numSelDescr);
	int BestNumIter;
	void SelectDescriptors();
	int GetBestNumDescr() const;
	void SelectDescriptors(int nSel);
	void SetMaxDescrInModel(int maxDescr);
	void RunMLR();
	virtual void DescriptorAddedCallback(int, int);
	void CalcSigmaPrediction();
	void CalcSigmaValidation();
	void CalcSigmaTraining();
	void CalcCorrCoef();
	void InitArrays();
	FSMLR(int nComp_, int nDescr_, double *Y_, double **X_, UINT8 *cSet_, int mdd_, double def_, double dnf_,
	  double shrinkage_);
	virtual ~FSMLR();

  
	// Function values that should be fitted
  double *y;

  // Initial values of the current property
  double *y_init;

  // Predicted values 
  double *y_pred;

  // Boolean array showing whether a given point (actually chemical structure)
  // is used for building statistical model
  BOOL *bPointUsedForTraining;

  // The number of points (actually chemical structures) in training set
  // that can be used in building regression model
  int nCompUsedForTraining;

  // Boolean array showing whether a given point (actually chemical structure)
  // is used for validating statistical model
  BOOL *bPointUsedForValidating;

  // The number of points (actually chemical structures) in training set
  // that can be used in validating regression model
  int nCompUsedForValidating;

  // Boolean array showing whether a given point (actually chemical structure)
  // is used for prediction
  BOOL *bPointUsedForPrediction;

  // The number of points (actually chemical structures) in training set
  // that can be used for prediction
  int nCompUsedForPrediction;

  // Boolean array showing whether a given descriptor is selected
  BOOL *bDescrSelected;

  // The number of selected descriptors
  int nDescrSelected;

  // Regression coefficients
  double *b;

  // Regression free coefficients
  double *b0;

  // The number of descriptors in MLR model
  int nDescrInModel;

  // Shows whether a given descriptor is used in regression model
  BOOL *bUsedInModel;

  // Maximal number of descriptors in regression model
  int nMaxDescrInModel;

  // Correlation coefficient
  double CorrCoef;

  // Sigma for training set
  double SigmaTraining;

  // Sigma for validation set
  double SigmaValidation;

  // Sigma for prediction
  double SigmaPrediction;

  // "Best" number of descriptors providing the lowest validation error
  int BestNumDescr;

  // Correlation coefficient for the "best" number of descriptors
  double BestCorrCoef;

  // Sigma for training set for the "best" number of descriptors
  double BestSigmaTraining;

  // Sigma for validation set for the "best" number of descriptors
  double BestSigmaValidation;

  // Sigma for prediction set for the "best" number of descriptors
  double BestSigmaPrediction;

	// Array of descriptor IDs included in model
	int *nDescrIDInModel;

	// Show whether descriptors may be reused in model building
	BOOL bDescriptorReuse;

	// Elongation factor for descriptors reuse
	double EF;

	// Reduced correlation coefficients
	double *c;

	// Reduced intercept
	double c0;

	// Array of descriptor IDs included in reduced model
	int *nDescrIDInReducedModel;
	
	// Maximal difference between a current number of descriptors and
	// the best number of descriptors (stop criteria, default=20)
	int MaxDeltaDescr;

	int best_descr;
	double best_b1;

	void SaveModel ();
	
	double shrinkage;
  
    // The number of unique descriptors
	int numUniqueDescr;

  private:
	// Current property number
  int uCurrentProperty, uSelProperty; 	
};

#endif // !defined(AFX_FSMLR1_H__13D8C287_E073_408F_92BC_D07F7E44A71D__INCLUDED_)
