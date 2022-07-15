// FSMLR1.cpp: implementation of the FSMLR class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>

#include "FSMLR1.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FSMLR::FSMLR(int nComp_, int nDescr_, double *Y_, double **X_, UINT8 *cSet_, int mdd_, double def_, double dnf_,
  double shrinkage_)
{
  nComp = nComp_;
  nDescr = nDescr_;
  Y = Y_;
  X = X_;
  cSet = cSet_;
  bDescriptorReuse = TRUE;
  EF = def_;
  F5 = dnf_;
  MaxDeltaDescr = mdd_;
  shrinkage = shrinkage_;
  y = new double[nComp];
  y_init = new double[nComp];
  y_pred = new double[nComp];
  bPointUsedForTraining = new BOOL[nComp];
  bPointUsedForValidating = new BOOL[nComp];
  bPointUsedForPrediction = new BOOL[nComp];
  bDescrSelected = new BOOL[nDescr];
  b = new double[int(nDescr * EF)+10];
  b0 = new double[int(nDescr * EF)+10];
  c0 = 0;
  bUsedInModel = new BOOL[nDescr];
  nDescrIDInModel = new int[int(nDescr * EF)+10];
  c = new double[nDescr];
  nDescrIDInReducedModel = new int[nDescr];
  uCurrentProperty = 0;
  uSelProperty = 0;
  InitArrays();
}

FSMLR::~FSMLR()
{
  delete [] bUsedInModel;
  delete [] y;
  delete [] y_init;
  delete [] y_pred;
  delete [] bPointUsedForTraining;
  delete [] bPointUsedForValidating;
  delete [] bPointUsedForPrediction;
  delete [] bDescrSelected;
  delete [] b;
  delete [] b0;
  delete [] nDescrIDInModel;
  delete [] c;
  delete [] nDescrIDInReducedModel;
}

void FSMLR::InitArrays()
{
  int i;

  // Init y
  for (i = 0; i < nComp;  i++)
  {
    y_init[i] = y[i] = Y[i];
    y_pred[i] = 0.0;
  }

  // Init bPointUsedForTraining
  nCompUsedForTraining = 0;
  for (i = 0; i < nComp; i++)
	if (cSet[i] == cTraining)
	{
	  bPointUsedForTraining[i] = TRUE;
	  nCompUsedForTraining++;
	}
	else
      bPointUsedForTraining[i] = FALSE;

  // Init bPointUsedForValidating
  nCompUsedForValidating = 0;
  for (i = 0; i < nComp; i++)
	if (cSet[i] == cValidation)
	{
	  bPointUsedForValidating[i] = TRUE;
	  nCompUsedForValidating++;
	}
	else
    bPointUsedForValidating[i] = FALSE;

  // Init bPointUsedForPrediction
  nCompUsedForPrediction = 0;
  for (i = 0; i < nComp; i++)
	if (cSet[i] == cPrediction)
	{
	  bPointUsedForPrediction[i] = TRUE;
	  nCompUsedForPrediction++;
	}
	else
      bPointUsedForPrediction[i] = FALSE;

  // Init bDescrSelected
  nDescrSelected = nDescr;
  for (i = 0; i < nDescr; i++)
    bDescrSelected[i] = TRUE;

  // Init b
  for (i = 0; i < nDescr * EF; i++)
    b[i] = 0.0;

  // Init nDescrInModel
  nDescrInModel = 0;

  // Init bUsedInModel;
  for (i = 0; i < nDescr; i++)
    bUsedInModel[i] = FALSE;

  // Init nMaxDescrInModel
  nMaxDescrInModel = min(nDescrSelected, int(nComp / F5) - 1);

}

void FSMLR::CalcCorrCoef()
{
  double SumX, SumX2, SumXY, SumY, SumY2;
  double AveX, AveX2, AveXY, AveY, AveY2;
  double SigmaX, SigmaY;
  int i;

  SumX = SumX2 = SumXY = SumY = SumY2 = 0.0f;
  for (i = 0; i < nComp; i++)
    if (bPointUsedForTraining[i])
	{
	  SumX  += y_pred[i];
	  SumX2 += y_pred[i] * y_pred[i];
	  SumXY += y_pred[i] * y_init[i];
	  SumY  += y_init[i];
	  SumY2 += y_init[i] * y_init[i];
	}
  AveX  = SumX  / nCompUsedForTraining;
  AveX2 = SumX2 / nCompUsedForTraining;
  AveXY = SumXY / nCompUsedForTraining;
  AveY  = SumY  / nCompUsedForTraining;
  AveY2 = SumY2 / nCompUsedForTraining;
  if (AveX2 - AveX * AveX > 0.0)
    SigmaX = (double)sqrt(AveX2 - AveX * AveX);
  else
    SigmaX = 0.0;
  if (AveY2 - AveY * AveY > 0.0)
    SigmaY = (double)sqrt(AveY2 - AveY * AveY);
  else
    SigmaY = 0.0;
  if (SigmaX > 0.00001 && SigmaY > 0.00001)
    CorrCoef = (AveXY - AveX * AveY) / (SigmaX * SigmaY);
  else
	  CorrCoef = 0.0;

}

void FSMLR::CalcSigmaTraining()
{
  double SumX, SumX2, AveX, AveX2;
  int i;

  SumX = SumX2 = 0.0;
  for (i = 0; i < nComp; i++)
    if (bPointUsedForTraining[i])
	{
	  SumX2 += (y_pred[i] - y_init[i]) * (y_pred[i] - y_init[i]);
	}
  AveX2 = SumX2 / nCompUsedForTraining;
  SigmaTraining = sqrt(AveX2);
}

void FSMLR::CalcSigmaValidation()
{
  double SumX, SumX2, AveX, AveX2;
  int i;

  SumX2 = 0.0f;
  for (i = 0; i < nComp; i++)
    if (bPointUsedForValidating[i])
	{
	  SumX2 += (y_pred[i] - y_init[i]) * (y_pred[i] - y_init[i]);
	}
  AveX2 = SumX2 / nCompUsedForValidating;
  SigmaValidation = sqrt(AveX2);

}

void FSMLR::CalcSigmaPrediction()
{
  double SumX, SumX2, AveX, AveX2;
  int i;

  SumX2 = 0.0f;
  for (i = 0; i < nComp; i++)
    if (bPointUsedForPrediction[i])
	{
	  SumX2 += (y_pred[i] - y_init[i]) * (y_pred[i] - y_init[i]);
	}
  AveX2 = SumX2 / nCompUsedForPrediction;
  SigmaPrediction = sqrt(AveX2);

}

void FSMLR::DescriptorAddedCallback(int iter, int ides)
{
  if (nCompUsedForValidating && nCompUsedForPrediction) {
    if (bDescriptorReuse)
      printf("%4d(%d)   %.4f   %.4f   %.4f   %.4f  %8.4f  %8.4f  descr_%d\n", 
             iter, ides, CorrCoef,
             SigmaTraining, SigmaValidation, SigmaPrediction,
             best_b1, b[0], best_descr);
    else
      printf("%4d   %.4f   %.4f   %.4f   %.4f  %8.4f  %8.4f  descr_%d\n", 
		     iter, CorrCoef,
             SigmaTraining, SigmaValidation, SigmaPrediction,
             best_b1, b[0], best_descr);
  } else if (nCompUsedForValidating) {
    if (bDescriptorReuse)
      printf("%4d(%d)   %.4f   %.4f   %.4f  -  %8.4f  %8.4f  descr_%d\n", 
             iter, ides, CorrCoef, 
             SigmaTraining, SigmaValidation,
             best_b1, b[0], best_descr);
    else
      printf("%4d   %.4f   %.4f   %.4f  -  %8.4f  %8.4f  descr_%d\n", 
             iter, CorrCoef, 
             SigmaTraining, SigmaValidation,
             best_b1, b[0], best_descr);
  } else {
    if (bDescriptorReuse)
      printf("%4d(%d)   %.4f   %.4f   -  -  %8.4f  %8.4f  descr_%d\n", 
             iter, ides, CorrCoef, 
             SigmaTraining, 
             best_b1, b[0], best_descr);
    else
      printf("%4d   %.4f   %.4f   -  -  %8.4f  %8.4f  descr_%d\n", 
             iter, CorrCoef, 
             SigmaTraining, 
             best_b1, b[0], best_descr);
  }
}

void FSMLR::RunMLR()
{
  int i, j, k;
	int ti; // "true value of i"
  double SumX, SumX2, SumXY, SumY, SumY2;
  double AveXY, AveY, AveY2;
  double *AveX, *AveX2, *SigmaX, SigmaY;
  double curr_corr_coef, curr_R2, best_R2;
  double curr_b0, best_b0;
  double curr_b1;
  double best_SigmaY;

	InitArrays();

  AveX = new double[nDescr];
  AveX2 = new double[nDescr];
  SigmaX = new double[nDescr];

  // Compute SigmaX[]
  for (j = 0; j < nDescr; j++)
    if (bDescrSelected[j])
	{
	  SumX = SumX2 = 0.0;
	  for (k = 0; k < nComp; k++)
	    if (bPointUsedForTraining[k])
		{
		  double x = X[k][j];
		  SumX += x;
		  SumX2 += x * x;
		}
	  AveX[j] = SumX / nCompUsedForTraining;
	  AveX2[j] = SumX2 / nCompUsedForTraining;
	  if (AveX2[j] - AveX[j] * AveX[j] > 0.0f)
	    SigmaX[j] = (double)sqrt(AveX2[j] - AveX[j] * AveX[j]);
	  else
		SigmaX[j] = 0.0;
	}
  
  BestNumDescr = 0;
  // Main loop
//		for (i = 0, ti = 0; i < (nMaxDescrInModel * EF); i++)
		for (i = 0, ti = 0; i - BestNumIter < MaxDeltaDescr; i++)
  {
    // Compute SigmaY
		SumY = SumY2 = 0.0f;
    for (k = 0; k < nComp; k++)
      if (bPointUsedForTraining[k])
			{
	      SumY += y[k];
		    SumY2 += y[k] * y[k];
			}
	  AveY = SumY / nCompUsedForTraining;
	  AveY2 = SumY2 / nCompUsedForTraining;
	  if (AveY2 - AveY * AveY > 0.000001)
	    SigmaY = (double)sqrt(AveY2 - AveY * AveY);
	  else
	    SigmaY = 0.0;

	  // Find best correlation with y[]
		best_descr = -1;
	  best_R2 = -1000.0;
	  best_b1 = best_b0 = 0;
		// Loop over candidates
    for (j = 0; j < nDescr; j++)
			if (bDescrSelected[j] && (!bDescriptorReuse)?(!bUsedInModel[j]):true)
	  {
	    // Evaluate simple linear regression
			SumXY = 0.0;
		  for (k = 0; k < nComp; k++)
		    if (bPointUsedForTraining[k])
		      SumXY += X[k][j] * y[k];
		  AveXY = SumXY / nCompUsedForTraining;
		  if (fabs(SigmaX[j]) > 0.00001 && fabs(SigmaY) > 0.00001)
			{
		    curr_corr_coef = (AveXY - AveX[j] * AveY) / (SigmaX[j] * SigmaY);
		    curr_R2 = curr_corr_coef * curr_corr_coef;
		    curr_b1 = (curr_corr_coef * SigmaY) / SigmaX[j];
		    curr_b0 = AveY - curr_b1 * AveX[j];
			}
		  else
		    curr_corr_coef = curr_R2 = curr_b1 = curr_b0 = 0.0;
			// Is this regression the best one?
		  if (curr_R2 > best_R2)
			{
		    best_descr = j;
		    best_R2 = curr_R2;
		    best_b1 = curr_b1 * shrinkage;
		    best_b0 = curr_b0 * shrinkage;
		    best_SigmaY = SigmaY;
			}
		} // j-loop

//		if (best_R2 < 0.01) continue;				
		if (!bUsedInModel[best_descr]) ti++;
//		if (ti > nMaxDescrInModel) break;
//		if (ti - BestNumDescr > MaxDeltaDescr) break;
		
		bUsedInModel[best_descr] = TRUE;
		nDescrIDInModel[i] = best_descr;
    for (k = 0; k < nComp; k++)
		{
	    double fData = X[k][best_descr];
	    y[k] -= (best_b0 + best_b1 * fData);
	    y_pred[k] += (best_b0 + best_b1 * fData);
		}
    nDescrInModel++;
    b[nDescrInModel] = best_b1;
    b[0] += best_b0;
	b0[nDescrInModel] = b[0];
    CalcCorrCoef();
	if (CorrCoef > 0.9999) break;
	if (1.0 * i / ti > EF) break;
    CalcSigmaTraining();
    CalcSigmaValidation();
    CalcSigmaPrediction();		
    if (i == 0 || SigmaValidation < BestSigmaValidation) {
	    BestNumDescr = ti;
		BestNumIter = i;
	    BestCorrCoef = CorrCoef;
	    BestSigmaTraining = SigmaTraining;
	    BestSigmaValidation = SigmaValidation;
	    BestSigmaPrediction = SigmaPrediction;
	}
    DescriptorAddedCallback(i + 1, ti);

  } // i
  
  SaveModel ();
 
  delete [] AveX;
  delete [] AveX2;
  delete [] SigmaX;

}

void FSMLR::SaveModel ()
{

 ReduceModel(BestNumIter+1);

}

void FSMLR::SetMaxDescrInModel(int maxDescr)
{
  nMaxDescrInModel = maxDescr;
}

void FSMLR::SelectDescriptors(int nSel)
{
	int i;
  for (i = 0; i < nDescr; i++)
    bDescrSelected[i] = FALSE;
  for (i = 0; i < nSel; i++)
    bDescrSelected[nDescrIDInModel[i]] = TRUE;
  nDescrSelected = 0;
  for (i = 0; i < nDescr; i++)
    if (bDescrSelected[i])
	  nDescrSelected++;
}

void FSMLR::SelectDescriptors()
{
  SelectDescriptors(GetBestNumDescr());
}

int FSMLR::GetBestNumDescr() const
{
  return BestNumDescr + 1;
}

void FSMLR::ReduceModel(int numSelDescr)
{
	int i, j;
	
	numUniqueDescr = 0;
	for (i = 0; i < numSelDescr; i++)
	{
		bool found = false;
		for (j = 0; j < numUniqueDescr; j++)
			if (nDescrIDInModel[i] == nDescrIDInReducedModel[j])
			{
				found = true;
				break;
			}
		if (found)
			c[j] += b[i + 1];
		else
		{
			nDescrIDInReducedModel[numUniqueDescr] = nDescrIDInModel[i];
			c[numUniqueDescr] = b[i + 1];
			numUniqueDescr++;
		}
	}
	c0 = b0[numSelDescr];
}
