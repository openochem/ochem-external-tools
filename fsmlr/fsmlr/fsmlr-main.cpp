#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "common.h"
#include "FSMLR1.h"

#define MOLNAMESIZE 64

UINT8 nfolds = 10;
UINT8 start = 0;

UINT8 nprop = 0; 
char *asnn_file_name = 0;
char *model_file_name = 0;
char *graphics_file_name = 0;
char **molnames = 0;

double *Y = 0;
double **X = 0;
UINT8 *cSet = 0;

UINT32 nComp = 0;
UINT32 nTrain = 0;
UINT32 nTest = 0;
UINT32 nDescr = 0;
UINT32 nSelDescr = 0;

FSMLR *fsmlr = 0;

void free_memory()
{
  UINT32 i;

  free(Y);
  for (i = 0; i < nComp; i++) free(X[i]);
  free(X);
  free(cSet);
}

int read_asnn_file() {
  const UINT32 BUFFERSIZE = 100000;
  char buffer[BUFFERSIZE];
  char *token;
  char seps[] = " \t";
  int i,i1,i2;
  
  FILE *fd = fopen(asnn_file_name, "r");
  if (!fd) return 1;
  
  /* Count the number of compounds */
  fscanf(fd, "%d %d\n", &i1, &i2);nTrain=i1;nTest=i2;
  nComp = nTrain + nTest;
  printf("DBG> nComp=%d nTrain=%d nTest=%d\n", (int)nComp, (int)nTrain, (int)nTest);
  
  /* Allocate memory for arrays */
  Y = (double *)malloc(nComp * sizeof(double));
  X = (double **)malloc(nComp * sizeof(double *));
  cSet = (UINT8 *)malloc(nComp * sizeof(UINT8));
  for (i = 0; i < nTrain; i++) cSet[i] = cTraining;
  for (i = nTrain; i < nComp; i++) cSet[i] = cPrediction;
  for (i = start; i < nTrain; i += nfolds) cSet[i] = cValidation;
  molnames = (char **)malloc(nComp * sizeof(char *));
  for (i = 0; i < nComp; i++) {
    molnames[i] = (char *)malloc(MOLNAMESIZE * sizeof(char));
	strcpy(molnames[i], "");
  }
  
  long pos = ftell(fd);

  /* Read the first line with data and find the number of descriptors */
  fgets(buffer, BUFFERSIZE, fd);
  token = strtok(buffer, seps);
  token = strtok(NULL, seps);
  while (token != NULL) {
    nDescr++;
    token = strtok(NULL, seps);
  }
  nDescr--;
  printf("DBG> nDescr=%d\n", (int)nDescr);
  
  fseek(fd, 0, SEEK_SET);
  fscanf(fd, "%d %d\n", &i1, &i2);nTrain=i1;nTest=i2;
  
  /* Read and parse data lines */
  for (int i = 0; i < nComp; i++) {
    char s[100];
	fscanf(fd, "%s", s);
	strcpy(molnames[i], s);
    X[i] = (double *)malloc(nDescr * sizeof(double));
    for (int id=0; id < nDescr; id++) {
	  fscanf(fd, "%s", s);
      X[i][id] = atof(s);	  
	}
	fscanf(fd, "%s\n", s);
	Y[i] = atof(s);
  }
  
  fclose(fd);
  
  return 0;
}

int write_model_file() 
{
  int i;

  FILE *fd = fopen(model_file_name, "w");
  if (!fd) return 1;
  
  fprintf(fd, "<MODEL>\n");
  fprintf(fd, "  <VARIABLES>\n");
  for (i = 0; i < fsmlr->BestNumDescr; i++) {
    fprintf(fd, "    <VARIABLE NO=\"%d\">%d</VARIABLE>\n", fsmlr->nDescrIDInReducedModel[i], i+1);
  }
  fprintf(fd, "  </VARIABLES>\n");
  fprintf(fd, "  <TERMS>\n");
  fprintf(fd, "    <TERM VALUE=\"%g\">0</TERM>\n", fsmlr->c0);
  for (i = 0; i < fsmlr->BestNumDescr; i++) {
    fprintf(fd, "    <TERM VALUE=\"%g\">%d</TERM>\n", fsmlr->c[i], i+1);
  }
  fprintf(fd, "  </TERMS>\n");
  fprintf(fd, "</MODEL>\n");
  
  fclose(fd);
  return 0;
}

int write_exp_pred_file(const char *file_name)
{
  int i, j;

  FILE *fd = fopen(file_name, "w");
  if (!fd) return 1;

  for (i = 0; i < nComp; i++) {
    double ypred = fsmlr->c0;
    for (j = 0; j < fsmlr->BestNumDescr; j++) {
      ypred += (X[i][fsmlr->nDescrIDInReducedModel[j]] * fsmlr->c[j]);
    }
    fprintf(fd, "%s\t%g\t%g\n", molnames[i], Y[i], ypred);
  }  
  
  fclose(fd);
  return 0;
}

int main(int argc, char *argv[])
{
  int i;
  int mdd = 20;
  double def = 20.0;
  double dnf = 1.0;
  double shrinkage = 1.0;

  printf("================================================\n");
  printf("Fast Stagewise Multivariate Linear Regression (FSMLR)\n");
  printf("Special edition for OCHEM\n");
  printf("Version 1.5\n");
  printf("================================================\n");
  
  if (argc < 10) {
	  printf("Usage: fsmlr ");
	  printf("-in <input_file_in_asnn_format>");
	  printf("-model <model_file_in_xml_format>");
	  printf("[-mdd <max_delta_descr>] ");
	  printf("[-nfolds <nfolds>] ");
	  printf("[-start <start>] ");
//	  printf("[-def <descr_extra_factor>] ");
//	  printf("[-dnf <descr_num_factor>] ");
	  printf("[-shr <shrinkage>] ");
	  printf("[-graph <graphics_file_name>] ");
	  printf("\n");
  }
  
  for (i = 1; i < argc; i++) {
	if (!strcmp("-mdd", argv[i]))
	  mdd = atoi(argv[++i]);
	else if (!strcmp("-in", argv[i]))
	  asnn_file_name = argv[++i];
	else if (!strcmp("-model", argv[i]))
	  model_file_name = argv[++i];
	else if (!strcmp("-nfolds", argv[i]))
	  nfolds = atoi(argv[++i]);
	else if (!strcmp("-start", argv[i]))
	  start = atoi(argv[++i]);
//	else if (!strcmp("-def", argv[i]))
//	  def = atof(argv[++i]);
//	else if (!strcmp("-dnf", argv[i]))
//	  dnf = atof(argv[++i]);
	else if (!strcmp("-shr", argv[i]))
	  shrinkage = atof(argv[++i]);
	else if (!strcmp("-graph", argv[i]))
	  graphics_file_name = argv[++i];
  }
  
  if (!asnn_file_name) {
    fprintf(stderr, "Unknown input file name\n");
	return 1;
  }	

  printf("Reading input asnn-file %s...\n", asnn_file_name);
  if (read_asnn_file()) return 1;

  fsmlr = new FSMLR(nComp, nDescr, Y, X, cSet, mdd, def, dnf, shrinkage);

  fsmlr->RunMLR();

  nSelDescr = fsmlr->BestNumDescr;
  fsmlr->ReduceModel(nSelDescr);
  fsmlr->SelectDescriptors(nSelDescr);

  printf("================= REGRESSION EQUATION =================\n");
  printf("\t %g \t Intercept\n", fsmlr->c0);
  if (fsmlr->BestNumDescr)
    for (i = 0; i < fsmlr->BestNumDescr; i++)
      printf("%6d \t %g \t descr_%d\n", i + 1, fsmlr->c[i], fsmlr->nDescrIDInReducedModel[i]);
	
  printf("====================== STATISTICS =====================\n");
  printf("Number of iterations = %d\n", fsmlr->BestNumIter+1);
  printf("Number of descriptors = %d\n", fsmlr->BestNumDescr);
  printf("Correlation coefficient = %.4f\n", fsmlr->BestCorrCoef);
  printf("RMS error on the training set = %.4f\n", fsmlr->BestSigmaTraining);
  if (fsmlr->nCompUsedForValidating)
    printf("RMS error on the internal validation set = %.4f\n", fsmlr->BestSigmaValidation);
  if (fsmlr->nCompUsedForPrediction)
    printf("RMS error on the external prediction set = %.4f\n", fsmlr->BestSigmaPrediction);

  write_model_file();
  if (graphics_file_name) write_exp_pred_file(graphics_file_name);
  write_exp_pred_file("ochem");

  delete fsmlr;

  free_memory();

  return 0;
}
