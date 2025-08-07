* Encoding: UTF-8.
* Yannick Apedo
* Module 1/mProject1

*Rename October dataset variable to merge with February dataset

Rename variables PFT_Situps = Oct_PFT_Situps.
Rename variables PFT_Pushups = Oct_PFT_Pushups.
Rename variables PFT_Runtime = Oct_PFT_Runtime.
Rename variables PFT_Date = Oct_PFT_Date.
Rename variables PSI_Score = Oct_PSI_Score. 
SORT CASES by ID(A).

DATASET ACTIVATE DataSet3.
SORT CASES BY ID.
DATASET ACTIVATE DataSet2.
SORT CASES BY ID.
DATASET ACTIVATE DataSet3.
MATCH FILES /FILE=*
  /FILE='DataSet2'
  /BY ID.
EXECUTE.

SAVE OUTFILE='C:\Users\yannicks\Desktop\Summer 2025\CHIP690-294 (SS2)\PFT_Merged feb and oct.sav'
  /COMPRESSED.

DATASET ACTIVATE DataSet2.

ADD FILES /FILE=*
  /FILE='DataSet7'
  /RENAME (Age Gender PFT_Date PFT_Pushups PFT_Runtime PFT_Situps PSI_Score Race Rank Tenure 
    Unit=d0 d1 d2 d3 d4 d5 d6 d7 d8 d9 d10)
  /DROP=d0 d1 d2 d3 d4 d5 d6 d7 d8 d9 d10.
EXECUTE.

DATASET

DATASET ACTIVATE DataSet7.
COMPUTE PFT_Total_Score=PFT_Situps + PFT_Pushups + PFT_Runtime.
VARIABLE LABELS PFT_Total_Score 'Compute PFT_Total_Score=PFT_Situps + PFT_Pushups+ PFT_Runtime'.
EXECUTE.

COMPUTE AveScore = mean(PFT_Situps, PFT_Pushups, PFT_Runtime).

RECODE PSI_Score (1=10)(2=9)(3=8)(4=7)(5=6)(6=5)(7=4)(8=3)(9=2)(10=1) into PSI_Score_recoded.
EXECUTE.

COMPUTE PSI_Score_recoded = 11 - PSI_Score.
EXECUTE.

RECODE PSI_Score (1=10)(2=9)(3=8)(4=7)(5=6)(6=5)(7=4)(8=3)(9=2)(10=1) INTO
    PSI_Score_recoded.
VARIABLE LABELS PSI_Score_recoded 'recodede PSI scores'. 
EXECUTE.


CROSSTABS
  /TABLES=PSI_Score BY PSI_Score_recoded
  /FORMAT=AVALUE TABLES
  /CELLS=COUNT
  /COUNT ROUND CELL.


VARIABLE LABELS PSI_Score 'non recoded PSI score that was collected in feb'.

USE ALL.
COMPUTE filter_$=(Age > 20 AND Tenure > 30).
VARIABLE LABELS filter_$ 'Age > 20 AND Tenure > 30 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.
_SLINE OFF.

USE ALL.
