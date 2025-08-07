* Encoding: UTF-8.
* Parallel Analysis program.
*O'Connor, B. P. (2000). SPSS and SAS programs for determining the number of components using parallel analysis and Velicer's MAP test. Behavior Research Methods, Instrumentation, and Computers, 32, 396-402.
*http://flash.lakeheadu.ca/~boconno2/nfactors.html

set mxloops=9000 printback=off width=80  seed = 1953125.
matrix.

* enter your specifications here, probably just ncases and nvars.
compute ncases   = 265. 
compute nvars    = 16.
compute ndatsets = 1000.
compute percent  = 95.

* Specify the desired kind of parallel analysis, where:
  1 = principal components analysis
  2 = principal axis/common factor analysis.
compute kind = 1 .


* principal components analysis.
do if (kind = 1).
compute evals = make(nvars,ndatsets,-9999).
compute nm1 = 1 / (ncases-1).
loop #nds = 1 to ndatsets.
compute x = sqrt(2 * (ln(uniform(ncases,nvars)) * -1) ) &*
            cos(6.283185 * uniform(ncases,nvars) ).
compute vcv = nm1 * (sscp(x) - ((t(csum(x))*csum(x))/ncases)).
compute d = inv(mdiag(sqrt(diag(vcv)))).
compute evals(:,#nds) = eval(d * vcv * d).
end loop.
end if.

* principal axis / common factor analysis with SMCs on the diagonal.
do if (kind = 2).
compute evals = make(nvars,ndatsets,-9999).
compute nm1 = 1 / (ncases-1).
loop #nds = 1 to ndatsets.
compute x = sqrt(2 * (ln(uniform(ncases,nvars)) * -1) ) &*
            cos(6.283185 * uniform(ncases,nvars) ).
compute vcv = nm1 * (sscp(x) - ((t(csum(x))*csum(x))/ncases)).
compute d = inv(mdiag(sqrt(diag(vcv)))).
compute r = d * vcv * d.
compute smc = 1 - (1 &/ diag(inv(r)) ).
call setdiag(r,smc).
compute evals(:,#nds) = eval(r).
end loop.
end if.

* identifying the eigenvalues corresponding to the desired percentile.
compute num = rnd((percent*ndatsets)/100).
compute results = { t(1:nvars), t(1:nvars), t(1:nvars) }.
loop #root = 1 to nvars.
compute ranks = rnkorder(evals(#root,:)).
loop #col = 1 to ndatsets.
do if (ranks(1,#col) = num).
compute results(#root,3) = evals(#root,#col).
break.
end if.
end loop.
end loop.
compute results(:,2) = rsum(evals) / ndatsets.

print /title="PARALLEL ANALYSIS:".
do if   (kind = 1).
print /title="Principal Components".
else if (kind = 2).
print /title="Principal Axis / Common Factor Analysis".
print /title="Compare the random data eigenvalues below to the".
print /title="real-data eigenvalues that are obtained from a".
print /title="Common Factor Analysis in which the # of factors".
print /title="extracted equals the # of variables/items, and the".
print /title="number of iterations is fixed at zero;".
print /title="To obtain these real-data values using SPSS, see the".
print /title="sample commands at the end of the parallel-sps program.".
end if.
compute specifs = {ncases; nvars; ndatsets; percent}.
print specifs /title="Specifications for this Run:"
 /rlabels="Ncases" "Nvars" "Ndatsets" "Percent".
print results /title="Random Data Eigenvalues"
 /clabels="Root" "Means" "Prcntyle"  /format "f12.6".

end matrix.




* Commands for obtaining the necessary real-data eigenvalues for
  principal axis / common factor analysis using SPSS;
  make sure to insert valid filenames/locations, and
  remove the '*' from the first columns.
* corr var1 to var20 / matrix out ('filename') / missing = listwise.
* matrix.
* MGET /type= corr /file='filename' .
* compute smc = 1 - (1 &/ diag(inv(cr)) ).
* call setdiag(cr,smc).
* compute evals = eval(cr).
* print { t(1:nrow(cr)) , evals }
 /title="Raw Data Eigenvalues"
 /clabels="Root" "Eigen."  /format "f12.6".
* end matrix.
*PCA to get the number of ffactors.
DATASET ACTIVATE DataSet1.
FACTOR
  /VARIABLES pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /MISSING LISTWISE 
  /ANALYSIS pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /PRINT INITIAL EXTRACTION
  /PLOT EIGEN
  /CRITERIA MINEIGEN(1) ITERATE(25)
  /EXTRACTION PC
  /ROTATION NOROTATE
  /METHOD=CORRELATION.
*Factor analysis_ML_2 factors_promax rotation.
FACTOR
  /VARIABLES pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /MISSING LISTWISE 
  /ANALYSIS pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /PRINT INITIAL EXTRACTION ROTATION
  /PLOT EIGEN
  /CRITERIA FACTORS(2) ITERATE(25)
  /EXTRACTION ML
  /CRITERIA ITERATE(25)
  /ROTATION PROMAX(4).
*3 factors_promax.
FACTOR
  /VARIABLES pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /MISSING LISTWISE 
  /ANALYSIS pers1 pers2 pers3 pers4 pers5 pers6 pers7 pers8 eff9 eff10 eff11 eff12 eff13 eff14 
    eff15 eff16
  /PRINT INITIAL EXTRACTION ROTATION
  /PLOT EIGEN
  /CRITERIA FACTORS(5) ITERATE(25)
  /EXTRACTION ML
  /CRITERIA ITERATE(25)
  /ROTATION PROMAX(4).
