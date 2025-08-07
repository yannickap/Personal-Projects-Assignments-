* Encoding: UTF-8.
                           
/* PROCESS version 5.0 */.        
/* Written by Andrew F Hayes */.
/* www.afhayes.com */.
/* www.processmacro.org */.
/* Copyright 2013-2025 by Andrew F Hayes */.
/* Documented in http://www.guilford.com/p/hayes3 and supplements */.
/* This version of PROCESS works only on SPSS 26 or later */.

/* THIS CODE SHOULD BE DISTRIBUTED ONLY THROUGH PROCESSMACRO.ORG */.
/* AND NOT POSTED OR STORED ON ANY WEB PAGE OR SERVER WITHOUT PERMISSION */.
/* THIS SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND */.
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF */.
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT */.
/* IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, */.
/* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT */.
/* OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE */.
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE */.
/* USE OF THIS SOFTWARE IMPLIES AGREEMENT WITH THESE TERMS */.

set printback=off.
define !bcboot3 (databcbt = !charend ('/')/estmte = !charend ('/') !default(9999)).
compute temp = !databcbt.
compute temp(GRADE(!databcbt)) = !databcbt.
compute badlo = 0.
compute badhi = 0.
do if (!estmte <> 9999).
  compute pv=csum(temp < !estmte)/nrow(temp).
  compute ppv = pv.
  do if (pv > .5).
    compute ppv = 1-pv.
  end if.
  /* compute y5=sqrt(-2*ln(ppv)) */.
  /* compute xp=y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0) */.
  compute xp=-idf.normal(ppv, 0,1).
  do if (pv <= .5).
    compute xp = -xp.
  end if.
  compute cilow=rnd(nrow(temp)*(cdfnorm(2*xp-xp2))).
  compute cihigh=trunc(nrow(temp)*(cdfnorm(2*xp+(xp2))))+1.
  do if (cilow < 1).
    compute cilow = 1.
     compute booterr=1.
     compute badlo = 1.
  end if.
  do if (cihigh > nrow(temp)).
    compute cihigh = boot.
    compute booterr=1.
     compute badhi = 1.
  end if.
  compute llcit=temp(cilow,1).
  compute ulcit=temp(cihigh,1).
  do if (badlo = 1 and llcit <> priorlo).
    compute badend={badend, llcit}.
    compute priorlo = llcit.
  end if.
  do if (badhi = 1 and ulcit <> priorhi).
    compute badend={badend, ulcit}.
    compute priorhi = ulcit.
  end if.
end if.
do if (!estmte = 9999).
   compute llcit=temp(cilow,1).
   compute ulcit=temp(cihigh,1).
end if.
compute bootse=(cssq(temp)-((csum(temp)&**2)/nrow(temp)))/(nrow(temp)-1).
do if (bootse < 0).
  compute bootse=0.
end if.
compute bootse=t(sqrt(bootse)).
!enddefine.

define !llrtest (lm=!charend('/')).
compute btemphld=b.
compute llrdat=make(nrow(x),nrow(!lm)-csum(!lm),-999).
compute llrdf=ncol(x)-ncol(llrdat).
compute llrcnt=0.
loop llri=1 to nrow(!lm).
  do if (!lm(llri,1)=0).
    compute llrcnt=llrcnt+1.
    compute llrdat(:,llrcnt)=x(:,llri).
  end if.
end loop.
!modeles y=y/x=llrdat/type=2/full=0.
compute b=btemphld.
compute fresult={(LL2-basemod),llrdf,1-chicdf((LL2-basemod),llrdf)}.
!enddefine.

define !outform (outtodo=!charend('/')/outbig=!charend('/')).        
  compute resultm2=make(nrow(!outtodo),!outbig,99999).
  do if (ncol(!outtodo) <= !outbig).        
    compute resultm2(1:nrow(!outtodo),1:ncol(!outtodo))=!outtodo.
    compute resultm={resultm;resultm2}.
  end if.
  do if (ncol(!outtodo) > !outbig).        
    compute resultmt=make(nrow(resultm),ncol(!outtodo),99999).
    compute resultmt(1:nrow(resultm),1:ncol(resultm))=resultm.
    compute resultm=resultmt.
    compute resultm2=make(nrow(!outtodo),ncol(resultm),99999).
    compute resultm2(1:nrow(!outtodo),1:ncol(!outtodo))=!outtodo.
    compute resultm={resultm;resultm2}.
    compute maxresm=ncol(resultm).
  end if.        
!enddefine.

DEFINE !describ (descdatf=!charend('/')/type=!charend('/') !default(0)).
compute desctmp=make(8,ncol(!descdatf),-999).
/* compute desctmp=make((8-(4*!type)),ncol(!descdatf),-999) */.
/* mean, sd, min, max, 16th, 50th, 84th, dich toggle */.
compute modvals=-999.
loop jd=1 to ncol(!descdatf).
  compute descdat=!descdatf(:,jd).
  /* get the mean, sd, minimum, and maximum */.
  compute desctmp(1,jd) = csum(descdat)/nrow(descdat).
  compute desctmp(2,jd) = (nrow(descdat)*sscp(descdat))-(t(csum(descdat))*(csum(descdat))).
  compute desctmp(2,jd) = sqrt(desctmp(2,jd)/(nrow(descdat)*(nrow(descdat)-1))).
  compute desctmp(3,jd)=cmin(descdat).
  compute desctmp(4,jd)=cmax(descdat).
  /* check to see if the variable is a constant */.
  do if ((desctmp(3,jd)=desctmp(4,jd)) and novar=0).
    compute errcode(errs,1)=15.
    compute errs=errs+1.
    compute criterr=1.
    compute novar=1.
  end if.
  /* check and set the dichotomous toggle */.
  compute tmp=((descdat(:,1)=desctmp(3,jd))+(descdat(:,1)=desctmp(4,jd))).
  compute desctmp(8,jd)=(csum(tmp)=nrow(tmp)).
  do if (desctmp(8,jd)=1).
    compute modvals={desctmp(3,jd);desctmp(4,jd)}.
    compute mnotev=0.
    compute minwarn=0.
    compute maxwarn=0.
  end if.
  do if (!type=0).
     compute minwarn=0.
     compute maxwarn=0.
    /* here we calculate the percentiles of the distribution */.
    do if (quantile=1 and desctmp(8,jd) <> 1).
      compute tmp = descdat.
      compute tmp(GRADE(descdat),:) = descdat.
      compute descdat = tmp.
      release tmp.
      compute decval={.16;.5;.84}.
      loop kd=1 to 3.
        compute low=trunc(decval(kd,1)*(nrow(descdat)+1)).
        compute lowdec=decval(kd,1)*(nrow(descdat)+1)-low.
        compute value=descdat(low,1)+(descdat((low+1),1)-descdat(low,1))*lowdec.
        compute desctmp((4+kd),jd)=value.
      end loop.
      compute mnotev=1.
      compute modvals=desctmp(5:7,:).
    end if.
    do if (quantile <> 1 and desctmp(8,jd) <> 1).
      compute desctmp(5,jd)=desctmp(1,jd)-desctmp(2,jd).
      compute desctmp(6,jd)=desctmp(1,jd).
      compute desctmp(7,jd)=desctmp(1,jd)+desctmp(2,jd).
      compute modvals=desctmp(5:7,:).
      compute mnotev=2.
      do if (modvals(1,jd) < desctmp(3,jd)).
        compute modvals(1,jd)=desctmp(3,jd).
        compute minwarn=1.
      end if.
      do if (modvals(3,jd) > desctmp(4,jd)).
        compute modvals(3,jd)=desctmp(4,jd).
        compute maxwarn=1.
      end if.
    end if.
  end if.
end loop.
!ENDDEFINE.

DEFINE !ftest3 (lm=!charend('/')/bcoef=!charend('/')/cv=!charend('/') !default(0)/chr=!charend('/') !default(0)/
  brsq=!charend('/') !default(0)/skip=!charend('/') !default(0)).
compute lmat2=!lm.
do if (!skip=0).
  compute lmat2 = mdiag(!lm).
  compute lmat3=make(nrow(lmat2),1,0).
  loop flp=1 to ncol(lmat2).
    do if (csum(lmat2(:,flp))=1).
      compute lmat3={lmat3,lmat2(:,flp)}.
    end if.
  end loop.
  compute lmat2=lmat3(:,2:ncol(lmat3)).
end if.
compute fresult={99999,99999,99999,99999}.
do if (ncol(lmat2) > rank(!cv)).
  compute intncd=1.
end if.
do if (ncol(lmat2) <= rank(!cv)).
  compute fratio = (t(t(lmat2)*!bcoef)*inv(t(lmat2)*!cv*lmat2)*((t(lmat2)*!bcoef)))/ncol(lmat2).
  compute pfr = 1-fcdf(fratio,ncol(lmat2),dfres).
  compute fresult={fratio,ncol(lmat2),dfres,pfr}.
end if.
do if (i = (nms+nys) and (ydich=1)).
  compute fratio = (t(t(lmat2)*!bcoef)*inv(t(lmat2)*!cv*lmat2)*((t(lmat2)*!bcoef)))/ncol(lmat2).
  compute fratio=fratio*ncol(lmat2).
  compute pfr=1-chicdf(fratio,ncol(lmat2)).
  compute fresult={fratio,ncol(lmat2),pfr}.
end if.
do if (!chr=1).
  compute lmat3=1-rsum(lmat2).
  compute xfm=make(n,csum(lmat3),0).
  compute flpc=1.
  loop flp=1 to nrow(lmat3).
    do if (lmat3(flp,1)=1).
      compute xfm(:,flpc)=x(:,flp).
      compute flpc=flpc+1.
    end if.
  end loop.
  compute bfm=inv(t(xfm)*xfm)*t(xfm)*y.
  compute residy=y-(xfm*bfm).
  compute sstotal=(y-(csum(y)/n)).
  compute sstotal=csum(sstotal&*sstotal).
  compute ssresid=csum(residy&*residy).
  compute rsqch=!brsq-((sstotal-ssresid)/sstotal).
    compute fresult={rsqch,fresult}.
  release xfm,flpc, residy, ssresid, bfm.
end if.
!ENDDEFINE.

DEFINE !hcest3 (x=!charend('/')/resid=!charend('/')/hc=!charend('/')/mse=!charend('/')).
  compute n1=nrow(!x).
  compute invXtX = inv(t(!x)*!x).
  compute varb = !mse*invXtX.
  compute k3 = ncol(!x).
  compute xhc=0.
  /* cluster robust varb */.
  do if (!hc=6).
    compute xcrub={clusdat(:,3),!x,!resid}.
    compute temp = xcrub.
    compute temp(GRADE(xcrub(:,1)),:) = xcrub.
    compute xcrub = temp.
    compute crres=xcrub(:,ncol(xcrub)).
    compute xcrub=xcrub(:,2:(ncol(xcrub)-1)).
    release temp.
    compute meat=make(ncol(xcrub),ncol(xcrub),0).
    loop clbt=1 to nclus.
        compute xmeat=xcrub(current(clbt,3):current(clbt,4),:).
        compute meate=crres(current(clbt,3):current(clbt,4),1).
        compute meat=meat+(t(xmeat)*meate*t(meate)*xmeat).
     end loop.
     compute varb=(invXtX*meat*invXtX)*((nclus/(nclus-1))*((n1-1))/(n1-ncol(xcrub))).
     release xcrub.
     release crres.
  end if.
  do if (!hc < 5).
    compute xhc=!x.
    compute hat = xhc(:,1).
    loop i3=1 to nrow(xhc).
      compute hat(i3,1)= xhc(i3,:)*invXtX*t(xhc(i3,:)). 
    end loop.
    do if (!hc = 0 or !hc =1).
      loop i3 = 1 to k3.
        compute xhc(:,i3)=xhc(:,i3)&*!resid.
      end  loop.
    end if.
    do if (!hc =3 or !hc=2).
      loop i3=1 to k3.
        compute xhc(:,i3) = (!resid&/(1-hat)&**(1/(4-!hc)))&*xhc(:,i3).
      end loop.
    end if.
    do if (!hc = 4).
      compute hcmn=make(n1,2,4).
      compute hcmn(:,2)=(n1*hat)/k3.
      loop i3= 1 to k3.
        compute xhc(:,i3) = (!resid&/(1-hat)&**(rmin(hcmn)/2))&*xhc(:,i3).
      end loop.
    end if.
    compute varb=(invXtX*t(xhc)*xhc*invXtX).
    do if (!hc=1).
      compute varb=(n1/(n1-ncol(!x)))&*varb.
    end if.
  end if.
  compute hclab={"se(HC0)","se(HC1)","se(HC2)","se(HC3)","se(HC4)","se","se(CR)"}.
  compute hclab=hclab(1,(!hc+1)).
  compute hcflab={"F(HC0)","F(HC1)","F(HC2)","F(HC3)","F(HC4)","F", "F(CR)"}.
  compute hcflab=hcflab(1,(!hc+1)).
  release xhc.
!ENDDEFINE.

DEFINE !longchk (variab=!charend('/')).
/* check for any long variable names */.
  compute save4sv={"CASENUM","DV","PRED","RESID", "D_RESID","STRESID","TRESID","H","MAHAL","COOK","DMSRES","DRSQ","DSKEW",t(dfbetasC),t(dfindsC)}.
  !let !toomany=0.
  !do !i !in (!variab).
    !do !j = 1 !to !length(!i).
      !if ((!j > 8) !and (!toomany = 0)) !then.
        !let !toomany = 1.
        compute toomany=1.
      !ifend.
    !doend.
      do if (savediag=1).
        !let !itsh=!upcase(!quote(!i)).
        loop sv4chk=1 to ncol(save4sv).
          do if (!itsh=save4sv(1,sv4chk)).
            do if (sv4match=0).
              compute sv4match=1.    
              compute notecode(notes,1) = 39.
              compute notes = notes + 1. 
            end if.
          end if.
        end loop.
      end if.
  !doend.
   do if (toomany=1 and longname=0).
    compute criterr=1.
    do if (tooman=0).
      compute tooman=1.     
      compute errcode(errs,1) = 61.
      compute errs = errs + 1. 
    end if.
  end if.
!ENDDEFINE.

DEFINE !drinef3 (decpnt=!charend('/') !default(F10.4)).
do if (criterr = 0 and nms > 0).
  compute paths=paths(:,2:ncol(paths)).
  compute pathsw=pathsw(:,2:ncol(pathsw)).
  compute pathsz=pathsz(:,2:ncol(pathsz)).
  compute pathswz=pathswz(:,2:ncol(pathswz)).
  compute pathsmod=pathsw+pathsz+pathswz.
  compute pathsdv=pathsdv(:,2:ncol(pathsdv)).
  compute pathsfoc=pathsfoc(:,2:ncol(pathsfoc)).
  compute pathtype=pathtype(:,2:ncol(pathtype)).
  compute anymod=(rsum(pathsmod) > 0).
  compute obscoeff=obscoeff(1,2:ncol(obscoeff)).
  do if (outscree=1).
    do if (xmint=1).
      print/title = "************************* COUNTERFACTUALLY DEFINED ************************".      
      else.
      print/title = " ".  
    end if.
    do if (dototal = 0 and alttotal=0).
      print/title = "****************** DIRECT AND INDIRECT EFFECTS OF X ON Y ******************"/space=0.
    end if.
    do if (alttotal=1).
      print/title = "************** TOTAL, DIRECT, AND INDIRECT EFFECTS OF X ON Y **************"/space=0.
    end if.
  end if.
  do if (dototal = 1).
    do if (outscree=1).
      print/title = "************** TOTAL, DIRECT, AND INDIRECT EFFECTS OF X ON Y **************"/space=0.   
    end if.
    compute totefflb(1,1)="Effect".
    /* do if (xmint <>  1) */.
      do if (effsize=1).
        compute toteffsz=toteff(:,1)/ysd.
        do if (xdich=1 or mcx > 0).
          compute totefflb={totefflb,"c_ps"}.
        end if.
        do if (xdich = 0 and mcx = 0).
          compute toteffsz=toteffsz*xsd.
          compute totefflb={totefflb,"c_cs"}.
        end if.
        compute toteff={toteff,toteffsz}.
      end if.
      !outform outtodo=toteff/outbig=maxresm.
      do if (nxvls > 1).
        !outform outtodo=totomni/outbig=maxresm.
        compute clabtmp={"R2-chng", hcflab, "df1","df2","p"}.
        do if (outscree=1).
          print toteff/title="Relative total effects of X on Y"/rnames=toteffl2/cnames=totefflb/format=!decpnt.
          do if (eivdo=1).
             compute clabtmp={hcflab, "df1","df2","p"}.
          end if. 
          print totomni/title="Omnibus test of total effect of X on Y"/cnames=clabtmp/format=!decpnt.
          print/title= "----------"/space=0.
        end if.
      else.
        do if (outscree=1).
          print toteff/title="Total effect of X on Y"/cnames=totefflb/format=!decpnt.
        end if.
      end if.
    /* end if */.
  end if.
  compute moddir=wcmat(nrow(bcmat),1)+zcmat(nrow(bcmat),1).
  do if (xmint=1).
    compute moddir=1.
  end if.
  /* do if (bcmat(nrow(bcmat),1)=1 and xmint <> 1) */.
  do if (bcmat(nrow(bcmat),1)=1).
    do if (ydich=1).
      compute direfflb(:,(ncol(direfflb)-5):ncol(direfflb))={"Effect","se","Z","p","LLCI","ULCI"}.
    end if. 
    do if (moddir=0 or xmint=1).
      compute direfflb(1,1)="Effect".
    end if.
    do if (effsize=1 and moddir=0 and anymod = 0).
      compute direffsz=direff(:,1)/ysd.
      do if (xdich=1 or mcx > 0).
        compute direfflb={direfflb,"c'_ps"}.
      end if.
      do if (xdich = 0 and mcx = 0).
        compute direffsz=direffsz*xsd.
        compute direfflb={direfflb,"c'_cs"}.
      end if.
      compute direff={direff,direffsz}.
    end if.
    !outform outtodo=direff/outbig=maxresm.
    do if (xmint=1).
      !outform outtodo=codireff/outbig=maxresm.
    end if.
    do if (moddir = 0 and nxvls=1).
      do if (outscree=1).
        print direff/title="Direct effect of X on Y"/cnames=direfflb/format=!decpnt.
      end if.
    end if.
    do if (moddir = 0 and nxvls>1). 
      !outform outtodo=diromni/outbig=maxresm.
      do if (outscree=1).
        print direff/title="Relative direct effects of X on Y"/rnames=direffl2/cnames=direfflb/format=!decpnt.  
        do if (ydich = 0).  
          compute clabtmp={"R2-chng", hcflab, "df1","df2","p"}.
          do if (eivdo=1).
             compute clabtmp={hcflab, "df1","df2","p"}.
          end if. 
          print diromni/title="Omnibus test of direct effect of X on Y:"/cnames=clabtmp/format=!decpnt.
        end if.
        do if (ydich = 1).  
          compute clabtmp={"Chi-sq", "df", "p"}.
          print diromni/title="Omnibus likelihood ratio test of direct effect of X on Y"/cnames=clabtmp/format=!decpnt.
        end if.
        print/title= "----------"/space=0.
      end if.
    end if.
    do if (moddir > 0 and nxvls=1).
      do if (xmint=0 and outscree=1).
        print direff/title="Conditional direct effects of X on Y"/cnames=direfflb/format=!decpnt.
      end if.
      do if (xmint=1).
        do if (outscree=1).
          print direff/title="(Pure) Natural direct effect of X on Y"/cnames=direfflb/format=!decpnt.
          print codireff/title="Controlled direct effect of X on Y"/cnames=direfflb/format=!decpnt.
          print/title= "----------"/space=0.
        end if.
        compute obnatdfx=direff(:,1).
      end if.         
    end if.
    compute direffl4=direffl2.
    do if (moddir > 0 and nxvls>1).
      compute direffl2={" "}.
      loop i = 1 to nxvls.
        loop j = 1 to (nrow(direff)/nxvls).
          compute direffl2={direffl2;xcatlab(i,1)}.
        end loop.
      end loop.
      compute direffl2=direffl2(2:nrow(direffl2),1).
      do if (xmint=0 and outscree=1).
        print direff/title="Relative conditional direct effects of X on Y"/cnames=direfflb/rnames=direffl2/format=!decpnt.
      end if.
      do if (xmint=1).
        do if (outscree=1).
          print direff/title="Relative (pure) natural direct effects of X on Y"/cnames=direfflb/rnames=direffl2/format=!decpnt.
          print codireff/title="Relative controlled direct effects of X on Y"/cnames=direfflb/rnames=direffl2/format=!decpnt.
          print/title= "----------"/space=0.
        end if.
        compute obnatdfx=direff(:,1).
      end if.
    end if.
    compute direffl2=direffl4.    
  end if.
  do if (bcmat(nrow(bcmat),1)=0 and xmint <> 1 and outscree=1).
      print/title="The direct effect of X on Y is fixed to zero.".
  end if.

  /* Here is the start of the indirect effects */.

  do if (nms = 1).
    compute indmark={2}.
    compute indsets={1,2}.
  end if.
  do if (nms = 2).
  compute indmark={2,2,3}.
  compute indsets={1,4,2,5,1,3,5}.
    compute thetam=1.
  end if.
  do if (nms = 3).
    compute indmark={2,2,2,3,3,3,4}.
    compute indsets={1,7,2,8,4,9,1,3,8,1,5,9,2,6,9,1,3,6,9}.
    compute thetam={1,2,3}.
  end if.
  do if (nms = 4).
    compute indmark={2,2,2,2,3,3,3,3,3,3,4,4,4,4,5}.
    compute indsets={1,11,2,12,4,13,7,14,1,3,12,1,5,13,1,8,14,2,6,13,2,9,14,4,10,14,1,3,6,13,1,3,9,14,1,5,10,14,2,6,10,14,1,3,6,10,14}.
    compute thetam={1,2,4,3,5,6}.
  end if.
  do if (nms = 5).
    compute indmark={2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,6}.
    compute indsets={1,16,2,17,4,18,7,19,11,20,1,3,17,1,5,18,1,8,19,1,12,20,2,6,18,2,9,19,2,13,20,4,10,19,4,14,20,7,15,20,1,3,6,18,
                                    1,3,9,19,1,3,13,20,1,5,10,19,1,5,14,20,1,8,15,20,2,6,10,19,2,6,14,20,2,9,15,20,4,10,15,20,1,3,6,10,19,
                                     1,3,6,14,20,1,3,9,15,20,1,5,10,15,20,2,6,10,15,20,1,3,6,10,15,20}.
    compute thetam={1,2,5,3,6,8,4,7,9,10}.
  end if.
  do if (nms = 6).
    compute indmark={2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,7}.
    compute indsets={1,22,2,23,4,24,7,25,11,26,16,27,1,3,23,1,5,24,1,8,25,1,12,26,1,17,27,2,6,24,2,9,25,2,13,26,2,18,27,4,10,25,4,14,26,
                                    4,19,27,7,15,26,7,20,27,11,21,27,1,3,6,24,1,3,9,25,1,3,13,26,1,3,18,27,1,5,10,25,1,5,14,26,1,5,19,27,1,8,15,26,1,8,20,27,
                                    1,12,21,27,2,6,10,25,2,6,14,26,2,6,19,27,2,9,15,26,2,9,20,27,2,13,21,27,4,10,15,26,4,10,20,27,4,14,21,27,7,15,21,27,
                                    1,3,6,10,25,1,3,6,14,26,1,3,6,19,27,1,3,9,15,26,1,3,9,20,27,1,3,13,21,27,1,5,10,15,26,1,5,10,20,27,1,5,14,21,27,
                                    1,8,15,21,27,2,6,10,15,26,2,6,10,20,27,2,6,14,21,27,2,9,15,21,27,4,10,15,21,27,1,3,6,10,15,26,1,3,6,10,20,27,
                                     1,3,6,14,21,27,1,3,9,15,21,27,1,5,10,15,21,27,2,6,10,15,21,27,1,3,6,10,15,21,27}.
    compute thetam={1,2,6,3,7,10,4,8,11,13,5,9,12,14,15}.
  end if.
  do if (nms = 7).
    compute indmark={2,2,2,2,2,2,2}.
    compute indsets={1,29,2,30,4,31,7,32,11,33,16,34,22,35}.
  end if.
  do if (nms = 8).
    compute indmark={2,2,2,2,2,2,2,2}.
    compute indsets={1,37,2,38,4,39,7,40,11,41,16,42,22,43,29,44}.
  end if.
  do if (nms = 9).
    compute indmark={2,2,2,2,2,2,2,2,2}.
    compute indsets={1,46,2,47,4,48,7,49,11,50,16,51,22,52,29,53,37,54}.
  end if.
  do if (nms=10).
    compute indmark={2,2,2,2,2,2,2,2,2,2}.
    compute indsets={1,56,2,57,4,58,7,59,11,60,16,61,22,62,29,63,37,64,46,65}.
  end if.
  compute indmake=make(ncol(indmark),(nms+2),0).
  compute indmod=make(ncol(indmark),1,999).
  compute indmmm=make(ncol(indmark),1,0).
  compute indmmmt=make(ncol(indmark),1,0).
  compute start=1.
  compute end=0.
  compute nindfx=0.
  compute indlocs=make(nrow(thetaxmb),ncol(paths),999).
  compute indkey=make(ncol(indmark),1+((rmax(indmark)*2)+1),"     ").
  compute c1=1.
  compute c2=1.
  compute c3=1.
  loop i = 1 to ncol(paths).
    do if (pathtype(1,i)=1).
      compute indlocs(:,i)=thetaxmb(:,c1).
      compute c1=c1+1.
    end if.
    do if (pathtype(1,i)=3).
      compute indlocs(:,i)=thetamyb(:,c2).
      compute c2=c2+1.
    end if.
    do if (pathtype(1,i)=2 and nms < 7 and serial=1).
      compute indlocs(:,i)=thetammb(:,thetam(1,c3)).
      compute c3=c3+1.
    end if.
  end loop.
  /* print indlocs */.
  
  loop i = 1 to ncol(indlocs).
    compute c1=2.
    loop j = 2 to nrow(indlocs).
      do if (indlocs(j,i) <> 0).
        compute indlocs(c1,i)=indlocs(j,i).
        compute c1=c1+1.
      end if.
    end loop.
    compute indlocs(1,i)=c1-2.
  end loop.
  compute indlocs=indlocs(1:rmax((indlocs(1,:))+1),:).
  loop i = 1 to ncol(indmark).
    compute numget=indmark(1,i).
    compute end=end+numget.
    compute gotcha=indsets(1,start:end).
    compute start=end+1.
    compute ok=1.
    compute temp=0.
    compute repoman=make(4,1,0).
      loop j = 1 to ncol(gotcha).
        do if paths(1,gotcha(1,j))=0.
          compute ok=0.
        end if.
        do if (pathsmod(1,gotcha(1,j)) > 0).
          compute temp=1.
          compute temp2={pathsw(1,gotcha(1,j));pathsz(1,gotcha(1,j));pathswz(1,gotcha(1,j));0}.
          do if (temp2(1,1)=1 and temp2(2,1)=1 and temp2(3,1)=0).
            compute temp2(4,1)=1.
          end if.
          compute repoman=repoman+temp2.
        end if.
      end loop.
      compute temp=0.  
      compute tempmmm=0.
      compute typemmm=0.
      do if ((repoman(1,1) > 0) and (repoman(2,1) = 0)).
        compute temp=1.
        do if (repoman(1,1)=1).
          compute tempmmm=1. 
        end if.
        do if (repoman(1,1) > 1 and (wdich=1 or mcw > 0)).
          compute tempmmm=12. 
          compute typemmm=mcw.
          do if (wdich=1).
            compute typemmm=1.
          end if.
        end if.
        do if (repoman(1,1) > 1 and (wdich=0 and mcw = 0)).
          compute tempmmm=101. 
        end if.
      end if.
      do if ((repoman(1,1) = 0) and (repoman(2,1) > 0)).
        compute temp=2.
        do if (repoman(2,1)=1).
          compute tempmmm=2.
        end if.
        do if (repoman(2,1) > 1 and (zdich = 1 or mcz > 0)).
          compute tempmmm=22.
          compute typemmm=mcz.
          do if (zdich=1).
            compute typemmm=1.
          end if.
        end if.
        do if (repoman(2,1) > 1 and (zdich = 0 and mcw = 0)).
          compute tempmmm=102.
        end if.
      end if.
      do if (repoman(1,1)>0 and repoman(2,1)>0).
        compute temp=3.
        do if (repoman(1,1)=1 and repoman(2,1)=1).
          do if (repoman(4,1)=1).
            compute tempmmm=31.
          end if.
          do if (repoman(3,1)=1).
            compute tempmmm=41.
          end if.
        end if.
      end if.
      do if (repoman(1,1)=1 and repoman(2,1)=1 and repoman(3,1)=0 and repoman(4,1)=0).
        compute tempmmm=51.
      end if.
      do if (ok = 1).
        compute nindfx=nindfx+1.
        compute indmake(nindfx,1)=numget.
        compute indmod(nindfx,1)=temp.
        compute indmmm(nindfx,1)=tempmmm.
        compute indmmmt(nindfx,1)=typemmm.
        compute indmake(nindfx,2:(1+numget))=gotcha.
        compute indkey(nindfx,1)=xnames.
          loop j = 1 to numget.
            compute indkey(nindfx,(j*2+1))=pathsdv(1,gotcha(1,j)).
            compute indkey(nindfx,(j*2))="   ->   ".
          end loop.
      end if.
  end loop.
  compute indkey=indkey(1:nindfx,1:((cmax(indmake(:,1))*2)+1)).
  compute indmake=indmake(1:nindfx,1:(cmax(indmake(:,1))+1)).
  compute indmod=indmod(1:nrow(indmake),1).
  compute indmmm=indmmm(1:nrow(indmake),1).
  compute indmmmt=indmmmt(1:nrow(indmake),1).
  compute ncpairs = (((nindfx)*(nindfx-1))/2).
  
  do if ((contrast = 1 or contrast = 2 or contrast=3) and (ncpairs > 105)).
    compute contrast=0.
    compute notecode(notes,1) = 13.
    compute notes = notes + 1.
  end if.
  do if (contrast = 4).
    do if (ncol(contvec) <> nindfx).
      compute contrast=0.
      compute notecode(notes,1) = 14.
      compute notes = notes + 1.
    end if.
  end if.

  /* This is for models with no moderator */.
  do if (anymod=0).
    do if (nms = 1 and contrast > 0).
      compute contrast=0.
    end if.
    compute efloop=(((1-(effsize=0))*2)+1)-(((mcx>0 or xdich=1))*(1-(effsize=0))).
    do if (diagnose=1 or savediag=1).
      compute alldfbs=alldfbs(:,2:ncol(alldfbs)).
    end if.
    compute dfindfx=make(n,1,999).
    loop kk=1 to efloop.
      do if (boot = 0).
        compute bootres=obscoeff.
        do if (kk=1).
          compute totbtvec=make(1,nxvls,0).
        end if.
        compute bootdir=obsdirfx.
        compute indtab=999.
        compute inddiff=999.
        compute bootysd=ysd.
        compute bootxsd=xsd.
      end if.
    
      do if (boot > 0).
        compute bootres={obscoeff;bootres}.
        do if (kk=1).
          compute totbtvec=make(nrow(bootres),nxvls,0).
        end if.
        compute indtab=make(1,4,999).
        compute inddiff=make(nrow(bootres),1,999).
      end if.
      compute indtotal=make(nrow(bootres),1,0).
      loop i = 1 to nrow(indmake).
        loop j = 1 to nxvls.
          compute indtemp=make(nrow(bootres),1,1).
          compute dfabtemp=make(n,1,1).
          loop k = 1 to indmake(i,1).
            compute jtemp=1.
            do if (j > 1 and k=1).
              compute jtemp=j.
            end if.
            compute indtemp=indtemp&*bootres(:,pathsfoc(jtemp,indmake(i,(k+1)))).
            do if (kk=1).
              do if (ydich=0 and ((diagnose=1) or (savediag=1))).
                compute dfabtemp=dfabtemp&*(bootres(1,(pathsfoc(jtemp,indmake(i,(k+1)))))-alldfbs(:,pathsfoc(jtemp,indmake(i,(k+1))))).
              end if.
            end if.
          end loop.
          do if (kk=1 and ydich=0 and ((diagnose=1) or (savediag=1))).
            compute dfindfx={dfindfx,indtemp(1,1)-dfabtemp}.
          end if.
          do if (kk = 2).
            compute indtemp=indtemp/bootysd.
          end if.
          do if (kk = 3).
            compute indtemp=(bootxsd&*indtemp)/bootysd.
          end if.
          do if (contrast <> 0).
            compute inddiff={inddiff,indtemp}.
          end if.
          do if (nxvls=1).
            compute indtotal=indtotal+indtemp.
          end if.
          compute indeff=indtemp(1,1).
          do if (kk=1).
            compute totbtvec(:,j)=totbtvec(:,j)+indtemp.
          end if.
          do if (boot > 0).
            !bcboot3 databcbt = indtemp(2:nrow(indtemp),1)/estmte=(indtemp(1,1)*bc)+(9999*(1-bc)).
            compute indeff={indeff,bootse,llcit,ulcit}.
          end if.
          do if (kk=1).
            compute indtabn={indtab;indeff}.
          end if.
          compute indtab={indtab;indeff}.
        end loop.
      end loop.
      compute indtab=indtab(2:nrow(indtab),:).
      do if (kk=1).
        compute indtabn=indtabn(2:nrow(indtabn),:).
        do if (ydich=0 and ((diagnose=1) or (savediag=1))).
          compute wadfie=make((ncol(dfindfx)-1),2,0).
          loop wadf1ie=1 to (ncol(dfindfx)-1).
            loop wadf2ie=1 to nrow(dfindfx).
              do if ((abs(dfindfx(wadf2ie,(wadf1ie+1)))) > abs(wadfie(wadf1ie,2))).
                compute wadfie(wadf1ie,2)=dfindfx(wadf2ie,(wadf1ie+1)).
                compute wadfie(wadf1ie,1)=rownum(wadf2ie,1).
              end if.                     
            end loop.
          end loop.
        end if.
      end if.
      compute rowlbs=indlbl(1:nrow(indtab),1).
      do if (mc > 0).
        compute inddiff=make(mc,1,-999).
        compute indtab2=make(nrow(indtab),4,-999).
        compute indtab2(:,1)=indtab.
        compute indtab=indtab2.
        release indtab2.
        compute mcct=0.
        compute indtotal=make(mc,1,0).
        do if (kk = 1).
          compute x1 = sqrt(-2*ln(uniform(mc,nrow(mcsopath))))&*cos((2*3.14159265358979)*uniform(mc,nrow(mcsopath))).
          compute x1=x1*chol(indcov).
          loop ii=1 to nrow(x1).
            compute x1(ii,:)=x1(ii,:)+t(mcsopath).
          end loop.
        end if.
        loop ii=1 to nms.
          compute tmpb=x1(:,((nms*nxvls)+ii)).
          compute tmpb2=tmpb.
          do if (nxvls > 1).
            loop jj=1 to (nxvls-1).
              compute tmpb2={tmpb2,tmpb}.
            end loop.
          end if.
          compute indtemp=x1(:,(((ii-1)*nxvls)+1):(ii*nxvls))&*tmpb2.
          loop jj=1 to ncol(indtemp).
            do if (kk = 2).
              compute indtemp(:,jj)=indtemp(:,jj)/ysd.
            end if.
            do if (kk = 3).
              compute indtemp(:,jj)=(xsd*indtemp(:,jj))/ysd.
            end if.            
            !bcboot3 databcbt = indtemp(:,jj).           
            compute mcct=mcct+1.
            compute indtab(mcct,2:4)={bootse,llcit,ulcit}.
          end loop.
          do if (nxvls=1).
            compute indtotal=indtotal+indtemp.
            do if (contrast <> 0).
              compute inddiff={inddiff,indtemp}.
            end if.
          end if.
        end loop.
        release indtemp,tmpb.
      end if.
      
      do if (normal = 1 and sobelok=1).
        compute sobelmat=indtab(:,1).
        compute sobelmat={sobelmat,(sobelmat/2),sobelmat,sobelmat}.
        loop ii=1 to nms.
          compute se2b=indcov(((nms*nxvls)+ii),((nms*nxvls)+ii)).
          compute bpath2=mcsopath(((nms*nxvls)+ii),1)&**2.
          compute se2a=diag(indcov((((ii-1)*nxvls)+1):(ii*nxvls),(((ii-1)*nxvls)+1):(ii*nxvls))).
          compute apath2=mcsopath((((ii-1)*nxvls)+1):(ii*nxvls) ,1)&**2.
          compute sesobel=sqrt(apath2*se2b+bpath2*se2a+se2a*se2b).
          compute sobelmat((((ii-1)*nxvls)+1):(ii*nxvls),2)=sesobel.
        end loop.
        release se2b,bpath2,se2a,apath2,sesobel,ii.
        compute sobelmat(:,3)=sobelmat(:,1)&/sobelmat(:,2).
        compute sobelmat(:,4) = 2*(1-cdfnorm(abs(sobelmat(:,3)))).
      end if.
      do if (serial = 0).
        compute rowlbs=t(mnames).
      end if.
      do if (nxvls=1 and nms > 1).
        compute rowlbs={"TOTAL";rowlbs}.
        compute indtemp=indtotal(1,1).
        do if (boot > 0 and nxvls=1).
          !bcboot3 databcbt = indtotal(2:nrow(indtotal),1)/estmte=(indtotal(1,1)*bc)+(9999*(1-bc)).
          compute indtemp={indtemp, bootse,llcit,ulcit}.
        end if.
        do if (mc > 0).
          compute obtmc=indtab(:,1).
          compute indtemp=csum(obtmc).          
          !bcboot3 databcbt = indtotal(:,1).
          compute indtemp={indtemp, bootse,llcit,ulcit}.
        end if.
        compute indtab={indtemp;indtab}.
      end if.
      compute bootlbs={"Effect", "BootSE","BootLLCI","BootULCI"}.
      do if (mc > 0).
        compute bootlbs={"Effect", "MC SE","MC LLCI","MC ULCI"}.
      end if.      
      do if (nxvls = 1).
        do if (contrast <> 0).
          compute inddiff=inddiff(:,2:ncol(inddiff)).
          do if (mc > 0).
            compute inddiff={t(obtmc);inddiff}.
          end if.
          do if (contrast = 4).
            compute inddifft=inddiff*t(contvec).
            compute indtemp=inddifft(1,1).
            do if (boot > 0 or mc > 0).
              do if (mc > 0).
                !bcboot3 databcbt = inddifft(2:nrow(inddifft),1).
              end if.
              do if (boot > 0).
                !bcboot3 databcbt = inddifft(2:nrow(inddifft),1)/estmte=(inddifft(1,1)*bc)+(9999*(1-bc)).
              end if.
              compute indtemp={indtemp, bootse,llcit,ulcit}.
            end if.
            compute indtab={indtab;indtemp}.
          end if.
          do if (contrast = 1 or contrast = 2 or contrast=3).
            compute conkey=make(1,4," ").
            loop i = 1 to ncol(inddiff)-1.
              loop j = (i+1) to ncol(inddiff).           
                compute inddifft=inddiff(:,i)-inddiff(:,j).
                do if (contrast=2).
                  compute inddifft=abs(inddiff(:,i))-abs(inddiff(:,j)).
                end if.
                do if (contrast=3).
                  compute inddifft=inddiff(:,i)+inddiff(:,j).
                end if.
                compute indtemp=inddifft(1,1).
                do if (contrast=1 or contrast=2).
                  compute conkeyt={" ", rowlbs((i+1),1)," minus  ",rowlbs((j+1),1)}.
                end if.
                do if (contrast=3).
                  compute conkeyt={" ", rowlbs((i+1),1),"  plus  ",rowlbs((j+1),1)}.
                end if.
                compute conkey={conkey;conkeyt}.
                do if (boot > 0 or mc > 0).
                  do if (mc > 0).
                    !bcboot3 databcbt = inddifft(2:nrow(inddifft),1).
                  end if.
                  do if (boot > 0).
                    !bcboot3 databcbt = inddifft(2:nrow(inddifft),1)/estmte=(inddifft(1,1)*bc)+(9999*(1-bc)).
                  end if.
                  compute indtemp={indtemp, bootse,llcit,ulcit}.
                end if.
                compute indtab={indtab;indtemp}.
              end loop.
            end loop.
          end if.
          release inddiff.
          compute contlbs=cntname(1:(((nindfx)*(nindfx-1))/2),1).
          compute rowlbs={rowlbs;contlbs}.
        end if.    
        do if ((kk=1) or (kk=2 and ((xdich=1) or (mcx > 0))) or (kk=3)).
          !outform outtodo=indtab/outbig=maxresm.
        end if.
        do if (kk=1 and ydich=0 and diagnose=1).
          !outform outtodo=wadfie/outbig=maxresm.
        end if.
        
        do if (outscree=1). 
          do if (kk=1).
            print indtab/title = "Indirect effect(s) of X on Y:"/rnames=rowlbs/cnames=bootlbs/format=!decpnt.
            compute dfielabs=rowlbs((2-(nms=1)):((2-(nms=1))+nrow(indmake)-1),1).
            
            do if (ydich=0 and diagnose=1).
              print wadfie/title="Cases with greatest influence on indirect effect(s):"/format=!decpnt/clabels="casenum", "dfb_ie"/rnames=dfielabs.  
            end if.
          end if.
          do if (kk = 2 and (xdich=1 or mcx > 0)).
            print indtab/title = "Partially standardized (StandY) indirect effect(s) of X on Y:"/rnames=rowlbs/cnames=bootlbs/format=!decpnt.
          end if.
          do if (kk = 3).
            print indtab/title = "Completely standardized (StandYX) indirect effect(s) of X on Y:"/rnames=rowlbs/cnames=bootlbs/format=!decpnt.
          end if.
        end if.
        do if (normal=1 and sobelok=1 and kk=1).
          !outform outtodo=sobelmat/outbig=maxresm.
          compute sobellab={"Effect",hclab,"Z","p"}.
          compute sobelrlb=rowlbs.
          do if (nms > 1).
            compute sobelrlb=rowlbs(2:(1+nms),1).
          end if.
          do if (outscree=1).
            print sobelmat/title="   Normal theory test for indirect effect(s):"/cnames=sobellab/rnames=sobelrlb/format=!decpnt.
          end if.
        end if.
        do if (contrast <> 0).
          do if ((contrast=1 or contrast = 2 or contrast=3) and kk=efloop ).
            compute conkey=conkey(2:nrow(conkey),:).
            do if (outscree=1).
              print conkey/title = "Specific indirect effect contrast definition(s):"/rnames=contlbs/format=A8.
            end if.
          end if.
          do if (outscree=1).
            do if (contrast = 4 and kk=efloop).
              compute crowlbs=rowlbs(2:(nindfx+1),1).
              print contvec/title = "Specific indirect effect contrast weights:"/cnames=crowlbs/rlabels="(C1)"/format=!decpnt.
            end if.
            do if (contrast = 2 and kk=efloop).
              print/title = "Contrasts are differences between absolute values of indirect effects".
            end if.
          end if.
        end if.
        do if (serial = 1 and kk=efloop).
          compute rowlbst=rowlbs(2:nrow(rowlbs),1).
          do if (outscree=1).
            print indkey/title = "Indirect effect key:"/rnames=rowlbst/format = A8.
            compute dfielabs=rowlbst.
          end if.
        end if.
        else.
        do if (outscree=1).
          do if (kk = 1).
            compute dfielabs={" "}.
            print/title = "Relative indirect effects of X on Y".
          end if.
          do if (kk = 2).
            print/title = "Partially standardized (StandY) relative indirect effect(s) of X on Y:".
          end if.
          do if (kk = 3).
            print/title = "Completely standardized (StandYX) relative indirect effect(s) of X on Y:".
          end if.
        end if.
        
        loop i = 1 to nrow(indmake).
          compute indtabsm=indtab((((i-1)*nxvls)+1):(nxvls*i),:).
          do if (diagnose=1 and ydich=0).
            compute wadfies=wadfie((((i-1)*nxvls)+1):(nxvls*i),:).
          end if.
          compute indkeyt=indkey(i,:).
          do if (outscree=1).
            print indkeyt/title=" "/space=0/format=A8.
          end if.
          do if (bcmat(nrow(bcmat),1)=0)).
            compute direffl2=xcatlab(1:nxvls,1).
          end if.
          !outform outtodo=indtabsm/outbig=maxresm.
          do if (kk=1 and ydich=0 and diagnose=1).
            !outform outtodo=wadfies/outbig=maxresm.
          end if.
          do if (outscree=1).
            print indtabsm/title = " "/cnames=bootlbs/rnames=direffl2/format=!decpnt/space=0.
            do if (ydich=0 and kk=1 and diagnose=1).
              print wadfies/title="Cases with greatest influence on relative indirect effects:"/format=!decpnt/clabels="casenum", "dfb_ie"/rnames=direffl2.  
              print/title="---".
            end if.
          end if.
          do if (normal=1 and sobelok=1 and kk=1).
            compute sobelsm=sobelmat((((i-1)*nxvls)+1):(nxvls*i),:).
            compute sobellab={"Effect",hclab,"Z","p"}.
            !outform outtodo=sobelsm/outbig=maxresm.
            do if (outscree=1).
              print sobelsm/title="      Normal theory test for relative indirect effects:"/cnames=sobellab/rnames=direffl2/format=!decpnt.
            end if.
          end if.
        end loop.
      end if.
      do if (effsize = 1 and boot > 0).
        compute bootres=bootres(2:nrow(bootres),:).
      end if.
    end loop.
     do if (alttotal=1). 
       compute altcnms="Effect". 
       compute totbtvec=totbtvec+bootdir.
       compute alttotfx=t(totbtvec(1,:)).
       do if (boot > 0).
         compute alttotfx=make(ncol(totbtvec),4,0).
         compute alttotfx(:,1)=t(totbtvec(1,:)).
         loop cec=1 to ncol(totbtvec).
           !bcboot3 databcbt = totbtvec(2:nrow(totbtvec),cec)/estmte=(totbtvec(1,cec)*bc)+(9999*(1-bc)).
           compute alttotfx(cec,2:4)={bootse,llcit,ulcit}.
         end loop.
         compute altcnms={altcnms,"BootSE", "BootLLCI","BootULCI"}.
       end if.
        !outform outtodo=alttotfx/outbig=maxresm.
       do if (outscree=1).
         do if (nxvls > 1).
           print/title= "----------"/space=0.
           print alttotfx/title="Relative total effects of X on Y (sum of direct and indirect effects)"/rnames=direffl2/cnames=altcnms/format=!decpnt.
         end if.
         do if (nxvls = 1).
           print alttotfx/title="Total effect of X on Y (sum of direct and indirect effects)"/cnames=altcnms/format=!decpnt.
         end if.
       end if.
     end if.      
  end if.
  /* this is the end of the no moderators loop */.
  
  
  do if (anymod > 0).
   do if (boot = 0).
      compute bootres=obscoeff.
      compute indtab=999.
      /* compute inddiff=999 */.
    end if.
    do if (boot > 0).
      compute bootres={obscoeff;bootres}.
      compute indtab=make(1,4,999).
      /* compute inddiff=make(nrow(bootres),1,999) */.
    end if.
    do if (csum((indmod > 0))=nrow(indmod)).
      do if (outscree=1).
        do if (nxvls > 1).
          do if (xmint=0).
            print/title= "Relative conditional indirect effects of X on Y:".
          end if.
          do if (xmint=1).
            print/title= "Relative (total) natural indirect effects of X on Y:".
          end if.
        end if.
        do if (nxvls = 1).
          do if (xmint=0).
            print/title= "Conditional indirect effects of X on Y:".
          end if.
          do if (xmint=1).
            print/title= "(Total) Natural indirect effect(s) of X on Y:".
          end if.
        end if.
      end if.
    end if.
    do if (csum((indmod > 0)) < nrow(indmod)).
      do if (outscree=1).
        do if (nxvls > 1).
          print/title= "Relative conditional and unconditional indirect effects of X on Y:".
        end if.
        do if (nxvls = 1).
          print/title= "Conditional and unconditional indirect effects of X on Y:".
        end if.
      end if.
    end if.
    compute cftotfx=make(nrow(bootres),nxvls,0).
    loop i = 1 to nrow(indmake).
      /* I added this but not sure */.
      compute indtab=indtab(1,:)*0.
      do if (outscree=1).
        do if (xmint=0).
          print indkey(i,:)/title = "INDIRECT EFFECT:"/format=A8.
        end if.
        do if (xmint=1).
          print indkey(i,:)/title = " "/format=A8/space=0.
        end if.
      end if.
      do if (indmod(i,1)=0).
        loop j = 1 to nxvls.
          compute indtemp=make(nrow(bootres),1,1).
          loop k = 1 to indmake(i,1).
            compute jtemp=1.
            do if (j > 1 and k=1).
              compute jtemp=j.
            end if.
            compute indtemp=indtemp&*bootres(:,pathsfoc(jtemp,indmake(i,(k+1)))).
          end loop.
          compute indeff=indtemp(1,1).
          do if (boot > 0).          
            !bcboot3 databcbt = indtemp(2:nrow(indtemp),1)/estmte=(indtemp(1,1)*bc)+(9999*(1-bc)).
            compute indeff={indeff,bootse,llcit,ulcit}.
          end if.
          compute indtab={indtab;indeff}.
        end loop.
        !outform outtodo=indtab/outbig=maxresm.
        do if (outscree=1).
          do if (nxvls > 1).
            compute indefflb=xcatlab(1:nxvls,1).
            print indtab(2:nrow(indtab),:)/title = " "/clabels="Effect", "BootSE",  "BootLLCI",  "BootULCI"/rnames=indefflb/format=!decpnt/space=0.
          end if.
          do if (nxvls = 1).
            print indtab(2:nrow(indtab),:)/title = " "/clabels="Effect", "BootSE",  "BootLLCI",  "BootULCI"/format=!decpnt/space=0.
          end if.
        end if.
      end if.
      /* end of unmoderated */.
      /* start of moderated */.
      do if (indmod(i,1)>0).
        do if (indmod(i,1)=1).
          compute indmodva=wmodvals.
          compute indprova=wprobval.
          compute condlbs={wnames}.
          compute printw=1.
        else if (indmod(i,1)=2).
          compute indmodva=zmodvals.
          compute indprova=zprobval.
          compute condlbs={znames}.
          compute printz=1.
        else if (indmod(i,1)=3).
          compute cntmp=1. 
          compute printz=1.
          compute printw=1.
          compute indmodva=make((nrow(wmodvals)*nrow(zmodvals)),2,999).
          loop k7 = 1 to nrow(wmodvals).
            loop k8 = 1 to nrow(zmodvals).
              compute indmodva(cntmp,:)={wmodvals(k7,1),zmodvals(k8,1)}.
              compute cntmp=cntmp+1.
            end loop.
          end loop.
          compute condlbs={wnames,znames}.
        end if.
        compute condres=make(nrow(indmodva),1,999).
        do if (boot > 0).
          compute condres=make(nrow(indmodva),4,999).
        end if.
        compute condres={indmodva,condres}.
        /* Here is where the computations start */.
       loop k4 = 1 to nxvls.
         compute imm3=make(nrow(bootres),1,1).
         compute imm4=make(nrow(bootres),1,1).
         compute indcontr=0.  
         do if (indmod(i,1)=3).
           compute tihsw=wprobval.
           compute tihsz=zprobval.
         end if.
         loop k1=1 to nrow(indmodva).
           compute tucker2=make(nrow(bootres),1,1).
           compute imm2=make(nrow(bootres),1,1).
           compute wfirst=0.
           compute zfirst=0.     
           compute immset=0.
           loop k2=1 to indmake(i,1).
             compute colnumb=indmake(i,(k2+1)).
             do if (k2=1).
                compute wbb=make(nrow(bootres),(nwvls*nxvls),0).
                compute zbb=make(nrow(bootres),(nzvls*nxvls),0).
                compute wzbb=make(nrow(bootres),(nwvls*nzvls*nxvls),0).
             end if.
             do if (k2<>1).
                compute wbb=make(nrow(bootres),nwvls,0).
                compute zbb=make(nrow(bootres),nzvls,0).
                compute wzbb=make(nrow(bootres),(nwvls*nzvls),0).
             end if.
             compute cnt=1.
             compute tihs=indlocs(2:((indlocs(1,colnumb))+1),colnumb).
             do if (k2 = 1).
               compute focbb=tihs(1:nxvls,1).
               compute focbb=bootres(:,focbb).
               do if (indmmm(i,1)>0).
                 compute imm=focbb(:,k4).
                 compute condbb=imm.
               end if.
               compute focaddon=make(1,nxvls,0).
               compute focaddon(1,k4)=1.
               compute cnt=cnt+nxvls.
               compute placeh=nxvls.
               do if (indmod(i,1)=1).
                 compute tihsz=make(nrow(wprobval),(nzvls*nxvls),0).
                 compute tihswz=make(nrow(wprobval),(nwvls*nzvls*nxvls),0).
                 do if (pathsw(1,colnumb)=1).
                   compute temp=make(nrow(wprobval),(nxvls*nwvls),0).
                   loop k5 = 1 to nrow(wprobval).
                     loop k6=1 to nwvls.
                       compute temp(k5, (((k4-1)*nwvls)+k6))=wprobval(k5,k6).
                     end loop.
                   end loop.
                   compute indprova={temp,tihsz,tihswz}.
                 else.
                   compute indprova={wprobval,tihsz,tihswz}.
                 end if.
               end if.
               do if (indmod(i,1)=2).
                 compute tihsw=make(nrow(zprobval),(nwvls*nxvls),0).
                 compute tihswz=make(nrow(zprobval),(nwvls*nzvls*nxvls),0).
                 do if (pathsz(1,colnumb)=1).
                   compute temp=make(nrow(zprobval),(nxvls*nzvls),0).
                   loop k5 = 1 to nrow(zprobval).
                     loop k6 =1 to nzvls.
                       compute temp(k5,(((k4-1)*nzvls)+k6))=zprobval(k5,k6).
                     end loop.
                   end loop.
                   compute indprova={tihsw,temp,tihswz}.
                 else.
                   compute indprova={tihsw,zprobval,tihswz}.
                 end if.
               end if.
               do if (indmod(i,1)=3).
                 compute indprova=make((nrow(wprobval)*nrow(zprobval)),((ncol(wprobval)*nxvls)+(ncol(zprobval)*nxvls)+(nwvls*nzvls*nxvls)),0).
                 compute cntemp=1.           
                 loop k7=1 to nrow(wprobval).
                   loop k8 =1 to nrow(zprobval).
                     compute temp=wprobval(k7,:)*focaddon(1,k4).
                     compute indprova(cntemp,(((k4-1)*nwvls)+1):(k4*(nwvls)))=temp.
                     compute temp=zprobval(k8,:)*focaddon(1,k4).
                     compute indprova(cntemp,  ((((k4-1)*nzvls)+1)+(nxvls*nwvls)) : ((((k4-1)*nzvls)+1)+(nxvls*nwvls)+(nzvls-1)))=temp.
                     compute cntemp=cntemp+1.
                   end loop.
                 end loop.
                 do if (pathsz(1,colnumb)=0).
                   compute temp=make(nrow(indprova),(ncol(zprobval)*nxvls),0).
                   compute indprova(:,((ncol(wprobval)*nxvls)+1):((ncol(wprobval)+ncol(zprobval))*nxvls))=temp.
                 end if.
                 do if (pathsw(1,colnumb)=0).
                   compute temp=make(nrow(indprova),(ncol(wprobval)*nxvls),0).
                   compute indprova(:,1:(ncol(wprobval)*nxvls))=temp.   
                 end if.
                 do if (pathswz(1,colnumb)=1).
                    compute cntemp=(ncol(wprobval)*nxvls)+(ncol(zprobval)*nxvls)+((k4-1)*ncol(wprobval)*ncol(zprobval))+1.
                    loop k7=1 to ncol(wprobval).
                      loop k8=1 to ncol(zprobval).
                        compute indprova(:,cntemp)=indprova(:,((ncol(wprobval)*(k4-1))+k7))&*indprova(:,((((k4-1)*ncol(zprobval))+k8)+(nxvls*ncol(wprobval)))).
                        compute cntemp=cntemp+1.
                      end loop.
                    end loop.
                 end if.
               end if.
             end if.
             do if (k2 > 1).
               compute focbb=tihs(1,1).
               compute focbb=bootres(:,focbb).
               do if (indmmm(i,1)>0).
                 compute imm=focbb(:,1).
                 compute condbb=imm.
               end if.
               compute focaddon=1.
               compute cnt=cnt+1.
               compute placeh=1.
               do if (indmod(i,1)=1).
                 compute tihsz=make(nrow(wprobval),nzvls,0).
                 compute tihswz=make(nrow(wprobval),(nwvls*nzvls),0).
                 compute indprova={wprobval,tihsz,tihswz}.
               end if.
               do if (indmod(i,1)=2).
                 compute tihsw=make(nrow(zprobval),nwvls,0).
                 compute tihswz=make(nrow(zprobval),(nwvls*nzvls),0).
                 compute indprova={tihsw,zprobval,tihswz}.
               end if.
               do if (indmod(i,1)=3).
                 compute indprova=make((nrow(wprobval)*nrow(zprobval)),((ncol(wprobval)+ncol(zprobval))+(nwvls*nzvls)),0).
                 compute cntemp=1.           
                 loop k7=1 to nrow(wprobval).
                   loop k8 =1 to nrow(zprobval).
                     compute indprova(cntemp,1:(ncol(wprobval)+ncol(zprobval)))={wprobval(k7,:),zprobval(k8,:)}.
                     compute cntemp=cntemp+1.
                   end loop.
                 end loop.
                 do if (pathsz(1,colnumb)=0).
                   compute temp=make(nrow(indprova),ncol(zprobval),0).
                   compute indprova(:,(ncol(wprobval)+1):(ncol(wprobval)+ncol(zprobval)))=temp.
                 end if.
                 do if (pathsw(1,colnumb)=0).
                   compute temp=make(nrow(indprova),ncol(wprobval),0).
                   compute indprova(:,1:ncol(wprobval))=temp.   
                 end if.
                 do if (pathswz(1,colnumb)=1).
                    compute cntemp=ncol(wprobval)+ncol(zprobval)+1.
                    loop k7=1 to ncol(wprobval).
                      loop k8=1 to ncol(zprobval).
                        compute indprova(:,cntemp)=indprova(:,k7)&*indprova(:,(ncol(wprobval)+k8)).
                        compute cntemp=cntemp+1.
                      end loop.
                    end loop.
                 end if.
               end if.
             end if.  
             do if (pathsw(1,colnumb)) = 1.
               compute wbb=tihs(cnt:(cnt+(placeh*nwvls)-1),1).
               compute wbb=bootres(:,wbb).
               compute immlbs2=wcatlab(1:nwvls,1).
               do if (zfirst=0).
                 compute wfirst=1.
               end if.               
               do if (indmmm(i,1)=1 or indmmm(i,1)=31 or indmmm(i,1)=51).
                 compute imm=wbb(:,1).
                 loop k7 = 1 to nwvls.
                   compute imm={imm,wbb(:,(((k4-1)*nwvls*(k2=1))+k7))}.
                 end loop.
                 compute imm=imm(:,2:ncol(imm)).               
               end if.
               do if (indmmm(i,1)=41 or indmmm(i,1)=51).
                 compute condbb=make(nrow(bootres),1,0).
                 loop k7 = 1 to nwvls.
                   compute condbb={condbb,wbb(:,(((k4-1)*nwvls*(k2=1))+k7))}.
                 end loop.
                 compute condbb=condbb(:,2:ncol(condbb)).
               end if.
               compute cnt=cnt+(placeh*nwvls).
             end if.
             do if (pathsz(1,colnumb)) = 1.
               compute zbb=tihs(cnt:(cnt+(placeh*nzvls)-1),1).
               compute zbb=bootres(:,zbb).
               do if (wfirst=0).
                 compute zfirst=1.
               end if.
               do if (indmmm(i,1) <> 31). 
                 compute immlbs2=zcatlab(1:nzvls,1).
               end if.
               do if (indmmm(i,1)=2 or indmmm(i,1)=31 or indmmm(i,1)=51).
                 do if (indmmm(i,1)=2).
                   compute imm=zbb(:,1).
                 end if.
                 loop k7 = 1 to nzvls.
                   compute imm={imm,zbb(:,(((k4-1)*nzvls*(k2=1))+k7))}.
                 end loop.
                 do if (indmmm(i,1)=2 or indmmm(i,1)=51).
                   compute imm=imm(:,2:ncol(imm)).  
                   do if (indmmm(i,1)=51).
                     compute condbb={condbb,imm}.
                   end if.
                 end if.         
               end if.
               compute cnt=cnt+(placeh*nzvls).
             end if.
             do if (pathswz(1,colnumb)) = 1.
               compute wzbb=tihs(cnt:(cnt+(placeh*nwvls*nzvls)-1),1).
               compute wzbb=bootres(:,wzbb).
               do if (indmmm(i,1)=41).
                 compute imm=wzbb(:,1).              
                 loop k7=1 to nwvls*nzvls.
                   compute imm={imm,wzbb(:,(((k4-1)*nzvls*nwvls*(k2=1))+k7))}.
                 end loop.
               end if.
               do if (indmmm(i,1)=41).
                 compute imm=imm(:,2:ncol(imm)).
                 compute condbb={condbb,imm(:,(ncol(imm)-(nwvls*nzvls)+1):ncol(imm))}. 
               end if.
               compute cnt=cnt+(placeh*nzvls*nwvls).
             end if.
             compute indprobe={focaddon,indprova(k1,:)}.
             compute tucker={focbb,wbb,zbb,wzbb}.
             loop k3=1 to ncol(indprobe).
               compute tucker(:,k3)=tucker(:,k3)*indprobe(1,k3).
             end loop.
             compute tucker2=tucker2&*rsum(tucker).
             do if (indmmm(i,1) = 1 or indmmm(i,1)=2 or indmmm(i,1)=31 or indmmm(i,1)=41 or indmmm(i,1)=51).
               do if (immset=1).
                 do if (ncol(imm2)=1 and ncol(imm) = 1).
                   compute imm2=imm2&*imm.
                 end if.
                 do if (indmmm(i,1)=41 or indmmm(i,1)=51).
                   do if ((ncol(condbb2) > 1) and (ncol(condbb)>1)).
                     compute condbb2t=make(nrow(condbb2),(ncol(condbb2)*ncol(condbb)),-999999).
                     compute k9=1.
                     do if (wfirst=1).
                       loop k7=1 to ncol(condbb2).
                         loop k8 = 1 to ncol(condbb).
                           compute condbb2t(:,k9)=condbb2(:,k7)&*condbb(:,k8).
                           compute k9=k9+1.
                         end loop.
                       end loop.
                     end if.
                     do if (zfirst=1).
                       loop k7=1 to ncol(condbb).
                         loop k8 = 1 to ncol(condbb2).
                           /* compute imm2t(:,k9)=condbb(:,k7)&*condbb2(:,k8) */.
                           compute condbb2t(:,k9)=condbb(:,k7)&*condbb2(:,k8).
                           compute k9=k9+1.
                         end loop.
                       end loop.
                     end if.
                     compute condbb2=condbb2t.
                     release condbb2t.
                   end if.
                   do if ((ncol(condbb2) > 1) and (ncol(condbb)=1)).
                     loop k7 = 1 to ncol(condbb2).
                       compute condbb2(:,k7)=condbb2(:,k7)&*condbb.
                     end loop.
                   end if.
                   do if ((ncol(condbb2) = 1) and (ncol(condbb)>1)).
                     loop k7 = 1 to ncol(condbb).
                       compute condbb(:,k7)=condbb2&*condbb(:,k7).
                     end loop.
                     compute condbb2=condbb.
                   end if.                 
                 end if.
                 do if (ncol(imm2) <> 1 and ncol(imm) <> 1).
                   compute imm2t=make(nrow(imm2),(ncol(imm2)*ncol(imm)),-999999).
                   compute k9=1.
                   do if (wfirst=1).
                     loop k7=1 to ncol(imm2).
                       loop k8 = 1 to ncol(imm).
                           compute imm2t(:,k9)=imm2(:,k7)&*imm(:,k8).
                           compute k9=k9+1.
                       end loop.
                     end loop.
                   end if.
                   do if (zfirst=1).
                     loop k7=1 to ncol(imm).
                       loop k8 = 1 to ncol(imm2).
                           compute imm2t(:,k9)=imm(:,k7)&*imm2(:,k8).
                           compute k9=k9+1.
                       end loop.
                     end loop.
                   end if.
                   compute imm2=imm2t.
                   release imm2t.
                 end if.
                 do if ((ncol(imm2) > 1) and (ncol(imm)=1)).
                   loop k7=1 to ncol(imm2).
                     compute imm2(:,k7)=imm2(:,k7)&*imm.
                   end loop.
                 end if.
                 do if ((ncol(imm2) = 1) and (ncol(imm) > 1)).
                   loop k7=1 to ncol(imm).
                     compute imm(:,k7)=imm2&*imm(:,k7).
                   end loop.
                   compute imm2=imm.
                 end if.
               end if.
               do if (immset=0).
                 compute imm2=imm.
                 do if (indmmm(i,1)=41 or indmmm(i,1)=51).
                   compute condbb2=condbb.
                 end if.
                 compute immset=1.
               end if.
             end if.
             /* that is it for mmm loop */.
           end loop.
           /* end of looping through paths: k2 */.
           compute indtemp=tucker2(1,1).
           do if (indmmm(i,1)=12 or indmmm(i,1)=22).
             compute imm3={imm3,tucker2}.
             do if (k1=nrow(indmodva)).
               compute imm3=imm3(:,2:ncol(imm3)).
               compute immstop=ncol(imm3).
               loop k8=2 to immstop.
                 do if (indmmmt(i,1)=1).
                   compute imm3={imm3,(imm3(:,k8)-imm3(:,1))}.
                 end if.
                 do if (indmmmt(i,1)=2).
                   compute imm3={imm3,(imm3(:,k8)-imm3(:,(k8-1)))}.
                 end if.
                 do if (indmmmt(i,1)=3).
                   compute imm3={imm3,((rsum(imm3(:,(k8:immstop)))/(immstop-k8+1))-imm3(:,(k8-1)))}.
                 end if.
                 do if (indmmmt(i,1)=4).
                   compute imm3={imm3,(imm3(:,k8)-(rsum(imm3(:,1:immstop))/immstop))}.
                 end if.
               end loop.
               do if (indmmmt(i,1)<5).
                 compute imm2=imm3(:,(immstop+1):ncol(imm3)).
                 release imm3.
               end if.
             end if.
           end if.      
          /* -1 used to be 0 */.
           do if (indmmm(i,1)>-1 and (contrast = 1 or contrast = 2 or contrast=3)).
             compute imm4={imm4,tucker2}.
             do if (k1=nrow(indmodva) and k1 > 1).
               compute imm4=imm4(:,2:ncol(imm4)).
               compute immstop=ncol(imm4).
               compute condcont=make((immstop*(immstop-1)/2),6,-999).              
               loop k8 = 1 to (immstop-1).
                 loop k9 = (k8+1) to immstop.                   
                   do if (contrast=1).
                     compute imm4={imm4,(imm4(:,k9)-imm4(:,k8))}.
                   end if.
                   do if (contrast=2).
                     compute imm4={imm4,(abs(imm4(:,k9))-abs(imm4(:,k8)))}.
                   end if.
                   do if (contrast=3).
                     compute imm4={imm4,(imm4(:,k9)+imm4(:,k8))}.
                   end if.
                   compute condcont((ncol(imm4)-immstop),1)=imm4(1,k9).
                   compute condcont((ncol(imm4)-immstop),2)=imm4(1,k8).
                 end loop.
               end loop.
               compute imm4=imm4(:,(immstop+1):ncol(imm4)).   
               loop k8=1 to ncol(imm4).
                 compute condcont(k8,3)=imm4(1,k8).
                 do if (boot > 0).
                   !bcboot3 databcbt = imm4(2:nrow(imm4),k8)/estmte=(imm4(1,k8)*bc)+(9999*(1-bc)).
                   compute condcont(k8,4:6)={bootse,llcit,ulcit}.  
                 end if.
               end loop.
               do if (boot=0).
                 compute condcont=condcont(:,1:3).
               end if.
               compute indcontr=1.
             end if.
           end if.
           do if (xmint=1).
             do if (k1=(k4+1)).
               compute cftotfx(:,k4)=cftotfx(:,k4)+(tucker2*xscaling).
             end if.
           end if.
           do if (boot > 0).
             compute tucker2=tucker2*xscaling.
             !bcboot3 databcbt = tucker2(2:nrow(tucker2),1)/estmte=(tucker2(1,1)*bc)+(9999*(1-bc)).
             compute indtemp={indtemp*xscaling, bootse,llcit,ulcit}.  
           end if.
           do if (boot=0).
             compute indtemp=indtemp*xscaling.
           end if.
           compute condres(k1,(ncol(indmodva)+1):ncol(condres))=indtemp.
         end loop.
         /* end of looping through indirect effects: k1 */.
         /* Here is where the computations end */.
         do if (xmint=0).
           !outform outtodo=condres/outbig=maxresm.
         end if.      
          compute condlbs={condlbs,"Effect"}.
          do if (boot > 0).
            compute condlbs={condlbs,"BootSE",  "BootLLCI",  "BootULCI"}.
          end if.          
          do if (xmint = 1).
            do if (k4 = 1).
              compute natindfx=condres((2+(k4-1)),2:ncol(condres)).
            end if.
            do if (k4 > 1).
              compute natindfx={natindfx;condres((2+(k4-1)),2:ncol(condres))}.
            end if.
            do if (k4=nxvls).
              !outform outtodo=natindfx/outbig=maxresm.
            end if.
          end if.          
          do if (outscree=1 and xmint=0).
            do if (nxvls=1).
              print condres/title=" "/cnames=condlbs/format=!decpnt/space=0.
            else.
              compute condrlb=make(nrow(condres),1,xcatlab(k4,1)).
              print condres/title=" "/cnames=condlbs/rnames=condrlb/format=!decpnt/space=0.
            end if.
          end if.
          do if (outscree=1 and xmint=1 and k4=nxvls).
            do if (nxvls=1).
              print natindfx/title=" "/clabels="Effect", "BootSE", "BootLLCI", "BootULCI"/format=!decpnt/space=0.
            else.
              compute condrlb=make(nrow(condres),1,xcatlab(k4,1)).
              print natindfx/title=" "/clabels="Effect", "BootSE", "BootLLCI", "BootULCI"/rnames=direffl2/format=!decpnt/space=0.
            end if.
          end if.
          compute dichadj=0.
          compute immcat=0.
          do if ((indmmm(i,1)>0) and (xmint=0)).
            do if (indmmm(i,1)=1 or indmmm(i,1)=12 or indmmm(i,1)=31).
              do if (wdich=1 and mcw=0).
                do if (indmmm(i,1) <> 12).
                  compute imm2(:,1)=imm2(:,1)*(wmax-wmin).
                end if.
                do if (indmmm(i,1) <> 31).
                  compute dichadj=1.
                end if.
              end if.           
              do if ((mcw = 1 or mcw = 2) and indmmm(i,1) <> 31)).
                compute immcat=1.
              end if.
            end if.
            do if (indmmm(i,1)=2 or indmmm(i,1)=22 or indmmm(i,1)=31).
              do if (zdich=1 and mcz=0).
                do if (indmmm(i,1) = 31).
                  compute imm2(:,(nwvls+1):ncol(imm2))=imm2(:,(nwvls+1):ncol(imm2))*(zmax-zmin).
                end if.
                do if (indmmm(i,1) = 2).
                  compute imm2(:,1)=imm2(:,1)*(zmax-zmin).
                end if.
                do if (indmmm(i,1) <> 31).
                  compute dichadj=1.
                end if.
              end if.
              do if ((mcz = 1 or mcz = 2) and indmmm(i,1) <> 31)).
                  compute immcat=1.
              end if.
            end if.
            compute immtemp2=t(imm2(1,:)).
            compute immtemp=immtemp2.
            compute immlbs={"Index"}.
            do if (boot > 0).
              compute immtemp=make(1,3,0).
              loop k7=1 to ncol(imm2).
                !bcboot3 databcbt = imm2(2:nrow(imm2),k7)/estmte=(imm2(1,k7)*bc)+(9999*(1-bc)).
                compute temp={bootse,llcit,ulcit}.
                compute immtemp={immtemp;temp}.
              end loop.  
              compute immtemp=immtemp(2:nrow(immtemp),:).
              compute immtemp={immtemp2,immtemp}.
              compute immlbs={immlbs,"BootSE",  "BootLLCI",  "BootULCI"}.
            end if.
            !outform outtodo=immtemp/outbig=maxresm.
            do if (dichadj=0 and immcat=0 and indmmmt(i,1)<>5 and indmmm(i,1) < 100 ).
              do if (outscree=1).
                do if (indmmm(i,1) < 30).
                  print immtemp/title="      Index of moderated mediation:"/cnames=immlbs/rnames=immlbs2/format=!decpnt.
                end if.
                do if (indmmm(i,1) = 31).
                  compute immlbs2={immlbs2;zcatlab(1:nzvls,1)}.
                  print immtemp/title="      Indices of partial moderated mediation:"/cnames=immlbs/rnames=immlbs2/format=!decpnt.
                end if.
              end if.
             /* work on this */.
              do if ((nzvls = 1) and (nwvls = 1)).

              do if (indmmm(i,1) = 41 or indmmm(i,1)=51).
                loop k7=1 to nwvls.
                  compute immlbs2=zcatlab(1:nzvls,1).
                  compute immtemp2=immtemp((((k7-1)*nzvls)+1):(((k7-1)*nzvls)+nzvls),:).       
                  do if (outscree=1).       
                    do if (nwvls > 1).
                      compute primodv={"        ", wcatlab(k7,1)}.
                      print primodv/title="      Primary moderator:"/format=A8.
                    end if.
                    do if (nzvls=1).
                      print immtemp2/title="      Index of moderated moderated mediation"/cnames=immlbs/format=!decpnt.
                    else.
                      print immtemp2/title="      Indices of moderated moderated mediation"/cnames=immlbs/rnames=immlbs2/format=!decpnt.
                    end if.
                  end if.
                  compute cmmtemp=make(nrow(zprobval),4,0).        
                  loop k8=1 to nrow(zprobval).
                    compute condbb3=condbb2(:,((nwvls+1)+((k7-1)*nzvls)):((nwvls+1)+((k7-1)*nzvls)+(nzvls-1))).
                    do if (ncol(zprobval) > 1).
                      compute condbb3=condbb3*mdiag(zprobval(k8,:)).
                    else.
                      compute condbb3=condbb3&*zprobval(k8,:).
                    end if.
                    compute condbb3={condbb2(:,k7),condbb3}.
                    compute icmm=rsum(condbb3).
                    compute cmmtemp(k8,1)=icmm(1,1).
                    do if (boot > 0).
                      !bcboot3 databcbt = icmm(2:nrow(icmm),1)/estmte=(icmm(1,1)*bc)+(9999*(1-bc)).
                      compute cmmtemp(k8,2:4)={bootse,llcit,ulcit}.
                    end if.
                  end loop.
                  compute cmmtemp={zmodvals,cmmtemp}.
                  do if (boot=0).
                    compute cmmtemp=cmmtemp(:,1:2).
                  end if.                    
                  compute cmmlbs={znames,immlbs}.
                  !outform outtodo=cmmtemp/outbig=maxresm.
                  do if (outscree=1).
                    print cmmtemp/title="      Indices of conditional moderated mediation by W"/cnames=cmmlbs/format=!decpnt.
                  end if.
                end loop.
              end if.

              end if.

            end if.
            do if (dichadj=1 or immcat=1 and indmmm(i,1) < 30).
              do if (outscree=1).
                print immtemp/title="Index of moderated mediation (difference between conditional indirect effects):"/cnames=immlbs/rnames=immlbs2/format=!decpnt.
              end if.
            end if.
          end if.
          do if (indcontr=1).
               !outform outtodo=condcont/outbig=maxresm.
               compute condctlb={"Effect1","Effect2","Contrast","BootSE", "BootLLCI","BootULCI"}.
               do if (outscree=1).
                 print condcont/title=" Pairwise contrasts between conditional indirect effects (Effect1 minus Effect2)"/cnames=condctlb/format=!decpnt.
               end if.
          end if.
          do if (outscree=1).
            /* print/title= "---"/space=0 */.
          end if.
        end loop.
      end if.
      /* end of moderated */.
    end loop.
    
    do if (alttotal=1).
      compute altcnms="Effect". 
      compute alttotfx=t(cftotfx(1,:))+obnatdfx.
      do if (boot > 0).
        compute alttotfx=make(ncol(cftotfx),4,0).
        compute alttotfx(:,1)=t(cftotfx(1,:))+obnatdfx.
        compute natdirbt(1,:)=t(obnatdfx).
        compute cftotfx=cftotfx+natdirbt.
        loop cec=1 to ncol(cftotfx).
          !bcboot3 databcbt = cftotfx(2:nrow(cftotfx),cec)/estmte=(cftotfx(1,cec)*bc)+(9999*(1-bc)).
          compute alttotfx(cec,2:4)={bootse,llcit,ulcit}.
        end loop.
        compute altcnms={altcnms,"BootSE", "BootLLCI","BootULCI"}.
      end if.
      !outform outtodo=alttotfx/outbig=maxresm.
      do if (outscree=1).
        print/title="----------".
        do if (nxvls > 1).
          print alttotfx/title="Relative total effects of X on Y (sum of direct and indirect effects)"/cnames=altcnms/rnames=direffl2/format=!decpnt.
        end if.
        do if (nxvls = 1).
          print alttotfx/title="Total effect of X on Y (sum of direct and indirect effects)"/cnames=altcnms/format=!decpnt.
        end if.
      end if.
    end if.
  end if.
  /* This is the end of the moderated loop */.
end if.
!ENDDEFINE.


DEFINE !modeles (y=!charend('/')/x=!charend('/')/type=!charend('/')/full=!charend('/')).
    compute eiverr=0.
    do if (!type=1).
      compute n1=nrow(!x).    
      compute n=n1.
      compute invtxx=inv(t(!x)*!x).
      compute hatmat=invtxx*t(!x).
      do if (eivdo=0).
        compute b = hatmat*!y. 
      end if.    
      do if (eivdo=1).
        compute relcorv=1.
        do if (eiv > 9).
          compute relcorv=(1-(ncol(!x)/nrow(!x))).
        end if.
        compute varx = ((nrow(!x)*cssq(!x))-(csum(!x)&**2))/((nrow(!x)-1)*nrow(!x)).
        compute varx=varx*(nrow(!x)-1)/nrow(!x).     
        compute varxr=nrow(!x)*(1-releiv)&*varx.        
        compute aiev=(t(!x)*!x)-(relcorv)*mdiag(varxr).
        compute itsh=inv(aiev).
        compute hatval2=make(n1,1,0).
        /* do if (eiv=5 or eiv=6) */.
        /* loop i3= 1 to n1 */.
        /*    compute hatval2(i3,1)=!x(i3,:)*itsh*t(!x(i3,:)) */.
        /* end loop */.     
        /* end if */.     
        do if ((eiv=2) or (eiv=3) or (eiv=12) or (eiv=13)).
          loop i3= 1 to n1.
            compute hatval2(i3,1)=!x(i3,:)*invtxx*t(!x(i3,:)).
          end loop.    
        end if.
        compute b=itsh*t(!x)*!y.    
        compute yeiv=!y.                  
        do if (!full=0).
          compute mseb=((t(!y)*!y)-(t(b)*aiev*b))/(n-ncol(!x)).
          compute sstotalb = t(!y-(csum(!y)/n1))*(!y-(csum(!y)/n1)).
          compute r2b=1-((mseb*(n1-ncol(!x)))/sstotalb). 
          do if (r2b < 0 or mseb < 0)).
            compute eiverr=1.
          end if.  
        end if. 
      end if.
      compute modres=b.
        
      do if (!full=1).
        do if (eivdo=1).
          compute mse=((t(!y)*!y)-(t(b)*aiev*b))/(n1-ncol(!x)).
          compute sstotal = t(!y-(csum(!y)/n1))*(!y-(csum(!y)/n1)).
          compute r2=1-((mse*(n1-ncol(!x)))/sstotal).
          do if (r2 < 0 or mse < 0)).
            compute eiverr=1.
            compute errcode(errs,1)=85.
            compute errs=errs+1.
            compute criterr=1.
          end if.
          do if ((eiv=5 or eiv=15)) and eiverr=0).
            compute varb=mse*inv(aiev)*t(!x)*!x*inv(aiev).
          end if.  
          do if (((eiv=4) or (eiv=14)) and eiverr=0).
            compute varxr2=mdiag(varxr).
            compute eivtmp=(t(b)*varxr2)&*(t(b)*varxr2).
            compute eivtmp=mdiag(eivtmp).
            compute varb=((mse)*inv(aiev))+((inv(aiev)/nrow(!x)) *((varxr2*mse)+(varxr2*b*t(b)*varxr2)+2*eivtmp)*inv(aiev)).
          end if.
          do if ((eiv=0 or eiv=1 or eiv=2 or eiv=3 or eiv=10 or eiv=11 or eiv=12 or eiv=13) and eiverr=0).
            compute resideiv=yeiv-!x*b.
            compute hb1eiv=!x.
            compute hb2eiv=!x.
            loop ieiv = 1 to ncol(!x).
              do if (eiv=3 or eiv=0 or eiv=1 or eiv=13 or eiv=10 or eiv=11).
                compute hb1eiv(:,ieiv)=hb1eiv(:,ieiv)&*(resideiv/(1-hatval2)).
              end if.
              do if (eiv=2 or eiv=12).
                compute hb1eiv(:,ieiv)=hb1eiv(:,ieiv)&*(resideiv/sqrt(1-hatval2)).
              end if.
              compute xmeaneiv=csum(!x(:,ieiv))/n1.
              compute xmeaneiv=make(n1,1,xmeaneiv).
              compute hb2eiv(:,ieiv)=(hb2eiv(:,ieiv)-xmeaneiv)&**2.
            end loop.
            compute hbeiv=(1-releiv)&*t(b).
            compute hbeiv=mdiag(hbeiv).
            compute hb3eiv=make(n1,ncol(!x),1).
            compute hb3eiv=hb3eiv*hbeiv.
            compute Heiv=hb1eiv+(hb2eiv&*hb3eiv).
            compute varb=inv(aiev)*t(Heiv)*Heiv*inv(aiev).
            do if ((eiv=1) or (eiv=11)).
              compute varb=varb*(n1/(n1-ncol(!x))).
            end if.
            do if (eiv < 4 and eivhcnt=0).
              compute notecode(notes,1) = 4.
              compute notes = notes + 1.
              compute eivhcnt=1.
            end if.
            release hb1eiv, hb2eiv, hb3eiv.
          end if.   
          /* compute hclab="se" */.
          /* compute hcflab="F" */.          
          compute hclab={"se(HC0)","se(HC1)","se(HC2)","se(HC3)","se","se"}.
          compute hclab=hclab(1,(eiv+1)).
          compute hcflab={"F(HC0)","F(HC1)","F(HC2)","F(HC3)","F","F"}.
          compute hcflab=hcflab(1,(eiv+1)).
          
          compute adjr2=0.  
          compute dfres=(n1-(ncol(!x)))*(1-robustse)+((nclus-1)*robustse).     
        end if.
        
        /* start no eiv1 */.
        do if (eivdo=0).        
          compute ypredv=!x*b.
          compute resid=!y-ypredv.                  
          /* compute xdattmp=!x */.
          compute hatval=make(n1,1,0).
          loop i3= 1 to n1.
            compute hatval(i3,1)=!x(i3,:)*invtxx*t(!x(i3,:)).
          end loop.                      
          compute yh=(ypredv-hatval&*!y)&/(1-hatval).
          compute ex=resid&/(1-hatval)).          
          compute modresdf=make(n1,ncol(!x),0).                    
          loop i3 = 1 to ncol(!x).
            compute modresdf(:,i3)=t(hatmat(i3,:)&*t(ex)).
          end  loop.
          compute dfres=(n1-(ncol(!x)))*(1-robustse)+((nclus-1)*robustse).
          compute sstotal = t(!y-(csum(!y)/n1))*(!y-(csum(!y)/n1)).
          compute ssresid = csum((resid)&**2).
          compute mse=ssresid/(n1-ncol(!x)).
          compute skewres=(n1*csum(resid&**3))/((n1-2)*ssresid*sqrt(ssresid/(n1-1))).
          compute kurtres=(((n1+1)*n1*(n1-1))/((n1-2)*(n1-3)))*(csum(resid&**4)/(ssresid*ssresid)).  
          compute kurtres=kurtres-(3*((n1-1)*(n1-1))/((n1-2)*(n1-3))).
          compute skewse=sqrt(((6*n1)*(n1-1))/((n1-2)*(n1+1)*(n1+3))).
          compute kurtse=sqrt(   (24*n1*(n1-1)*(n1-1))/((n1-3)*(n1-2)*(n1+3)*(n1+5))).
          compute skewres={skewres,kurtres;skewse,kurtse}.
          compute ssreg=sstotal-ssresid.      
          compute sumtable={ssreg,(ncol(!x)-1),(ssreg/(ncol(!x)-1));ssresid,(n1-(ncol(!x))),(ssresid/(n1-(ncol(!x))));sstotal,(n1-1),(sstotal/(n1-1))}.
          compute r2 = (sstotal-ssresid)/sstotal.
          do if (r2 < 0).
            compute r2=0.
          end if.
          compute adjr2 = 1-((1-r2)*(n1-1)/(n1-ncol(!x))).
          compute mndevy=(!y-(csum(!y)/n1)).
          compute msresch= ((ex&*resid)-mse)/(n1-ncol(!x)-1). 
          compute skewwo=make(n1,1,9999).
          /* compute scresnoi=make(n1,1,0) */.
          /* loop wadf2=1 to nrow(modresdf) */.
          /*   compute cresnoi=(!y-(!x*(b-t(modresdf(wadf2,:)))))&**3 */.
          /*  compute scresnoi(wadf2,1)=(n1-1)*(csum(cresnoi)-cresnoi(wadf2,1)) */.      
          /* end loop */.        
          /* compute denomsk=(n1-3)*(ssresid-(ex&*resid))&*sqrt((ssresid-(ex&*resid))/(n1-2)) */.     
          /* compute skewwo=skewres(1,1)-(scresnoi&/denomsk) */.    
          /* end new here */.
          compute rchngwo=r2-(1-((ssresid-(ex&*resid))&/(sstotal-((n1*(mndevy&*mndevy))/(n1-1))))).
          compute stri=resid&/sqrt((1-hatval)*mse).
          compute tricheck=(((n1-ncol(!x))-1)/((n1-ncol(!x))-(stri&*stri))).
          compute tricheck=csum(tricheck < 0).
          do if (tricheck <> 0).
            do if (triprob=0).
              compute notecode(notes,1) = 40.
              compute notes = notes + 1.
              compute triprob=1.      
            end if.        
          end if.
          do if (tricheck = 0).
            compute tri=stri&*sqrt(((n1-ncol(!x))-1)/((n1-ncol(!x))-(stri&*stri))).
          end if.
          compute mahal=(n1*hatval-1)*((n1-1)/n1).
          compute cook=(stri&*stri)&*(hatval&/((1-hatval)*(ncol(!x)))).
          do if (triprob=0).
            compute regdiag2={ypredv,resid,ex,stri,tri,hatval,mahal,cook,msresch,rchngwo,skewwo,modresdf}. 
          end if. 
          do if (diagnose=1 or savediag=1).
            /* find the most extreme dfbeta values */. 
            compute wadf=make(ncol(modresdf),2,0).
            loop wadf1=1 to ncol(modresdf).
              loop wadf2=1 to nrow(modresdf).
                do if ((abs(modresdf(wadf2,wadf1))) > abs(wadf(wadf1,2))).
                  compute wadf(wadf1,2)=modresdf(wadf2,wadf1).
                  compute wadf(wadf1,1)=rownum(wadf2,1).
                end if.               
              end loop.
            end loop.      
            compute alldfbs={alldfbs,modresdf}.
            do if (tricheck=0).  
              compute ptri = 2*(1-tcdf(abs(tri), ((n1-ncol(!x))-1))).
              compute bp=n1*ptri.
              compute bpm=bp-1.
              compute bp=bp-(bpm&*(bp > 1)).
              compute bp = {abs(tri),bp,rownum}.
              compute temp = bp.
              compute temp(GRADE(bp(:,1)),:) = bp.
              compute bp = temp.
              compute bp=bp(n1,:).   
              compute sqruncr=ssresid/n1.   
              compute qprime=((resid&*resid)-sqruncr)/sqruncr.
              compute bptest=0.5*(t(qprime)*!x*hatmat*qprime).
              compute bprobust=bptest*(2/((t(qprime)*qprime)/n1)).
              compute bpresult={bptest,(ncol(!x)-1), (1-chicdf(bptest,(ncol(!x)-1)));bprobust,(ncol(!x)-1),(1-chicdf(bprobust,(ncol(!x)-1)))}.
            end if.        
          end if.    
        
          /* compute msreg=ssreg/ncol(!x)) */.
          !hcest3 x=x/resid=resid/hc=hc/mse=mse.
        
        end if.
        /* end no eiv1 */.
        
        do if (eiverr=0).
          compute seb=sqrt(diag(varb)).
          compute trat = b&/seb.
          compute p = 2*(1-tcdf(abs(trat), (dfres))).
          /* compute tval =  sqrt(dfres* (exp((dfres-(5/6))*((xp2/(dfres-(2/3)+(.11/dfres)))*(xp2/(dfres-(2/3)+(.11/dfres)))))-1)) */.
          compute tval=-idf.t(alpha2,dfres).
          compute modres={modres,seb,trat,p}.
          compute modres={modres,(b-tval&*seb),(b+tval&*seb)}.
          compute modresl={"coeff",hclab,"t","p","LLCI","ULCI"}.
          compute lmat = ident(ncol(!x)).
          compute lmat = lmat(:,2:ncol(lmat)).
          /*compute ife=1-sqrt(1-r2) */.
          /* new here */.
          /* compute modsum={sqrt(r2),r2,adjr2,99999,99999,sqrt(mse),ife} */.
          compute modsum={sqrt(r2),r2,adjr2,99999,99999,sqrt(mse)}.
          compute modsuml={"R","R-sq","Adj R-sq",hcflab,"p", "SEest","IFE"}.      
          /* end new here */.
          do if (ncol(lmat) > rank(varb)).
            compute intncd=1.
          end if.
          do if (ncol(lmat) <= rank(varb)).
          compute fratio = (t(t(lmat)*b)*inv(t(lmat)*varb*lmat)*((t(lmat)*b)))/(ncol(!x)-1).
          compute pfr = 1-fcdf(fratio,(ncol(!x)-1),dfres).
          /* new here */.
          /* compute modsum={sqrt(r2),r2,adjr2,fratio,pfr, sqrt(mse),ife} */.
         compute modsum={sqrt(r2),r2,adjr2,fratio,pfr, sqrt(mse)}.
          /* end new here */.
          end if. 
        end if.
        
        /* start no eiv2 */.
        do if (eivdo=0).                                    
          compute sehomo=sqrt(diag(mse*invtxx)).
          compute srrat=b&/sehomo.
          compute sign=b&/abs(b).
          compute semir=srrat*sqrt((1-r2)/(n1-ncol(!x))).      
          compute partialr=sqrt((srrat&*srrat)&/((srrat&*srrat)+(n1-ncol(!x))))&*sign.     
          compute iv5={!y,yh}.      
          compute numdesc = (n*sscp(iv5)-(t(csum(iv5))*(csum(iv5)))).
          compute dendesc = sqrt(t((n*cssq(iv5))-(csum(iv5)&**2))*((n*cssq(iv5))-(csum(iv5)&**2))).
          compute corall = numdesc&/dendesc. 
          compute lvout1=corall(2,1).
          compute browne=0.
          do if (adjr2 >= 0).
            compute rho4=(adjr2*adjr2)-((2*(ncol(!x)-1)*(1-adjr2)*(1-adjr2))/((n1-1)*(n1-(ncol(!x)-1)+1))).
            do if (rho4 >=0).
              compute browne=sqrt(((n1-(ncol(!x)-1)-3)*rho4+adjr2)/((n1-2*(ncol(!x)-1)-2)*adjr2+(ncol(!x)-1))).
            end if.
          end if.
          compute newsse=ssresid&*(n1-(ncol(!x)-1)-2)&/((tri&*tri)+n1-(ncol(!x)-1)-2).
          compute mycv=csum(!y)/n1.
          compute dleave=(mycv-!y)&/(n1-1).
          compute yd=!y-mycv.
          compute newtss=sstotal+n1&*(dleave&*dleave)-((!y-mycv-dleave)&*(!y-mycv-dleave)).
          compute newvhaty = (newtss&/n1)&*(1-newsse&/newtss).
          compute zh=(yh-dleave-mycv)&/sqrt(newvhaty).                    
          compute temptxd={y,x(:,2:ncol(x)),yh,zh}.
          /* compute ones=make(n1,1,1) */.
          compute sdvec = t(sqrt(((n*cssq(temptxd))-(csum(temptxd)&**2))/((n-1)*n))).
          compute numdesc = (n*sscp(temptxd)-(transpos(csum(temptxd))*(csum(temptxd)))).
          compute dendesc = sqrt(transpos((n*cssq(temptxd))-(csum(temptxd)&**2))*((n*cssq(temptxd))-(csum(temptxd)&**2))).
          compute corall = numdesc&/dendesc.
          compute corbxx=corall(2:ncol(x),2:ncol(x)).
          compute tolstat=make((ncol(x)-1),1,9999).
          loop tolcom=1 to ncol(x)-1.
            compute riichoose=ident(ncol(x)-1).
            compute riichoose(tolcom,tolcom)=0.
            compute riix=riichoose*corbxx*t(riichoose).
            compute riix(tolcom,tolcom)=1.
            compute riiy=corbxx(:,tolcom)&*diag(riichoose).
            compute tolstat(tolcom,1)=1-(t(riiy)*inv(riix)*riiy).
          end loop.
          compute sdvec=sdvec(1:(nrow(sdvec)-2),1).
          compute tildeb=(modres(:,1)&*sdvec)/sdvec(1,1). 
          compute tildeby=modres(:,1)/sdvec(1,1).
          compute tildebx=modres(:,1)&*sdvec.
          compute crossr={browne,corall((nrow(corall)-1),1),corall(nrow(corall),1)}. 
          compute zppout={corall(2:(nrow(corall)-2),1),semir(2:nrow(semir),1),partialr(2:nrow(partialr),1),tildeb(2:nrow(tildeb),1),tildeby(2:nrow(tildeby),1),tildebx(2:nrow(tildebx),1)}.   
          compute zppout2={zppout(:,2)&**2,zppout(:,3)&**2,(zppout(:,2)&**2)/(1-r2)}.   
          /* new here */.
          /* compute ipfe=1-(mse/sqrt((sstotal*(zppout2(:,1))+ssresid)/(n1-ncol(x)+1))) */.
          /* compute zppout2={zppout2,ipfe} */.           
          /* end new here */.                  
        end if.
        /* end no eiv2 */. 
        
        do if (ssquares=0 and eiverr=0).
          compute modsum={sqrt(r2),r2,mse,99999,99999,99999,99999}.
          compute modsuml={"R","R-sq","MSE",hcflab,"df1","df2", "p"}.
          do if (ncol(lmat) <= rank(varb)).
            compute modsum={sqrt(r2),r2,mse,fratio,(ncol(!x)-1),dfres,pfr}.
            do if (eivdo=1).
              compute modsuml=modsuml(1,2:ncol(modsuml)).
              compute modsum=modsum(1,2:ncol(modsum)).              
            end if.
          end if.
        end if.
        
        do if (diagnose=1 and tricheck=0).
          compute maxnmin={cmin(ypredv),cmax(ypredv);cmin(resid),cmax(resid);cmin(tri),cmax(tri)}.
          compute trmin=cmin(tri).
          do if ((bp(1,1)=abs(trmin)) and (trmin < 0)).
            compute bp(1,1)=-bp(1,1).
          end if.             
        end if.
        
      end if.
   end if.
  /* for logistic Y model */.
   do if (!type = 2 or !type=3).
     compute xlp=!x.
     compute ylp=!y.
     compute pt2 = make(nrow(ylp),1,(csum(ylp)/nrow(ylp))).
     do if (!type=2).
       compute LL3 = ylp&*ln(pt2)+(1-ylp)&*ln(1-pt2).
     end if.  
     compute LL3 = -2*csum(LL3).
     compute bt1 = make(ncol(xlp),1,0).
     compute LL1 = 0.
     compute pt1 = make(nrow(ylp),1,0.5).
     compute pt1lp=pt1.
     loop jjj = 1 to iterate.
       compute xlptmp=t(xlp).
       compute vecprb=pt1lp&*(1-pt1lp).     
       loop kkk=1 to ncol(xlp).
         compute xlptmp(kkk,:)=xlptmp(kkk,:)&*t(vecprb).
       end loop.
       compute b = bt1+inv(xlptmp*xlp)*t(xlp)*(ylp-pt1lp).
       do if (!type=2).
         compute xlpb=xlp*b.
         compute xlpbt=(xlpb > -709.7).
         compute xlpb709=(1-xlpbt)*(-709.7).
         compute xlpb=(xlpb&*xlpbt)+xlpb709.
         compute pt1lp = 1/(1+exp(-(xlpb))).
       end if.
       compute itprob = csum((pt1lp < .00000001) or (pt1lp > .9999999)).
       do if (itprob > 0).
         loop kkk = 1 to nrow(pt1lp).
           do if (pt1lp(kkk,1) > .9999999).
             compute pt1lp(kkk,1) = .9999999.
           end if.
           do if (pt1lp(kkk,1) < .00000001).
             compute pt1lp(kkk,1) = .00000001.
           end if.
         end loop.
         compute itprob = 0.
       end if.
       do if (itprob = 0).
         do if (!type=2).         
           compute LL = ylp&*ln(pt1lp)+(1-ylp)&*ln(1-pt1lp).
         end if.
         compute LL2 = -2*csum(ll).
       end if.
       do if (abs(LL1-LL2) < converge).
         do if (!full=1).
           compute xlptmp=t(xlp).
           compute vecprb=pt1lp&*(1-pt1lp).
           loop kkk=1 to ncol(xlp).
             compute xlptmp(kkk,:)=xlptmp(kkk,:)&*t(vecprb).
           end loop.          
           compute varb = inv(xlptmp*xlp).
           compute seb = sqrt(diag(varb)).
           release xlptmp.
         end if.
         break.
       end if.
       compute bt1 = b.
       compute LL1 = LL2.
     end loop.
     compute modres=b.
     do if (jjj > iterate).
       compute itprob = 2.
       do if (booting=0).
          compute iterrmod=1.
       end if.
       do if (booting=1).
         compute bootiter=1.
       end if.
       do if (itprobtg=0).   
         compute itprobtg=1.
         compute errcode(errs,1) = 47.
         compute errs = errs + 1.
         do if (booting = 0 and !full=1).
           compute vt1 = mdiag(pt1lp&*(1-pt1lp)).
           compute varb = inv(t(xlp)*vt1*xlp).
           compute seb = sqrt(diag(varb)).
         end if.
       end if.
     end if.
     do if (!full=1).
       compute trat = b&/seb.
       compute dfres=nrow(xlp).     
       compute p = 2*(1-cdfnorm(abs(trat))).
       compute modres={modres,seb,trat,p}. 
       compute modres={modres,(b-xp2&*seb),(b+xp2&*seb)}.
       compute pvchi=1-chicdf((LL3-LL2),(nrow(modres)-1)).
       compute mcF = (LL3-LL2)/LL3.
       compute cox = 1-exp(-(LL3-LL2)/nrow(xlp)).
       compute nagel = cox/(1-exp(-(LL3)/nrow(xlp))).
       compute modsum={LL2,(LL3-LL2),(nrow(modres)-1),pvchi, mcF,cox,nagel}.
       compute modsuml={"-2LL","ModelLL", "df", "p", "McFadden", "CoxSnell", "Nagelkrk"}.
       compute modresl={"coeff","se","Z","p","LLCI","ULCI"}.
     end if.
   end if.
!ENDDEFINE.

DEFINE !makdumm (dd=!charend('/')/method=!charend('/')/custcodv=!charend('/') !default (0)/custcode=!charend('/') !default("-999")).
  compute dd=!dd.
  compute temp = dd.
  compute temp(GRADE(dd(:,2)),:) = dd.
  compute dd = temp.
  compute dummy = design(dd(:,2)).
  compute nvls = ncol(dummy).
  compute nnvls = csum(dummy).
  compute mnvls = cmin(t(nnvls)).
  compute conmat1=1.
  do if (mnvls < 2).
    compute errcode(errs,1) = 5.
    compute errs = errs + 1.
    compute criterr = 1.
  end if.
  do if (nvls > 9).
    compute errcode(errs,1) = 4.
    compute errs = errs+1.
    compute criterr = 1.
  end if.
  do if (criterr = 0).
    compute dumok = 1.
    compute nnvls=make(nvls,1,0).
    compute nnvls(1,1)=dd(1,2).
    compute temp = 2.
    loop i = 2 to n.
      do if (dd(i,2) <> nnvls((temp-1),1)).
         compute nnvls(temp,1)=dd(i,2).
         compute temp = temp+1.
      end if.
    end loop.
    do if (!method > 0).
      compute x = dummy(:,2:ncol(dummy)).
      compute nx = ncol(x).
      compute minus1 = make(1,ncol(x),-1).
      compute xdes=make((nx+1),3,0).
      compute xdes(1,1)=dd(1,2).
      compute xdes(1,2)=1.
      compute temp = 2.
      loop k = 2 to n.
        do if (dd(k,2) <> dd((k-1),2)).
          compute xdes(temp,2) = k.
          compute xdes(temp,1) = dd(k,2).
          compute xdes((temp-1),3) = k-1.
          compute temp=temp+1.
        end if.
      end loop.
      compute xdes((temp-1),3)=n.
      compute xdes = {xdes, (xdes(:,3)-xdes(:,2)+1)}.
      do if (!method = 4).
        loop k = 1 to n.
          do if (rsum(x(k,:)) = 0).
            compute x(k,:) = minus1.
          end if.
        end loop.
      end if.
      do if (!method = 2 or !method = 3 or !method=5).
         loop k = 1 to n.
           do if (rsum(x(k,:)) > 0).
               loop i = 1 to ncol(x).
                   do if (x(k,i) = 0).
                       compute x(k,i) = 1.
                   else.
                        break.
                   end if.
                end loop.
           end if.
         end loop.
         do if (!method = 3).
           compute conmat1={-8,1,1,1,1,1,1,1,1;
                                              0,-7,1,1,1,1,1,1,1;
                                              0,0,-6,1,1,1,1,1,1;
                                              0,0,0,-5,1,1,1,1,1;
                                              0,0,0,0,-4,1,1,1,1;
                                              0,0,0,0,0,-3,1,1,1;
                                              0,0,0,0,0,0,-2,1,1;
                                              0,0,0,0,0,0,0,-1,1}.
           loop i = 1 to 8.
             compute conmat1(i,:)=conmat1(i,:)/(10-i).
           end loop.
           compute conmat1=t(conmat1((10-nvls):8,(10-nvls):9)).
           loop k=1 to n.
              compute x(k,:)=conmat1((rsum(x(k,:))+1),:). 
           end loop.
         end if.
      end if.
      do if (!method = 5).
        compute custcode={!custcode}.
        do if (ncol(custcode) <> (nvls*(nvls-1))).
          compute errcode(errs,1) = (37+!custcodv).
          compute errs = errs + 1.
          compute criterr = 1.
        end if.
        do if (ncol(custcode) = (nvls*(nvls-1))).
          compute conmat1=make(nvls,(nvls-1),0).
          compute cnt=1.
          loop i = 1 to nvls.
            loop k = 1 to (nvls-1).
              compute conmat1(i,k)=custcode(1,cnt).
              compute cnt=cnt+1.
            end loop.
          end loop.
          loop k=1 to n.
            compute x(k,:)=conmat1((rsum(x(k,:))+1),:). 
          end loop.      
        end if.
      end if.
      compute xskip = 1.
      compute dummat = make((nx+1),nx,0).
      compute dummat((2:nrow(dummat)),:)=ident(nx).
      do if (!method = 4).
         compute dummat(1,:) = minus1.
      end if.
      do if (!method = 2).
         loop i = 2 to nrow(dummat).
            loop j = 1 to (i-1).
               compute dummat(i,j) = 1.
            end loop.
          end loop.
      end if.
      do if (!method = 3).
        compute dummat=conmat1.
      end if.
      do if (!method = 5 and criterr=0).
        compute dummat=conmat1.
      end if.
      compute dummat={nnvls, dummat}.
      compute x={dd(:,1),x}.
      compute temp = x.
      compute temp(GRADE(x(:,1)),:) = x.
      compute x = temp.
      release conmat1,temp,dd,xskip,xdes,dummy.
    end if.
  end if.
!ENDDEFINE.


DEFINE !probe3 (decpnt=!charend('/') !default(F10.4)).
do if (criterr=0).
  compute threeway=0.
  compute didprint=0.
  compute didsome=0.
  compute sigintct=0.
  /* could this loop be changed to <= i? */.
  /* start R */.
  loop jmed =1 to (nms+1).  
    compute hasw=0.
    compute hasz=0.
    compute jnok=0.
    compute nm1vls=0.
    compute nm2vls=0.
    compute panelgrp=0.
    compute graphixs={"WITH", outnames(1,i), "BY"}.
    compute focpred4={" "}.
    compute intprint=0.
    compute modcat=0.
   /* start A */.
    do if (jmed <= i).
      do if ((jmed = 1) and ((i+1) = nrow(bcmat))).
        compute pathscnt=pathscnt+1.
      else.
        compute paths={paths,bcmat((i+1),jmed)}.
        compute pathsw={pathsw,wcmat((i+1),jmed)}.
        compute pathsz={pathsz,zcmat((i+1),jmed)}.
        compute pathswz={pathswz,wzcmat((i+1),jmed)}.
        compute temp=fochigh(:,pathscnt)&*bootloc(:,i).
        compute pathsfoc={pathsfoc,pathsfoc(:,1)}.
        do if (jmed=1).
          compute pathtype={pathtype,1}.
        end if.
        do if ((i+1)=nrow(bcmat)).
          compute pathtype={pathtype,3}.
        end if.
        do if (jmed > 1) and ((i+1) < nrow(bcmat)).
          compute pathtype={pathtype,2}.
        end if.
        do if (jmed=1 and nxvls > 1 and (bcmat((i+1),jmed)=1)).
          compute pathsfoc(:,(pathscn2+1))=temp(2:(nxvls+1),1).
        end if.
        do if ((jmed > 1) or (jmed=1 and nxvls=1)).
          compute temp=cmax(temp).
          compute pathsfoc(1,(pathscn2+1))=temp.
        end if.
        compute pathscnt=pathscnt+1.
        compute pathscn2=pathscn2+1.
        do if (i <= nms).
          compute pathsdv={pathsdv,mnames(1,i)}.
        end if.
        do if (i > nms).
          compute pathsdv={pathsdv,ynames}.
        end if.
      end if.
      compute coeffcol=coeffcol+1.
      compute probettt=coeffs(1:nrow(b),coeffcol).
      do if (jmed=1 and (bcmat((i+1),jmed)=1)).
        compute omni=make(nrow(probettt),nxvls,0).
        compute omnitmp=ident(nxvls).
        compute omni(2:(1+nxvls),:)=omnitmp.
      end if.
      do if (csum(probettt)>0).
        compute probvarb=make(csum(probettt),csum(probettt),999).
        compute probcoef=make(csum(probettt),1,999).
        compute coefflp2=1.
        loop coefflp=1 to nrow(probettt).
          do if (probettt(coefflp,1)=1).
            compute probcoef(coefflp2,1)=b(coefflp,1).
            compute coefflp2=coefflp2+1.
          end if.
        end loop.
        compute coefflp=0.
        compute coefflp2=0.
        loop iclp=1 to nrow(probettt).
          do if probettt(iclp,1)=1.
            compute coefflp=coefflp+1.
            compute coefflp2=coefflp.
            compute probvarb(coefflp,coefflp) = varb(iclp,iclp).
            do if (iclp < nrow(probettt)).
              loop jclp=(iclp+1) to nrow(probettt).       
                do if (probettt(jclp,1)=1).
                  compute coefflp2=coefflp2+1.
                  compute probvarb(coefflp,coefflp2)=varb(iclp, jclp).
                  compute probvarb(coefflp2,coefflp)=varb(iclp, jclp).
                end if.
              end loop.
            end if.
          end if.
        end loop. 
      end if.
    end if.
    /* end A */.

    compute xprobval=xmodvals. 
    do if ((nxvls > 1 or mcx > 0) and nspl < 1).
      compute xprobval=dummatx(:,2:ncol(dummatx)).
    end if.
    /* start B */.
    do if ((wcmat((i+1),jmed)=1) and (zcmat((i+1),jmed)=0)).
      compute numplps=1.
      compute modvals=wmodvals.
      compute probeval=wmodvals. 
      compute wheremv1=wherexw.
      compute nm1vls=nwvls.
      compute lpstsp={1,1}.
      compute modcat=0.
      compute jnmod=wtmp.
      compute jnmodlab=wnames.
      compute jnok=1.
      compute jnmin=wmin.
      compute jnmax=wmax.
      compute wherejn1=2.
      do if (jmed=1).
        compute wherejn3=wherexw(1,i).
        do if (nxvls > 1).
          compute jnok=0.
        end if.
      end if.
      do if (jmed > 1).
        compute wherejn1=wherem((jmed-1),i).
        compute wherejn3=wheremw(((2*jmed)-3),i).
      end if.
      do if (nwvls > 1).
        /* compute probeval=dummatw(:,2:ncol(dummatw)) */.
        compute probeval=wprobval.
        compute lpstsp(1,2)=ncol(probeval).
        compute modcat=1.
        compute jnok=0.
      end if.
      do if (wdich = 1).
        compute modcat=1.
        compute jnok=0.
      end if.
      compute problabs=wnames.
      compute focpred3={wnames,"(W)"}. 
      do if (xmint=1).
        compute focpred3={wnames,"(X)"}.       
      end if.
      compute hasw=1.
      compute modgrph=wnames.
      compute intprint=1.
      compute sigintct=sigintct+1.
      compute printpbe=intpb(sigintct,1).
    end if.
    /* end B */.
    /* start C */.
    do if ((wcmat((i+1),jmed)=0) and (zcmat((i+1),jmed)=1)).
      compute numplps=1.
      compute modvals=zmodvals.
      compute probeval=zmodvals.
      compute wheremv1=wherexz.
      compute nm1vls=nzvls.
      compute lpstsp={1,1}.
      compute jnok=1.
      compute jnmod=ztmp.
      compute jnmin=zmin.
      compute jnmax=zmax.
      compute jnmodlab=znames.
      compute wherejn1=2.
      do if (jmed=1).
        compute wherejn3=wherexz(1,i).
        do if (nxvls > 1).
          compute jnok=0.
        end if.
      end if.
      do if (jmed > 1).
        compute wherejn1=wherem((jmed-1),i).
        compute wherejn3=wheremz(((2*jmed)-3),i).
      end if.
      do if (nzvls > 1).
        /* compute probeval=dummatz(:,2:ncol(dummatz)) */.
        compute probeval=zprobval.
        compute lpstsp(1,2)=ncol(probeval).
        compute modcat=1.
        compute jnok=0.
      end if.
      do if (zdich = 1).
        compute modcat=1.
        compute jnok=0.
      end if.
      compute problabs=znames.
      compute focpred3={znames,"(Z)"}. 
      compute modgrph=znames.
      compute hasz=1.
      compute intprint=1.
      compute sigintct=sigintct+1.
      compute printpbe=intpb(sigintct,1).
    end if.
    /* end C */.
    /* start D */.
    do if ((wzcmat((i+1),jmed)=1) or ((wcmat((i+1),jmed)=1) and (zcmat((i+1),jmed)=1))).
      compute numplps=2.
      compute probecnt=1.
      compute intprint=1.
      do if (wzcmat((i+1),jmed)=1).
        compute sigintct=sigintct+1.
        compute printpbe=intpb(sigintct,1).
      else.
        compute sigintct=sigintct+2.    
        compute printpbe=cmin(intpb((sigintct-1):sigintct,1)).
      end if.
      compute panelgrp=1.
      compute hasw=1.
      compute hasz=1.
      compute panelcde={"/PANEL","ROWVAR=",znames,"."}.
      compute modgrph=wnames.
      compute lpstsp={1,1;1,1}.
      compute wheremv1=wherexw.
      compute nm1vls=nwvls.
      compute wheremv2=wherexz.
      compute nm2vls=nzvls.
      /* compute wprobval=wmodvals */.
      compute jnok=0.
      do if (wzcmat((i+1),jmed)=1).
        compute jnok=1.
      end if.
      do if (jmed > 1).
        compute mprobval=mmodvals.
      end if.
      do if (jmed=1).
        do if (nxvls > 1).
          compute jnok=0.
        end if.
      end if.
      do if (nwvls > 1).
        /* compute wprobval=dummatw(:,2:ncol(dummatw)) */.
        compute lpstsp(1,2)=ncol(wprobval).
        compute modcat=1.
        compute jnok=0.
      end if.
      do if (zdich=1).
        compute modcat=1.
        compute jnok=0.
      end if.
      /* compute zprobval=zmodvals */.
      compute lpstsp(2,1)=lpstsp(1,2)+1.
      compute lpstsp(2,2)=lpstsp(1,2)+1.
      do if (nzvls > 1).
        /* compute zprobval=dummatz(:,2:ncol(dummatz)) */.
        compute lpstsp(2,1)=lpstsp(1,2)+1.
        compute lpstsp(2,2)=lpstsp(1,2)+ncol(zprobval).
        compute jnok=0.
      end if.
      do if (zdich=1).
        compute jnok=0.
      end if.
      compute omni3=make(nrow(b),(nxvls*nwvls),0).
      do if (jmed > 1).
        compute omni3=make(nrow(b),nwvls,0).
      end if.
      compute focpred3={wnames,"(W)"}.
      do if (xmint=1).
        compute focpred3={wnames,"(X)"}.      
      end if.
      compute focpred4={"    ", "Mod var:", znames, "(Z)"}. 
      compute modvals=make((nrow(wmodvals)*nrow(zmodvals)),2,0).
      compute probeval=make((nrow(wmodvals)*nrow(zmodvals)),(ncol(wprobval)+ncol(zprobval)),0).
      loop probei= 1 to nrow(wmodvals).
        loop probej =1 to nrow(zmodvals).
          compute modvals(probecnt,1)=wmodvals(probei,1).
          compute probeval(probecnt,1:nwvls)=wprobval(probei,:).
          compute modvals(probecnt,2)=zmodvals(probej,1).
          compute probeval(probecnt,(nwvls+1):(nwvls+nzvls))=zprobval(probej,:).
          compute probecnt=probecnt+1.
        end loop.
      end loop.
      do if (wzcmat((i+1),jmed)=1).
        compute numplps=numplps+1.
        compute probprod=make(1,(ncol(wprobval)*ncol(zprobval)),0).
        compute lpstsp2={1,1}.
        compute lpstsp={lpstsp;lpstsp2}.
        compute lpstsp(3,1)=lpstsp(2,2)+1.
        compute lpstsp(3,2)=lpstsp(2,2)+ncol(probprod).
        compute jnmod=ztmp.
        compute jnmin=zmin.
        compute jnmax=zmax.
        compute jnmodlab=znames.
        do if (jmed = 1).
          compute wherejn1=wherexw(1,i).
          compute wherejn3=wherexwz(1,i).
        end if.
        do if (jmed > 1).
          /* compute wherejn1=wheremw((jmed-1),i) */.
          compute wherejn1=wheremw(((2*jmed)-3),i).
          compute wherejn3=wheremwz(((2*jmed)-3),i).
        end if.
        loop probei = 1 to nrow(wmodvals).
          loop probej = 1 to nrow(zmodvals).
            compute probtemp=1.
            loop probek = 1 to ncol(wprobval).
              compute probtemp={probtemp,(wprobval(probei,probek)&*zprobval(probej,:))}.
            end loop.
            compute probprod={probprod;probtemp(1,2:ncol(probtemp))}.
          end loop.
        end loop.
        compute probprod=probprod(2:nrow(probprod),:).
        compute probeval={probeval,probprod}.
      end if.
      compute problabs={wnames,znames}.
      release probecnt, probei, probej.
    end if.
    /* end D */.
    /* Start E */.
    do if (intprint=1).
      compute focpred={"   Focal", "predict:"}.
      do if (jmed=1).
        compute focpred={focpred,xnames,"(X)"}.
        compute focplotv=xmodvals.
      end if.
      do if (jmed >1).
        do if (nms > 1).
          compute focpred={focpred,mnames(1,(jmed-1)), medlb2(1,(jmed-1))}.
        end if.
        do if (nms = 1).
          compute focpred={focpred,mnames(1,(jmed-1)), "(M)"}.
        end if.
        compute focplotv=mmodvals(:,(jmed-1)).
      end if.
      compute focpred2={"    ", "Mod var:",focpred3}.
      compute focpred={focpred;focpred2}.
      do if (ncol(focpred4) > 1).
        compute focpred={focpred;focpred4}.
        compute focpred4={" "}.
      end if.
      release focpred2,focpred3.
      do if (outscree=1).
        do if ((plot = 1 or plot = 2) or (printpbe <= intprobe)).
          print focpred/title="----------"/format=A8/space=0.
        end if.
      end if.
      compute foctmp=make(nrow(modvals),1,1).
      compute probexpl=1.
      compute probeva2={foctmp,probeval}.
      do if (jmed=1 and nxs > 0 and mcx > 0).
        compute probexpl=nxvls.
      end if. 
      compute foctmp=make(nrow(modvals),1,1).
      compute modvals3=make(1,(6+ncol(problabs)),0).
      compute probrown=make(nrow(probeval),1,0).
      compute jtmp=1.
      loop probei = 1 to nrow(probeval).
        compute probrown(probei,1)=jtmp.
        compute jtmp=jtmp+nxvls.
      end loop.
      release jtmp.
      compute probrow=999.
      compute modvarl=problabs.
       /* start F */. 
      do if (plot = 1 or plot = 2 or nxvls > 1).
        compute plotvals=make((nrow(modvals)*nrow(focplotv)),(ncol(modvals)+1),999).
        loop ploti=1 to nrow(modvals).
          loop plotj=1 to nrow(focplotv).
            compute plotvals((((ploti-1)*nrow(focplotv))+plotj),2:ncol(plotvals))=modvals(ploti,:).
            compute plotvals((((ploti-1)*nrow(focplotv))+plotj),1)=focplotv(plotj,1).
          end loop.
        end loop.
        compute focpredn=3.
        do if (jmed=1).
          do if (nxvls > 1).
            compute focpredn=(nxvls+1).
          end if.
          do if (nxvls=1 and xdich=1).
            compute focpredn=2.
          end if.
        end if.
        compute meanmat=mdiag(means).
        compute onesmat=make(nrow(meanmat),(nrow(probeval)*focpredn),1).
        compute probeplt=t(mdiag(means)*onesmat).
        /* start G */.
        do if (jmed=1).
          do if ((wcmat((i+1),1)=1) or (zcmat((i+1),1)=1)).
            compute plotcnt=1.
            compute iloops=nwpval*nzpval.
            compute plotmx=nxpval*nzpval.
            do if ((wcmat((i+1),1)=1) and (zcmat((i+1),1)=0)).
              compute iloops=nwpval.
              compute plotmx=nxpval.
            end if.
            do if ((wcmat((i+1),1)=0) and (zcmat((i+1),1)=1)).
              compute iloops=nzpval.
              compute plotmx=nxpval.
            end if.
            compute xestvals=make((nxpval*iloops),ncol(xprobval),-999).
            do if (wcmat((i+1),1)=1).
              compute westvals=make(nrow(xestvals),ncol(wprobval),-999).
            end if.
            do if (zcmat((i+1),1)=1).
              compute zestvals=make(nrow(xestvals),ncol(zprobval),-999).
            end if.
            loop ploti=1 to iloops.
              loop plotj=1 to nxpval.
                compute xestvals(plotcnt,:)=xprobval(plotj,:).
                compute plotcnt=plotcnt+1.
              end loop.
            end loop.
            compute plotcnt=1.
            compute plotcnt1=1.
            compute plotcnt2=1.
            compute plotcntz=1.
            loop ploti = 1 to (iloops*nxpval).
              do if (wcmat((i+1),1)=1).
                compute westvals(ploti,:)=wprobval(plotcnt1,:).
              end if.        
              do if ((wcmat((i+1),1)=0) and (zcmat((i+1),1)=1)).
                compute zestvals(ploti,:)=zprobval(plotcnt1,:).
              end if.   
              do if ((wcmat((i+1),1)=1) and (zcmat((i+1),1)=1)).
                compute zestvals(ploti,:)=zprobval(plotcnt2,:).
                compute plotcntz=plotcntz+1.
              end if.
              compute plotcnt=plotcnt+1.
              do if (plotcnt > plotmx).
                compute plotcnt=1.
                compute plotcnt1=plotcnt1+1.
              end if.
              do if (plotcntz > nxpval).
                compute plotcnt2=plotcnt2+1.
                compute plotcntz=1.
                do if (plotcnt2 > nzpval).
                  compute plotcnt2=1. 
                end if.
              end if.
            end loop.
            compute probeplt(:,2:(1+(ncol(xestvals))))=xestvals. 
            do if (wcmat((i+1),1)=1).
              compute probeplt(:,(wherew(1,i)):(wherew(2,i)))=westvals.
            end if.    
            do if (zcmat((i+1),1)=1).
              compute probeplt(:,(wherez(1,i)):(wherez(2,i)))=zestvals.
            end if.  
          end if.
        end if.
        /* end G */.
        /* start H */.
        do if (jmed > 1).
          do if ((wcmat((i+1),jmed)=1) or (zcmat((i+1),jmed)=1)).
            compute plotcnt=1.
            compute iloops=nwpval*nzpval.
            compute plotmx=3*nzpval.
            do if ((wcmat((i+1),jmed)=1) and (zcmat((i+1),jmed)=0)).
              compute iloops=nwpval.
              compute plotmx=3.
             end if.
             do if ((wcmat((i+1),jmed)=0) and (zcmat((i+1),jmed)=1)).
               compute iloops=nzpval.
               compute plotmx=3.
             end if.
             compute mestvals=make((3*iloops),1,-999).
             do if (wcmat((i+1),jmed)=1).
               compute westvals=make(nrow(mestvals),ncol(wprobval),-999).
             end if.
             do if (zcmat((i+1),jmed)=1).
               compute zestvals=make(nrow(mestvals),ncol(zprobval),-999).
             end if.
             loop ploti=1 to iloops.
               loop plotj=1 to 3.
                 compute mestvals(plotcnt,:)=mprobval(plotj,(jmed-1)).
                 compute plotcnt=plotcnt+1.
               end loop.
             end loop.
             compute plotcnt=1.
             compute plotcnt1=1.
             compute plotcnt2=1.
             compute plotcntz=1.
             loop ploti = 1 to (iloops*3).
               do if (wcmat((i+1),jmed)=1).
                 compute westvals(ploti,:)=wprobval(plotcnt1,:).
               end if.        
               do if ((wcmat((i+1),jmed)=0) and (zcmat((i+1),jmed)=1)).
                 compute zestvals(ploti,:)=zprobval(plotcnt1,:).
               end if.   
               do if ((wcmat((i+1),jmed)=1) and (zcmat((i+1),jmed)=1)).
                 compute zestvals(ploti,:)=zprobval(plotcnt2,:).
                 compute plotcntz=plotcntz+1.
               end if.
               compute plotcnt=plotcnt+1.
               do if (plotcnt > plotmx).
                 compute plotcnt=1.
                 compute plotcnt1=plotcnt1+1.
               end if.
               do if (plotcntz > 3).
                 compute plotcnt2=plotcnt2+1.
                 compute plotcntz=1.
                 do if (plotcnt2 > nzpval).
                   compute plotcnt2=1. 
                 end if.
               end if.
             end loop.
             compute probeplt(:,wherem((jmed-1),i))=mestvals.
             do if (wcmat((i+1),jmed)=1).
               do if (model <> 74).
                 compute probeplt(:,(wherew(1,i)):(wherew(2,i)))=westvals.
               end if.
               do if (model = 74).
                 compute probeplt(:,(wherex(1,i)):(wherex(2,i)))=westvals.
               end if.
             end if.    
             do if (zcmat((i+1),jmed)=1).
               compute probeplt(:,(wherez(1,i)):(wherez(2,i)))=zestvals.
             end if.  
          end if.
        end if.
        /* end H */.
        /* Here I am doing the multiplications to produce data for the plot */.
        compute prodloop = 1.
        do if (jmed=1).
          compute prodloop=ncol(xestvals).
        end if.
        do if (wcmat((i+1),jmed))=1.
          compute plotcnt=0.
          loop ploti = 1 to prodloop.
            loop plotj = 1 to ncol(westvals).
              do if (jmed=1).
                compute probeplt(:,(wherexw(1,i)+plotcnt))=xestvals(:,ploti)&*westvals(:,plotj).
              end if.
              do if (jmed > 1).
                compute probeplt(:,(wheremw(((jmed*2)-3) ,i)+plotcnt))=mestvals(:,ploti)&*westvals(:,plotj).
              end if.
              compute plotcnt=plotcnt+1.
            end loop.
          end loop. 
        end if.
        do if (zcmat((i+1),jmed))=1.
          compute plotcnt=0.
          loop ploti = 1 to prodloop.
            loop plotj = 1 to ncol(zestvals).   
              do if (jmed = 1).
                compute probeplt(:,(wherexz(1,i)+plotcnt))=xestvals(:,ploti)&*zestvals(:,plotj).
              end if.
              do if (jmed > 1).
                compute probeplt(:,(wheremz(((jmed*2)-3),i)+plotcnt))=mestvals(:,ploti)&*zestvals(:,plotj).
              end if.
              compute plotcnt=plotcnt+1.
            end loop.
          end loop. 
        end if.
        do if (wzcmat((i+1),jmed))=1.
          compute plotcnt=0.
          compute threeway=1.
          loop ploti = 1 to ncol(westvals).
            loop plotj = 1 to ncol(zestvals).
              compute probeplt(:,(wherewz(1,i)+plotcnt))=westvals(:,ploti)&*zestvals(:,plotj).
              compute plotcnt=plotcnt+1.
            end loop.
          end loop. 
          compute plotcnt=0.
          loop plotk = 1 to prodloop.
            loop ploti = 1 to ncol(westvals).
              loop plotj = 1 to ncol(zestvals).
                do if (jmed = 1).
                  compute probeplt(:,(wherexwz(1,i)+plotcnt))=xestvals(:,plotk)&*westvals(:,ploti)&*zestvals(:,plotj).
                end if.
                do if (jmed > 1).
                  compute probeplt(:,(wheremwz(((jmed*2)-3),i)+plotcnt))=mestvals(:,plotk)&*westvals(:,ploti)&*zestvals(:,plotj).
                end if.
                compute plotcnt=plotcnt+1.
              end loop.
            end loop. 
          end loop.
        end if.
        /* here is where we add holding constant products as needed */.
        /* start I */.
        loop newplp=1 to i.
          do if (newplp <> jmed).
            do if (wcmat((i+1),newplp))=1.
              compute prodloop=1.
              do if (newplp=1).
                compute prodloop=nxvls.
              end if.
              compute plotcnt=0.
              loop ploti = 1 to prodloop.
                loop plotj = 1 to nwvls.
                  do if (newplp = 1).
                    compute probeplt(:,(wherexw(1,i)+plotcnt))=probeplt(:,(1+ploti))&*probeplt(:,(wherew(1,i)+plotj-1)).
                  end if.
                  do if (newplp > 1).
                    do if (model <> 74).
                      compute probeplt(:,(wheremw(((newplp*2)-3) ,i)+plotcnt))=probeplt(:,wherem((newplp-1),i))&*probeplt(:,(wherew(1,i)+plotj-1)).
                    end if.
                    do if (model = 74).
                      compute probeplt(:,(wheremw(((newplp*2)-3) ,i)+plotcnt))=probeplt(:,wherem((newplp-1),i))&*probeplt(:,(wherex(1,i)+plotj-1)).
                    end if.
                  end if.
                  compute plotcnt=plotcnt+1.
                end loop.
              end loop. 
            end if.  
            do if (zcmat((i+1),newplp))=1.
              compute prodloop=1.
              do if (newplp=1).
                compute prodloop=nxvls.
              end if.
              compute plotcnt=0.
              loop ploti = 1 to prodloop.
                loop plotj = 1 to nzvls.
                  do if (newplp = 1).
                    compute probeplt(:,(wherexz(1,i)+plotcnt))=probeplt(:,(1+ploti))&*probeplt(:,(wherez(1,i)+plotj-1)).
                  end if.
                  do if (newplp > 1).
                    compute probeplt(:,(wheremz(((newplp*2)-3),i)+plotcnt))=probeplt(:,wherem((newplp-1),i))&*probeplt(:,(wherez(1,i)+plotj-1)).
                  end if.
                  compute plotcnt=plotcnt+1.
                end loop.
              end loop. 
            end if.
            do if (wzcmat((i+1),newplp))=1.
              compute plotcnt=0.
              do if (threeway=0).
                loop ploti = 1 to nwvls.
                  loop plotj = 1 to nzvls.
                    compute probeplt(:,(wherewz(1,i)+plotcnt))=probeplt(:,(wherew(1,i)+ploti-1))&*probeplt(:,(wherez(1,i)+plotj-1)).
                    compute plotcnt=plotcnt+1.
                  end loop.
                end loop.
              end if. 
              compute prodloop=1.
              do if (newplp=1).
                compute prodloop=nxvls.
              end if.
              compute plotcnt=0.
              loop plotk = 1 to prodloop.
                loop ploti = 1 to nwvls.
                  loop plotj = 1 to nzvls.
                    do if (newplp = 1).
                      compute probeplt(:,(wherexwz(1,i)+plotcnt))=probeplt(:,(1+plotk))&*probeplt(:,(wherew(1,i)+ploti-1))&*probeplt(:,(wherez(1,i)+plotj-1)).
                    end if.
                    do if (newplp > 1).
                      compute probeplt(:,(wheremwz(((newplp*2)-3),i)+plotcnt))=probeplt(:,wherem((newplp-1),i))&*probeplt(:,(wherew(1,i)+ploti-1))&*probeplt(:,(wherez(1,i)+plotj-1)).
                    end if.
                    compute plotcnt=plotcnt+1.
                  end loop.
                end loop. 
              end loop.
            end if.      
          end if.
        end loop.
        /* END I */.
        compute predvals=probeplt*b.
        do if (debug <> 0).
          print probeplt.
          print b.
        end if.
        
        do if (i = nms+nys) and (ydich=1).
          compute predvalt=(predvals < 709.7).
          compute prevalt7=(1-predvalt)*(709.7).
          compute predvals=(predvals&*predvalt)+prevalt7.
          compute expyhat=exp(predvals)&/(1+exp(predvals)).
        end if. 
        compute sepred=make(nrow(plotvals),3,999).
        loop sei=1 to nrow(plotvals).
          compute ask=probeplt(sei,:).       
          compute sepred(sei,1)=sqrt(ask*varb*t(ask)).
          do if ((i < nms+nys) or (ydich=0)).
              compute sepred(sei,2)=predvals(sei,1)-tval*sepred(sei,1).
              compute sepred(sei,3)=predvals(sei,1)+tval*sepred(sei,1).
          end if.
          do if ((i = nms+nys) and (ydich=1)).
            compute sepred(sei,2)=predvals(sei,1)-xp2*sepred(sei,1).
            compute sepred(sei,3)=predvals(sei,1)+xp2*sepred(sei,1).
          end if.
        end loop.
        compute prevloc=ncol(plotvals)+1.
        compute probeplt={plotvals,predvals}.
        do if (plot = 2).
          compute probeplt={probeplt,sepred}.
        end if.
        do if ((i = nms+nys) and (ydich=1)).
          compute probeplt={probeplt,expyhat}.
        end if.
        compute didsome=0.
      end if.
      /* END F */.
      /* here is the loop that is printing the conditional effects */.
      /* this does conditional two way interactions */.
      /* start J */.
      do if ((wzcmat((i+1),jmed)=1) and (printpbe <= intprobe)).
        do if (jmed=1).
          compute omnilp2=nxvls*nwvls.
          compute omnitmp=ident(omnilp2).
          compute omni3(wherexw(1,i):wherexw(2,i),:)=omnitmp.
        end if.
        do if (jmed>1).
          compute omnilp2=nwvls.
          compute omnitmp=ident(omnilp2).
          compute omni3(wheremw(((jmed*2)-3),i):wheremw(((jmed*2)-2),i),:)=omnitmp.
        end if.
        compute omnif=make(1,4,0).
        do if ((i = nms+nys) and (ydich=1)).
          compute omnif=make(1,3,0).
        end if.
        compute condeff3=0.
        loop omnilp1=1 to nrow(zprobval).
          loop omnilp=1 to (omnilp2).
            do if (jmed=1).
              compute omni3((wherexwz(1,i)+((omnilp-1)*nzvls)):(wherexwz(1,i)+((omnilp-1)*nzvls)+(nzvls-1)),omnilp)=t(zprobval(omnilp1,:)).
            end if.
            do if (jmed > 1).
              compute omni3((wheremwz(((jmed*2)-3),i)+((omnilp-1)*nzvls)):(wheremwz(((jmed*2)-3),i)+((omnilp-1)*nzvls)+(nzvls-1)),omnilp)=t(zprobval(omnilp1,:)).
            end if.
          end loop.
          compute condeff=t(omni3)*b.
          compute condeff3={condeff3;condeff}.
          !ftest3 lm=omni3/bcoef=b/cv=varb/skip=1.
         compute omnif={omnif;fresult}.
        end loop.
        compute omnif=omnif(2:nrow(omnif),:).
        compute clabtmp=znames.
        compute condeff3=condeff3(2:nrow(condeff3),:).
        do if ((nxvls*nwvls)=1).
          compute omnif={condeff3,omnif}.
          compute clabtmp={clabtmp,"Effect"}.
        end if. 
        compute omnif={zmodvals,omnif}.
        do if ((i < nms+nys) or (ydich=0)).
          compute clabtmp={clabtmp,hcflab,"df1","df2","p"}.
        end if.
        do if ((i = (nms+nys)) and (ydich=1)).
          compute clabtmp={clabtmp,"Chi-sq","df","p"}.
        end if.
        !outform outtodo=omnif/outbig=maxresm.
        do if (outscree=1).
          do if (jmed=1).
            print omnif/title="Test of conditional X*W interaction at value(s) of Z:"/cnames=clabtmp/format=!decpnt.
          end if.
          do if (jmed>1).
            print omnif/title="Test of conditional M*W interaction at value(s) of Z:"/cnames=clabtmp/format=!decpnt.
          end if.
        end if.
        release omni3.
      end if.
      /* end J */.
      /* start O */.
      loop probei = 1 to probexpl.
        do if (probexpl > 1).
          compute foctmp=make(nrow(modvals),probexpl,0).
          compute foctmp(:,probei)=foctmp(:,probei)+1.
          compute probtemp=make(nrow(modvals),1,0).
          loop probem = 1 to numplps.
            loop probek = 1 to nxvls.
              loop probej=lpstsp(probem,1) to lpstsp(probem,2).
                compute probtemp={probtemp,foctmp(:,probek)&*probeval(:,probej)}.
              end loop.
            end loop.
          end loop.
          compute probeva2=probtemp(:,2:ncol(probtemp)).
          compute probeva2={foctmp,probeva2}.
        end if.            
        compute probres=probeva2*probcoef.
        compute probrese=sqrt(diag(probeva2*probvarb*t(probeva2))). 
        compute tratio = probres&/probrese.
        compute p = 2*(1-tcdf(abs(tratio), dfres)).
        do if (ydich=1 and i = (nms+nys)).
          compute p = 2*(1-cdfnorm(abs(tratio))).
        end if.
        compute modvals2={modvals,probres,probrese,tratio, p}.
        do if ((i < nms+nys) or (ydich=0)).
          compute modvals2={modvals2,(probres-tval&*probrese),(probres+tval&*probrese)}.   
          compute problabs={problabs,"Effect",hclab,"t", "p", "LLCI", "ULCI"}.
        end if.
        do if (ydich=1 and i = (nms+nys)).
          compute modvals2={modvals2,(probres-xp2&*probrese),(probres+xp2&*probrese)}.
          compute problabs={problabs,"Effect","se","Z", "p", "LLCI", "ULCI"}.
        end if.
        /* start L */.
        do if (probexpl > 1 and (printpbe <= intprobe)).
          do if (hasz = 1).
            compute printz=1.
          end if.
          do if (hasw=1).
            compute printw=1.
          end if.
          compute probrlab=make(nrow(modvals),1,xcatlab(probei,1)).
          compute modvals3={modvals3; modvals2}.
          compute probrow={probrow;probrown}.
          compute probrown=probrown+1.
          do if (probei=probexpl).
            compute xproblab=xcatlab(1:nxvls,1).
            compute probrow=probrow(2:nrow(probrow),1).
            compute modvals3=modvals3(2:nrow(modvals3),:).
            compute temp=modvals3.
            compute temp(GRADE(probrow(:,1)),:)=modvals3.
            compute modvals3=temp.
            compute start2=1.
            compute problabs=problabs(1,(1+(ncol(modvarl))):ncol(problabs)).
            compute pstart=1.
            /* start K */.
            loop probek= 1 to nrow(probeval).
              compute endstart=start2+(nxvls-1).
              compute temp=modvals3(start2:endstart,(1+ncol(modvarl)):ncol(modvals3)).
              compute temp2=t(modvals3(start2:start2,1:ncol(modvarl))).
              compute trnames=t(modvarl).
              do if (outscree=1).
                do if (probek > 1).
                  print/title="----------"/space=0.
                else.
                  print/title = "Conditional effects of the focal predictor at values of the moderator(s):".
                  do if ((jmed=1) and (i = (nms+nys)) and (nms > 0)).
                    do if (nxvls = 1).
                      print/title = "(These are also the conditional direct effects of X on Y)"/space=0.
                    else.
                      print/title = "(These are also the relative conditional direct effects of X on Y)"/space=0.
                    end if.
                  end if.
                  print.
                end if.
                print temp2/title = "Moderator value(s):"/rnames=trnames/format=!decpnt/space=0.
                print temp/title = " "/cnames=problabs/rnames=xproblab/format=!decpnt/space=0.
              end if.
              !outform outtodo=temp2/outbig=maxresm.
              !outform outtodo=temp/outbig=maxresm.
              compute start2=start2+nxvls.
              compute didsome=1.
              do if (jmed=1).
                compute mod1val=probeval(probek,1:nm1vls).
                loop omnilp=1 to nxvls.
                  compute omni((wheremv1(1,i)+((omnilp-1)*nm1vls)):(wheremv1(1,i)+((omnilp-1)*nm1vls)+(nm1vls-1)),omnilp)=t(mod1val).
                  do if (nm1vls < ncol(probeval)).
                    compute mod2val=probeval(probek,(nm1vls+1):(nm1vls+nm2vls)).
                    compute omni((wheremv2(1,i)+((omnilp-1)*nm2vls)):(wheremv2(1,i)+((omnilp-1)*nm2vls)+(nm2vls-1)),omnilp)=t(mod2val).
                    do if ((nm1vls+nm2vls) < ncol(probeval)).
                      compute intlen=nm1vls*nm2vls.
                      compute modintvl=probeval(probek,(nm1vls+nm2vls+1):ncol(probeval)).
                      compute omni((wherexwz(1,i)+((omnilp-1)*intlen)):(wherexwz(1,i)+((omnilp-1)*intlen)+(intlen-1)),omnilp)=t(modintvl).
                    end if.
                  end if.
                end loop.
                /* compute condeff=t(omni)*b */. 
               !ftest3 lm=omni/bcoef=b/cv=varb/skip=1.
               !outform outtodo=fresult/outbig=maxresm.
               do if (outscree=1).
                 do if (i < (nms + nys) or (ydich = 0)).
                   compute clabtmp={hcflab,"df1","df2","p"}.
                   print fresult/title="Test of equality of conditional means"/cnames=clabtmp/format=!decpnt.
                 end if.
                 do if (i = (nms + nys) and (ydich = 1)).
                   compute clabtmp={"Chi-sq","df","p"}.
                   print fresult/title="Test of equality of conditional logits or probabilities"/cnames=clabtmp/format=!decpnt.
                 end if.
               end if.
                compute probetmp=probeplt(pstart:(pstart+nxvls),1).
                compute probetmp={probetmp,probeplt(pstart:(pstart+nxvls),prevloc:ncol(probeplt))}.
                compute pstart=pstart+(nxvls+1).
                !outform outtodo=probetmp/outbig=maxresm.
                do if (i < ((nms + nys)) or (ydich = 0)).   
                  compute clabtmp={xnames, outnames(1,i), hclab, "LLCI", "ULCI"}.
                  do if (outscree=1).
                    print probetmp/title = "Estimated conditional means being compared:"/cnames=clabtmp/format=!decpnt.
                  end if.
                end if.
                do if ((i = (nms+nys)) and (ydich = 1)).
                   compute clabtmp={xnames, outnames(1,i), "prob"}.
                  compute probetm2=probetmp(:,1:2).
                  compute probetm2={probetm2,probetmp(:,ncol(probetmp))}.
                  do if (outscree=1).
                    print probetm2/title = "Estimated conditional logits and probabilities:"/cnames=clabtmp/format=!decpnt.
                  end if.
                end if.
              end if.
            end loop.
            /* end K */.
            release probrow, start2, endstart, temp, temp2.
          end if.
        end if.
        /* end L */.
        /* start N */.
        do if (probexpl = 1 and (printpbe <= intprobe)).
          !outform outtodo=modvals2/outbig=maxresm.
          do if (outscree=1).
            print/title = "Conditional effects of the focal predictor at values of the moderator(s):".
            print modvals2/cnames=problabs/title = " "/space=0/format=!decpnt.
          end if.
          compute didsome=1.
          do if (hasz = 1).
            compute printz=1.
          end if.
          do if (hasw=1).
            compute printw=1.
          end if.
          /* start M */.
          do if (jn = 1 and jnok=1).
            do if (criterr = 0).
              /* compute dfres=n-nrow(b) */.
              compute roots=99999.                           
              /* compute jncrit =(dfres* (exp((dfres-(5/6))*((xp2/(dfres-(2/3)+(.11/dfres)))*(xp2/(dfres-(2/3)+(.11/dfres)))))-1)) */.
              compute jncrit=(-idf.t(alpha2,dfres))**2.
             do if (i = (nms+nys) and (ydich=1)).
               compute jncrit=xp2*xp2.
             end if.
              compute jnb1=b(wherejn1,1).
              compute jnb3=b(wherejn3,1).
              compute jnsb1=varb(wherejn1,wherejn1).
              compute jnsb3=varb(wherejn3,wherejn3).
              compute jnsb1b3=varb(wherejn1,wherejn3).
              compute ajn =(jncrit*jnsb3)-(jnb3*jnb3).
              compute bjn = 2*((jncrit*jnsb1b3)-(jnb1*jnb3)).
              compute cjn = (jncrit*jnsb1)-(jnb1*jnb1).
              compute radarg = (bjn*bjn)-(4*ajn*cjn).
              compute den = 2*ajn.
              compute nrts = 0.
              do if (radarg >= 0 and den <> 0).
                compute x21 = (-bjn+sqrt(radarg))/den.
                compute x22 = (-bjn-sqrt(radarg))/den.
                compute roots = 0.
                do if (x21 >= jnmin and x21 <= jnmax).
                  compute nrts = 1.
                  compute roots = {roots; x21}.
                end if.
                do if (x22 >= jnmin and x22 <= jnmax).
                  compute nrts = nrts + 1.
                  compute roots = {roots; x22}.
                end if.
                compute roots={roots,make(nrow(roots),2,0)}.
              end if.
              do if (nrts > 0).                           
                compute roots = roots(2:nrow(roots),1:3).
                compute roots(1,2)=(csum(jnmod < roots(1,1))/n)*100.
                compute roots(1,3)=(csum(jnmod > roots(1,1))/n)*100.
                do if (nrow(roots)=2).
                  compute roots(2,2)=(csum(jnmod < roots(2,1))/n)*100.
                  compute roots(2,3)=(csum(jnmod > roots(2,1))/n)*100.
                end if.
                do if (outscree=1).
                  print roots/title = "Moderator value(s) defining Johnson-Neyman significance region(s):"/clabels = "Value", "% below", "% above"/format !decpnt.
                end if.
                do if (nrts=1).
                  compute tmprts=make(1,3,99999).
                  compute roots={roots;tmprts}.
                end if.                
              end if.
              do if (nrts = 0).
                compute roots=make(2,1,99999).
                do if (outscree=1).
                  print/title = "There are no statistical significance transition points within the observed".
                  print/title = "range of the moderator found using the Johnson-Neyman method."/space=0.
                end if.
              end if.
              !outform outtodo=roots/outbig=maxresm.
              compute jnvals=make(23,7,0).
              loop jni= 0 to (21-(nrts)).
                compute jnvals((jni+1),1)=jnmin+(jni*((jnmax-jnmin)/(21-nrts))).
              end loop.
              do if (nrts > 0).
                loop jni = 1 to nrts.
                  loop jnj = 2 to (nrow(jnvals)-1).
                    do if ((roots(jni,1) > jnvals((jnj-1),1)) and (roots(jni,1) < jnvals(jnj,1))).                      
                      compute jnvals((jnj+1):(21+jni),1)=jnvals(jnj:(20+jni),1).
                      compute jnvals(jnj,1)=roots(jni,1).
                    end if.
                  end loop.
                end loop.
              end if. 
              compute jnvals=jnvals(1:22,:).
              loop jni = 1 to nrow(jnvals).
                compute jnvals(jni,2)=jnb1+jnb3*jnvals(jni,1).
                compute jnvals(jni,3)=sqrt(jnsb1+2*jnvals(jni,1)*jnsb1b3+(jnvals(jni,1)*jnvals(jni,1))*jnsb3).
                compute jnvals(jni,4)=jnvals(jni,2)/jnvals(jni,3).
                compute jnvals(jni,5)=2*(1-tcdf(abs(jnvals(jni,4)), dfres)).
                compute jnvals(jni,6)=jnvals(jni,2)-sqrt(jncrit)*jnvals(jni,3).
                compute jnvals(jni,7)=jnvals(jni,2)+sqrt(jncrit)*jnvals(jni,3).
                do if ((i = nms + nys) and (ydich=1)).
                  compute jnvals(jni,5)=2*(1-cdfnorm(abs(jnvals(jni,4)))).
                  compute jnvals(jni,6)=jnvals(jni,2)-xp2*jnvals(jni,3).
                  compute jnvals(jni,7)=jnvals(jni,2)+xp2*jnvals(jni,3).
                end if.
              end loop. 
              !outform outtodo=jnvals/outbig=maxresm.
              do if ((i < nms+nys) or (ydich=0)).
                compute jnclbs={jnmodlab,"Effect",hclab,"t", "p", "LLCI", "ULCI"}.
              end if.
              do if ((i = nms + nys) and (ydich=1)).
                compute jnclbs={jnmodlab,"Effect","se","Z", "p", "LLCI", "ULCI"}.
              end if.
              do if (((wcmat((i+1),jmed)=1) or (zcmat((i+1),jmed)=1)) and (wzcmat((i+1),jmed)=0)).
                do if (outscree=1).
                  print jnvals/title = "Conditional effect of focal predictor at values of the moderator:"/cnames =jnclbs/format = !decpnt.
                end if.
              end if.
              do if (outscree=1).
                do if ((jmed = 1) and (wzcmat((i+1),jmed)=1)).
                  print jnvals/title = "Conditional X*W interaction at values of the moderator Z:"/cnames =jnclbs/format = !decpnt.
                end if.
                do if ((jmed > 1) and (wzcmat((i+1),jmed)=1)).
                  print jnvals/title = "Conditional M*W interaction at values of the moderator Z:"/cnames =jnclbs/format = !decpnt.
                end if.
              end if.
            end if.
          end if.
          /* END M */.
        end if.
        /* END N */.
        do if ((i = (nms+nys)) and (jmed=1) and (bcmat(nrow(bcmat),1)=1)).
          do if (probei=1).
            compute direfflb=problabs.
            compute direff=modvals2.
          end if.
          do if (probei>1).
            compute direff={direff;modvals2}.
            /* compute direfflb=problabs */.
          end if.  
        end if.
        compute intprint=0.
        /* This does the contrast for conditional effects */.
        do if ((jmed=1) and (i=1) and nms=0) and modcok=1).
          compute contvec2=make(2,1,1).
          compute contvec2={contvec2,wcontval,zcontval}.
          do if (wzcmat((i+1),jmed)=1).
            loop conti= 1 to ncol(wcontval).
              loop contj = 1 to ncol(zcontval).
                compute contvec2={contvec2,wcontval(:,conti)&*zcontval(:,contj)}.
              end loop.
            end loop.
          end if.
          compute conteff=contvec2*probcoef.
          compute contdiff=contvec2(1,:)-contvec2(2,:).
          compute contse=sqrt(contdiff*probvarb*t(contdiff)).
          compute conteffd=conteff(1,1)-conteff(2,1).
          compute contvec={contvec,conteff}.
          compute contvecm=contvec.
          !outform outtodo=contvecm/outbig=maxresm.
          do if (outscree=1).
            print/title="Contrast between conditional effects of X:".
            print contvec/title=" "/rlabels="Effect1:","Effect2:"/cnames=problabs/format = !decpnt/space=0.
          end if.
          do if (ydich=0). 
            compute p=2*(1-tcdf(abs(conteffd/contse), dfres)).
            compute contvec={conteffd,contse,conteffd/contse, p}.
            compute contvec={contvec,(conteffd-(tval*contse))}.
            compute contvec={contvec,(conteffd+(tval*contse))}.
            compute contlabs={"Contrast", hclab, "t", "p", "LLCI", "ULCI"}.
          end if.
          do if (ydich=1). 
            compute p=2*(1-cdfnorm(abs(conteffd/contse))).
            compute contvec={conteffd,contse,conteffd/contse, p}.
            compute contvec={contvec,(conteffd-(xp2*contse))}.
            compute contvec={contvec,(conteffd+(xp2*contse))}.
            compute contlabs={"Contrast", "se", "Z", "p", "LLCI", "ULCI"}.
          end if.
          !outform outtodo=contvec/outbig=maxresm.
          do if (outscree=1).
            print contvec/title="Test of Effect1 minus Effect2"/format=!decpnt/cnames=contlabs.
          end if.
        end if.
      end loop.
      /* end O */.
      do if (plot = 1 or plot = 2).
        compute datalabs={t(focpred(:,3)),outnames(1,i)}.
        do if (plot = 1).
          compute datalabs={datalabs}.
        end if.
        do if (plot = 2).
          compute datalabs={datalabs,"se", "LLCI", "ULCI"}.
        end if.
        do if ((i = nms+nys) and (ydich=1)).
          compute datalabs={datalabs,"prob"}.
        end if.
        !outform outtodo=probeplt/outbig=maxresm.
        compute datalabs={datalabs,"."}.
        do if (outscree=1).
          print/title = "Data for visualizing the conditional effect of the focal predictor:".
          print/title = "Paste text below into a SPSS syntax window and execute to produce plot."/space=0.
          compute dumb = {" ", " ", " ", " ", " ", " ", " "}.
          print datalabs/title = "DATA LIST FREE/"/format=A10. 
          print probeplt/title = "BEGIN DATA."/format=!decpnt/space=0.
          print/title = "END DATA."/space=0.
        end if.
        compute focgrph=datalabs(1,1).
        compute graphix={focgrph,graphixs,modgrph}.
        do if (((xdich=1) or (nxvls > 1)) and ((modcat=0) and (focgrph = xnames))).
          compute graphix={modgrph,graphixs,focgrph}.     
        end if.
        do if (panelgrp = 0).
          compute graphix={graphix,"."}.
        else.
          compute graphix={graphix,panelcde}.
        end if.
        do if (outscree=1).
          print graphix/title = "GRAPH/SCATTERPLOT="/format=A8/space=0.
        end if.
        do if (i = (nms + nys) and ydich = 1).
          compute graphixd=graphix.
          compute graphixd(1,3)="prob". 
          do if (outscree=1).
            print graphixd/title = "GRAPH/SCATTERPLOT="/format=A8/space=0.
          end if.
        end if.
      end if.  
    end if.
    /* end E */.
  end loop.
  /* end R */.
  release jmed, intprint,didprint.
end if.
!ENDDEFINE.

DEFINE PROCESS (vars=!charend('/') !default(xxxxx)/y=!charend('/') !default(xxxxx)/x=!charend('/') !default(xxxxx)
   /m=!charend('/') !default(xxxxx)/cov=!charend('/') !default(xxxxx)/z=!charend('/') !default(xxxxx)/w=!charend('/') !default(xxxxx)
   /v=!charend('/') !default(xxxxx)/q=!charend('/') !default(xxxxx)/total=!charend('/') !default(0)/varorder=!charend('/') !default(0)
   /mcw=!charend('/') !default(0)/mcx=!charend('/') !default(0)/mcz=!charend('/') !default(0)/normal=!charend('/') !default(0)
   /quantile=!charend('/') !default(999)/hc=!charend('/') !default(5)/hc3=!charend('/') !default(0)/moments=!charend('/') !default(0) 
   /vmodval=!charend('/') !default(999)/qmodval=!charend('/') !default(999)/percent=!charend('/') !default(1)/subsets=!charend('/') !default(0)
   /cluster=!charend('/') !default(xxxxx)/conf=!charend('/') !default(95)/boot=!charend('/') !default(5000)/ssquares=!charend('/') !default(9)
   /seed=!charend('/') !default(random)/bmatrix=!charend('/') !default("-999")/jn=!charend('/') !default(0)/debug=!charend('/') !default(0)
   /wmatrix=!charend('/') !default("-999")/zmatrix=!charend('/') !default("-999")/wzmatrix=!charend('/') !default("-999")/dominate=!charend('/') !default(0)
   /cmatrix=!charend('/') !default("-999")/model = !charend('/') !default(999)/decimals=!charend('/') !default(F10.4)/exclude=!charend('/') !default("0")
   /matrices=!charend('/') !default (0)/covmy=!charend('/') !default(0)/center=!charend('/') !default(0)/ws=!charend('/') !default(0)
   /iterate = !charend('/') !default(100)/converge = !charend('/') !default(0.00001)/stand=!charend('/') !default(0)/zpp=!charend('/') !default(0)
   /contrast=!charend('/') !default("999")/effsize=!charend('/') !default(0)/mc=!charend('/') !default(0)/mdichok=!charend('/') !default(0)
   /save=!charend('/') !default(0)/maxboot=!charend('/') !default(0)/covcoeff=!charend('/') !default(0)/activate=!charend('/') !default(0)
   /wmodval=!charend('/') !default("999")/intprobe=!charend('/') !default(.1)/zmodval=!charend('/') !default("999")/plot=!charend('/') !default(0)
    /modelbt = !charend('/') !default(0)/xcatcode=!charend('/')!default("-999")/xmtest=!charend('/') !default(0)/listmiss=!charend('/') !default(0)
    /wcatcode=!charend('/') !default("-999")/zcatcode=!charend('/') !default("-999")/describe=!charend('/') !default(0)/coval=!charend('/') !default ("-999")/
    modelres=!charend('/') !default(0)/bc=!charend('/') !default(0)/longname=!charend('/') !default(0)/cdeval=!charend ('/') !default("-999")/settest=!charend('/') !default(0)
    /outscreen=!charend('/') !default(1)/linsum=!charend('/') !default("-999")/xmint=!charend('/') !default(0)/xrefval=!charend('/') !default("999")
    /crossv=!charend('/') !default(0)/spline=!charend('/') !default("-999")/diagnose=!charend('/') !default(0)/robustse=!charend('/') !default(0)
    /clusboot=!charend('/') !default (0)/display=!charend('/') !default("text")/eiv=!charend('/') !default (3)
    /relx=!charend('/') !default("-999")/relm=!charend('/') !default("-999")/relcov=!charend('/') !default("-999")).
preserve.
set width=255.
set printback = off.
set seed = !seed.
set mdisplay=!display.
set mxloop = 100000000.
MATRIX.
compute wnames="xxxxx".
compute znames="xxxxx".
compute triprob=0.
compute mcerpt=0.
compute wiscov=0.
compute ziscov=0.
compute tooman=0.
compute sv4match=0.
compute nclsv=0.
compute nclus=0.
compute tricheck=0.
compute errcode=make(100,1,0).
compute notecode=make(100,1,0).
compute model = trunc(!model).
compute iterate = abs(trunc(!iterate)).
compute converge = abs(!converge).
compute itprobtg=0.
compute v2tag=0.
compute ydich=0.
compute maxwwarn=0.
compute minwwarn=0.
compute maxzwarn=0.
compute minzwarn=0.
compute toomany=0.
compute wdich=0.
compute zdich=0.
compute wnotev=0.
compute znotev=0.
compute nxpval=1.
compute nwpval=1.
compute nzpval=1.
compute errs=1.
compute notes=1.
compute criterr=0.
compute novar=0.
compute adjust=0.
compute ncs=0.
compute serial=0.
compute sobelok=0.
compute hasw=0.
compute hasz=0.
compute printw=0.
compute printz=0.
compute relx={!relx}.
compute relm={!relm}.
compute relcov={!relcov}.
compute eiv={!eiv}.
compute eiv=trunc(eiv).
compute yrelx=0.
compute yrelm=0.
compute yrelcov=0.
compute relxphm=0.
compute eivhcnt=0.
do if ((eiv < 0) or (eiv > 5) or (eiv=4)).
 compute errcode(errs,1)=88.
 compute errs=errs+1.
 compute criterr=1.
end if.
do if (relx(1,1)=-999).
  compute relx=1.
  else.
  compute yrelx=1.
end if.
do if (relm(1,1)=-999).
  compute relm=1.
  else.
  compute yrelm=1.
end if.
do if (relcov(1,1)=-999).
  compute relcov=1.
  else.
  compute yrelcov=1.
end if.
compute eivbad=0.
compute eivdo=0.
compute fulldatr=0.
compute robustse=(!robustse=1).
compute clusboot=trunc(!clusboot).
compute diagnose=(!diagnose=1).
compute xmint=(!xmint=1).
compute settest=(!settest=1).
compute zpp=(!zpp=1).
compute setxcat=0.
compute subsets=(!subsets=1).
compute dominate=(!dominate=1).
compute crossv=(!crossv=1).
compute spl={!spline}.
compute nspl=ncol(spl).
do if (spl(1,1) = -999).
  compute nspl=0.
end if.
compute wmodcust=0.
compute zmodcust=0.
compute badboot=0.
compute eiverrc=0.
compute eiverc=0.
compute booting=0.
compute bootiter=0.
compute iterrmod=0.
compute cov = !quote(!cov).
compute varorder=(!varorder <> 0).
compute nws=0.
compute w=!quote(!w).
compute nzs=0.
compute z = !quote(!z).
compute nms=0.
compute m = !quote(!m).
compute nys=0.
compute y = !quote(!y).
compute nxs=0.
compute x = !quote(!x).
compute effsize=(!effsize=1).
compute stand=(!stand=1).
compute intprobe = !intprobe.
compute xrefvals={!xrefval}.
compute center=trunc(!center) .
compute xcontcf=0.
compute xscaling=1.
compute cdeval={!cdeval}.
compute cuscoval=0.

compute bmatrix={!bmatrix}.
do if (x <> "xxxxx" and  y <> "xxxxx" and m <> "xxxxx" and w = "xxxxx" and z = "xxxxx" and model=999 and bmatrix(1,1)=-999).
  compute model=4.
end if.

do if (model=74).
 compute errcode(errs,1)=7.
 compute errs=errs+1.
 compute criterr=1.
end if.
do if (xmint=1 and model <> 4).
 compute errcode(errs,1)=63.
 compute errs=errs+1.
 compute criterr=1.
end if.
do if (xmint=1 and model=4).
  compute w=x.
  compute model=74.
  compute intprobe=1.
  compute notecode(notes,1) = 32.
  compute notes = notes + 1.
  do if (effsize=1 or stand=1).
    compute notecode(notes,1) = 34.
    compute notes = notes + 1.
    compute stand=0.
    compute effsize=0.
  end if.
  do if (center <> 0).
    compute center=0.
     compute errcode(errs,1)=71.
     compute errs=errs+1.
     compute criterr=1.
  end if.
end if.

compute v = !quote(!v).
compute q = !quote(!q).
compute linsum={!linsum}.
compute nlinsum=ncol(linsum).
do if (linsum(1,1) = -999).
  compute nlinsum=0.
end if.
compute deleteme={!exclude}.
compute oldvars=!quote(!vars).
compute mcxok=0.
compute mcwok=0.
compute mczok=0.
compute xprod=0.
compute zprod=0.
compute wprod=0.
compute modcok=0.
compute alttotal=0.
compute hc3=trunc(!hc3).
compute jn=(!jn = 1).
compute listmiss=(!listmiss=1).
compute modelres=(!modelres=1).
compute outscree=(!outscreen=1).
compute activate=(!activate=1).
compute booterr=0.
compute normal=(!normal=1).
compute xmtest=(!xmtest=1).
compute describe=(!describe=1).
compute longname=(!longname=1).

/* here is where the variable check would go */.

do if (stand=1).
  compute effsize=1.
end if.
compute pstog=0.
compute sobelok=0.
compute normal=(!normal=1).
compute mdichok=(!mdichok=1).
compute contrast={!contrast}.
compute ncontr=ncol(contrast).
compute ncontrow=nrow(contrast).
do if (contrast(1,1) = 999).
  compute ncontr=1.
  compute contrast=0.
end if.
do if (ncontr = 1).
  compute contrast=trunc(contrast).
  do if (contrast > 3 or contrast < 0)).
    compute ncontr=1.
    compute contrast = 0.
  end if.
end if.
do if (ncontr > 1).
    compute contvec=contrast.
    compute contrast=4.
    do if (ncontrow > 1).
      compute contrast=0.
      compute modcok=1.
      compute wcontval=contvec(:,1).
      compute zcontval=contvec(:,2).
      do if ((ncontr <> 2) or (ncontrow <> 2)).
        compute notecode(notes,1) = 19.
        compute notes = notes + 1.
        compute modcok=0.
      end if.
    end if.
end if.
do if (xmint = 1 and contrast(1,1) <> 0) .
  compute contrast=0.
  compute notecode(notes,1) = 37.
  compute notes = notes + 1.
end if.
do if (varorder = 1).
  compute notecode(notes,1) = 21.
  compute notes = notes + 1.
end if.
do if (!vmodval <> 999 or !qmodval <> 999).
  compute notecode(notes,1) = 22.
  compute notes = notes + 1.
end if.
compute modelbt=(!modelbt=1).
compute cluster=!quote(!cluster).
do if ((cluster="xxxxx") and (robustse <> 0 or clusboot <> 0)).
  compute errcode(errs,1)=79.
  compute errs=errs+1.
  compute criterr=1.
end if.
compute matrices=(!matrices=1).
compute covcoeff=(!covcoeff=1).
compute covmy=trunc(!covmy).
do if (covmy < 0 or covmy > 2).
   compute covmy = 0.
end if.
compute boot = abs(trunc(!boot)).
compute bc=(!bc=1).
compute mc=abs(trunc(!mc)).
compute hc=trunc(!hc).
do if (intprobe < 0 or intprobe > 1).
  compute intprobe = .10.
end if.
compute plot=trunc(!plot).
do if (plot < 0 or plot > 2).
  compute plot=0.
end if.
compute total=(!total=1).
compute dototal=0.
compute saveboot = (!save = 1).
compute saveest=(!save = 2).
compute savediag=(!save=4).
do if (saveest=1).
  compute intprobe=1.
end if.
do if (hc >= 0 and hc < 5 and robustse=0).
  compute notecode(notes,1) = 4.
  compute notes = notes + 1.
end if.
do if (hc > 5 or hc < 0).
  compute hc=5.
end if.
compute mcw=trunc(!mcw).
compute mcz=trunc(!mcz).
compute mcx=trunc(!mcx).
do if (mcx > 0 and mcx < 3 and model = 74).
  compute mcw=mcx.
  compute xscaling=1.
end if.
do if (mcx > 2 and model = 74).
  compute errcode(errs,1)=65.
   compute errs=errs+1.
   compute criterr=1.
end if.
do if (model = 74 and normal=1).
  compute notecode(notes,1) = 33.
  compute normal=0.
  compute notes = notes + 1.
end if.
do if (mcx > 0 and contrast > 0).
  compute notecode(notes,1) = 28.
  compute notes = notes + 1.
  compute contrast=0.
end if.
compute nxvls=1.
compute nmvls=1.
compute nwvls=1.
compute nzvls=1.
compute paths=999.
compute pathsw=999.
compute pathsz=999.
compute pathswz=999.
compute pathsmod=999.
compute pathtype=999.
compute obscoeff=999.
compute pathsdv={" "}.
compute quantile=1.
do if (!quantile<>999).
  compute notecode(notes,1) = 23.
  compute notes = notes + 1.
end if.
compute moments=(!moments=1).
do if (moments=1).
  compute quantile=0.
end if.

compute wmatrix={!wmatrix}.
compute zmatrix={!zmatrix}.
compute wzmatrix={!wzmatrix}.
compute cmatrix={!cmatrix}.
compute xcatcode={!xcatcode}.
compute wcatcode={!wcatcode}.
compute zcatcode={!zcatcode}.
compute needed=0.
compute conf=!conf.
compute ssquares=!ssquares.
do if (x <> "xxxxx" and  y <> "xxxxx" and w = "xxxxx" and z = "xxxxx" and m = "xxxxx" and (model=999 or model=0)).
  compute model=0.
  compute total=0.
  comute dototal=0.
  do if (ssquares=9).
    compute ssquares=1.
  end if.
end if.
compute intncd=0.
do if (x <> "xxxxx" and  y <> "xxxxx" and w <> "xxxxx" and z = "xxxxx" and m = "xxxxx" and model=999).
  compute model=1.
end if.
compute bmatrix={!bmatrix}.
do if (x <> "xxxxx" and  y <> "xxxxx" and m <> "xxxxx" and w = "xxxxx" and z = "xxxxx" and model=999 and bmatrix(1,1)=-999).
  compute model=4.
end if.
compute ssquares=(ssquares=1).
do if (trunc(!conf) >= 100 or (trunc(!conf) <= 50)).
  compute conf = 95.
  compute notecode(notes,1)=2.
  compute notes=notes+1.
end if.
do if (model >= 0 and model < 4 and modelbt=0).
  compute boot=0.
  compute mc=0.
  compute bc=0.
end if.
/* mc takes priority over bootstrap */.
do if (mc > 0 and boot > 0).
  compute boot=0.
  compute bc=0.
  compute clusboot=0.
end if.
do if ((boot < 1000) and (mc = 0) and (boot > 0)).
  compute boot=5000.
end if.
do if ((mc < 1000) and (boot = 0) and (mc > 0)).
  compute mc=5000.
end if.
compute badend=0.
compute priorlo = -9999999.
compute priorhi = 9999999.
compute alpha2 = (1-(conf/100))/2.
compute cilm=alpha2*2.
/* compute y5=sqrt(-2*ln(alpha2)) */.
/* compute xp2=(y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0)).
compute xp2=-idf.normal(alpha2, 0,1).
compute medlb={"     M1:";"     M2:";"     M3:";"     M4:";"     M5:";"     M6:";"     M7:";"     M8:";"     M9:";"    M10:"}.
compute medlb2={"(M1)","(M2)","(M3)","(M4)","(M5)","(M6)","(M7)","(M8)","(M9)","(M10)"}.
compute xlb={"   X1  :";"   X2  :";"   X3  :";"   X4  :";"   X5  :";"   X6  :";"   X7  :";"   X8  :";"   X9  :"}.
compute highlbw={"M1*W"; "M2*W"; "M3*W"; "M4*W"; "M5*W"; "M6*W"; "M7*W"; "M8*W"; "M9*W"; "M10*W"}.
do if (xmint=1).
  compute highlbw={"M1*X"; "M2*X"; "M3*X"; "M4*X"; "M5*X"; "M6*X"; "M7*X"; "M8*X"; "M9*X"; "M10*X"}.
end if.
compute highlbz={"M1*Z"; "M2*Z"; "M3*Z"; "M4*Z"; "M5*Z"; "M6*Z"; "M7*Z"; "M8*Z"; "M9*Z";"M10*Z"}.
compute highlbwz={"M1*W*Z"; "M2*W*Z"; "M3*W*Z"; "M4*W*Z"; "M5*W*Z"; "M6*W*Z"; "M7*W*Z"; "M8*W*Z"; "M9*W*Z";"M10*W*Z"}.
compute highlbbt={"BOTH(M1)"; "BOTH(M2)"; "BOTH(M3)"; "BOTH(M4)"; "BOTH(M5)"; "BOTH(M6)"; "BOTH(M7)"; "BOTH(M8)"; "BOTH(M9)";"BTH(M10)"}.
compute highlbx={"M1*X"; "M2*X"; "M3*X"; "M4*X"; "M5*X"; "M6*X"; "M7*X"; "M8*X"; "M9*X"; "M10*X"}.

/* check to see if model number specified is acceptable for this version */.
compute skipwz=0.
compute validm={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1}.

do if (activate=1).
 compute errcode(errs,1)=60.
 compute errs=errs+1.
 compute criterr=1.
end if.

do if (criterr=0).
do if ((!ws=1)=1).
    compute errcode(errs,1)=42.
    compute errs=errs+1.
    compute criterr=1.
end if.
do if (model > 0 and model < 93).
  do if (validm(1,model)=0).
    compute errcode(errs,1)=6.
    compute errs=errs+1.
    compute criterr=1.
  end if.
  release validm.
end if.



do if ((model > 92 or model < 0) and model <> 999)).
  compute errcode(errs,1)=7.
  compute errs=errs+1.
  compute criterr=1.
end if.
do if (model = 999 and bmatrix(1,1)=-999).
  compute errcode(errs,1)=24.
  compute errs=errs+1.
  compute criterr=1.
end if.
do if (model <> 999 and bmatrix(1,1) <> -999).
  compute errcode(errs,1)=25.
  compute errs=errs+1.
  compute criterr=1.
end if.
do if ((model = 74 or (model > 0 and model < 4)) and ((wmatrix(1,1) <> -999) or (zmatrix(1,1)<>-999) or (wzmatrix(1,1)<>-999)))).
  compute errcode(errs,1)=41.
  compute errs=errs+1.
  compute criterr=1.
end if.
do if (hc3 <> 0).
  compute notecode(notes,1) = 5.
  compute notes = notes + 1.
  do if (hc3 = 1).
    compute hc=3.
  end if.
end if.
/* Do not allow V or Q */.
do if ((v <> "xxxxx") or (q <> "xxxxx")).
  compute errcode(errs,1)=14.
  compute errs=errs+1.
  compute errcode(errs,1)=48.
  compute errs=errs+1.
  compute criterr=1.
end if.
do if (oldvars <> "xxxxx").
  compute errcode(errs,1)=48.
  compute errs=errs+1.
  compute criterr=1.
end if.
/* Cluster option is disabled */.
do if (cluster <> "xxxxx").
  /* compute errcode(errs,1)=27 */.
  /* compute errs=errs+1 */.
  /* compute criterr=1 */.
end if.
/* Check to make sure an X and a Y are specified */.
do if ((y = "xxxxx") or (x = "xxxxx")).
  compute errcode(errs,1)=1.
  compute errs=errs+1.
  compute criterr=1.
end if.
/* Check to make sure an M is specified if needed */.
do if ((m = "xxxxx") and model > 3 and model <> 999).
  compute errcode(errs,1)=8.
  compute errs=errs+1.
  compute criterr=1.
end if.
end if.
/* start A */.

do if (criterr=0).  

  /* get fntmp/names = fnames/MISSING = 99999 */.
  /* release fntmp */. 
  /* print fnames/format=a8 */.
  
  compute intlab=make(200,1," ").
  compute dfinds=make(200,1," ").
  compute dfindsC=make(200,1," ").
  compute indlbl=make(200,1," ").
  compute cntname=make(200,1," ").
  !do !i = 1 !to 200.
    !let !v2=!concat(Int_,!unquote(!quote(!i))).
    !let !v3=!concat(dfie_,!unquote(!quote(!i))).
    !let !v4=!concat(DFIE_,!unquote(!quote(!i))).
    !let !v5=!concat(Ind,!unquote(!quote(!i))).
    !let !v6=!concat("(C",!unquote(!quote(!i)),")").
    compute intlab(!i,1)=!quote(!v2).
    compute dfinds(!i,1)=!quote(!v3).
    compute dfindsC(!i,1)=!quote(!v4).
    compute indlbl(!i,1)=!quote(!v5).
    compute cntname(!i,1)=!quote(!v6).
  !doend.
  compute dfbetas=make(201,1," ").
  compute dfbetasC=make(201,1," ").
  !do !i = 0 !to 200.
    !let !v2=!concat(dfb_,!unquote(!quote(!i))).
    !let !v3=!concat(DFB_,!unquote(!quote(!i))).
    compute dfbetas((!i+1),1)=!quote(!v2).
    compute dfbetasC((!i+1),1)=!quote(!v3).
  !doend.

  /* extract data and variable names */.
  get ytmp/variables = !y/names = ynames/MISSING = 99999.
  compute nys=ncol(ytmp).
  compute needed=nys.
  compute n=nrow(ytmp).
  compute varnames={ynames}.
  compute dat=ytmp.
  !longchk variab=!y.
  compute modelvar={!quote(!model);t(ynames)}.
  do if (!model=999).
    compute modelvar(1,1)="CUSTOM".
    do if (model=0).
        compute modelvar(1,1)="0".
    end if.
    do if (model=1).
        compute modelvar(1,1)="1".
    end if.
    do if ((model=4) or (model=74)).
      compute modelvar(1,1)="4".
    end if.
  end if.
  get xtmp/variables = !x/names = xnames/MISSING = 99999.
  compute nxs=ncol(xtmp).
  do if ((ncol(relx) <> nxs) and (yrelx=1) and (eivbad=0)).
   compute errcode(errs,1)=83.
   compute errs=errs+1.
   compute criterr=1.
   compute eivbad=1.
  end if.
  do if (yrelx=1 and eivbad=0).
    compute eivdo=1.
  end if.
  compute n=nrow(xtmp).
  compute nxset=nxs.
  compute xtraxs=0.
  compute xnames2=xnames.
  compute xfakecov=0.
  compute xmodelbs= "      X:".
  do if ((nxs > 1) and (model=0)).
    compute xotmp=xtmp(:,2:ncol(xtmp)).
    compute xfakecov=ncol(xotmp).
    compute ncs=xfakecov.
    compute xonames=xnames(1,2:ncol(xnames)).
    compute covnames=xonames.
    compute xnames=xnames(1,1).
    do if (yrelx=0).
      compute relx=make(1,nxs,1).
    end if.
    do if (eivbad=0).
      compute relxph=relx(1,2:ncol(relx)).
      compute relxphm=1.
      compute relx=relx(1,1).
    end if.
    compute nxs=1.
    compute xtralbs=make((ncol(xtmp)-1),1," ").
    compute xtmp=xtmp(:,1).
    compute xtraxs=1.
    compute xmodelbs={xmodelbs;xtralbs}.
  end if.
  compute eivrel=relx.
  compute needed=needed+nxs.
  compute varnames={varnames,xnames}.
  /* do if (xtraxs=1) */.
  /*  compute varnames={varnames,xonames} */.
  /* end if */.
  compute xcatlab=t(xnames).
  compute dat={dat,xtmp}.
  !longchk variab=!x.
  compute modelvar={modelvar;t(xnames2)}.
  /* do if (nxs = 1 and xtraxs=0) */.
    compute modelvlb={"  Model:";"      Y:";xmodelbs}.
  /* else */.
    /* compute modelvlb={"  Model:";"      Y:";xlb(1:nxs,1)} */.
  /* end if */.

  do if (m <> "xxxxx").
    get mtmp/variables = !m/names = mnames/MISSING = 99999.
    compute nms=ncol(mtmp).
    do if ((ncol(relm) <> nms) and (yrelm=1) and (eivbad=0)).
      compute errcode(errs,1)=83.
      compute errs=errs+1.
      compute criterr=1.
      compute eivbad=1.
    end if.
    do if (yrelm=0).
      compute relm=make(1,nms,1).
    end if.
    do if (yrelm=1 and eivbad=0).
      compute eivdo=1.
    end if.
    compute eivrel={eivrel,relm}.
    compute mprod=make(1,nms,0).
    compute n=nrow(mtmp).
    compute needed=needed+nms.
    compute varnames={varnames,mnames}.
    compute dat={dat,mtmp}.
    compute modelvar={modelvar;t(mnames)}.
    compute x2m=make(99,nms,0).
    compute m2y=make(99,nms,0).
    compute onem=make(nms,1,1).
    !longchk variab=!m.
    do if (nms > 1 and nms < 11).
      compute modelvlb={modelvlb;medlb(1:nms,1)}.
    else.
      compute modelvlb={modelvlb;"      M:"}.
    end if.
    do if (nms > 0 and model < 4).
      compute errcode(errs,1)=9.
      compute errs=errs+1.
      do if (model <> 0).
        compute errcode(errs,1)=48.
        compute errs=errs+1.
      end if.
      compute criterr=1.
    end if.
  end if.
  compute wlocatet=0.
  compute wlocate=0.
  do if (w <> "xxxxx").
    do if (xmint=0).
      get wtmp/variables = !w/names = wnames/MISSING = 99999.
      !longchk variab=!w.
    end if.
    do if (xmint=1).
      get wtmp/variables = !x/names = wnames/MISSING = 99999.
    end if.
    compute nws=ncol(wtmp).
    compute n=nrow(wtmp).
    compute varnames={varnames,wnames}.
    compute wlocate=ncol(varnames).
    do if (model=74).
      compute wlocatet=1.
      do if (xnames <> wnames).
        compute errcode(errs,1)=45.
        compute errs=errs+1.
        compute criterr=1.
      end if.
    end if.
    compute wcatlab=t(wnames).
    compute dat={dat,wtmp}.
    do if (xmint <> 1).
      compute modelvar={modelvar;t(wnames)}.
      compute modelvlb={modelvlb;"      W:"}.
    end if.
  end if.
  do if (z <> "xxxxx").
    get ztmp/variables = !z/names = znames/MISSING = 99999.
    compute nzs=ncol(ztmp).
    compute n=nrow(ztmp).
    !longchk variab=!z.
    compute varnames={varnames,znames}.
    compute zcatlab=t(znames).
    compute dat={dat,ztmp}.
    compute modelvar={modelvar;t(znames)}.
    compute modelvlb={modelvlb;"      Z:"}.
  end if.
  do if (xtraxs=1 and cov = "xxxxx").
    compute dat={dat,xotmp}.
    compute ctmp=xotmp.
    do if (eivbad=0).
      compute relcov=relxph.
      compute eivrel={eivrel,relcov}.
    end if.
    compute varnames={varnames,xonames}.
  end if. 
  do if (cov <> "xxxxx").
    get ctmp/variables = !cov/names = covnames/MISSING = 99999.
    compute ncs=ncol(ctmp)+xfakecov.
    do if ((ncol(relcov) <> ncol(ctmp)) and (yrelcov=1) and (eivbad=0)).
      compute errcode(errs,1)=83.
      compute errs=errs+1.
      compute criterr=1.
      compute eivbad=1.
    end if.
    do if (yrelcov=0).
      compute relcov=make(1,ncs,1).
      /* do if (model=0) */.
        do if relxphm=1.
          compute relcov(1,1:ncol(relxph))=relxph.
        end if.
        compute eivrel={eivrel,relcov}.
      /* end if  */.
    end if.
    do if (yrelcov=1 and eivbad=0).
      /* do if (model=0) */.
        do if (relxphm=1).
        compute relcov={relxph,relcov}.
        end if.
        compute eivrel={eivrel,relcov}.
      /* end if */.
      compute eivdo=1.
    end if.
    do if ((ncol(ctmp) > 0) and model=0).
      compute settest=1.
    end if.
    do if (xtraxs=1).
      compute ctmp={xotmp,ctmp}.
      compute covnames={xonames,covnames}.
    end if.
    compute n=nrow(ctmp).
    !longchk variab=!cov.
    compute varnames={varnames,covnames}.
    compute dat={dat,ctmp}.
  end if.
  compute clusok=0.
  do if (cluster <> "xxxxx").
    get clusdat/variables = !cluster/names = clunames/MISSING = 99999.
    compute nclsv=ncol(clusdat).
    compute dat={dat,clusdat}.
    !longchk variab=!cluster.
    do if (nclsv=1).
      compute modelvar={modelvar;clunames}.
      compute modelvlb={modelvlb;"Cluster:"}.
      compute clusok=1.
      do if (robustse=1).
        compute hc=6.
        compute notecode(notes,1) = 4.
        compute notes = notes + 1.
      end if.
      do if (robustse=0 and clusboot=0).
        compute errcode(errs,1)=78.
        compute errs=errs+1.
        compute criterr=1.
      end if.
    end if.
  end if.

 
  /* check for too many variables specified */.
  do if (nws > 1 or nzs > 1 or nys > 1 or (nxs > 1 and xtraxs=0) or nclsv > 1).
    compute errcode(errs,1)=3.
    compute errs=errs+1.
    compute criterr=1.
  end if.

  /* check for proper number of mediators in some models */.
  do if ((model = 80 or model = 81) and (nms < 3 or nms > 6)).
    compute errcode(errs,1)=32.
    compute errs=errs+1.
    compute criterr=1.    
  end if.
  do if (model = 82 and nms <> 4).
    compute errcode(errs,1)=33.
    compute errs=errs+1.
    compute criterr=1.    
  end if.
  do if (nms > 10).
    compute errcode(errs,1)=37.
    compute errs=errs+1.
    compute criterr=1.    
  end if.
  do if ((model = 6 or (model > 82 and model < 999)) and (nms < 2 or nms > 6)).
    compute errcode(errs,1)=34.
    compute errs=errs+1.
    compute criterr=1.    
  end if.
  /* check for redundant variable names */.
  compute match=0.
  compute match2=0.
  compute mcwzcov=0.
  loop i = 1 to (ncol(varnames)-1).
    loop j = (i+1) to ncol(varnames).
      do if (varnames(i)=varnames(j)).
        do if (i < (nxs+nms+nys+1)).
          compute match2=match2+1.
        end if.
        do if (model=0).
           compute match2=match2+1.
        end if.
        do if (wlocatet=1 and i=2 and j=wlocate).
          compute match2=match2-1.
        end if. 
        do if ((wnames=znames) and (nws > 0 or nzs > 0))).
          compute match2=match2+1.
        end if.
        do if (i < (ncol(varnames)-ncs+1)) and j > (ncol(varnames)-ncs)).              
          do if ((varnames(j)=wnames) and mcw=0)).
            compute match=0.
            compute wiscov=(j-(ncol(varnames)-ncs)).
          end if.
          do if ((varnames(j)=wnames) and mcw <>0)).
            compute mcwzcov=1.
          end if.
          do if ((varnames(j)=znames) and mcz=0).
            compute match=0.
            compute ziscov=(j-(ncol(varnames)-ncs)).
          end if.
          do if ((varnames(j)=znames) and mcz<>0)).
             compute mcwzcov=1.
          end if.
        end if.
      end if.
    end loop.
  end loop.
  do if (match2>0 or match=1).
    compute errcode(errs,1)=2.
    compute errs=errs+1.
    compute criterr=1.
  end if.
  do if (mcwzcov=1).
    compute errcode(errs,1)=50.
    compute errs=errs+1.
    compute criterr=1.
  end if.  
    /* listwise deletion */.
  compute ninit=nrow(dat).
  compute rownum=make(ninit,1,0).
  compute rownumd=make(ninit,1,0).
  loop i = 1 to ninit.
    compute rownum(i,1)=i.
    loop k = 1 to ncol(deleteme).
      do if (deleteme(1,k)=i).
        compute rownumd(i,1)=1.
      end if.
    end loop.
  end loop.
  compute dat={rownum,dat}.
  compute j=1.
  compute missrow=0.
  compute delrow=0.
  compute nmiss=0.
  compute delident=0.
  loop i = 1 to ninit.
    compute delskip=0.
    do if ((rownumd(i,1)=1) and (delskip=0)).
      compute delskip=1.
      compute delident=1.
    end if.
    do if ((rsum(dat(i,2:ncol(dat))=99999)<>0) and delskip=0).
      compute missrow={missrow;dat(i,1)}.
      compute nmiss=nmiss+1.
      compute delskip=1.
    end if.
    do if (delskip=0).
      compute dat(j,:)=dat(i,:).
      compute j=j+1.
    end if.
  end loop.
  do if (j < 5).  
    compute errcode(errs,1)=62.
    compute errs=errs+1.
    compute criterr=1.
  end if.
  do if (delident=1).  
    compute notecode(notes,1)=38.
    compute notes=notes+1.
  end if.
  compute jjjj=j.
  /* startit */.
  
  compute bcmat=make(needed,needed,0).
  compute wcmat=make(needed,needed,0).
  compute zcmat=make(needed,needed,0).
  compute wzcmat=make(needed,needed,0).
  compute wsum=0.
  compute zsum=0.
  compute wzsum=0.
  
  /* DEFINE MODEL MATRICES FOR CANNED MODELS */.
  do if (criterr = 0 and model <> 999).  
    /* X->MW, X->MZ, X->MWZ, M->YW, M->YZ, M->YWZ, X->YW, X->YZ, X->YWZ */.
    compute modelmat=
    {1,0,0,0,0,0,0,1,0,0;2,0,0,0,0,0,0,1,1,0;3,0,0,0,0,0,0,1,1,1;4,0,0,0,0,0,0,0,0,0;
    5,0,0,0,0,0,0,1,0,0;6,0,0,0,0,0,0,0,0,0;7,1,0,0,0,0,0,0,0,0;8,1,0,0,0,0,0,1,0,0;
    9,1,1,0,0,0,0,0,0,0;10,1,1,0,0,0,0,1,1,0;11,1,1,1,0,0,0,0,0,0;12,1,1,1,0,0,0,1,1,1;
    13,1,1,1,0,0,0,1,0,0;14,0,0,0,1,0,0,0,0,0;15,0,0,0,1,0,0,1,0,0;16,0,0,0,1,1,0,0,0,0;
    17,0,0,0,1,1,0,1,1,0;18,0,0,0,1,1,1,0,0,0;19,0,0,0,1,1,1,1,1,1;20,0,0,0,1,1,1,1,0,0;
    21,1,0,0,0,1,0,0,0,0;22,1,0,0,0,1,0,1,0,0;23,0,0,0,0,0,0,0,0,0;24,0,0,0,0,0,0,0,0,0;
    25,0,0,0,0,0,0,0,0,0;26,0,0,0,0,0,0,0,0,0;27,0,0,0,0,0,0,0,0,0;28,1,0,0,0,1,0,0,1,0;
    29,1,0,0,0,1,0,1,1,0;30,0,0,0,0,0,0,0,0,0;31,0,0,0,0,0,0,0,0,0;32,0,0,0,0,0,0,0,0,0;
    33,0,0,0,0,0,0,0,0,0;34,0,0,0,0,0,0,0,0,0;35,0,0,0,0,0,0,0,0,0;36,0,0,0,0,0,0,0,0,0;
    37,0,0,0,0,0,0,0,0,0;38,0,0,0,0,0,0,0,0,0;39,0,0,0,0,0,0,0,0,0;40,0,0,0,0,0,0,0,0,0;
    41,0,0,0,0,0,0,0,0,0;42,0,0,0,0,0,0,0,0,0;43,0,0,0,0,0,0,0,0,0;44,0,0,0,0,0,0,0,0,0;
    45,0,0,0,0,0,0,0,0,0;46,0,0,0,0,0,0,0,0,0;47,0,0,0,0,0,0,0,0,0;48,0,0,0,0,0,0,0,0,0;
    49,0,0,0,0,0,0,0,0,0;50,0,0,0,0,0,0,0,0,0;51,0,0,0,0,0,0,0,0,0;52,0,0,0,0,0,0,0,0,0;
    53,0,0,0,0,0,0,0,0,0;54,0,0,0,0,0,0,0,0,0;55,0,0,0,0,0,0,0,0,0;56,0,0,0,0,0,0,0,0,0;
    57,0,0,0,0,0,0,0,0,0;58,1,0,0,1,0,0,0,0,0;59,1,0,0,1,0,0,1,0,0;60,1,1,0,1,0,0,0,0,0;
    61,1,1,0,1,0,0,1,0,0;62,1,1,0,1,0,0,0,1,0;63,1,1,0,1,0,0,1,1,0;64,1,0,0,1,1,0,0,0,0;
    65,1,0,0,1,1,0,1,0,0;66,1,0,0,1,1,0,0,1,0;67,1,0,0,1,1,0,1,1,0;68,1,1,1,1,0,0,0,0,0;
    69,1,1,1,1,0,0,1,1,1;70,1,0,0,1,1,1,0,0,0;71,1,0,0,1,1,1,1,1,1;72,1,1,1,1,1,1,0,0,0;
    73,1,1,1,1,1,1,1,1,1;74,0,0,0,1,0,0,0,0,0;75,1,1,0,1,1,0,0,0,0;76,1,1,0,1,1,0,1,1,0;
    77,0,0,0,0,0,0,0,0,0;78,0,0,0,0,0,0,0,0,0;79,0,0,0,0,0,0,0,0,0;80,0,0,0,0,0,0,0,0,0;
    81,0,0,0,0,0,0,0,0,0;82,0,0,0,0,0,0,0,0,0;83,1,0,0,0,0,0,0,0,0;84,1,0,0,0,0,0,0,0,0;
    85,1,0,0,0,0,0,1,0,0;86,1,0,0,0,0,0,1,0,0;87,0,0,0,1,0,0,0,0,0;88,0,0,0,1,0,0,0,0,0;
    89,0,0,0,1,0,0,1,0,0;90,0,0,0,1,0,0,1,0,0;91,0,0,0,0,0,0,0,0,0;92,1,0,0,1,0,0,1,0,0}.
    do if (model > 0).
      compute tmp=modelmat(model,2:ncol(modelmat)).
    end if.
    do if (model=0).
      compute tmp={0,0,0,0,0,0,0,0,0}.
    end if.
    do if (model < 4).
      compute bcmat((nxs+1),1)=1.
    end if.
    do if ((model > 3) and (model <> 6)).
      compute bcmat((nxs+1):(nxs+nms),1)=onem.
      compute bcmat(nrow(bcmat),(nxs+1):(nxs+nms))=t(onem).
      compute bcmat(nrow(bcmat),1)=1.
    end if.
    /* models 6, 80, 81, and 82 are serial models */.
    do if ((model = 6) or (model > 82 and model < 93)).
      loop j = 2 to nrow(bcmat).
        loop i = 1 to (j-1).
          compute bcmat(j,i)=1.
        end loop.
      end loop.
    end if.
    do if (model = 80).
      loop i = 1 to nms.
        compute bcmat((nrow(bcmat)-1),i)=1.
      end loop.
    end if.
    do if (model = 81).
      loop j = 3 to nrow(bcmat).
        compute bcmat(j,2)=1.
      end loop.
    end if.
    do if (model = 82).
      compute bcmat(3,2)=1.
      compute bcmat(5,4)=1.
    end if.
    /* set up moderator matrices */.
    do if (tmp(1,1)=1).
      compute wcmat((nxs+1):(nxs+nms),1)=onem.
      compute wprod=1.
      compute xprod=1.
        do if (model = 83 or model = 86).
         compute onemsx=onem.
         loop i = 1 to (nms-1).
           compute onemsx(i+1,1)=0.
         end loop.
         compute wcmat((nxs+1):(nxs+nms),1)=onemsx.
      end if.
    end if.
    do if (tmp(1,4)=1).
      compute wcmat(nrow(wcmat),(nxs+1):(nxs+nms))=t(onem).
      compute wprod=1.
      do if (model = 87 or model = 90).
        compute onemsx=onem.
        loop i = 1 to (nms-1).
          compute onemsx(i,1)=0.
        end loop.
        compute wcmat(nrow(wcmat),(nxs+1):(nxs+nms))=t(onemsx).
      end if.
    end if.
    do if (tmp(1,7)=1).
      compute wcmat(nrow(wcmat),1)=1.
      compute wprod=1.
      compute xprod=1.
    end if.
    do if (tmp(1,2)=1).
      compute zcmat((nxs+1):(nxs+nms),1)=onem.
      compute zprod=1.
      compute xprod=1.
    end if.
    do if (tmp(1,5)=1).
      compute zcmat(nrow(zcmat),(nxs+1):(nxs+nms))=t(onem).
      compute zprod=1.
    end if.
    do if (tmp(1,8)=1).
      compute zcmat(nrow(zcmat),1)=1.
      compute zprod=1.
      compute xprod=1.
    end if.
    do if (tmp(1,3)=1).
      compute wzcmat((nxs+1):(nxs+nms),1)=onem.
      compute xprod=1.
      compute wprod=1.
      compute zprod=1.
    end if.
    do if (tmp(1,6)=1).
      compute wzcmat(nrow(wzcmat),(nxs+1):(nxs+nms))=t(onem).
      compute zprod=1.
      compute wprod=1.
    end if.
    do if (tmp(1,9)=1).
      compute wzcmat(nrow(wzcmat),1)=1.
      compute xprod=1.
      compute wprod=1.
      compute zprod=1.
    end if.
    do if (model = 91 or model = 92).
      loop j = 1 to (nms-1).
        loop i = 1 to j.
          compute wcmat((nxs+1+j),(nxs+i))=1.
        end loop.
      end loop.
    end if.
    do if (nms < 0).
      loop i = 1 to nms.
        /* compute tmp=wcmat(nrow(wcmat),(1+i))+zcmat(nrow(zcmat),(1+i))+wzcmat(nrow(wzcmat),(1+i)) */.
        compute tmp=csum(wcmat(:,(1+i)))+csum(zcmat(:,(1+i)))+csum(wzcmat(:,(1+i))).
        compute mprod(1,i)=(tmp>0).
      end loop.
    end if.
  end if.  

  /* DEFINE MATRICES FOR CUSTOM MODELS AND DO SOME ERROR CHECKING */.
  do if (criterr=0).
    compute needed=needed*(needed-1)/2.
    compute nopath=0.
    do if (bmatrix(1,1) <> -999).
      compute tmp=1.
      do if ((ncol(bmatrix) <> needed) or (csum(rsum(bmatrix))=0)).
        compute errcode(errs,1)=16.
        compute errs=errs+1.
        compute criterr=1.
      else.
        loop i = 2 to nrow(bcmat).
          loop j = 1 to (i-1).
            compute bcmat(i,j)=1-(bmatrix(1,tmp) = 0).
            compute tmp=tmp+1.
          end loop.
        end loop.
      end if.
     /* check to make sure X affects something */.
      do if ((csum(bcmat(:,1))=0) and criterr=0).
        compute errcode(errs,1)=22.
        compute errs=errs+1.
        compute criterr=1.
      end if.
      /* check to make sure Y is affected by something */.
      do if ((rsum(bcmat(nrow(bcmat),:))=0) and criterr=0).
        compute errcode(errs,1)=23.
        compute errs=errs+1.
        compute criterr=1.
      end if.
      /* check for dangling mediators */.
      compute dm=0.
      do if (nms > 0).
        loop i = 1 to nms.
          do if (((rsum(bcmat((nxs+i),:)) = 0) or (csum(bcmat(:,(nxs+i))) = 0)) and (dm=0) and (criterr=0)).
            compute errcode(errs,1)=26.
            compute errs=errs+1.
            compute criterr=1.
            compute dm=1.
          end if.
        end loop.
      end if.
      release dm.
    end if.
  end if.



/* start b */.
do if (criterr=0).
  do if (wmatrix(1,1) <> -999).
    compute tmp=1.
    do if (ncol(wmatrix) <> needed).
      compute errcode(errs,1)=17.
      compute errs=errs+1.
      compute criterr=1.
    else.
      compute modelvar(1,1)="CUSTOM".
      loop i = 2 to nrow(wcmat).
        loop j = 1 to (i-1).
          compute wcmat(i,j)=1-(wmatrix(1,tmp) = 0).
          /* dont allow to specify moderation of a path that doesnt exist */.
          do if ((wcmat(i,j)=1) and (bcmat(i,j)=0) and (nopath=0)).
            compute errcode(errs,1)=20.
            compute errs=errs+1.
            compute criterr=1.
            compute nopath=1.
          end if.
          compute tmp=tmp+1.
        end loop.
      end loop.
    end if.
  end if.
  do if (zmatrix(1,1) <> -999).
    compute tmp=1.
    do if (ncol(zmatrix) <> needed).
      compute errcode(errs,1)=18.
      compute errs=errs+1.
      compute criterr=1.
    else.
      compute modelvar(1,1)="CUSTOM".
      /* dont allow for a ZMODEL without W being used somewhere in the model */.
      do if (csum(rsum(wcmat))=0 and model=999).
        compute errcode(errs,1)=21.
        compute errs=errs+1.
        compute criterr=1.          
      end if.
      loop i = 2 to nrow(zcmat).
        loop j = 1 to (i-1).
          compute zcmat(i,j)=1-(zmatrix(1,tmp) = 0).
          /* dont allow to specify moderation of a path that doesnt exist */.
          do if ((zcmat(i,j)=1) and (bcmat(i,j)=0) and (nopath=0)).
            compute errcode(errs,1)=20.
            compute errs=errs+1.
            compute criterr=1.
            compute nopath=1.
          end if.
          compute tmp=tmp+1.
        end loop.
      end loop.
    end if.
  end if.
  compute tmp=1.

  do if (wzmatrix(1,1) <> -999).
    do if (ncol(wzmatrix) <> needed).
      compute errcode(errs,1)=19.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    compute modelvar(1,1)="CUSTOM".
  end if.
  do if (criterr=0).
    loop i = 2 to nrow(wzcmat).
      loop j = 1 to (i-1).
        /* set corresponding elements in W and Z for three way interaction */.
        do if (wzmatrix(1,1) <> -999).
          compute wzcmat(i,j)=1-(wzmatrix(1,tmp) = 0).
        end if.
        do if (wzcmat(i,j)=1).
          compute wcmat(i,j)=1.
          compute zcmat(i,j)=1.
        end if.
        /* dont allow to specify moderation of a path that doesnt exist */.
        do if ((wzcmat(i,j)=1) and (bcmat(i,j)=0) and (nopath=0)).
          compute errcode(errs,1)=20.
          compute errs=errs+1.
          compute criterr=1.
          compute nopath=1.
        end if.
        compute tmp=tmp+1.
      end loop.
    end loop.
  end if.
end if.
/* end B */.
    

  do if (criterr=0).
    compute xprod=csum(wcmat(:,1))+csum(zcmat(:,1))+csum(wzcmat(:,1)).
    compute xprod=(xprod > 0).
    compute wsum=csum(rsum(wcmat)).
    compute wprod=(wsum > 0).
    compute zsum=csum(rsum(zcmat)).
    compute zprod=(zsum > 0).
    do if (nms > 0).
      loop i = 1 to nms.
        /* compute tmp=wcmat(nrow(wcmat),(1+i))+zcmat(nrow(zcmat),(1+i))+wzcmat(nrow(wzcmat),(1+i)) */.
        compute tmp=csum(wcmat(:,(1+i)))+csum(zcmat(:,(1+i)))+csum(wzcmat(:,(1+i))).
        compute mprod(1,i)=(tmp>0).
      end loop.
    end if.
    do if ((wsum > 0) and (w = "xxxxx")).
      compute errcode(errs,1)=11.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((wsum = 0) and (w <> "xxxxx")).
      compute errcode(errs,1)=10.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((zsum > 0) and (z = "xxxxx")).
      compute errcode(errs,1)=13.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((zsum = 0) and (z <> "xxxxx")).
      compute errcode(errs,1)=12.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((zsum > 0) and (wsum = 0)).
      compute errcode(errs,1)=35.
      compute errs=errs+1.
      compute criterr=1.
    end if.
  end if.
  
 
  
  do if (criterr=0).
    compute j=jjjj.
    compute rownum=dat(1:(j-1),1).
    do if (nrow(missrow) > 1).
      compute missrow=t(missrow(2:nrow(missrow),1)).
      compute notecode(notes,1) = 29.
      compute notes = notes + 1.
    end if.  
    compute dat=dat(1:(j-1),2:ncol(dat)).
    compute n=nrow(dat).  
    do if (clusok=1).
      compute clusdat=dat(:,ncol(dat)).
      compute clusdat={rownum,clusdat,clusdat}.
      loop i = 1 to nrow(clusdat).
        compute clusdat(i,1)=i.
      end loop.
      compute temp = clusdat.
      compute temp(GRADE(clusdat(:,2)),:) = clusdat.
      compute clusdat = temp.
      loop i = 1 to nrow(clusdat).
        compute clusdat(i,3)=i.
      end loop.
      compute current={clusdat(1,2),1,1,1}.
      compute currentc=1.
      loop i = 2 to nrow(clusdat).
        do if (clusdat(i,2)=current(currentc,1)).
          compute current(currentc,2)=current(currentc,2)+1.
          compute current(currentc,4)=i.
        else.
          compute current={current;clusdat(i,2),1,i,i}.
          compute currentc=currentc+1.
        end if.        
      end loop.
      compute nclus=nrow(current).
      compute clusinfo={nrow(current);csum(current(:,2))/nclus;cmin(current(:,2));cmax(current(:,2))}.
      compute temp(GRADE(clusdat(:,1)),:) = clusdat.
      compute clusdat = temp.
      compute dat=dat(:,1:(ncol(dat)-clusok)).
    end if.
    
  /* reextract the data columns after deletion */.
  compute ytmp=dat(:,1:nys).
  !describ descdatf=ytmp/type=1.
  compute ysd=desctmp(2,:).
  compute ovsd=ysd.
  do if (desctmp(8,1)=1).
    compute ydich=1.
    do if (robustse=1).
      compute robustse=0.
      compute errcode(errs,1) = 81.
      compute errs = errs + 1.
      compute criterr=1.
    end if.
    do if (ydich=1 and eivdo=1).
      compute eivdo=0.
      compute errcode(errs,1) = 86.
      compute errs = errs + 1.
      compute criterr=1.
    end if.
    do if (total=1).
      compute total=0.
      compute notecode(notes,1) = 24.
      compute notes = notes + 1.
    end if.
    do if (effsize=1).
      compute effsize=0.
      compute notecode(notes,1) = 25.
      compute notes = notes + 1.
    end if.
    do if (model=74).
      compute errcode(errs,1) = 72.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    compute omx = cmax(ytmp).
    compute omn = cmin(ytmp).
    compute ytmp = (ytmp = omx).
    compute dat(:,1:nys)=(dat(:,1:nys)=omx).
    compute rcd = {omn, 0; omx, 1}.
  end if.
  compute xtmp=dat(:,(nys+1):(nys+nxs)).
  !describ descdatf=xtmp/type=(1-xprod).
  compute xsd=desctmp(2,:).
  compute xmodvals=modvals.
  compute xdich=desctmp(8,1).
  compute xmx=cmax(xtmp).
  compute xmn=cmin(xtmp).
  do if ((mcx > 0) and (xrefvals(1,1) <> 999) and (xmint=1) and (model= 74)).
    compute notecode(notes,1) = 36.
    compute notes = notes + 1.
  end if.
  do if (mcx=0).
    do if (ncol(xrefvals)>2 and model=74 and xmint=1).
      compute errcode(errs,1)=67.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((model=74) and (xmint=1)).
      do if ((xrefvals(1,1)=999) and nxvls=1 and xdich=0).
        compute errcode(errs,1)=66.
        compute errs=errs+1.
        compute criterr=1.
      end if.
      do if ((xrefvals(1,1)=999) and (xdich=1)).
        compute xrefvals={xmn,xmx}.
        compute xscaling=xrefvals(1,2)-xrefvals(1,1).
      end if.
      do if (ncol(xrefvals)=1 and xrefvals(1,1) <> 999).
        do if (xdich=0).
          compute xrefvals={xrefvals,(xrefvals(1,1)+1)}.
          compute xscaling=xrefvals(1,2)-xrefvals(1,1).
        end if.
        do if (xdich=1).
          do if ((xrefvals(1,1) <> xmx)  and (xrefvals(1,1) <> xmn)).
            compute errcode(errs,1)=70.
            compute errs=errs+1.
            compute criterr=1.
          end if.
          do if (xrefvals(1,1) = xmx).
             compute xrefvals={xrefvals,xmn}.
             compute xscaling=xrefvals(1,2)-xrefvals(1,1).
          end if.
          do if (xrefvals(1,1) = xmn).
             compute xrefvals={xrefvals,xmx}.
             compute xscaling=xrefvals(1,2)-xrefvals(1,1).
          end if.
        end if.
      end if.
      do if (ncol(xrefvals)=2).
        compute xscaling=xrefvals(1,2)-xrefvals(1,1).
        do if (xdich=1).
           compute xreferr=1.
           do if (((xrefvals(1,1) = xmx) and (xrefvals(1,2)=xmn)) or  ((xrefvals(1,1) = xmn) and (xrefvals(1,2)=xmx))).
             compute xreferr=0.
           end if.
           do if (xreferr=1).
             compute errcode(errs,1)=70.
             compute errs=errs+1.
             compute criterr=1.
           end if.         
        end if.
      end if.
    end if.
  end if.
  
  do if (xmint=1 and model=74 and mcx=0).
    compute xmodvals=t(xrefvals).
    compute xcontcf=1.
  end if.
  compute nxpval=nrow(xmodvals).
  compute xprobval=xmodvals.
  do if (xdich =1 and mcx > 0).
    compute mcx=0.
    compute errcode(errs,1) = 52.
    compute errs = errs + 1.
    compute criterr = 1.
  end if.
  do if (nms > 0).
    compute mtmp=dat(:,(nys+nxs+1):(nys+nxs+nms)).
    !describ descdatf=mtmp.
    compute ovsd={desctmp(2,:),ysd}.    
    compute medmeans=cdeval.
    do if ((cdeval(1,1) <> -999) and (ncol(medmeans) <> nms) and (model=74)).
      compute errcode(errs,1)=64.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((cdeval(1,1)=-999) and (model=74)).
      compute medmeans=desctmp(1,:).
    end if.
    do if ((cdeval(1,1)<> -999) and (model=74) and (criterr=0)).
      compute notecode(notes,1) = 31.
      compute notes = notes + 1.
    end if.
    do if ((rsum(desctmp(8,:))>0) and (mdichok <> 1)).
      compute errcode(errs,1)=43.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    compute mmodvals=modvals.
    compute mprobval=mmodvals.
  end if.
  do if (nws > 0).
    compute wtmp=dat(:,(nys+nxs+nms+1):(nys+nxs+nms+nws)).
    !describ descdatf=wtmp/type=(1-wprod).
    compute wmodvals=modvals.
    compute wdich=desctmp(8,1).
    do if (wdich =1 and mcw > 0).
      compute mcw=0.
      compute errcode(errs,1) = 52.
      compute errs = errs + 1.
      compute criterr = 1.
    end if.
    compute wmin=desctmp(3,1).
    compute wmax=desctmp(4,1).
    compute minwwarn=minwarn.
    compute maxwwarn=maxwarn.
    compute wnotev=mnotev.
    compute wmodval={!wmodval}.
    do if (xmint=1 and model=74 and mcx=0).
      compute wmodval=xrefvals.
    end if.
    compute nwcontr=ncol(wmodval).
    do if (wmodval(1,1) <> 999).
      compute wmodvals=wmodval(1,1).
      compute wmodcust=1.
      do if (nwcontr > 1).
        compute wmodvals=t(wmodval).
      end if.
      compute minwwarn=0.
      compute maxwwarn=0. 
      compute wnotev=0.
    end if.
    compute wprobval=wmodvals.
    compute nwpval=nrow(wmodvals).
  end if.
  do if (nzs > 0).
    compute ztmp=dat(:,(nys+nxs+nms+nws+1):(nys+nxs+nms+nws+nzs)).
    !describ descdatf=ztmp/type=(1-zprod).
    compute zmodvals=modvals.
    compute zdich=desctmp(8,1).
    do if (zdich =1 and mcz > 0).
      compute mcz=0.
      compute errcode(errs,1) = 52.
      compute errs = errs + 1.
      compute criterr = 1.
    end if.
    compute zmin=desctmp(3,1).
    compute zmax=desctmp(4,1).
    compute minzwarn=minwarn.
    compute maxzwarn=maxwarn.
    compute znotev=mnotev.
    compute zmodval={!zmodval}.
    compute nzcontr=ncol(zmodval).
    do if (zmodval(1,1) <> 999).
      compute zmodvals=zmodval(1,1).
      compute zmodcust=1.
      do if (nzcontr > 1).
        compute zmodvals=t(zmodval).
      end if.
      compute minzwarn=0.
      compute maxzwarn=0. 
      compute znotev=0.
    end if.  
    compute zprobval=zmodvals.
    compute nzpval=nrow(zmodvals).
  end if.
  do if (ncs > 0).
    compute ctmp=dat(:,(nys+nxs+nms+nws+nzs+1):(nys+nxs+nms+nws+nzs+ncs)).
    !describ descdatf=ctmp/type=1.
    compute covmeans=desctmp(1,:).
    compute coval={!coval}.
    do if ((coval(1,1) <> -999) and (ncol(coval) <> ncs) and (model=74)).
      compute errcode(errs,1)=69.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if ((coval(1,1)<> -999) and (criterr=0) and (model=74)).
      compute notecode(notes,1) = 35.
      compute notes = notes + 1.
      compute cuscoval=1.
    end if.    
  end if.
  compute n=nrow(ytmp).
  compute ones=make(n,1,1).
  /* create codes for categorical variables */.
  do if (nws > 0 and mcw > 0).
    compute tmp={rownum,wtmp(:,1)}.
    !makdumm dd=tmp/method=mcw/custcodv=2/custcode=!wcatcode.
    compute wmodvals=nnvls.
    compute nwpval=nrow(wmodvals).
    do if (criterr=0).
      compute minwwarn=0.
      compute maxwwarn=0.
      compute wnotev=0.
      compute wtmp=x(:,2:ncol(x)).
      compute wcatlab={"W1";"W2";"W3";"W4";"W5";"W6";"W7";"W8";"W9"}.
      do if (xmint=1).
        compute wcatlab={"X1";"X2";"X3";"X4";"X5";"X6";"X7";"X8";"X9"}.
      end if.
      compute nwvls=nvls-1.
      compute mcwok=1.
      compute dummatw=dummat.
      compute wprobval=dummatw(:,2:ncol(dummatw)).
      do if (modcok=1).
        compute wcontval=make(2,ncol(wprobval),-999).
        compute temp=0.
        loop i = 1 to 2.
          loop j = 1 to nrow(dummatw).
            do if (contvec(i,1)=dummatw(j,1)).
              compute wcontval(i,:)=wprobval(j,:).
              compute temp=temp+1.
            end if.
          end loop.
        end loop.
        do if (temp < 2).
          compute notecode(notes,1) = 20.
          compute notes = notes + 1.
          compute modcok=0.
        end if.
      end if.
      do if ((wmodval(1,1) <> 999) and (xmint <> 1)).
        compute notecode(notes,1) = 9.
        compute notes = notes + 1.
      end if.
      release tmp, dummat.
    end if.
  end if.
  do if (nzs > 0 and mcz > 0).
    compute tmp={rownum,ztmp(:,1)}.
    !makdumm dd=tmp/method=mcz/custcodv=3/custcode=!zcatcode.
    compute zmodvals=nnvls.
    compute nzpval=nrow(zmodvals).
    do if (criterr=0).
      compute minzwarn=0.
      compute maxzwarn=0.
      compute znotev=0.
      compute ztmp=x(:,2:ncol(x)).
      compute zcatlab={"Z1";"Z2";"Z3";"Z4";"Z5";"Z6";"Z7";"Z8";"Z9"}.
      compute nzvls=nvls-1.
      compute mczok=1.
      compute dummatz=dummat.
      compute zprobval=dummatz(:,2:ncol(dummatz)).
      do if (modcok=1).
        compute zcontval=make(2,ncol(zprobval),-999).
        compute temp=0.
        loop i = 1 to 2.
          loop j = 1 to nrow(dummatz).
            do if (contvec(i,2)=dummatz(j,1)).
              compute zcontval(i,:)=zprobval(j,:).
              compute temp=temp+1.
            end if.
          end loop.
        end loop.
        do if (temp < 2).
          compute notecode(notes,1) = 20.
          compute notes = notes + 1.
          compute modcok=0.
        end if.
      end if.
      do if (zmodval(1,1) <> 999).
        compute notecode(notes,1) = 10.
        compute notes = notes + 1.
      end if.
      release tmp, dummat.
    end if.
  end if.
  do if (nxs > 0 and mcx > 0).
    compute tmp={rownum,xtmp(:,1)}.
    !makdumm dd=tmp/method=mcx/custcodv=1/custcode=!xcatcode.
    do if (criterr=0).
      compute xtmp=x(:,2:ncol(x)).
      compute xcatlab={"X1";"X2";"X3";"X4";"X5";"X6";"X7";"X8";"X9"}.
      compute nxvls=nvls-1.
      compute setxcat=nxvls-1.
      compute xdich=(nvls=2).
      compute mcxok=1.
      compute dummatx=dummat.
      compute xmodvals=dummatx(:,1).
      compute nxpval=nrow(xmodvals).
      release tmp, dummat.
    end if.
  end if.
  /* making spline data */.
  do if (nspl > 0).
    compute spld3=make(n,nspl,0).
    compute spld2=make(n,nspl,0).
    compute xsplmin=cmin(xtmp).
    compute mcheck=0.
    compute ocheck=0.
    loop i = 1 to nspl.
      compute spld3(:,i)=xtmp > spl(1,i).
      compute spld2(:,i)=xtmp - spl(1,i).
      compute spld=spld3&*spld2.
      compute mcheck=mcheck+(spl(1,i)=xsplmin).
      do if (i > 1).
        do if (spl(1,i) <= spl(1,(i-1))).
          compute ocheck=1.
        end if.
      end if.
    end loop.
    compute spld4=csum(spld3).
    compute fff=rsum((csum(spld3) < 2) > 0).
    do if (ocheck=0).
      loop i = 1 to nspl.
        do if (i > 1).
          do if (spld4(1,i) >= ((spld4(1,i-1))-1)).
            compute fff=1.
          end if.
        end if.
      end  loop.
    end if.    
    do if (fff > 0).
       compute errs=errs+1.
       compute errcode(errs,1)=76.
       compute criterr=1.
     end if.
     do if (mcheck > 0).
       compute errs=errs+1.
       compute errcode(errs,1)=75.
       compute criterr=1.
     end if.
     do if (ocheck > 0).
       compute errs=errs+1.
       compute errcode(errs,1)=74.
       compute criterr=1.
     end if.
     do if (nspl > 8).
       compute errs=errs+1.
       compute errcode(errs,1)=73.
       compute criterr=1.
     end if.
    release spld2, spld3.
    compute xtmp={xtmp,spld}.        
    do if (nspl < 9).
      compute jlab={"Joint1";"Joint2";"Joint3";"Joint4";"Joint5";"Joint6";"Joint7";"Joint8"}.
      compute jlab=jlab(1:nspl,1).
      compute xcatlab={xnames;jlab}.
      compute nmspl=jlab.
      compute xnmspl=xnames.
      compute nxvls=nspl+1.
      compute setxcat=nxvls-1.
    end if.   
  end if.
end if. 
/* end A */.

end if.
/* endit */.

/* DEFINE AND CHECK COVARIATES MATRIX */.
do if (ncs > 0). 
  compute ccmat=make((nms+nys),ncs,1).
  compute ccmatoff=ccmat.
  do if (covmy=1).
    compute ccmat(nrow(ccmat),:)=make(1,ncs,0).
  end if.
  do if (covmy=2).
    compute ccmat(1:nms,:)=make(nms,ncs,0).
  end if.
  do if (cmatrix(1,1) <> -999).
    do if (ncol(cmatrix) <> ((nms+nys)*ncs)).
      compute errcode(errs,1)=29.
      compute errs=errs+1.
      compute criterr=1.
    end if.
    do if (criterr = 0).
      compute tmp=1.
      loop i = 1 to (nms+nys).
        loop j = 1 to ncs.
          compute ccmat(i,j)=1-(cmatrix(1,tmp) = 0).
          compute tmp=tmp+1.
        end loop.
      end loop.
      do if (rsum((csum(ccmat)=0)) <> 0).
        compute errcode(errs,1)=30.
        compute errs=errs+1.
        compute criterr=1.
      end if.
    end if.
    do if (covmy <> 0).
      compute notecode(notes,1)=1.
      compute notes=notes+1.
    end if.
  end if.
  do if (xmint=1 and (covmy > 0 or (csum(rsum(ccmat)) <> ((nms+nys)*ncs)))).
    compute errcode(errs,1)=68.
    compute errs=errs+1.
    compute criterr=1.
  end if.      
end if. 



do if (criterr=0).
  compute anymod2=csum(rsum(wcmat(:,1)+zcmat(:,1)+wzcmat(:,1))).
  compute anymod3=csum(rsum(wcmat+zcmat+wzcmat)).
  do if (((anymod3 > 0) or (nspl > 0)) and (yrelx = 1 or yrelm = 1 or yrelcov = 1)).
    compute errcode(errs,1)=82.
    compute errs=errs+1.
    compute criterr=1.
    compute eivbad=1.
  end if.
  do if ((anymod2 > 0) and nspl > 0).
    compute errcode(errs,1)=77.
    compute errs=errs+1.
    compute criterr=1.
  end if.
  do if ((savediag=1 and ydich=1) or (diagnose=1 and ydich=1)).
     compute notecode(notes,1) = 27.
     compute notes = notes + 1.
     compute savediag=0.
     compute diagnose=0.
  end if.
end if.




do if (criterr=0 and nms > 1).
  compute serchk=bcmat(2:(nrow(bcmat)-1),2:ncol(bcmat)).
  do if (csum(rsum(serchk))) > 0.
    compute serial=1.
    do if (nms > 6).
      compute errcode(errs,1)=36.
      compute errs=errs+1.
      compute criterr=1.
    end if.
  end if.
end if.


/* MEAN CENTER IF NEEDED */.
do if (center > 0 and criterr=0).
  compute centvar={" "}.
  do if (criterr=0).
    do if ((center = 1) or (center = 2 and wdich = 0)).
    do if (wprod=1 and mcwok=0 and nwpval > 0).
      loop i = 1 to nws.
        compute wtmp(:,i)=wtmp(:,i)-(csum(wtmp(:,i))/n).
        compute centvar={centvar,wnames(1,i)}.
      end loop.
      !describ descdatf=wtmp/type=wmodcust.
      compute wmin=desctmp(3,1).
      compute wmax=desctmp(4,1).
      do if (wmodcust=0).
        compute wmodvals=modvals.
        compute wprobval=wmodvals.
      end if.
    end if.
    end if.    
    do if ((center = 1) or (center = 2 and zdich = 0)). 
    do if (zprod=1 and mczok=0 and nzpval > 0).
      loop i = 1 to nzs.
        compute ztmp(:,i)=ztmp(:,i)-(csum(ztmp(:,i))/n).
        compute centvar={centvar,znames(1,i)}.
      end loop.
      !describ descdatf=ztmp/type=zmodcust.
      compute zmin=desctmp(3,1).
      compute zmax=desctmp(4,1).
      do if (zmodcust=0).
        compute zmodvals=modvals.
        compute zprobval=zmodvals.
      end if.
    end if.
    end if.
    do if ((center = 1) or (center = 2 and xdich = 0)).
    do if (xprod=1 and mcxok=0).
      loop i = 1 to nxs.
        compute xtmp(:,i)=xtmp(:,i)-(csum(xtmp(:,i))/n).
        compute centvar={centvar,xnames(1,i)}.
      end loop.
      !describ descdatf=xtmp.
        compute xmodvals=modvals.
        compute xprobval=xmodvals.
    end if.
    end if.
    do if (nms > 0).
      loop i = 1 to nms.
        do if (mprod(1,i)=1).
          compute mtmp(:,i)=mtmp(:,i)-(csum(mtmp(:,i))/n).
          compute centvar={centvar,mnames(1,i)}.
        end if.
      end loop.
      !describ descdatf=mtmp.
      do if ((cdeval(1,1)=-999) and (model=74)).
        compute medmeans=desctmp(1,:).
      end if.
    compute mmodvals=modvals.
    compute mprobval=mmodvals.
    end if.
  end if.
  do if (ncol(centvar) > 1).
    /* compute centvar=centvar(1,2:ncol(centvar)) */.
    compute notecode(notes,1)=3.
    compute notes=notes+1.
  end if.
end if.


/* start D */.
/* CONSTRUCT THE DATA MATRICES FOR EACH OF THE MODELS */.
do if (criterr=0).
  /* The i loop is the dependent variable in the model matrices */.
  /* The j loop is the predictor variables in the model matrices */.
  compute wsum=rsum(csum(wcmat)).
  compute zsum=rsum(csum(zcmat)).
  compute wzsum=rsum(csum(wzcmat)).
  compute nump=make(1,(nys+nms),-999).
  compue numint=make(1,(nys+nms),0).
  /* DV */.
  compute datcount=1.
  compute xtmpuse=0.
  compute wtmpuse=0.
  compute ztmpuse=0.
  compute xwtmpus=0.
  compute xztmpus=0.
  compute wztmpus=0.
  compute xwztmpu=0.
  compute xtmploc=-999.
  compute wtmploc=-999.
  compute xwtmplo=-999.
  compute ztmploc=-999.
  compute xztmplo=-999.
  compute wztmplo=-999.
  compute xwztmplo=-999.
  compute vlabs={" "}.
  do if (ncs > 0).
    compute ctmpuse=make(1,ncs,0).
  end if.
  do if (nms > 0).
    compute mtmpuse=make(1,nms,0).
    compute mwtmpus=make(1,nms,0).
    compute mztmpus=make(1,nms,0).
    compute mwztmpu=make(1,nms,0).
    compute mtmploc=make(1,nms,0).
    compute mwtmplo=make(nwvls,nms,-999).
    compute mztmplo=make(nzvls,nms,-999).
    compute mwztmplo=make((nwvls*nzvls),nms,-999).
  end if.
  do if (ncs > 0).
    compute ctmploc=make(1,ncs,0).
  end if.
  compute fulldat=make(n,1,1).
  compute fulldatn="xyz".
  compute datindx=make(1000,(nms+nys),-999).
  compute wherew=make(2,(nms+nys),-999).
  compute wherex=make(2,(nms+nys),-999).
  compute wherez=make(2,(nms+nys),-999).
  compute wherexw=make(2,(nms+nys),-999).
  compute wherexz=make(2,(nms+nys),-999).
  compute wherewz=make(2,(nms+nys),-999).
  compute wherexwz=make(2,(nms+nys),-999).
  do if (nms > 0).
    compute wherem=make(nms,(nms+nys),-999).
    compute wheremw = make(nms*2,(nms+nys),-999).
    compute wheremz = make(nms*2,(nms+nys),-999).
    compute wheremwz = make(nms*2,(nms+nys),-999).
    compute intbrdmw=make(nwvls,nms,"xxxx").
    compute intbrdmz=make(nzvls,nms,"xxxx").
    compute intbdmwz=make((nwvls*nzvls),nms,"xxxx").
  end if.
  compute wzhigh=make(1000,(((nms+1)*(nms+2))/2),0).
  compute whigh=make(1000,(((nms+1)*(nms+2))/2),0).
  compute zhigh=make(1000,(((nms+1)*(nms+2))/2),0).
  compute fochigh=make(1000,(((nms+1)*(nms+2))/2),0).
  compute xcoefloc={1;2;3;4;5;6;7;8;9}.
  compute intkey = {" ", " ", " ", " ", " ", " ", " "}.
  compute wzhighct=0.
  compute whighct=0.
  compute zhighct=0.   
  compute foccnt=0.
  compute cntmp2=1.
  loop i = 2 to nrow(bcmat). 
    compute wdid=0. 
    compute zdid=0.
    compute wzdid=0.
    compute cntmp=1.
    compute start=1.
    do if (i < nrow(bcmat)).
      compute outv=mtmp(:,(i-1)).
      compute modlabel={mnames(1,(i-1));"constant"}.
    end if.
    do if (i = nrow(bcmat)).
      compute outv=ytmp.
      compute modlabel={ynames;"constant"}.
    end if.
    /* The j loop is the mediator */.
    /* MED */.
    /* do if (nms > 0) */.
    loop j = 1 to (i-1).
      compute foccnt=foccnt+1.
      do if (j = 1 and bcmat(i,j)=1).
        compute outv={outv,xtmp}.
        compute modlabel={modlabel;xcatlab(1:nxvls,1)}.
        do if (xtmpuse=0).
          compute fulldat={fulldat,xtmp}.
          compute relxmc=make(1,nxvls,relx).
          compute fulldatr={fulldatr,relxmc}.
          compute fulldatn={fulldatn,t(xcatlab(1:nxvls,1))}.
          compute xtmpuse=1.
          loop k4=datcount to (datcount+(nxvls-1)).
            compute xtmploc={xtmploc;k4}.
          end loop.
          compute xtmploc=xtmploc(2:nrow(xtmploc),1).
          compute datcount=datcount+nxvls.
        end if.
        compute datindx(start:(start+nrow(xtmploc)-1),(i-1))=xtmploc.
        compute wherex(1,(i-1))=start+1.
        compute wherex(2,(i-1))=start+nrow(xtmploc)-1+1.
        do if (model = 74).
        end if.
        compute onebl=make(nrow(xtmploc),1,1).
        compute fochigh((start+1):(start+nrow(xtmploc)),foccnt)=onebl.
        compute start=start+nrow(xtmploc).
      end if.
      do if (j > 1 and bcmat(i,j)=1).
        compute outv={outv,mtmp(:,(j-1))}.
        compute modlabel={modlabel;mnames(1,(j-1))}.
        do if (mtmpuse(1,(j-1))=0).
          compute fulldat={fulldat,mtmp(:,(j-1))}.
          compute fulldatr={fulldatr,relm(1,(j-1))}.
          compute fulldatn={fulldatn,mnames(1,(j-1))}.
          compute mtmpuse(1,(j-1))=1.
          compute mtmploc(1,(j-1))=datcount.
          compute datcount=datcount+1.
        end if.
        compute datindx(start:(start+nrow(mtmploc)-1),(i-1))=mtmploc(1,(j-1)).
        compute wherem((j-1),(i-1))=start+1.
        compute onebl=make(nrow(mtmploc(1,j-1)),1,1).
        compute ttt=nrow(mtmploc(1,(j-1)))+start-1.
        compute fochigh((start+1):(start+nrow(mtmploc(1,(j-1)))),foccnt)=onebl.
        compute start=start+nrow(mtmploc(1,(j-1))).
      end if.
    end loop.
    /* end if */.
    /* END MED */.
    /* This starts the W loop */.
    /* W */.
    do if (wsum > 0).
      loop j = 1 to (i-1).
      compute whighct=whighct+1.
        do if (j = 1 and wcmat(i,j)=1).
          do if (wdid=0).
            compute outv={outv,wtmp}.
            do if (ncs > 0 and wiscov > 0).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            compute modlabel={modlabel;wcatlab(1:nwvls,1)}.
            compute wdid=1.
            do if (wtmpuse=0).
              compute fulldat={fulldat,wtmp}.
              compute fulldatn={fulldatn,t(wcatlab(1:nwvls,1))}.
              do if (ncs > 0 and wiscov > 0).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              compute wtmpuse=1.
              loop k4=datcount to (datcount+(nwvls-1)).
                compute wtmploc={wtmploc;k4}.
              end loop.
              compute wtmploc=wtmploc(2:nrow(wtmploc),1).
              compute datcount=datcount+nwvls.
            end if.
          end if.
          compute datindx(start:(start+nrow(wtmploc)-1),(i-1))=wtmploc.
          compute wherew(1,(i-1))=start+1.
          compute wherew(2,(i-1))=start+nrow(wtmploc)-1+1.
          compute start=start+nrow(wtmploc).
          compute intlbrd="xxx".
          compute cntmp1=1.
          loop k1=1 to nxvls.
            loop k2 = 1 to nwvls.
              compute outv={outv,(xtmp(:,k1)&*wtmp(:,k2))}.
              do if (ncs > 0 and wiscov > 0).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              do if (xwtmpus=0).
                compute modlabel={modlabel;intlab(cntmp2,1)}.  
                compute intkey={intkey;intlab(cntmp2,1),":",xcatlab(k1,1),"x",wcatlab(k2,1)," "," "}.
                compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                compute cntmp=cntmp+1.
                compute cntmp2=cntmp2+1.
              end if.    
              do if (xwtmpus=1).
                compute modlabel={modlabel;intbrdxw(cntmp1,1)}.  
                compute intkey={intkey;intbrdxw(cntmp1,1),":",xcatlab(k1,1),"x",wcatlab(k2,1)," "," "}.
                compute cntmp1=cntmp1+1.
                compute cntmp=cntmp+1.
              end if.
            end loop.
          end loop.
          do if (xwtmpus=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-(nxvls*nwvls)+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
            compute intbrdxw=t(intlbrd(1,2:ncol(intlbrd))).
            compute xwtmpus=1.
            do if (ncs > 0 and wiscov > 0).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            loop k4=datcount to (datcount+((nwvls*nxvls)-1)).
              compute xwtmplo={xwtmplo;k4}.
            end loop.
            compute xwtmplo=xwtmplo(2:nrow(xwtmplo),1).
            compute datcount=datcount+(nxvls*nwvls).
          end if.
          compute datindx(start:(start+nrow(xwtmplo)-1),(i-1))=xwtmplo.
          compute wherexw(1,(i-1))=start+1.
          compute wherexw(2,(i-1))=start+nrow(xwtmplo)-1+1.
          compute onebl=make(nrow(xwtmplo),1,1).
          compute whigh((start+1):(start+nrow(xwtmplo)),whighct)=onebl.
          compute start=start+nrow(xwtmplo).
        end if.
        do if (j > 1 and wcmat(i,j)=1).
          do if (wdid=0 and model <> 74).
            compute outv={outv,wtmp}.
            do if (ncs > 0 and wiscov > 0).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            compute modlabel={modlabel;wcatlab(1:nwvls,1)}.
            compute wdid=1.
            do if (wtmpuse=0).
              compute fulldat={fulldat,wtmp}.
              compute fulldatn={fulldatn,t(wcatlab(1:nwvls,1))}.
              do if (ncs > 0 and wiscov > 0).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              compute wtmpuse=1.
              loop k4=datcount to (datcount+(nwvls-1)).
                compute wtmploc={wtmploc;k4}.
              end loop.
              compute wtmploc=wtmploc(2:nrow(wtmploc),1).
              compute datcount=datcount+nwvls.
            end if.
            compute datindx(start:(start+nrow(wtmploc)-1),(i-1))=wtmploc.
            compute wherew(1,(i-1))=start+1.
            compute wherew(2,(i-1))=start+nrow(wtmploc)-1+1.
            compute start=start+nrow(wtmploc).
          end if.         
          compute intlbrd="xxx".
          compute cntmp1=1.
          loop k2 = 1 to nwvls.
            compute outv={outv,(mtmp(:,(j-1))&*wtmp(:,k2))}.
            do if (ncs > 0 and wiscov > 0).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            do if (mwtmpus(1,(j-1))=0).
              compute modlabel={modlabel;intlab(cntmp2,1)}.
              compute intkey={intkey;intlab(cntmp2,1),":", mnames(1,(j-1)),"x",wcatlab(k2,1)," "," "}.
              compute intlbrd={intlbrd,intlab(cntmp2,1)}.
              compute cntmp=cntmp+1.
              compute cntmp2=cntmp2+1.
            end if.
            do if (mwtmpus(1,(j-1))=1).
                compute modlabel={modlabel;intbrdmw(cntmp1,(j-1))}.
                compute intkey={intkey;intbrdmw(cntmp1,(j-1)),":", mnames(1,(j-1)),"x",wcatlab(k2,1)," "," "}.
                compute cntmp1=cntmp1+1.
                compute cntmp=cntmp+1.
            end if.
          end loop.
          do if (mwtmpus(1,(j-1))=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-nwvls+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
            compute intbrdmw(1:(ncol(intlbrd)-1),(j-1))=t(intlbrd(1,2:ncol(intlbrd))).
            do if (ncs > 0 and wiscov > 0).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            compute mwtmpus(1,(j-1))=1.
            compute mw22=-999.
            loop k4=datcount to (datcount+(nwvls-1)).
              compute mw22={mw22;k4}.
            end loop.
            compute mwtmplo(:,(j-1))=mw22(2:nrow(mw22),1).
            compute datcount=datcount+nwvls.
          end if.
          compute datindx(start:(start+nrow(mwtmplo)-1),(i-1))=mwtmplo(:,(j-1)).
          compute wheremw(((2*j)-3),(i-1))=start+1.
          compute wheremw(((2*j)-2),(i-1))=start+nrow(mwtmplo)-1+1.
          compute onebl=make(nrow(mwtmplo),1,1).
          compute whigh((start+1):(start+nrow(mwtmplo)),whighct)=onebl.
          compute start=start+nrow(mwtmplo).
        end if.
      end loop.
    end if.
    /* END W */.
    /* This starts the Z loop */.
    /* Z */.
    do if (zsum > 0).
      loop j = 1 to (i-1).
        compute zhighct=zhighct+1.
        do if (j = 1 and zcmat(i,j)=1).
          do if (zdid=0).
            compute outv={outv,ztmp}.
            do if (ncs > 0 and ziscov > 0).
              compute ccmatoff((i-1),ziscov)=0.
            end if.
            compute modlabel={modlabel;zcatlab(1:nzvls,1)}.
            compute zdid=1.
            do if (ztmpuse=0).
              compute fulldat={fulldat,ztmp}.
              compute fulldatn={fulldatn,t(zcatlab(1:nzvls,1))}.
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
              compute ztmpuse=1.
              loop k4=datcount to (datcount+(nzvls-1)).
                compute ztmploc={ztmploc;k4}.
              end loop.
              compute ztmploc=ztmploc(2:nrow(ztmploc),1).
              compute datcount=datcount+nzvls.
            end if.
          end if.
          compute datindx(start:(start+nrow(ztmploc)-1),(i-1))=ztmploc.
          compute wherez(1,(i-1))=start+1.
          compute wherez(2,(i-1))=start+nrow(ztmploc)-1+1.
          compute start=start+nrow(ztmploc).
          comput intlbrd="xxx".
          compute cntmp1=1.
          loop k1=1 to nxvls.
            loop k2 = 1 to nzvls.
              compute outv={outv,(xtmp(:,k1)&*ztmp(:,k2))}.
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
              do if (xztmpus=0).
                compute modlabel={modlabel;intlab(cntmp2,1)}.
                compute intkey={intkey;intlab(cntmp2,1),":",xcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                compute cntmp=cntmp+1.
                compute cntmp2=cntmp2+1.
              end if.
              do if (xztmpus=1).
                compute modlabel={modlabel;intbrdxz(cntmp1,1)}.
                compute intkey={intkey;intbrdxz(cntmp1,1),":",xcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                compute cntmp1=cntmp1+1.
                compute cntmp=cntmp+1.
              end if.
            end loop.
          end loop.
          do if (xztmpus=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-(nxvls*nzvls)+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
            compute intbrdxz=t(intlbrd(1,2:ncol(intlbrd))).
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
            compute xztmpus=1.
            loop k4=datcount to (datcount+((nzvls*nxvls)-1)).
              compute xztmplo={xztmplo;k4}.
            end loop.
            compute xztmplo=xztmplo(2:nrow(xztmplo),1).
            compute datcount=datcount+(nxvls*nzvls).
          end if.
          compute datindx(start:(start+nrow(xztmplo)-1),(i-1))=xztmplo.
          compute wherexz(1,(i-1))=start+1.
          compute wherexz(2,(i-1))=start+nrow(xztmplo)-1+1.
          compute onebl=make(nrow(xztmplo),1,1).
          compute zhigh((start+1):(start+nrow(xztmplo)),zhighct)=onebl.
          compute start=start+nrow(xztmplo).
        end if.
        do if (j > 1 and zcmat(i,j)=1).
          do if (zdid=0).
            compute outv={outv,ztmp}.
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
            compute modlabel={modlabel;zcatlab(1:nzvls,1)}.
            compute zdid=1.
            do if (ztmpuse=0).
              compute fulldat={fulldat,ztmp}.
              compute fulldatn={fulldatn,t(zcatlab(1:nzvls,1))}.
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
              compute ztmpuse=1.
              loop k4=datcount to (datcount+(nzvls-1)).
                compute ztmploc={ztmploc;k4}.
              end loop.
              compute ztmploc=ztmploc(2:nrow(ztmploc),1).
              compute datcount=datcount+nzvls.
            end if.
            compute datindx(start:(start+nrow(ztmploc)-1),(i-1))=ztmploc.
            compute wherez(1,(i-1))=start+1.
            compute wherez(2,(i-1))=start+nrow(ztmploc)-1+1.
            compute start=start+nrow(ztmploc).
          end if.    
          compute intlbrd="xxx".
          compute cntmp1=1.
          loop k2 = 1 to nzvls.
            compute outv={outv,(mtmp(:,(j-1))&*ztmp(:,k2))}.
              do if (ncs > 0 and ziscov > 0).
                compute ccmatoff((i-1),ziscov)=0.
              end if.
            do if (mztmpus(1,(j-1))=0).
              compute modlabel={modlabel;intlab(cntmp2,1)}.
              compute intkey={intkey;intlab(cntmp2,1),":", mnames(1,(j-1)),"x",zcatlab(k2,1)," "," "}.
              compute intlbrd={intlbrd,intlab(cntmp2,1)}.
              compute cntmp=cntmp+1.
              compute cntmp2=cntmp2+1.
              end if.
            do if (mztmpus(1,(j-1))=1).
              compute modlabel={modlabel;intbrdmz(cntmp1,(j-1))}.
              compute intkey={intkey;intbrdmz(cntmp1,(j-1)),":", mnames(1,(j-1)),"x",zcatlab(k2,1)," "," "}.
              compute cntmp1=cntmp1+1.
              compute cntmp=cntmp+1.
            end if.
          end loop.
          do if (mztmpus(1,(j-1))=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-nzvls+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
             compute intbrdmz(1:(ncol(intlbrd)-1),(j-1))=t(intlbrd(1,2:ncol(intlbrd))).
            do if (ncs > 0 and ziscov > 0).
              compute ccmatoff((i-1),ziscov)=0.
            end if.
            compute mztmpus(1,(j-1))=1.
            compute mz22=-999.
            loop k4=datcount to (datcount+(nzvls-1)).
              compute mz22={mz22;k4}.
            end loop.
            compute mztmplo(:,(j-1))=mz22(2:nrow(mz22),1).
            compute datcount=datcount+nzvls.
          end if.
          compute datindx(start:(start+nrow(mztmplo)-1),(i-1))=mztmplo(:,(j-1)).
          compute wheremz(((2*j)-3),(i-1))=start+1.
          compute wheremz(((2*j)-2),(i-1))=start+nrow(mztmplo)-1+1.
          compute onebl=make(nrow(mztmplo),1,1).
          compute zhigh((start+1):(start+nrow(mztmplo)),zhighct)=onebl.
          compute start=start+nrow(mztmplo).
        end if.
      end loop.
    end if.
    /* END Z */.
    /* This starts the WZ loop */.
    /* WZ */.
    do if (wzsum > 0).  
      loop j = 1 to (i-1).
        compute wzhighct=wzhighct+1.
        do if (j = 1 and wzcmat(i,j)=1).
          do if (wzdid=0).
            compute intlbrd="xxx".
            compute cntmp1=1.
            loop k1=1 to nwvls.
              loop k2 = 1 to nzvls.
                compute outv={outv,(wtmp(:,k1)&*ztmp(:,k2))}. 
                 do if (ncs > 0 and (ziscov > 0)).
                   compute ccmatoff((i-1),ziscov)=0.
                 end if.  
                 do if (ncs > 0 and (wiscov > 0)).
                   compute ccmatoff((i-1),wiscov)=0.
                 end if.
                do if (wztmpus=0).
                  compute modlabel={modlabel;intlab(cntmp2,1)}.
                  compute intkey={intkey;intlab(cntmp2,1),":",wcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                  compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                  compute cntmp=cntmp+1.
                  compute cntmp2=cntmp2+1.
                end if.
                do if (wztmpus=1).
                  compute modlabel={modlabel;intbrdwz(cntmp1,1)}.
                  compute intkey={intkey;intbrdwz(cntmp1,1),":",wcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                  compute cntmp1=cntmp1+1.
                  compute cntmp=cntmp+1.
                end if.
              end loop.
            end loop.
            do if (wztmpus=0).
              compute fulldat={fulldat,outv(:,(ncol(outv)-(nwvls*nzvls)+1):ncol(outv))}.
              compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
              compute intbrdwz=t(intlbrd(1,2:ncol(intlbrd))).
              do if (ncs > 0 and (ziscov > 0)).
                compute ccmatoff((i-1),ziscov)=0.
              end if.  
              do if (ncs > 0 and (wiscov > 0)).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              compute wztmpus=1.
              loop k4=datcount to (datcount+((nwvls*nzvls)-1)).
                compute wztmplo={wztmplo;k4}.
              end loop.
              compute wztmplo=wztmplo(2:nrow(wztmplo),1).
              compute datcount=datcount+(nzvls*nwvls).
            end if.
            compute wzdid=1.         
          end if.
          compute datindx(start:(start+nrow(wztmplo)-1),(i-1))=wztmplo.
          compute wherewz(1,(i-1))=start+1.
          compute wherewz(2,(i-1))=start+nrow(wztmplo)-1+1.
          compute start=start+nrow(wztmplo).
          compute intlbrd="xxx".
          compute cntmp1=1.
          loop k1=1 to nxvls.
            loop k2=1 to nwvls.
              loop k3=1 to nzvls.
                compute outv={outv,(xtmp(:,k1)&*wtmp(:,k2)&*ztmp(:,k3))}. 
                do if (ncs > 0 and (ziscov > 0)).
                 compute ccmatoff((i-1),ziscov)=0.
                end if.  
                do if (ncs > 0 and (wiscov > 0)).
                 compute ccmatoff((i-1),wiscov)=0.
                end if.    
                do if (xwztmpu=0).
                  compute modlabel={modlabel;intlab(cntmp2,1)}.
                  compute intkey={intkey;intlab(cntmp2,1),":",xcatlab(k1,1),"x",wcatlab(k2,1),"x",zcatlab(k3,1)}.
                  compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                  compute cntmp=cntmp+1. 
                  compute cntmp2=cntmp2+1.
                end if.
                do if (xwztmpu=1).
                  compute modlabel={modlabel;intbdxwz(cntmp1,1)}.
                  compute intkey={intkey;intbdxwz(cntmp1,1),":",xcatlab(k1,1),"x",wcatlab(k2,1),"x",zcatlab(k3,1)}.
                  compute cntmp1=cntmp1+1. 
                  compute cntmp=cntmp+1.
                end if.                
              end loop.
            end loop.
          end loop.
          do if (xwztmpu=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-(nxvls*nwvls*nzvls)+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
            compute intbdxwz=t(intlbrd(1,2:ncol(intlbrd))).
            do if (ncs > 0 and (ziscov > 0)).
              compute ccmatoff((i-1),ziscov)=0.
            end if.  
            do if (ncs > 0 and (wiscov > 0)).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            compute xwztmpu=1.
            loop k4=datcount to (datcount+((nzvls*nxvls*nwvls)-1)).
              compute xwztmplo={xwztmplo;k4}.
            end loop.
            compute xwztmplo=xwztmplo(2:nrow(xwztmplo),1).
            compute datcount=datcount+(nxvls*nzvls*nwvls).
          end if.
          compute datindx(start:(start+nrow(xwztmplo)-1),(i-1))=xwztmplo.
          compute wherexwz(1,(i-1))=start+1.
          compute wherexwz(2,(i-1))=start+nrow(xwztmplo)-1+1.
          compute onebl=make(nrow(xwztmplo),1,1).
          compute wzhigh((start+1):(start+nrow(xwztmplo)),wzhighct)=onebl.
          compute start=start+nrow(xwztmplo).
        end if.
        do if (j > 1 and wzcmat(i,j)=1).
          do if (wzdid=0).
            compute intlbrd="xxx".
            compute cntmp1=1.
            loop k1=1 to nwvls.
              loop k2 = 1 to nzvls.
                compute outv={outv,(wtmp(:,k1)&*ztmp(:,k2))}.
                do if (ncs > 0 and (ziscov > 0)).
                  compute ccmatoff((i-1),ziscov)=0.
                end if.  
                do if (ncs > 0 and (wiscov > 0)).
                  compute ccmatoff((i-1),wiscov)=0.
                end if.
                do if (wztmpus=0).
                  compute modlabel={modlabel;intlab(cntmp2,1)}.
                  compute intkey={intkey;intlab(cntmp2,1),":",wcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                  compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                  compute cntmp=cntmp+1.
                  compute cntmp2=cntmp2+1.
                end if.
                do if (wztmpus=1).
                  compute modlabel={modlabel;intbrdwz(cntmp1,1)}.
                  compute intkey={intkey;intbrdwz(cntmp1,1),":",wcatlab(k1,1),"x",zcatlab(k2,1)," "," "}.
                  compute cntmp1=cntmp1+1.
                  compute cntmp=cntmp+1.
                end if.
              end loop.
            end loop.
            do if (wztmpus=0).
              compute fulldat={fulldat,outv(:,(ncol(outv)-(nwvls*nzvls)+1):ncol(outv))}.
              compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
              compute intbrdwz=t(intlbrd(2:ncol(intlbrd))).
              do if (ncs > 0 and (ziscov > 0)).
                compute ccmatoff((i-1),ziscov)=0.
              end if.  
              do if (ncs > 0 and (wiscov > 0)).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              compute wztmpus=1.
              loop k4=datcount to (datcount+((nwvls*nzvls)-1)).
                compute wztmplo={wztmplo;k4}.
              end loop.
              compute wztmplo=wztmplo(2:nrow(wztmplo),1).
              compute datcount=datcount+(nzvls*nwvls).
            end if.
            compute wzdid=1.
            compute datindx(start:(start+nrow(wztmplo)-1),(i-1))=wztmplo.
            compute wherewz(1,(i-1))=start+1.
            compute wherewz(2,(i-1))=start+nrow(wztmplo)-1+1.
            compute start=start+nrow(wztmplo).
          end if.
          compute intlbrd="xxx".
          compute cntmp1=1.
          loop k1 = 1 to nwvls.
            loop k2 = 1 to nzvls.
              compute outv={outv,(mtmp(:,(j-1))&*wtmp(:,k1)&*ztmp(:,k2))}.
              do if (ncs > 0 and (ziscov > 0)).
                compute ccmatoff((i-1),ziscov)=0.
              end if.  
              do if (ncs > 0 and (wiscov > 0)).
                compute ccmatoff((i-1),wiscov)=0.
              end if.
              do if (mwztmpu(1,(j-1))=0).
                compute modlabel={modlabel;intlab(cntmp2,1)}.
                compute intkey={intkey;intlab(cntmp2,1),":",mnames(1,(j-1)),"x",wcatlab(k1,1),"x",zcatlab(k2,1)}.
                compute intlbrd={intlbrd,intlab(cntmp2,1)}.
                compute cntmp=cntmp+1.
                compute cntmp2=cntmp2+1.
              end if.
              do if (mwztmpu(1,(j-1))=1).
                compute modlabel={modlabel;intbdmwz(cntmp1,(j-1))}.
                compute intkey={intkey;intbdmwz(cntmp1,(j-1)),":",mnames(1,(j-1)),"x",wcatlab(k1,1),"x",zcatlab(k2,1)}.
                compute cntmp1=cntmp1+1.
                compute cntmp=cntmp+1.
              end if.
            end loop.
          end loop.
          do if (mwztmpu(1,(j-1))=0).
            compute fulldat={fulldat,outv(:,(ncol(outv)-(nwvls*nzvls)+1):ncol(outv))}.
            compute fulldatn={fulldatn,intlbrd(2:ncol(intlbrd))}.
            compute intbdmwz(1:(ncol(intlbrd)-1),(j-1))=t(intlbrd(1,2:ncol(intlbrd))).
            do if (ncs > 0 and (ziscov > 0)).
              compute ccmatoff((i-1),ziscov)=0.
            end if.  
            do if (ncs > 0 and (wiscov > 0)).
              compute ccmatoff((i-1),wiscov)=0.
            end if.
            compute mwztmpu(1,(j-1))=1.
            compute mz22=-999.
            loop k4=datcount to (datcount+(nwvls*nzvls)-1).
              compute mz22={mz22;k4}.
            end loop.
            compute mwztmplo(:,(j-1))=mz22(2:nrow(mz22),1).
            compute datcount=datcount+(nwvls*nzvls).
          end if.
          compute datindx(start:(start+nrow(mwztmplo)-1),(i-1))=mwztmplo(:,(j-1)).
          compute wheremwz(((2*j)-3),(i-1))=start+1.
          compute wheremwz(((2*j)-2),(i-1))=start+nrow(mwztmplo)-1+1.
          compute onebl=make(nrow(mwztmplo),1,1).
          compute wzhigh((start+1):(start+nrow(mwztmplo)),wzhighct)=onebl.
          compute start=start+nrow(mwztmplo).
        end if.
      end loop.
    end if.
    /* END WZ */.
    /* This starts the COV loop */.
    /* COV */.
    do if (ncs > 0).
      compute ccmat=ccmat&*ccmatoff.
      loop j = 1 to ncs.
        do if (ccmat((i-1),j))=1.
           do if (j=wiscov).
             compute ctmp(:,j)=wtmp.
           end if.
           do if (j=ziscov).
             compute ctmp(:,j)=ztmp.
           end if.
           compute outv={outv,ctmp(:,j)}.
           compute modlabel={modlabel;covnames(1,j)}.
           do if (ctmpuse(1,j)=0).
             compute fulldat={fulldat,ctmp(:,j)}.
             compute fulldatr={fulldatr,relcov(1,j)}.
             compute fulldatn={fulldatn,covnames(1,j)}.
             compute ctmpuse(1,j)=1.
             compute ctmploc(1,j)=datcount.
             compute datcount=datcount+1.
           end if.
           compute datindx(start:(start+nrow(ctmploc)-1),(i-1))=ctmploc(1,j).
           compute start=start+nrow(ctmploc(1,j)).
        end if.
      end loop.
    end if.
    /* END COV */.
    compute wdid=0.
    compute zdid=0.
    compute wzdid=0.
    compute vlabs={vlabs;modlabel(2:nrow(modlabel),1)}.
    compute numint(1,(i-1))=cntmp-1.
    compute nump(1,(i-1))=nrow(modlabel)-1.
  end loop.
  /* END DV*/.
  



  release datcount, xtmpuse, wtmpuse, ztmpuse, xwtmpus, xztmpus, wztmpus, xwztmpu.
  release xtmploc, wtmploc, xwtmplo, ztmploc, xztmplo, wztmplo, xwztmplo, foccnt.
  do if (modcok=1 and ((nms > 0) or (zcmat(2,1) <> 1) or (mcx <> 0))).
    compute notecode(notes,1) = 19.
    compute notes = notes + 1.
    compute modcok=0.
  end if.
  do if ((serial = 1 or (rsum(numint)>0) or nms=0) and mc > 0).
    compute notecode(notes,1) = 15.
    compute notes = notes + 1.
    compute boot=mc.
    compute mc=0.
  end if.
  do if (boot <> 0 or mc <> 0).
    compute bootsz=boot.
    do if (mc > 0).
      compute bootsz=mc.
      compute saveboot=0.
    end if.
    loop.
      compute cilow = rnd(bootsz*(1-(conf/100))/2).
      compute cihigh = trunc((bootsz*(conf/100)+(bootsz*(1-(conf/100))/2)))+1.
      do if (cilow < 1 or cihigh > bootsz).
        compute bootsz=trunc((bootsz+1000)/1000)*1000.
        compute adjust = 1.
      end if.
    end loop if (cilow gt 0 and cihigh le bootsz).
    do if (boot > 0).
      compute boot=bootsz.
    end if.
    do if (mc > 0).
      compute mc=bootsz.
    end if.
    do if (adjust = 1 and boot > 0).
      compute notecode(notes,1) = 8.
      compute notes = notes + 1.
    end if.
    do if (adjust = 1 and mc > 0).
      compute notecode(notes,1) = 16.
      compute notes = notes + 1.
    end if.
  end if.
  compute maxboot = trunc(2*boot).
  do if (!maxboot > maxboot).
    compute maxboot=trunc(!maxboot).
  end if.
  do if (nms > 0).
    release mtmpuse, mwtmpus, mwztmpu, mtmploc, mwtmplo, mztmplo, mwztmplo.
  end if.
  /* release wdid, zdid, wzdid, start,modlabel */.
  release wdid, zdid, wzdid, start.
  compute vlabs=vlabs(2:nrow(vlabs),1).
  do if (rsum(numint) > 0).
    compute intkey=intkey(2:nrow(intkey),:).
  end if.
  compute fulldat=fulldat(:,2:ncol(fulldat)).
  compute fulldatr=fulldatr(1,2:ncol(fulldatr)).
  compute fochigh=fochigh(1:rmax(nump),:).
  compute whigh=whigh(1:rmax(nump),:).
  compute zhigh=zhigh(1:rmax(nump),:).
  compute wzhigh=wzhigh(1:rmax(nump),:).
  compute coeffs=fochigh+whigh+zhigh+wzhigh.
  compute bootloc=make(rmax(nump),ncol(nump),0).
  /* here i am deriving the locations in boot file needed for indirect effects */.
  do if (nms > 0).
    compute cntmp=1.
    loop i = 1 to ncol(nump).
      loop j = 1 to nump(1,i).
        compute bootloc(j,i)=cntmp.
        compute cntmp=cntmp+1.
      end loop.
    end loop.
    compute fochighb=make(nrow(fochigh),ncol(fochigh),0).
    compute whighb=fochighb.
    compute zhighb=fochighb.
    compute wzhighb=fochighb.
    compute thetaxmb=make(nrow(fochighb),nms,0).
    compute thetaxyb=make(nrow(fochighb),1,0).
    compute pathsfoc=make(nxvls,1,0).
    compute cntmp=1.
    loop i = 1 to (nms+nys).
      loop j = 1 to i.
        compute fochighb(:,cntmp)=fochigh(:,cntmp)&*bootloc(:,i).
        compute whighb(:,cntmp)=whigh(:,cntmp)&*bootloc(:,i).
        compute zhighb(:,cntmp)=zhigh(:,cntmp)&*bootloc(:,i).
        compute wzhighb(:,cntmp)=wzhigh(:,cntmp)&*bootloc(:,i).
        compute coeffsb=fochighb+whighb+zhighb+wzhighb.
        do if ((i < (nms+nys)) and (j = 1)).
          compute thetaxmb(:,i)=coeffsb(:,cntmp).
        end if.
        do if ((i = (nms+nys)) and (j = 1)).
          compute thetaxyb(:,1)=coeffsb(:,cntmp).
        end if.
        compute cntmp=cntmp+1.
      end loop.
    end loop.
    compute thetamyb=coeffsb(:,(ncol(coeffsb)-nms+1):ncol(coeffsb)).
    do if (serial = 1).
      compute thetammb=make(nrow(coeffsb),((nms*(nms-1))/2),0).
    end if.
    compute cntmp=1.
    do if (nms > 1 and serial = 1).
      loop i = 1 to (nms-1).
        compute start=((i+2)*(i+1))/2.
        loop j = 2 to (nms-i+1).
          compute thetammb(:,cntmp)=coeffsb(:,start).
          compute start=start+j+i-1.
          compute cntmp=cntmp+1.
        end loop.
      end loop. 
    end if.
  end if.
  do if ((total = 1) and (rsum(numint)=0 or (xmint=1))).
    compute dototal=1.
    do if ((csum(bcmat(:,1)) <> (nms+nys)) or (rsum(bcmat(nrow(bcmat),:)) <> (nms+nys))).
      compute dototal=0.
      compute alttotal=1.
      compute notecode(notes,1) = 12.
      compute notes = notes + 1.
    end if.
    do if (ncs > 0).
      do if ((csum(rsum(ccmat))) < (nrow(ccmat)*ncol(ccmat))).
        compute dototal=0.
        compute alttotal=1.
        compute notecode(notes,1) = 11.
        compute notes = notes + 1.
      end if.
    end if.
    do if (model=74).
      do if (xdich=0 and nxvls=1).
        compute dototal=0.
        compute alttotal=1.
      end if.
      do if ((xdich=1 or nxvls > 1) and (ncs > 0) and model=74).
        compute dototal=0.
        compute alttotal=1.
      end if.
    end if.
  end if.
  
do if (eivbad=0).
  compute reltmp={relx,relm,relcov}.
  compute reliabck=(rsum(reltmp > 1)+rsum(reltmp <=0)).
  do if (reliabck <> 0).
    compute errcode(errs,1)=84.
    compute errs=errs+1.
    compute criterr=1.
    compute eivbad=1.
  end if.
end if.  
  
do if (eivdo=1).
  /* what to do about spline models is a question */.
  do if (nxvls > 1).
    compute mcxrel=make(1,(nxvls-1),relx(1,1)).
    compute relx={mcxrel,relx}.
    compute eivrel={mcxrel,eivrel}.
  end if.
  /* compute relx={1,relx} */.
  /* compute eivrel={1,eivrel} */.
  compute ssquares=0.
  compute crossv=0.
  compute stand=0.
  compute effsize=0.
  compute subsets=0.
  compute dominate=0.
  compute diagnose=0.
  do if (robustse=1 or hc <> 5).
    compute robustse=0.
    compute hc=5.
    compute errcode(errs,1)=87.
    compute errs=errs+1.
    compute criterr=1.
  end if.
  compute modelres=0.
  /* compute linsum=-999 */.
  /* compute nlinsum=0 */.
  compute savediag=0.
end if.
  
  
end if.
/* END D */.

do if (criterr=0 and ncs > 0).
  do if (rsum((csum(ccmat)=0)) <> 0).
    compute errcode(errs,1)=51.
    compute errs=errs+1.
    compute criterr=1.
  end if.
end if.
compute debug=!debug.

do if (outscree=1).
  print/title = "***************** PROCESS Procedure for SPSS Version 5.0 *****************".
  print/title = "          Written by Andrew F. Hayes, Ph.D.       www.afhayes.com".  
  print/title = "    Documentation available in Hayes (2022). www.guilford.com/p/hayes3"/space=0.
end if.
do if (criterr=0).
  compute modresid=make(n,1,99999).
  compute alldfbs=modresid.
  do if (stand=1 and ydich=1).
    compute stand=0.
  end if.
  compute anymod2=csum(rsum(wcmat+zcmat+wzcmat)).
  /* do if (anymod2 > 0 and stand = 1) */.
  /*  compute notecode(notes,1) = 27 */.
  /*  compute notes = notes + 1 */.
  /*  compute stand=0 */.
  /* end if */.
  compute funny=1.
  do if (outscree=1).
    print modelvar/title = "**************************************************************************"/format = A8/rnames=modelvlb.
    do if (ncs > 0).
      compute okcovpr=1.
      compute covname9=covnames.
      do if (ncs=xfakecov)).
        compute okcovpr=0.
      end if.
      do if (ncs > xfakecov).
        compute covname9=covname9(1,(xfakecov+1):ncol(covname9)).
      end if.
      do if (okcovpr=1).
        print covname9/title="Covariates:"/format=A8.
      end if.
    end if.
    print n/title="Sample"/rlabel="Size:".
    do if ((!quote(!seed) <> "random")).
      compute seedt=!quote(!seed).
      print seedt/title="Custom"/format=A12/rlabel = "Seed:".
    end if.
    do if (nspl > 0).
      print spl/title = "Location of spline joints:"/rnames=xnmspl/cnames=nmspl/format=!decimals.
    end if.
    do if (clusok=1).
      print clusinfo/title="Number of clusters and cluster size:"/rlabels="Number","Mean","Min","Max"/format=!decimals.
    end if.
  end if.
  compute maxresm=9.
  compute resultm=make(1,maxresm,99999).
  do if (describe=1).
    compute means=csum(dat)/n.  
    compute sdvecdes = sqrt(((n*cssq(dat))-(csum(dat)&**2))/((n-1)*n)).
    compute numdesc = (n*sscp(dat)-(transpos(csum(dat))*(csum(dat)))).
    compute dendesc = sqrt(transpos((n*cssq(dat))-(csum(dat)&**2))*((n*cssq(dat))-(csum(dat)&**2))).
    compute corall = numdesc&/dendesc.
    compute means={means;sdvecdes;cmin(dat);cmax(dat)}.
    !outform outtodo=means/outbig=maxresm.
    !outform outtodo=corall/outbig=maxresm.
    do if (outscree=1).
      print means/title="Variable descriptive statistics"/cnames=varnames/rlabels="Mean","SD","Min","Max"/format=!decimals.
      print corall/title="Variable intercorrelations (Pearson r)"/cnames=varnames/rnames=varnames/format=!decimals.
    end if.
  end if.

  
  do if (outscree=1).
    do if (mcxok=1).
      compute labtmp={xnames,t(xcatlab(1:nxvls,1))}.
      print dummatx/title = "Coding of categorical X variable for analysis:"/cnames = labtmp/format = F6.3.
    end if.
    do if (mcwok=1 and xmint=0).
      compute labtmp={wnames,t(wcatlab(1:nwvls,1))}.
      print dummatw/title = "Coding of categorical W variable for analysis:"/cnames = labtmp/format = F6.3.
    end if.
    do if (mczok=1).
      compute labtmp={znames,t(zcatlab(1:nzvls,1))}.
      print dummatz/title = "Coding of categorical Z variable for analysis:"/cnames = labtmp/format = F6.3.
    end if.
  end if.
end if.
/* cycle through the models */.
do if (criterr = 0).
  do if (savediag=1).
    compute bigny=rmax(cmax(datindx))+1.
    compute regdiag=make((nrow(fulldat)*(nms+nys)),(2+bigny+11+rmax(nump)),999999).
    compute regdiagl=make(1,ncol(regdiag),"xxxxx").
  end if.
  compute outnames=ynames.
  compute outvars=ytmp.
  do if (nms > 0).
    compute outnames={mnames,ynames}.
    compute outvars={mtmp,ytmp}.
    compute indcov=make(((nms*2)+(nms*(nxvls-1))),((nms*2)+(nms*(nxvls-1))),0).
    compute mcsopath=make(((nms*2)+(nms*(nxvls-1))),1,0).
  end if.
  compute labstart=1.
  compute intstart=1.
  compute start=1.
  compute coeffmat=make(1,6,0).
  compute conseq={"        "}.
  compute dfmat=0.
  compute coeffcol=0.
  compute pathscnt=1.
  compute pathscn2=1.
  /* START G LOOP */.
  loop i = 1 to (nms+nys).
  
  /* snip cut */.
    compute highf=make(1,5,0).
    compute highf2=highf.
    do if ((i = (nms+nys)) and (ydich=1)).
      compute highf=make(1,3,0).
      compute highf2=highf.
    end if.
    compute flabel={" "}.
    compute y=outvars(:,i).
    compute xindx=datindx(1:(nump(1,i)-1),i).   
    compute x = fulldat(:,xindx).
    do if (eivdo=1 and ydich=0).
      compute releiv=fulldatr(1,xindx).
      compute releiv={1,releiv}.
      compute x2eiv={y,x}.
      compute tolstat2=make((ncol(x2eiv)-1),1,9999).
      compute numdesc = (nrow(x2eiv)*sscp(x2eiv)-(transpos(csum(x2eiv))*(csum(x2eiv)))).
      compute dendesc = sqrt(transpos((nrow(x2eiv)*cssq(x2eiv))-(csum(x2eiv)&**2))*((n*cssq(x2eiv))-(csum(x2eiv)&**2))).
      compute corbxx = numdesc&/dendesc.
      compute tolstat2=make((ncol(x2eiv)-1),1,9999).
      loop tolcom=2 to ncol(x2eiv).
        compute riichoose=ident(ncol(x2eiv)).
        compute riichoose(tolcom,tolcom)=0.
        compute riichoose(1,1)=1.
        compute riix=riichoose*corbxx*t(riichoose).
        compute riix(tolcom,tolcom)=1. 
        compute riiy=corbxx(:,tolcom)&*diag(riichoose).
        compute tolstat2((tolcom-1),1)=(t(riiy)*inv(riix)*riiy).
      end loop.
      compute tolstat2=csum(t(releiv(1,2:ncol(releiv))) < tolstat2).
      do if (tolstat2 > 0).
        compute criterr=1.
        compute errcode(errs,1)=85.
        compute errs=errs+1.
      end if.  
      release x2eiv.
    end if.
    compute x={ones,x}.
    compute xsq=t(x)*x.
    compute exsq=eval(xsq).
    release xsq.
    compute zeroeig=csum(exsq <= 0.000000000002).
 
   /* snip */.
     
    compute means=csum(x)/n.
    compute vlabsm=vlabs(labstart:(labstart+(nump(1,i)-1)),1).
    /* type: 1 = ols */.
    /* full:  0 = just coefficients, 1 = everything */.
    /* START E */.
    do if (criterr=0).
      do if (ydich=0 or (i < (nms+nys))).
        !modeles y=y/x=x/type=1/full=1.
        /* compute dfmatt=make(nrow(modres),1,modsum(1,6)) */.
        do if (modelres=1).
          compute modresid={modresid,resid}.    
        end if.   
      end if.
      do if (ydich=1 and (i = (nms+nys))).
        !modeles y=y/x=x/type=2/full=1.      
        /* compute dfmatt=make(nrow(modres),1,-999) */.
      end if.
      
      do if (criterr=0).
    
    /* paste snip.  
    do if (outscree=1 and criterr=0).
      print/title = "**************************************************************************".
      do if (eivdo=1).
        print/title = "                     Errors-in-variables regression"/space=0.
      end if.
      do if (savediag=1 and nms > 0).
        print i/title="Diagnostics file equation #:"/format=F2.0.
        print/title=" ".
      end if.
      print outnames(1,i)/title = "OUTCOME VARIABLE:"/format = A8.
      do if (ydich=1 and (i = (nms+nys))).
        compute nmsd = {outnames(1,i), "Analysis"}.
        print rcd/title = "Coding of binary Y for logistic regression analysis:"/cnames = nmsd/format = F9.2.
      end if.
    end if.
    do if (zeroeig > 0).
      print / title = "SINGULAR OR NEAR SINGULAR DATA MATRIX.".
      compute criterr=1.
      compute errcode(errs,1)=31.
      compute errs=errs+1.
    end if.   
      
      !outform outtodo=modsum/outbig=maxresm.
      do if ((ydich=0) or (i < (nms+nys))).
        do if (crossv=1).
          !outform outtodo=crossr/outbig=maxresm.
        end if.
        do if (ssquares=1).
          !outform outtodo=sumtable/outbig=maxresm.
        end if.
      end if.
      !outform outtodo=modres/outbig=maxresm.
      compute obscoeff={obscoeff,t(b)}.
      do if (outscree=1).
        print modsum/title = "Model Summary"/cnames = modsuml/format=!decimals.
        do if ((ydich=0) or (i < (nms+nys))).
          do if (crossv=1).
            print crossr/clabels = "Browne" "LvOut1" "LvOut2"/title = "Shrunken R estimates"/format=!decimals.    
          end if.
          do if (ssquares=1).
            print sumtable/title=" "/rlabels="Regress","Residual","Total"/clabels="SS","df","MS"/format=!decimals/space=0. 
          end if.
        end if.
        print modres/title="Model"/rnames=vlabsm/cnames=modresl/format=!decimals.        
      end if.
      do if (savediag=1 and ydich=0 and triprob=0).        
        compute regdiag((((i-1)*n1)+1):(i*n1),(2+xindx))=fulldat(:,xindx).
        compute regdiag((((i-1)*n1)+1):(i*n1),(2+bigny))=y. 
        compute regdiag((((i-1)*n1)+1):(i*n1),(2+bigny+1):(2+bigny+11+nrow(modres)))=regdiag2.
        compute regdiag((((i-1)*n1)+1):(i*n1),2)=rownum.
        compute regdiag((((i-1)*n1)+1):(i*n1),1)=(ones+i-1).
        release regdiag2.
        do if (i=(nms+nys)).
          do if (i=1).
            compute regdiag=regdiag(:,2:ncol(regdiag)).
            compute regdiagl={"casenum",t(vlabsm(2:nrow(vlabsm),1)),outnames(1,i),"pred","resid","d_resid","stresid","tresid","h","mahal","cook","dmsres","drsq","dskew",t(dfbetas(1:nrow(modres),1))}.
          end if.
          do if (i > 1).
            compute regdiagl(1,(2+bigny):(2+bigny+11+rmax(nump)))={"dv","pred","resid","d_resid","stresid","tresid","h","mahal","cook","dmsres","drsq","dskew",t(dfbetas(1:rmax(nump),1))}.
            compute regdiagl(1,1:2)={"equation","casenum"}.
            loop ijp=1 to (nms+nys).
              compute xindxk=datindx(1:(nump(1,ijp)-1),ijp).
              compute regdiagl(1,(2+xindxk))=fulldatn(1,(xindxk+1)).
            end loop.
          end if.
        end if.
      end if.
      compute basemod=modsum(1,1).
      compute basemodx=basemod.
      do if (ydich=1 and (i = (nms+nys))).
         do if (outscree=1).
           print/title="These results are expressed in a log-odds metric.".
         end if.
         compute notecode(notes,1) = 26.
         compute notes = notes + 1.
      end if.
      compute coeffmat={coeffmat;modres}.
      compute conseqt=make(nrow(modres),1,outnames(1,i)).
      compute conseq={conseq;conseqt}.     
      compute labstart=labstart+nump(1,i).      
      do if (stand=1).
        compute predsd=make(nrow(modres),1,0).
        compute stdmod=modres(:,1)&/ovsd(1,i).
        loop jd=1 to ncol(x).
          compute descdat=x(:,jd).
          compute predsd(jd,1) = (nrow(descdat)*sscp(descdat))-(t(csum(descdat))*(csum(descdat))).
          compute predsd(jd,1) = sqrt(predsd(jd,1)/(nrow(descdat)*(nrow(descdat)-1))).
        end loop.  
        do if (wherex(1,i) <> -999 and ((nxvls > 1) or (xdich=1))).    
          compute sdmsone=make(nxvls,1,1).
          compute predsd(wherex(1,i):wherex(2,i),1)=sdmsone.    
          compute pstog=1.    
        end if.
        compute predsd(1,1)=1.
        compute stdmod=stdmod&*predsd.
        compute stdmod=stdmod(2:nrow(stdmod),1).
        compute sdvlabs=vlabsm(2:nrow(vlabsm),1).
        /* !outform outtodo=stdmod/outbig=maxresm */.
        /* do if (outscree=1) */.
        /*  print stdmod/title="Standardized coefficients"/clabels="coeff"/rnames=sdvlabs/format=!decimals */.
        /* end if */.
      end if.
      do if (stand=1 and   ((i < (nms+nys)) or ((i=(nms+nys)) and ydich=0))).      
        compute ivnames=vlabsm(2:nrow(vlabsm),1).
        !outform outtodo=zppout/outbig=maxresm.
        !outform outtodo=zppout2/outbig=maxresm.
        do if (outscree=1).
          print zppout/title = "Scale-free and standardized measures of association"/rnames=ivnames/clabels="r","sr","pr","standYX", "standY", "standX"/format=!decimals.
          print zppout2/title = " "/rnames=ivnames/clabels="eta-sq","p_eta-sq","f-sq","ipfe"/format=!decimals/space=0.
        end if.
      end if.
      do if (nms > 0 and serial = 0 and (rsum(numint) = 0) and (normal=1 or mc > 0)).
        do if (i < (nms+nys)).
          compute indcov((((i-1)*nxvls)+1):(i*nxvls),(((i-1)*nxvls)+1):(i*nxvls))=varb(2:(1+nxvls),2:(1+nxvls)).
          compute mcsopath((((i-1)*nxvls)+1):(i*nxvls) ,1)=modres(2:(1+nxvls),1).
        end if.
        do if (i = (nms+nys)).
          compute atm=ncol(wherem).
          compute indcov(((nms*nxvls)+1):nrow(mcsopath),((nms*nxvls)+1):nrow(mcsopath))=varb(wherem(1,atm):(wherem(1,atm)+nms-1),wherem(1,atm):(wherem(1,atm)+nms-1)).
          compute mcsopath(((nms*nxvls)+1):nrow(mcsopath),1)=modres(wherem(1,atm):(wherem(1,atm)+nms-1),1).
          compute sobelok=1.          
        end if.
      end if.    
      compute obsdirfx=make(1,nxvls,0).
      compute dirzes=make(1,nxvls,0).
      do if ((i = (nms+nys)) and (bcmat(nrow(bcmat),1)=1)).
        compute direff=modres(2:(1+nxvls),:).
        compute obsdirfx=t(direff(:,1)).
        compute direfflb=modresl.
        compute direffl2=vlabsm(2:(1+nxvls),:).
        compute lmat=make(nrow(b),1,0).
        compute lmat2=make(nxvls,1,1).
        compute lmat(2:(1+nxvls),1)=lmat2.
        do if (ydich <> 1).
          !ftest3 lm=lmat/bcoef=b/cv=varb/chr=1/brsq=r2.
          compute diromni=fresult(1,(1+eivdo):ncol(fresult)).
        end if.
        do if (ydich = 1 ).
          !llrtest lm=lmat.
          compute diromni=fresult.
        end if.
      end if.
      do if (numint(1,i) > 0).
        compute intkeym=intkey(intstart:(intstart+numint(1,i)-1),:).
        do if (outscree=1).
          print intkeym/title="Product terms key:"/format = A8.
        end if.
      end if.
      do if (covcoeff=1).
        do if (outscree=1).
          print varb/title="Covariance matrix of regression parameter estimates:"/rnames=vlabsm
        /cnames=vlabsm/format=!decimals.
        end if.
        !outform outtodo=varb/outbig=maxresm.
      end if.      
      do if ((model=0) and (settest=1)).
        compute lmat2=make(nrow(b),(nxset+setxcat),0).
        compute lmat2(2:(nxset+setxcat+1),:)=ident(nxset+setxcat).
        do if (ydich <> 1).
          !ftest3 lm=lmat2/bcoef=b/cv=varb/chr=(1-(eivdo))/brsq=r2/skip=1.          
          do if (outscree=1).
            compute fresltlb={"R2-chng",hcflab,"df1","df2","p"}.
            do if (eivdo=1).
              compute fresltlb=fresltlb(1,2:ncol(fresltlb)).
            end if.
            print fresult/title = "Hypothesis test for variables in X set:"/cnames=fresltlb/format !decimals.
          end if.
          do if (diagnose=1). 
            compute zzzx=x*lmat2.
            compute zzzx={ones,zzzx}.            
            compute bptestx=0.5*(t(qprime)*zzzx*inv(t(zzzx)*zzzx)*t(zzzx)*qprime).
            compute bprobusx=bptestx*(2/((t(qprime)*qprime)/nrow(zzzx))).
            compute bpresulx={bptestx,ncol(lmat2), (1-chicdf(bptestx,ncol(lmat2)));bprobusx,ncol(lmat2),(1-chicdf(bprobusx,ncol(lmat2)))}.
            release zzzx.
            do if (outscree=1).
              print bpresulx/title="Breusch-Pagan test of heteroskedasticity for variables in X set:"/rlabels="Normal", "Robust"/clabels="Chi-sq", "df", "p"/format=!decimals.
            end if.
          end if.
        end if.
        do if (ydich=1). 
            compute lmat2=rsum(lmat2).   
            !llrtest lm=lmat2.
            do if (outscree=1).
              print fresult/title="Likelihood ratio test for variables in X set:"/clabels="Chi-sq","df","p"/format !decimals.
            end if.
        end if.
        !outform outtodo=fresult/outbig=maxresm.   
        do if (ydich=0 and diagnose=1).
          !outform outtodo=bpresulx/outbig=maxresm.  
        end if. 
      end if.
      /* start X by M interaction */.
      do if (model <> 74 and xmtest=1 and nms > 0).
        compute r2tmp=r2.
        compute btmp=b.
        compute varbtmp=varb.
        compute dfrestmp=dfres.
        compute tvaltmp=tval.
        compute xmtst=make(nms,4,0).
        compute xmtstlbc={hcflab,"df1","df2","p"}.
        do if ((i=(nms+nys)) and (ydich=1)).
          compute xmtst=make(nms,3,0).
          compute xmtstlbc={"Chi-sq","df","p"}.
        end if.
        compute xmtstlb={" "}.
        compute xmtmat=x.
        compute numxint=0.
        loop xmints=2 to i.
          compute x=xmtmat.
          do if ((bcmat((i+1),xmints)=1) and (wzcmat((i+1),xmints) <> 1)). 
            do if (bcmat((i+1),1)=0).
              compute x={xmtmat,xtmp}.
              do if ((ydich=1) and (i=(nms+nys))).
                !modeles y=y/x=x/type=2/full=1.
                compute basemodx=LL2.
              end if.
            end if.
            loop xmtlp1=1 to nxvls.
              compute x={x,xtmp(:,xmtlp1)&*(mtmp(:,(xmints-1))-(csum(mtmp(:,(xmints-1)))/nrow(mtmp)))}.
            end loop.
            do if ((i < (nms+nys)) or (ydich=0)).
              !modeles y=y/x=x/type=1/full=1.
            end if.
            do if ((i = (nms+nys)) and (ydich=1)).
              !modeles y=y/x=x/type=2/full=1.
              compute chidfxm=basemodx-LL2.
            end if.
            compute lmat=make(nrow(b),nxvls,0).
            compute lmattmp=ident(nxvls).
            compute lmat((nrow(lmat)-nxvls+1):nrow(lmat),:)=lmattmp.
            !ftest3 lm=lmat/bcoef=b/cv=varb/skip=1.      
            compute numxint=numxint+1.
            compute xmtst(numxint,:)=fresult.
            do if ((i = (nms+nys)) and (ydich=1)).
              compute xmtst(numxint,1)=chidfxm.
              compute xmtst(numxint,3)=1-chicdf(chidfxm,nxvls).
            end if.
            compute xmtstlb={xmtstlb;highlbx((xmints-1),1)}.     
          end if.      
        end loop.
        compute x=xmtmat.
        release xmtmat.     
        do if (numxint>0).      
          compute xmtstlb=xmtstlb((2:(numxint+1)),:).
          compute xmtst=xmtst(1:numxint,:).
          !outform outtodo=xmtest/outbig=maxresm.
          do if (nms=1).
            compute xmtstlb=" ".
          end if.
          do if (outscree=1).
            do if ((i < (nms+nys)) or (ydich=0)).
              print xmtst/title="Test(s) of X by M interaction:"/rnames=xmtstlb/cnames=xmtstlbc/format=!decimals. 
            end if. 
            do if ((i = (nms+nys)) and (ydich=1)).
              print xmtst/title="Likelihood ratio test(s) of X by M interaction:"/rnames=xmtstlb/cnames=xmtstlbc/format=!decimals. 
            end if. 
          end if.
        end if.
        compute r2=r2tmp.
        compute b=btmp.
        compute varb=varbtmp.
        compute dfres=dfrestmp.
        compute tval=tvaltmp.
      end if.
     /* end X by M interaction */.
      
      /* here is where we do F tests */.
      /* START F */.
      do if (criterr = 0).
        compute jj=0.
        loop j = start to ((start+i)-1).
          compute dbint=0.
          compute lmat=whigh(1:nump(1,i),j).
          compute lmat2=wzhigh(1:nump(1,i),j).
          do if ((csum(lmat) > 0) and (csum(lmat2) = 0)).          
            do if ((i < (nms+nys)) or (ydich <> 1)).
            
              !ftest3 lm=lmat/bcoef=b/cv=varb/chr=1/brsq=r2.
              compute lmatdb=lmat.
              compute dbint=dbint+1.
            end if.
            do if ((ydich = 1) and (i = (nms+nys))).                      
              !llrtest lm=lmat.
              compute lmatdb=lmat.
              compute dbint=dbint+1.
            end if.
            compute highf={highf;fresult}.
            compute highf2={highf2;fresult}.
            do if (j = start).
              compute flabel={flabel;"X*W"}.
            end if.
            do if (j > start).
              do if (nms > 1).
                compute flabel={flabel;highlbw(jj,1)}.
              else if (nms = 1).
                do if (xmint=0).
                  compute flabel={flabel;"M*W"}.
                end if.
                do if (xmint=1).
                  compute flabel={flabel;"X*M"}.                
                end if.
              end if.
            end if.
          end if.      
          compute lmat=zhigh(1:nump(1,i),j).
          compute lmat2=wzhigh(1:nump(1,i),j).
          do if ((csum(lmat) > 0) and (csum(lmat2) = 0)).
            do if ((i < (nms+nys)) or (ydich <> 1)).
              !ftest3 lm=lmat/bcoef=b/cv=varb/chr=1/brsq=r2.
              compute dbint=dbint+1.             
            end if.
            do if ((ydich = 1) and (i = (nms+nys))).
              !llrtest lm=lmat.
              compute dbint=dbint+1.    
            end if.
            compute highf={highf;fresult}.   
            compute highf2={highf2;fresult}.                
            do if (j = start).
              compute flabel={flabel;"X*Z"}.
            end if.
            do if (j > start).
              do if (nms > 1).
                compute flabel={flabel;highlbz(jj,1)}.
              else if (nms = 1).
                compute flabel={flabel;"M*Z"}.
              end if.
            end if.
          end if.
          do if (dbint=2).
            compute lmatdb=lmatdb+lmat.
            do if ((ydich = 1) and (i = (nms+nys))).
              !llrtest lm=lmatdb.
            end if.     
            do if (ydich <> 1 or i < (nms+nys)).
              !ftest3 lm=lmatdb/bcoef=b/cv=varb/chr=1/brsq=r2.
            end if.
            compute dbint=0. 
            compute highf={highf;fresult}.
            do if (jj=0 and nms > 0).
              compute flabel={flabel;"BOTH(X)"}.  
            end if.
            do if (jj=0 and nms = 0).
              compute flabel={flabel;"BOTH"}.  
            end if.
            do if (jj>0 and nms = 1).
              compute flabel={flabel;"BOTH(M)"}.  
            end if.
            do if (nms > 1 and jj > 0). 
              compute flabel={flabel;highlbbt(jj,1)}.        
            end if.         
          end if. 
          compute lmat2=wzhigh(1:nump(1,i),j).
          do if (csum(lmat2) > 0).
            do if ((i < (nms+nys)) or (ydich <> 1)).
              !ftest3 lm=lmat2/bcoef=b/cv=varb/chr=1/brsq=r2.
            end if.
            do if ((ydich = 1) and (i = (nms+nys))).
              !llrtest lm=lmat2.
            end if.
            compute highf={highf;fresult}.
            compute highf2={highf2;fresult}.
            do if (j = start).
              compute flabel={flabel;"X*W*Z"}.
            end if.
            do if (j > start).
              do if (nms > 1).
                compute flabel={flabel;highlbwz(jj,1)}.
              else if (nms = 1).
                compute flabel={flabel;"M*W*Z"}.
              end if.
            end if.
          end if.
          compute jj=jj+1.
        end loop.
        release jj.
        compute start=start+i.
      end if.
      /* END F*/.
      do if (nrow(highf) > 1).
        compute highf=highf(2:nrow(highf),:).
        compute highf2=highf2(2:nrow(highf2),:).
        compute flabel=flabel(2:nrow(flabel),1).
        !outform outtodo=highf/outbig=maxresm.
        do if (outscree=1).
          do if ((i < nms+nys) or (ydich=0)).
             compute clabtmp={"R2-chng", hcflab,"df1","df2","p"}.
             print highf/format = !decimals/rnames=flabel/cnames=clabtmp/
            title = "Test(s) of highest order unconditional interaction(s):".
          end if.
          do if (ydich=1 and i=(nms+nys)).
            compute clabtmp={"Chi-sq", "df","p"}.
            print/title="Likelihood ratio test(s) of highest order".   
            print highf/format = !decimals/rnames=flabel/cnames=clabtmp/
            title = "unconditional interactions(s):"/space=0. 
          end if.      
        end if.  
        compute intpb=highf2(:,ncol(highf2)).
      end if.
      compute intstart=intstart+numint(1,i).
    end if.
    /* END E*/.    
    !PROBE3 decpnt=!decimals.
    do if (model = 74 and i <= nms).
      compute onetemp=make(nrow(xprobval),1,1).
      compute mestmt74={onetemp,xprobval}.
      do if (ncs > 0).
        compute ncovmdl=rsum(ccmat(i,:)).
        do if (ncovmdl > 0).
          compute cvmnc=make(nrow(mestmt74),ncovmdl,1).
          compute cvmnctmp=csum(x(:,(ncol(x)-ncovmdl+1):ncol(x)))/nrow(x). 
          do if (cuscoval > 0). 
            compute cvmnctmp=coval.
          end if.  
          loop mestlp=1 to ncovmdl.
            compute cvmnc(:,mestlp)=cvmnc(:,mestlp)*cvmnctmp(:,mestlp).
          end loop.  
          compute mestmt74={mestmt74,cvmnc}.
        end if.
      end if.
      compute mest74t=mestmt74*b.
      do if (i = 1).
        compute mest74=mest74t.
      end if.
      do if (i > 1).
        compute mest74={mest74,mest74t}.
      end if.
    end if.
    
    
    do if ((i = (nms+nys)) and (model >= 0 and model < 4) and (linsum(1,1) <> -999)).
      compute lhyprob=1.
      compute meansub=0.
      compute covsubok=0.
      /* snipped from here */.
      do if (nlinsum=nrow(b)).
        compute lhyprob=0.
      end if.
      do if ((nlinsum = (nrow(b)-ncs)) and ncs > 0 and model<>0).
        compute lhyprob=0.
        compute covsubok=1.
        /* pasted here */.
        compute covmean2=covmeans.
      end if.
      do if ((nlinsum=(nrow(b)-ncs+xfakecov)) and model=0 and lhyprob=1).
        compute covmean2=covmeans(1,(xfakecov+1):ncol(covmeans)).
        compute covsubok=1.
        compute lhyprob=0.
      end if.
      do if (lhyprob=0).      
      /* do if ( (nlinsum=nrow(b)) or ((nlinsum = (nrow(b)-ncs)) and model<>0) ) */.
        do if (covsubok=1).
          compute linsum={linsum,covmean2}.
          compute meansub=1.
        end if.
        compute hypest=linsum*b.
        compute sehypest=sqrt(linsum*varb*t(linsum)).
        do if (ydich=0).
          compute phypest=2*(1-tcdf(abs(hypest/sehypest), (dfres))).
          compute hypest={hypest,sehypest,(hypest/sehypest),phypest,(hypest-tval*sehypest),(hypest+tval*sehypest)}.
        end if.
        do if (ydich=1).
          compute phypest=2*(1-cdfnorm(abs(hypest/sehypest))).
          compute hypest={hypest,sehypest,(hypest/sehypest),phypest,(hypest-xp2*sehypest),(hypest+xp2*sehypest)}.
        end if.
        !outform outtodo=t(linsum)/outbig=maxresm.
        !outform outtodo=hypest/outbig=maxresm.
        do if (outscree=1).
          do if (ydich=0).
            compute hyplabs={"Estimate",hclab,"t","p","LLCI","ULCI"}.
          end if.
          do if (ydich=1).
            compute hyplabs={"Estimate","se","Z","p","LLCI","ULCI"}.
          end if.
          print/title="----------".
          print/title="Linear Combination Estimate and Hypothesis Test"/space=0.
          print t(linsum)/title= "Weight vector:"/clabels="weight"/rnames=vlabsm/format !decimals.
          print hypest/title=" "/cnames=hyplabs/format=!decimals/space=0.
          do if (meansub=1).
            print/title="Covariate weight(s) set to the sample mean.".
          end if.
        end if.
      end if.
      do if (lhyprob=1).   
        compute notecode(notes,1) = 30.
        compute notes = notes + 1. 
      end if.
    end if.
    
    do if ((subsets=1 or dominate=1) and model=0 and ydich=0).     
      compute numbx=ncol(x)-nxvls.
      compute subdomin=1.
      do if (numbx < 2 or numbx > 15).
        compute notecode(notes,1) = 40.
        compute notes = notes + 1. 
        compute subdomin=0.
        compute dominate=0.
      end if.
      do if (subdomin=1).
        compute rsqrmat = make(((2&**numbx)-1),numbx+2,0).
        compute dommat=make(numbx,numbx,0).
        compute desc5x=corall(2:(numbx+nxvls),2:(numbx+nxvls)).
        compute allcomp=(2&**(numbx-2)).
        compute domx = make(1,numbx,0).
        compute domx2=make(1,(numbx+(nxvls-1)),0).
        loop jss = 1 to ((2&**numbx)-1).
          compute domx(1,1)=domx(1,1)+1.
          loop iss = 1 to numbx.
            do if (domx(1,iss) = 2).
              compute domx(1,iss)=0.
              compute domx(1,(iss+1))=domx(1,(iss+1))+1.
            end if.
          end loop.
          compute rsqrmat(jss,1:(ncol(rsqrmat)-2))=domx.
          compute xdumsone=make(1,nxvls,0).
          do if (domx(1,1)=1).
            compute xdumsone=make(1,nxvls,1).
          end if.
          compute domx2={xdumsone,domx(1,2:ncol(domx))}.
          compute rii=(t(domx2)*domx2)&*desc5x.
          compute bin=make(1,(numbx+nxvls-1),0).
          loop ksss = 1 to (numbx+nxvls-1).
            compute rii(ksss,ksss)=1.
            compute bin(1,ksss)=2&**(ksss-1).
          end loop .
          compute riy=corall(1,2:(numbx+nxvls))&*domx2.
          compute rsqtmp=(riy*inv(rii)*t(riy)).
          compute rsqrmat(jss,numbx+1)=rsqtmp.
          compute rsqrmat(jss,numbx+2)=1-((1-rsqtmp)*(n1-1)/(n1-(rsum(domx)+(domx(1,1)*(nxvls-1))+1))).
          compute rsqrmat(jss,1:numbx)=domx.
        end loop.
        compute temp = rsqrmat. 
        compute temp(GRADE(rsqrmat(:,(ncol(rsqrmat)-0))),:) = rsqrmat.
        compute rsqrmats = temp. 
        compute ssvlabs={varnames(1,2:ncol(varnames)),"R-sq", "Adj R-sq"}.
        do if (subsets=1).
          !outform outtodo=rsqrmats/outbig=maxresm.
          do if (outscree=1).
            print rsqrmats/title="All subsets regression results"/cnames=ssvlabs/format=!decimals.
          end if.
        end if.
      end if.
    end if.
    do if (dominate = 1).
      loop idmt = 1 to (numbx-1).
        loop jdmt = (idmt+1) to numbx.
          compute critdiff=bin(1,jdmt)-bin(1,idmt).
          loop kdmt = 1 to nrow(rsqrmat).
            do if (rsqrmat(kdmt,idmt)=1 and rsqrmat(kdmt,jdmt)=0).
               compute tmps={rsqrmat(kdmt,(numbx+1)),rsqrmat((kdmt+critdiff),(numbx+1))}.
               do if (tmps(1,1) > tmps(1,2)).
                 compute dommat(idmt,jdmt)=dommat(idmt,jdmt)+1.
               end  if.
               do if (tmps(1,1) < tmps(1,2)).
                 compute dommat(jdmt,idmt)=dommat(jdmt,idmt)+1.
               end if.
            end if.
          end loop.
        end loop.
      end loop.
      compute dommat=dommat/allcomp.
      compute dovlabs=varnames(1,2:ncol(varnames)).
      !outform outtodo=dommat/outbig=maxresm.
      do if (outscree=1).
        /* print/title = "***************************************************************************" */.
        print dommat/title "Dominance matrix"/cnames=dovlabs/rnames=dovlabs/format = F4.3.
      end if.
    end if.
    do if (diagnose=1).
      do if (ydich=0).
         do if (tricheck=0).
           !outform outtodo=maxnmin/outbig=maxresm.
         end if.
         !outform outtodo=skewres/outbig=maxresm.
         do if (tricheck=0).
           !outform outtodo=bp/outbig=maxresm.
         end if.
         !outform outtodo=wadf/outbig=maxresm.
         compute tolstat={tolstat,(1/tolstat)}.    
         !outform outtodo=tolstat/outbig=maxresm.
         !outform outtodo=bpresult/outbig=maxresm.         
        do if (outscree=1).
          print/title= "-----------".
          print/title="Some regression diagnostics"/space=0.
          do if (tricheck=0).
            print maxnmin/title=" "/rlables="fitted","residual","t-resid"/clabels="Min.", "Max."/format=!decimals/space=0.
          end if.
          print skewres/title="Shape of residuals"/rlabels="Value","se"/clabels="Skewness","Kurtosis"/format=!decimals.
          do if (tricheck=0).
            print bp/title = "Bonferroni-corrected p for largest t-residual"/clabels "t-resid" "p-value" "casenum"/format=!decimals.  
          end if.
          print wadf/title = "Most influential observations"/rnames=vlabsm/clabels "casenum" "dfbeta"/format=!decimals.
          compute tolvnm=vlabsm(2:nrow(vlabsm),1).      
          print tolstat/title="Variable tolerance and VIF"/rnames=tolvnm/clabels="Tol.","VIF"/format=!decimals.
          print bpresult/title="Breusch-Pagan test of heteroskedasticity"/rlabels="Normal", "Robust"/clabels="Chi-sq", "df", "p"/format=!decimals.
        end if.
      end if.
      
    end if.
    
    end if.

  end loop if criterr=1.
  /* end G loop */.
  /* I added this do total condition just below because of R */.
  do if ((criterr=0)).
    compute lastb=b.
    compute lastcov=varb.
  end if.
  do if (criterr=0 and dototal = 1).
    compute x=xtmp.
    compute releiv={1,relx}.
    compute vlabsm={"constant";xcatlab(1:nxvls,1)}.
    do if (ncs > 0).
      compute x = {x,ctmp}.
      compute vlabsm={vlabsm;t(covnames)}.
      compute releiv={releiv,relcov}.
    end if.
    compute x = {ones,x}.
    !modeles y=y/x=x/type=1/full=1.
    
    do if (criterr=0).
    
      compute toteff=modres(2:(1+nxvls),:).
      compute nodotot=0.
      do if (xdich=1 and xmint=1 and model=74).
        compute toteff(:,1)=toteff(:,1)*xscaling.
        compute toteff(:,2)=toteff(:,2)*abs(xscaling).
        compute toteff(:,3)=toteff(:,3)*xscaling.
        compute toteff(:,5)=toteff(:,5)*xscaling.
        compute toteff(:,6)=toteff(:,6)*xscaling.
        compute citmp=toteff(:,5:6).
        do if (xscaling < 0).
          compute toteff(:,5)=citmp(:,2).
          compute toteff(:,6)=citmp(:,1).
        end if.
        compute nodotot=1.
      end if.
      compute totefflb=modresl.
      do if (outscree=1). 
        do if (nodotot=0).
          print/title = "************************** TOTAL EFFECT MODEL ****************************".
          do if (eivdo=1).
            print/title = "                     Errors-in-variables regression"/space=0.
          end if.
        end if.
        do if (nodotot=1).
          print/title = "**************************************************************************".
        end if.
        print outnames(1,ncol(outnames))/title = "OUTCOME VARIABLE:"/format = A8.
      end if.
      compute toteffl2=vlabsm(2:(1+nxvls),:).
      !outform outtodo=modsum/outbig=maxresm.   
      do if (outscree=1).
        print modsum/title = "Model Summary"/cnames = modsuml/format=!decimals.
      end if.
      do if ((ydich=0) or (i < (nms+nys))).
        do if (crossv=1).
          !outform outtodo=crossr/outbig=maxresm.
          do if (outscree=1).        
            print crossr/clabels = "Browne" "LvOut1" "LvOut2"/title = "Shrunken R estimates"/format=!decimals.    
          end if.
        end if.
        do if (ssquares=1).
          !outform outtodo=sumtable/outbig=maxresm.
          do if (outscree=1).
            print sumtable/title=" "/rlabels="Regress","Residual","Total"/clabels="SS","df","MS"/format=!decimals/space=0. 
          end if.
        end if.
      end if.
      !outform outtodo=modres/outbig=maxresm.
      do if (outscree=1).
        print modres/title="Model"/rnames=vlabsm/cnames=modresl/format=!decimals.              
      end if.
      compute lmat=make(nrow(b),1,0).
      compute lmat2=make(nxvls,1,1).
      compute lmat(2:(1+nxvls),1)=lmat2.
      do if (ydich <> 1).
        !ftest3 lm=lmat/bcoef=b/cv=varb/chr=1/brsq=r2.
         compute totomni=fresult(1,(1+eivdo):ncol(fresult)).
      end if.    
      do if (stand=1).
          compute predsd=make(nrow(modres),1,0).
          compute stdmod=modres(:,1)&/ovsd(1,ncol(ovsd)).
          loop jd=1 to ncol(x).
            compute descdat=x(:,jd).
            compute predsd(jd,1) = (nrow(descdat)*sscp(descdat))-(t(csum(descdat))*(csum(descdat))).
            compute predsd(jd,1) = sqrt(predsd(jd,1)/(nrow(descdat)*(nrow(descdat)-1))).
          end loop.  
          do if (wherex(1,ncol(wherex)) <> -999 and ((nxvls > 1) or (xdich=1))).    
            compute sdmsone=make(nxvls,1,1).
            compute predsd(wherex(1,ncol(wherex)):wherex(2,ncol(wherex)),1)=sdmsone.    
            compute pstog=1.    
          end if.
          compute predsd(1,1)=1.
          compute stdmod=stdmod&*predsd.
          compute stdmod=stdmod(2:nrow(stdmod),1).
          compute sdvlabs=vlabsm(2:nrow(vlabsm),1).
         /* !outform outtodo=stdmod/outbig=maxresm */.
          
          /* do if (outscree=1) */.
         /*   print stdmod/title="Standardized coefficients"/clabels="coeff"/rnames=sdvlabs/format=!decimals */.
          /* end if */.
        end if.
      do if (stand=1 and ydich=0).      
        compute ivnames=vlabsm(2:nrow(vlabsm),1).
        !outform outtodo=zppout/outbig=maxresm.
        !outform outtodo=zppout2/outbig=maxresm.
        do if (outscree=1).
          print zppout/title = "Scale-free and standardized measures of association"/rnames=ivnames/clabels="r","sr","pr","standYX","standY","standX"/format=!decimals.
          print zppout2/title = " "/rnames=ivnames/clabels="eta-sq","p_eta-sq","f-sq"/format=!decimals/space=0.
        end if.
      end if.    
      do if (covcoeff=1).
        do if (outscree=1).
          print varb/title="Covariance matrix of regression parameter estimates:"/rnames=vlabsm/cnames=vlabsm/format=!decimals.
        end if.
        !outform outtodo=varb/outbig=maxresm.
      end if.
      do if (diagnose=1).
        do if ((ydich=0) or (i < (nms+nys))).
          do if (tricheck=0).
             !outform outtodo=maxnmin/outbig=maxresm.
           end if.
             !outform outtodo=skewres/outbig=maxresm.
           do if (tricheck=0).
             !outform outtodo=bp/outbig=maxresm.
           end if.
           !outform outtodo=wadf/outbig=maxresm.
           compute tolstat={tolstat,(1/tolstat)}.    
           !outform outtodo=tolstat/outbig=maxresm.
           !outform outtodo=bpresult/outbig=maxresm.
          do if (outscree=1).
            print/title= "-----------".
            print/title="Some regression diagnostics"/space=0.
            do if (tricheck=0).
              print maxnmin/title=" "/rlables="fitted","residual","t-resid"/clabels="Min.", "Max."/format=!decimals/space=0.
            end if.
            print skewres/title="Shape of residuals"/rlabels="Value","se"/clabels="Skewness","Kurtosis"/format=!decimals.
            do if (tricheck = 0).
              print bp/title = "Bonferroni-corrected p for largest t-residual"/clabels "t-resid" "p-value" "Casenum"/format=!decimals.  
            end if.
            print wadf/title = "Most influential observations"/rnames=vlabsm/clabels "casenum" "dfbeta"/format=!decimals.
            compute tolvnm=vlabsm(2:nrow(vlabsm),1).      
            print tolstat/title="Variable tolerance and VIF"/rnames=tolvnm/clabels="Tol.","VIF"/format=!decimals.
            print bpresult/title="Breusch-Pagan test of heteroskedasticity"/rlabels="Normal", "Robust"/clabels="Chi-sq", "df", "p"/format=!decimals.
          end if.
        end if.
      end if.    
    end if.                    
  end if.
end if.
/* end cycle through the models */.

do if (criterr=0 and nms > 0 and ydich=0 and modelres=1).    
  compute modresid=modresid(:,2:ncol(modresid)).
  compute sdvecres = sqrt(((n*cssq(modresid))-(csum(modresid)&**2))/((n-1)*n)).
  compute numdesc = (n*sscp(modresid)-(transpos(csum(modresid))*(csum(modresid)))).
  compute dendesc = sqrt(transpos((n*cssq(modresid))-(csum(modresid)&**2))*((n*cssq(modresid))-(csum(modresid)&**2))).
  compute corall = numdesc&/dendesc.
  !outform outtodo=corall/outbig=maxresm.
  do if (outscree=1).
    print/title = "****************** CORRELATIONS BETWEEN MODEL RESIDUALS ******************".
    print corall/title=" "/format=!decimals/
      cnames=outnames/rnames=outnames/space=0.
  end if.
end if.

/* DO BOOTSTRAPPING */.
do if (criterr=0 and boot > 0).
  compute bootres=make(1,rsum(nump),-999).
  compute bootdir=obsdirfx.
  compute natdirbt=make(1,nxvls,-999).
  do if (effsize=1).
    compute bootysd=make(1,1,-999).
    compute bootxsd=make(1,1,-999).
  end if.
  compute badboot=0.
  compute goodboot=0.
  compute eiverrc=0.
  compute smallest=1.
  compute booting=1.
  compute singerc=0.
  compute v=make(n,1,0).
  compute ones2=ones.
  do if (clusok=1 and clusboot > 0).
    compute temp = clusdat.
    compute temp(GRADE(clusdat(:,3)),:) = clusdat.
    compute clusdat = temp.
  end if.
  loop j = 1 to maxboot.
    compute v4=999999.
    compute nobootx=1.
    compute modres2=999.
    compute eiverrct=0.
    compute singcnt=0.
    compute v=trunc(uniform(n,1)*n)+1.
    do if (clusok=1).      
      loop clbt=1 to nclus.
        do if (clusboot=1).
          compute isht=clusdat(current(clbt,3):current(clbt,4),1).
          compute v2=trunc(uniform(current(clbt,2),1)*current(clbt,2))+1.
          compute v(current(clbt,3):current(clbt,4),1)=isht(v2,1).
        end if.
        do if (clusboot > 1).
          compute v2=trunc(uniform(1,1)*nclus)+1.
          compute isht=clusdat(current(v2,3):current(v2,4),1).
          compute v4={v4;isht}.
        end if.
      end loop.
      do if (clusboot > 1).
        compute v=v4(2:nrow(v4),1).
        compute ones2=make(nrow(v),1,1).
      end if.
    end if.
    compute bad=0.
    loop i = 1 to (nms+nys).
      compute y=outvars(v,i).
      compute ynovar= (nrow(y)*sscp(y))-(t(csum(y))*(csum(y))).
      do if (ynovar <= 0.00000000001).
        compute bad=1.
      end if.
      compute xindx=datindx(1:(nump(1,i)-1),i).
      compute hello=0.
      compute x = fulldat(v,xindx).
      compute x={ones2,x}.
      compute xsq=t(x)*x.
      compute exsq=eval(xsq).
      release xsq.
      compute holymoly=cmin(exsq).
                                                         /*     .000000000002344350286  */.
      compute zeroeig=csum(exsq <= 0.000000000002).
      compute bad=bad + (zeroeig > 0).
      compute singcnt=(zeroeig > 0).
      /* !describ descdatf=y/type=1 */.
      /* compute bad=bad+ (desctmp(2,1) <= 0.00000000001) */.
      do if (bad = 0).
        do if (holymoly < smallest).
          compute smallest=holymoly.
        end if.
        /* print smallest/format=F16.15 */.
        do if (ydich=0 or (i < (nms+nys))). 
          do if (eivdo=1).
            compute releiv=fulldatr(1,xindx).
            compute releiv={1,releiv}.         
          end if.
          !modeles y=y/x=x/type=1/full=0.
          do if (eiverr=1).
            compute bad=1.
            compute eiverrct=1.
          end if.
        end if.
        do if (ydich=1 and (i = (nms+nys))).
          !modeles y=y/x=x/type=2/full=0.
        end if.
        do if (bad=0).
          compute modres2={modres2,t(modres)}.
        end if.
        do if ((i = (nms+nys)) and (bad=0)).
          do if (bcmat((i+1),1) = 1).
            compute bootdir={bootdir;t(modres(wherex(1,i):wherex(2,i),1))}.       
          end if.
          do if (bcmat((i+1),1) = 0).
            compute bootdir={bootdir;dirzes}.   
          end if.
        end if.             
        do if (model = 74 and i <= nms).
          compute onetemp=make(nrow(xprobval),1,1).
          compute mestmtb={onetemp,xprobval}.
          do if (ncs > 0).
            compute ncovmdlb=rsum(ccmat(i,:)).
            do if (ncovmdl > 0).
              compute cvmncb=make(nrow(mestmtb),ncovmdlb,1).
              compute cvmnctmp=csum(x(:,(ncol(x)-ncovmdlb+1):ncol(x)))/nrow(x).   
              do if (cuscoval > 0). 
                compute cvmnctmp=coval.
              end if.  
              loop mestlp=1 to ncovmdlb.
                compute cvmncb(:,mestlp)=cvmncb(:,mestlp)*cvmnctmp(:,mestlp).
              end loop.  
              compute mestmtb={mestmtb,cvmncb}.
            end if.
          end if.
          compute mestbt=mestmtb*modres.
          do if (i = 1).
            compute mestb=mestbt.
          end if.
          do if (i > 1).
            compute mestb={mestb,mestbt}.
          end if.
        end if.
        do if (model = 74 and i = (nms+nys)).
          compute xvalptmp=make(1,nxvls,0).
          compute mest74sp=make(1,nms,0).
          compute mest74sp={mest74sp;mestb}.
          compute xvalptmp={xvalptmp;ident(nxvls)}.
          compute dirfxcf=make(nrow(mestb),1,0).         
          loop cfloop1=1 to nrow(mestb).
            compute ndirfx=make(nrow(modres),1,0).
            compute ndirfx(2:(nxvls+1),1)=t(xvalptmp(cfloop1,:)).
            loop cfloop3=1 to nms.
              do if (mcx=1 or mcx=0).
                compute ndirfx(wheremw((1+((cfloop3-1)*2)),ncol(wheremw)):wheremw((2+((cfloop3-1)*2)),ncol(wheremw)))=t(xvalptmp(cfloop1,:)*mestb(1,cfloop3)).
              end if.
              do if (mcx=2).
                compute ndirfx(wheremw((1+((cfloop3-1)*2)),ncol(wheremw)):wheremw((2+((cfloop3-1)*2)),ncol(wheremw)))=t(xvalptmp(cfloop1,:)*mest74sp(cfloop1,cfloop3)).
              end if.
            end loop.
            do if (cfloop1 > 1).
              compute dirfxcf(cfloop1,1)=t(ndirfx)*modres*xscaling.
            end if.
          end loop.
          compute natdirbt={natdirbt;t(dirfxcf(2:nrow(dirfxcf),:))}.
        end if.                    
        do if ( (bcmat((i+1),1) = 1) and (nobootx=1) and (effsize=1)).
          compute nobootx=0.
          compute xsdtemp= (nrow(x)*sscp(x(:,2)))-(t(csum(x(:,2)))*(csum(x(:,2)))).
          compute xsdtemp= sqrt(xsdtemp/(nrow(x)*(nrow(x)-1))).
        end if.
      end if.
    end loop.
    do if (bad = 0).
      compute bootres={bootres;modres2(:,2:ncol(modres2))}.      
      do if (effsize=1).
        compute ysdtemp= (nrow(y)*sscp(y))-(t(csum(y))*(csum(y))).
        compute ysdtemp= sqrt(ysdtemp/(nrow(y)*(nrow(y)-1))).
        compute bootysd={bootysd;ysdtemp}.
        compute bootxsd={bootxsd;xsdtemp}.
      end if.
      compute goodboot=goodboot+1.
    end if.
    do if (bad <> 0).
      compute badboot=badboot+1.
      compute eiverc=eiverc+eiverrct.
      compute singerc=singerc+singcnt.
    end if.
  end loop if (goodboot = boot).
  /* print smallest/format=F16.15 */.
  compute bootres=bootres(2:nrow(bootres),:).
  do if (effsize=1).
    compute bootysd=bootysd(2:nrow(bootysd),:).
    do if (nrow(bootxsd) > 1).
      compute bootxsd=bootxsd(2:nrow(bootxsd),:).
    end if.
  end if.
  do if (goodboot < (boot)).
    compute boot=0.
    compute modelbt=0.
    compute notecode(notes,1) = 7.
    compute notes = notes + 1.
  end if.
  do if (boot > 0).
     do if (effsize=1).
      compute bootysd={ysd;bootysd}.
      compute bootxsd={xsd;bootxsd}.
    end if.
    do if (saveboot = 1).
      save bootres/outfile = *.
    end if.
    do if (modelbt=1).
      compute bootcim=make(ncol(bootres),5,-99999).
      compute bootcim(:,2) = t(csum(bootres)/nrow(bootres)).
      compute bootcim(:,1) = coeffmat(2:nrow(coeffmat),1).
      loop i = 1 to ncol(bootres).
        !bcboot3 databcbt = bootres(:,i)/estmte=(bootcim(i,1)*bc)+(9999*(1-bc)).
        compute bootcim(i,4:5)={llcit,ulcit}.
        compute bootcim(i,3)=bootse.
      end loop.
    end if.  
  end if.
  do if (badboot > 0).
    compute notecode(notes,1) = 6.
    compute notes = notes + 1.
  end if.
end if.

/* print direct and indirect effects */.
do if (criterr=0).
  release fulldat,x,y.  
end if.

do if (xmint=1 and criterr=0).
  compute xvalptmp=make(1,nxvls,0).
  compute mest74sp=make(1,nms,0).
  compute mest74sp={mest74sp;mest74}.
  compute xvalptmp={xvalptmp;ident(nxvls)}.
  loop kcfuhd=1 to 2.
    compute dirfxcf=make(nrow(mest74),6,0).
    loop cfloop1=1 to nrow(mest74).
      compute ndirfx=make(nrow(lastb),1,0).
      compute ndirfx(2:(nxvls+1),1)=t(xvalptmp(cfloop1,:)).
      loop cfloop3=1 to nms.
        do if (kcfuhd=1).
          do if (mcx=1 or mcx=0).
            compute ndirfx(wheremw((1+((cfloop3-1)*2)),ncol(wheremw)):wheremw((2+((cfloop3-1)*2)),ncol(wheremw)))=t(xvalptmp(cfloop1,:)*mest74(1,cfloop3)).
          end if.
          do if (mcx=2).
            compute ndirfx(wheremw((1+((cfloop3-1)*2)),ncol(wheremw)):wheremw((2+((cfloop3-1)*2)),ncol(wheremw)))=t(xvalptmp(cfloop1,:)*mest74sp(cfloop1,cfloop3)).
          end if.
        end if.
        do if (kcfuhd=2).
          compute ndirfx(wheremw((1+((cfloop3-1)*2)),ncol(wheremw)):wheremw((2+((cfloop3-1)*2)),ncol(wheremw)))=t(xvalptmp(cfloop1,:)*medmeans(1,cfloop3)).
        end if.
      end loop.
      do if (cfloop1 > 1).
        compute dirfxcf(cfloop1,1)=t(ndirfx)*lastb*xscaling.
        compute dirfxcf(cfloop1,2)=sqrt(diag(t(ndirfx)*lastcov*ndirfx))*abs(xscaling).
        compute dirfxcf(cfloop1,3)=dirfxcf(cfloop1,1)/dirfxcf(cfloop1,2).
        compute dirfxcf(cfloop1,4)=2*(1-tcdf(abs(dirfxcf(cfloop1,3)), dfres)).
        compute dirfxcf(cfloop1,5)=dirfxcf(cfloop1,1)-tval*dirfxcf(cfloop1,2).  
        compute dirfxcf(cfloop1,6)=dirfxcf(cfloop1,1)+tval*dirfxcf(cfloop1,2).
      end if.
    end loop.
    compute codireff=dirfxcf(2:nrow(dirfxcf),:).
    do if (kcfuhd=1).
      compute direff=codireff.
    end if.
  end loop.
end if.

!drinef3 decpnt=!decimals.

do if (criterr=0 and (saveest = 1 or saveboot=1)).
  compute coeffsav=coeffmat(2:nrow(coeffmat),:).  
  compute conseq=conseq(2:nrow(conseq),1).
  /* compute dfmat=dfmat(2:nrow(dfmat),1) */.
  /* compute dfmat=dfmat */.
  compute saven={"conseqnt","antecdnt","coeff","se","t","p","LLCI","ULCI","df"}.
  do if (ydich=1).
    compute saven={"conseqnt","antecdnt","coeff","se","t_or_Z","p","LLCI","ULCI","df"}.
  end if.
 /* do if (saveest=1) */.
      /* save {conseq,vlabs,coeffsav,dfmat}/outfile = /names=saven/strings=conseqnt antecdnt */.
 /* end if */.
  do if (saveboot=1 and boot > 0 ).
    compute colslab=make(300,1," ").
    !do !i = 1 !to 300.
      !let !v2=!concat("COL",!unquote(!quote(!i))).
      compute colslab(!i,1)=!quote(!v2).
    !doend.
    compute colslab=colslab(1:ncol(bootres),1).
    compute colslab={colslab,conseq,vlabs}.
    do if (outscree=1).
      print/title = "**************************************************************************".
      print/title="Bootstrap estimates were saved to a file"/space=0.
      print colslab/title="Map of column names to model coefficients:"/clabels=" ","Conseqnt","Antecdnt"/format=a8.
    end if.
  end if.
end if.


/* PRINT BOOTSTRAP RESULTS FOR MODEL PARAMETERS */.
do if (criterr=0 and boot > 0 and modelbt=1).
  compute labstart=1.
  do if (outscree=1).
    print/title = "*********** BOOTSTRAP RESULTS FOR REGRESSION MODEL PARAMETERS ************".
  end if.
  loop iboot = 1 to (nms+nys).
    do if (outscree=1).
      print outnames(1,iboot)/title = "OUTCOME VARIABLE:"/format = A8.
    end if.
    compute vlabsm=vlabs(labstart:(labstart+(nump(1,iboot)-1)),1).
    compute outnmtmp=bootcim(labstart:(labstart+(nump(1,iboot)-1)),:).
    !outform outtodo=outnmtmp/outbig=maxresm.

    do if (outscree=1).
      print bootcim(labstart:(labstart+(nump(1,iboot)-1)),:)/title=" "/rnames=vlabsm/
      clabels="Coeff" "BootMean" "BootSE" "BootLLCI" "BootULCI"/format=!decimals/space=0.
    end if.
    compute labstart=labstart+nump(1,iboot).
    do if (iboot < (nms+nys)).
      do if (outscree=1).
        print/title= "----------".
      end if.
    end if.
  end loop.
end if.

do if (criterr=0 and saveest=1).
  do if (boot > 0).
   compute erroutp={badboot, singerc, eiverc}.
   !outform outtodo=erroutp/outbig=maxresm.
  end if.
  compute resultm=resultm(2:nrow(resultm),:).
  compute bocaj=(resultm=99999).
  compute bocaj=csum(bocaj).
  compute bocaj=1-(bocaj=nrow(resultm)).
  compute j=1.
  loop i = 1 to ncol(resultm).
    do if (bocaj(1,i)=1). 
      compute j=j+1.
    end if.
  end loop.
  compute resultm=resultm(:,1:(j-1)).
  save resultm/outfile = *.
end if.

do if (criterr=0).
  /* do if (savediag=1 and ydich=0 and ((anymod2=0) or (model < 4))) */.     
  /* new here */.
  do if (savediag=1 and ydich=0 and triprob=0).  
    /* new here */.
    do if (nms > 0 and anymod2 = 0) .
      loop i = 1 to ((nms+nys)-1).
        compute dfindfx={dfindfx;dfindfx(1:n,:)}.
      end loop.
      compute regdiag={regdiag, dfindfx(:,2:ncol(dfindfx))}.
      compute regdiagl={regdiagl,t(dfinds(1:(ncol(dfindfx)-1),1))}. 
    end if.
    do if (sv4match=0).
      save regdiag/outfile=*/names=regdiagl.
        /* snip */.
        compute notecode(notes,1) = 41.
        compute notes = notes + 1.
    end if.
  end if.
end if.

/* PRINT MODEL MATRICES */.
do if (criterr = 0 and matrices=1 and outscree=1).
  print/title = "************************ MODEL DEFINITION MATRICES ************************".
  print/title = "FROM variables are columns, TO variables are rows.".
  compute temp2=make(nrow(bcmat),ncol(bcmat),"0").
  loop i = 2 to nrow(bcmat).
    loop j = 1 to (ncol(bcmat)-1).
      do if (bcmat(i,j)=1).
        compute temp2(i,j)="1".
      end if.
      do if (j >= i).
        compute temp2(i,j)=" ".
      end if.
    end loop.
  end loop.
  compute temp2=temp2(2:nrow(bcmat),(1:(ncol(bcmat)-1))).
  do if (nms > 0).
    compute cmatlabs={xnames,mnames}.
    compute rmatlabs={mnames,ynames}.
  end if.
  do if (nms = 0).
    compute cmatlabs={xnames}.
    compute rmatlabs={ynames}.
  end if.
  print temp2/title="BMATRIX: Paths freely estimated (1) and fixed to zero (0):"/cnames=cmatlabs/rnames=rmatlabs/format A3.
  compute z=0.
  do if (rsum(csum(wcmat))<>0).
    compute temp2=make(nrow(wcmat),ncol(wcmat),"0").
    loop i = 2 to nrow(wcmat).
      loop j = 1 to (ncol(wcmat)-1).
        do if (wcmat(i,j)=1).
          compute temp2(i,j)="1".
        end if.
        do if (j >= i).
          compute temp2(i,j)=" ".
        end if.
      end loop.
    end loop.
    compute temp2=temp2(2:nrow(wcmat),(1:(ncol(wcmat)-1))).
    print temp2/title="WMATRIX: Paths moderated (1) and not moderated (0) by W:"/cnames=cmatlabs/rnames=rmatlabs/format A3.
  end if.
  do if (rsum(csum(zcmat))<>0).
    compute temp2=make(nrow(zcmat),ncol(zcmat),"0").
    loop i = 2 to nrow(zcmat).
      loop j = 1 to (ncol(zcmat)-1).
        do if (zcmat(i,j)=1).
          compute temp2(i,j)="1".
        end if.
        do if (j >= i).
          compute temp2(i,j)=" ".
        end if.
      end loop.
    end loop.
    compute temp2=temp2(2:nrow(zcmat),(1:(ncol(zcmat)-1))).
    print temp2/title="ZMATRIX: Paths moderated (1) and not moderated (0) by Z:"/cnames=cmatlabs/rnames=rmatlabs/format a3.
  end if.
  do if (rsum(csum(wzcmat))<>0).
    compute temp2=make(nrow(wzcmat),ncol(wzcmat),"0").
    loop i = 2 to nrow(wzcmat).
      loop j = 1 to (ncol(wzcmat)-1).
        do if (wzcmat(i,j)=1).
          compute temp2(i,j)="1".
        end if.
        do if (j >= i).
          compute temp2(i,j)=" ".
        end if.
      end loop.
    end loop.
    compute temp2=temp2(2:nrow(wzcmat),(1:(ncol(wzcmat)-1))).
    print temp2/title="WZMATRIX: W moderated paths moderated (1) and not moderated (0) by Z:"/cnames=cmatlabs/rnames=rmatlabs/format a3.
  end if.
  do if (ncs > 0).
    print ccmat/title="CMATRIX: Covariates (columns) in (1) and not in (0) the models of M and Y (rows):"/rnames=rmatlabs/cnames=covnames.
  end if.
end if.


/* PRINT NOTES */.
do if (outscree=1).
  do if (activate=0).
    print/title = "*********************** ANALYSIS NOTES AND ERRORS ************************".
  end if.
end if.
do if (activate=1).
  print/title = "**************************************************************************".
end if.
do if (criterr=0).
  do if (outscree=1).
  print conf/title = "Level of confidence for all confidence intervals in output:"/format = F8.4.
  do if (boot > 0).
    do if (goodboot = boot and bc=0).
      do if (clusboot=0).
        print boot/title="Number of bootstrap samples for percentile bootstrap confidence intervals:".
      end if.
      do if (clusboot > 0).
        print boot/title="Number of bootstrap samples for (cluster) percentile bootstrap confidence intervals:".
      end if.
    end if.
    do if (goodboot = boot and bc=1).
      do if (clusboot=0).
        print boot/title="Number of bootstrap samples for bias-corrected bootstrap confidence intervals:".
      end if.
      do if (clusboot > 0).
        print boot/title="Number of bootstrap samples for (cluster) bias-corrected bootstrap confidence intervals:".
      end if.
    end if.
    do if (booterr = 1).
      compute badend = badend(1,2:ncol(badend)).
      print/title = "WARNING: Bootstrap CI endpoints below not trustworthy. Decrease confidence".
      print badend/title="or increase the number of bootstrap samples."/space=0/format = !decimals.
    end if.
  end if.
  do if (mc > 0).
    print mc/title="Number of samples for Monte Carlo confidence intervals:".
  end if.
  do if (wnotev > 0 and printw=1).
    do if (wnotev=1).
      print/title = "W values in conditional tables are the 16th, 50th, and 84th percentiles.". 
    else if (wnotev=2).
      do if (minwwarn=0 and maxwwarn=0).
        print/title = "W values in conditional tables are the mean and +/- SD from the mean.". 
      end if.
      do if (minwwarn=1).
        print/title = "W values in conditional tables are the minimum, the mean, and 1 SD above the mean.". 
      end if.
      do if (maxwwarn=1).
        print/title = "W values in conditional tables are 1 SD below the mean, the mean, and the maximum.". 
      end if.
    end if.
  end if.
  do if (znotev > 0 and printz=1).
    do if (znotev=1).
      print/title = "Z values in conditional tables are the 16th, 50th, and 84th percentiles.". 
    else if (znotev=2). 
      do if (minzwarn=0 and maxzwarn=0).
        print/title = "Z values in conditional tables are the mean and +/- SD from the mean.". 
      end if.
      do if (minzwarn=1).
        print/title = "Z values in conditional tables are the minimum, the mean, and 1 SD above the mean.". 
      end if.
      do if (maxzwarn=1).
        print/title = "Z values in conditional tables are 1 SD below the mean, the mean, and the maximum.". 
      end if.
    end if.
  end if.
  do if (minwwarn > 0).
    print/title = "NOTE: One SD below the mean is below the minimum observed in the data for W,".
    print/title = "      so the minimum measurement on W is used for conditioning instead."/space=0.
  end if.
  do if (maxwwarn > 0).
    print/title = "NOTE: One SD above the mean is above the maximum observed in the data for W,".
    print/title = "      so the maximum measurement for W is used for conditioning instead."/space=0.
  end if.
  do if (minzwarn > 0).
    print/title = "NOTE: One SD below the mean is below the minimum observed in the data for Z,".
    print/title = "      so the minimum measurement for Z is used for conditioning instead."/space=0.
  end if.
  do if (maxzwarn > 0).
    print/title = "NOTE: One SD above the mean is above the maximum observed in the data for Z,".
    print/title = "      so the maximum measurement for Z is used for conditioning instead."/space=0.
  end if.
  /* do if (pstog=1) */.
  /*  print/title= "NOTE: Standardized coefficients for dichotomous or multicategorical X are in" */.
  /*  print/title= "      partially standardized form."/space=0 */.
  /* end if */.
  loop i = 1 to 100.
    do if (notecode(i,1) = 32).
      print/title = "Direct, indirect, and total effects are counterfactually defined".
      do if (xcontcf=0).
        print/title = "assuming X by M interaction."/space=0. 
      end if.
      do if (xcontcf=1).
         print/title = "assuming X by M interaction and with the following reference (x_ref)"/space=0. 
        compute states=t(xrefvals).
        print states/title = "and counterfactual (x_cf) states for X:"/rlabels="x_ref" "x_cf"/format=!decimals/space=0.
      end if.
    end if.
    do if (notecode(i,1)=41).
      print/title="Regression diagnostic statistics were saved to a file.".
    end if.
    do if (notecode(i,1)=1).
      print/title = "NOTE: COVMY is ignored when using CMATRIX option.".
    end if.
    do if (notecode(i,1)=2).
      print/title = "NOTE: Confidence level restricted to between 50 and 99.9999%.  95% confidence is provided in output".
    end if.
    do if (notecode(i,1)=3).
      print centvar/title = "NOTE: The following variables were mean centered prior to analysis:"/format = a8.
    end if.
    do if (notecode(i,1) = 4).
      do if (robustse=0).
        print/title = "NOTE: A heteroskedasticity-consistent standard error and covariance matrix".
        print/title = "      estimator was used."/space=0.
      end if.
      do if (robustse=1).
        print/title = "NOTE: A cluster robust standard error and covariance matrix estimator was used.".
        print dfres/title = "NOTE: Degrees of freedom for t statistics and denominator of F tests is".
        do if (intncd=1).
          print/title="NOTE: Some F tests could not be completed due to an insufficient number of clusters.".
          print/title="      These are displayed in the output as 99999."/space=0.
        end if.
      end if.
    end if.
    do if (notecode(i,1) = 5).
      print/title = "NOTE: The HC3 option has been replaced with HC.  See the documentation.".
    end if.
    do if (notecode(i,1) = 6).
      print/title = "NOTE: Due to estimation problems, some bootstrap samples had to be replaced.".
      print badboot/title="      The number of times this happened was:"/space=0/format=F8.0.
    end if.
    do if (notecode(i,1) = 7).
      print/title = "NOTE: The bootstrapping was not completed due to problematic bootstrap samples.".
      print/title = "      Bootstrap confidence intervals are therefore suppressed."/space=0.
    end if.
    do if (notecode(i,1) = 8).
      print/title = "NOTE: The number of bootstrap samples was adjusted upward given your desired confidence.". 
    end if.
    do if (notecode(i,1) = 9).
      print/title = "NOTE: WMODVAL is ignored when W is specified as multicategorical.". 
    end if.
    do if (notecode(i,1) = 10).
      print/title = "NOTE: ZMODVAL is ignored when Z is specified as multicategorical.". 
    end if.
    do if (notecode(i,1) = 11).
      print/title = "NOTE: Total effect model generated only when all covariates are specified in all".
      print/title = "      models of M and Y."/space=0. 
    end if.
    do if (notecode(i,1) = 30).
      print/title = "NOTE: Your vector of linear hypothesis weights is of the wrong length for this model.".
    end if.
    do if (notecode(i,1) = 12).
      print/title = "NOTE: Total effect model generated only when X is freely estimated to affect each M".
      print/title = "      and both X and M are freely estimated to affect Y"/space=0. 
    end if.
    do if (notecode(i,1) = 13).
      print/title = "NOTE: There are too many pairwise contrasts to conduct with this model.". 
    end if.
    do if (notecode(i,1) = 40).
        print/title = "NOTE: All subsets regression and dominance analysis are not available with fewer".
        print/title = "      than 2 or more than 15 regressors."/space=0.
    end if.
    do if (notecode(i,1) = 14).
      print/title = "NOTE: The number of contrast weights must equal the number of indirect effects.". 
    end if.
    do if (notecode(i,1) = 15).
      print/title = "NOTE: Monte Carlo confidence intervals not available for this model.".
      print/title = "      Bootstrapping is used instead."/space=0.
    end if.
    do if (notecode(i,1) = 16).
      print/title = "NOTE: The number of Monte Carlo samples was adjusted upward given your desired confidence.". 
    end if.
    do if (notecode(i,1) = 19).
      print/title = "NOTE: Your contrast matrix is invalid or not applicable to this model.". 
    end if.
    do if (notecode(i,1) = 20).
      print/title = "NOTE: One of the groups specified by your contrast matrix does not exist in the data.". 
    end if.
    do if (notecode(i,1) = 21).
      print/title = "NOTE: The VARORDER option is not available in this release.". 
    end if.
    do if (notecode(i,1) = 22).
      print/title = "NOTE: The VMODVAL and QMODVAL options are not available in this release.". 
    end if.
    do if (notecode(i,1) = 23).
      print/title = "NOTE: The QUANTILE option is not available in this release.". 
    end if.
    do if (notecode(i,1) = 24).
      print/title = "NOTE: Total effect model not available with dichotomous Y.". 
    end if.
    do if (notecode(i,1) = 25).
      print/title = "NOTE: STAND/EFFSIZE options not available with dichotomous Y.". 
    end if.
    do if ((notecode(i,1) = 26) and nms > 0).
      print/title = "NOTE: Direct and indirect effects of X on Y are on a log-odds metric.". 
    end if.
    do if (notecode(i,1) = 27).
      print/title = "NOTE: Regression diagnostics are not available for models with a dichotmous Y. ". 
    end if.
    do if (notecode(i,1) = 28).
      print/title = "NOTE: The contrast option is not available with a multicategorical X.". 
    end if.
    do if (notecode(i,1) = 31).
      do if (nms > 1).
        print medmeans/title = "NOTE: Controlled direct effect(s) estimated at the following mediator values:"/cnames=mnames/format=!decimals. 
      end if.
      do if (nms=1).
        print medmeans/title = "NOTE: Controlled direct effect(s) estimated at the following mediator value:"/format=!decimals.      
      end if.
    end if.
    do if (notecode(i,1) = 33).
      print/title = "NOTE: Sobel test is not available with the XMINT option.". 
    end if.
    do if (notecode(i,1) = 37).
      print/title = "NOTE: The CONTRAST option is not available when using the XMINT option.". 
    end if.
    do if (notecode(i,1) = 34).
      print/title = "NOTE: Standardized effects are not available when using the XMINT option.". 
    end if.
    do if (notecode(i,1) = 36).
      print/title = "NOTE: The XREFVAL option is ignored when X is declared as multicategorical.". 
    end if.
    do if (notecode(i,1) = 38).
      print/title = "NOTE: One or more cases were deleted prior to analysis at your request.". 
    end if.
    do if (notecode(i,1) = 35).
       print coval/title = "NOTE: Counterfactual effects estimated at the following covariate values:"/cnames=covnames/format=!decimals. 
    end if.
    do if notecode(i,1)=39.
      print/title = "NOTE: Due to a naming conflict, regression diagnostics were not saved.".
    end if.  
    do if (notecode(i,1) = 40).
      print/title = "NOTE: Not all regression diagnostics could be computed. As a result, some".
      print/title = "      diagnostic tests are not provided and a diagnostics file is not saved."/space=0.
    end if.
    do if (notecode(i,1) = 29 and (listmiss=1)).
      compute a=missrow.
      compute conum=ncol(a).
      compute allgood=0.
      compute smremain=12.
      compute largesti=1.
      compute smallrow=0.
      do if (conum > 12).
        loop ii = 1 to 12.
          compute check=(conum/ii).
          do if (check = trunc(check)).
            compute check2=conum/ii.
            compute aok=ii.
            do if (aok > 2).
              compute allgood=1.
            end if.
          else.
            compute remain=conum-(ii*trunc(check)).
            do if (remain <= smremain).
              compute smremain=remain.
              compute largesti=ii.
              compute smallrow=trunc(conum/largesti).            
            end if.
          end if.
        end loop.
        compute atemp=a(1,1:(aok*check2)).
        compute atemp=reshape(atemp,check2,aok).
        do if (ncol(atemp) > 2).
          print/title="NOTE: Missing data resulted in the deletion of the following row(s) of".
          print atemp/title="      data from the analysis:"/format=F6.0/space=0.
        end if.
        do if (allgood=0).
          compute atemp=a(1,1:(smallrow*largesti)).
          compute atemp=reshape(atemp,smallrow,largesti).
          compute btemp=a(1,((largesti*smallrow)+1):conum).
          print/title="NOTE: Missing data resulted in the deletion of the following row(s) of".
          print atemp/title="      data from the analysis:"/format=F6.0/space=0.
          print btemp/title=" "/space=0/format=f6.0.
        end if.
      else if (conum <= 12).
        print/title="NOTE: Missing data resulted in the deletion of the following row(s) of".
        print a/title="      data from the analysis:"/format=F6.0/space=0.
      end if.
    end if.
  end loop.
  do if (eivdo=1).
    compute varneiv=varnames(1,2:ncol(varnames)).
    compute eivrel2=eivrel(1,nxvls:ncol(eivrel)).
    do if (ncol(eivrel2) < 9).
      print eivrel2/title="NOTE: This errors-in-variables analysis assumes the following reliabilities:"/cnames=varneiv/format=F8.4.
    else.
      compute eivrel2=t(eivrel2).
      print eivrel2/title="NOTE: This errors-in-variables analysis assumes the following reliabilities:"/rnames=varneiv/format=F8.4.
    end if.
  end if.
  do if (toomany=1).
    print/title="WARNING: Variables names longer than eight characters can produce incorrect output".
    print/title="when some variables in the data file have the same first eight characters. Shorter"/space=0.
    print/title="variable names are recommended. By using this output, you are accepting all risk"/space=0.
    print/title="and consequences of interpreting or reporting results that may be incorrect."/space=0.
  end if.
  end if.
end if.


/* PRINT ERRORS */.
loop i = 1 to 100.
  do if (errcode(i,1)=1).
    print/title = "ERROR: You must specify a Y and an X variable.".
  end if.
  do if (errcode(i,1)=2).
    print/title = "ERROR: X, M, or Y variable used more than once or W and Z are the same variable.".  
    do if (toomany = 1).
      print/title = "       This could be caused by the use of variables names longer"/space=0.
      print varnames/title = "       than eight characters. Here are the variables I see:"/space=0/format=A8.
    end if.
  end if.
  do if (errcode(i,1)=3).
    print/title = "ERROR: This model accepts only one variable for W, Y, X, Z, or Cluster".
  end if.
  do if (errcode(i,1)=4).
    print/title = "ERROR: A variable specified as multicategorical has more than nine categories.".
  end if.
  do if (errcode(i,1)=5).
    print/title = "ERROR: One of the categories contains only a single case.".
  end if.
  do if (errcode(i,1)=6).
    print/title = "ERROR: Invalid model number in this version of PROCESS.".
  end if.
  do if (errcode(i,1)=7).
    print/title = "ERROR: Invalid model number.".
  end if.
  do if (errcode(i,1)=8).
    print/title = "ERROR: You must specify an M variable for this model.".
  end if.
  do if (errcode(i,1)=9).
    print/title = "ERROR: You have specified an M variable in a model that does not use it.".
  end if.
  do if (errcode(i,1)=10).
    print/title = "ERROR: You have specified a W variable in a model that does not use it.".
  end if.
  do if (errcode(i,1)=11).
    print/title = "ERROR: You have not specified a W variable in a model that requires it.".
  end if.
  do if (errcode(i,1)=12).
    print/title = "ERROR: You have specified a Z variable in a model that does not use it.".
  end if.
  do if (errcode(i,1)=13).
    print/title = "ERROR: You have not specified a Z variable in a model that requires it.".
  end if.
  do if (errcode(i,1)=14).
    print/title = "ERROR: V and Q are not proper specifications in this release of PROCESS.".
    print/title = "       Moderators must be specified as W and/or Z."/space=0.
  end if.
  do if (errcode(i,1)=15).
    print/title = "ERROR: One of your model variables exhibits no variation (it is a constant).".
  end if.
  do if (errcode(i,1)=16).
    print/title = "ERROR: BMATRIX is not the correct length or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=17).
    print/title = "ERROR: WMATRIX is not the correct length or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=18).
    print/title = "ERROR: ZMATRIX is not the correct length or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=19).
    print/title = "ERROR: WZMATRIX is not the correct length or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=20).
    print/title = "ERROR: A path fixed at zero cannot be moderated.".
  end if.
  do if (errcode(i,1)=21).
    print/title = "ERROR: If only one moderator is specified, it must be specified as W.".
  end if.
  do if (errcode(i,1)=22).
    print/title = "ERROR: In BMATRIX, X must be specified to affect at least one variable.".
  end if.
  do if (errcode(i,1)=23).
    print/title = "ERROR: In BMATRIX, at least one variable must be specified to affect Y.".
  end if.
  do if (errcode(i,1)=24).
    print/title = "ERROR: You must specify a model number or a custom BMATRIX specification.".
  end if.
  do if (errcode(i,1)=25).
    print/title = "ERROR: BMATRIX cannot be used in conjunction with a model number.".
  end if.
  do if (errcode(i,1)=26).
    print/title = "ERROR: Your model has a dangling mediator (all Ms must affect and be affected).".
  end if.
  do if (errcode(i,1)=27).
    print/title = "ERROR: CLUSTER is not available on this release of PROCESS.".
  end if.
  do if (errcode(i,1)=29).
    print/title = "ERROR: CMATRIX is not the correct length or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=30).
    print/title = "ERROR: In CMATRIX, all covariates must be assigned to an M or a Y.".
  end if.
  do if (errcode(i,1)=31).
    print/title = "ERROR: A linear or near linear dependency (singularity) exists in the data.".
  end if.
  do if (errcode(i,1)=32).
    print/title = "ERROR: Models 80 and 81 require between 3 and 6 mediators.".
  end if.
  do if (errcode(i,1)=33).
    print/title = "ERROR: Model 82 requires 4 mediators.".
  end if.
  do if (errcode(i,1)=48 and v2tag=0).
    print/title = "ERROR: You are using outdated syntax. The current version of PROCESS".
    print/title = "       is documented in Appendices A and B of www.guilford.com/p/hayes3"/space=0.
    compute v2tag=1. 
  end if.
  do if (errcode(i,1)=34).
    print/title = "ERROR: This model number requires between 2 and 6 mediators.".
  end if.
  do if (errcode(i,1)=35).
    print/title = "ERROR: In a model with only one moderator, that moderator must be W.".
  end if.
  do if (errcode(i,1)=36).
    print/title = "ERROR: A serial mediation model cannot have more than 6 mediators.".
  end if.
  do if (errcode(i,1)=37).
    print/title = "ERROR: No more than 10 mediators are allowed in a PROCESS command.".
  end if.
  do if (errcode(i,1)=38).
    print/title = "ERROR: XCATCODE is not provided, not the correct length, or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=39).
    print/title = "ERROR: WCATCODE is not provided, not the correct length, or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=40).
    print/title = "ERROR: ZCATCODE is not provided, not the correct length, or is otherwise invalid.".
  end if.
  do if (errcode(i,1)=41).
    print/title = "ERROR: Models 1, 2, and 3 cannot be customized.".
  end if.
  do if (errcode(i,1)=42).
    print/title = "ERROR: WS option available only in PROCESS v2. Or use the MEMORE macro instead.".
    print/title = "       MEMORE can be downloaded from www.akmontoya.com."/space=0.
  end if.
  do if (errcode(i,1)=43).
    print/title = "ERROR: PROCESS does not allow dichotomous mediators.".
  end if.
  do if (errcode(i,1)=60).
    print/title = "PROCESS is now ready for use.".   
    print/title="Copyright 2013-2025 by Andrew F. Hayes. ALL RIGHTS RESERVED.". 
    print/title="This version of PROCESS requires SPSS version 26 or later"/space=0.
    print/title="Workshop schedule available at haskayne.ucalgary.ca/CCRAM"/space=0.
    print/title="In SPSS 29 and later, change default output font to Courier New for tidier"/space=0.
    print/title="output. More information about PROCESS at processmacro.org/faq.html."/space=0.
  end if. 
  do if (errcode(i,1)=50).
    print/title = "ERROR: A multicategorical moderator cannot be specified as a covariate.".
  end if.
  do if (errcode(i,1)=51).
    print/title = "ERROR: A variable you specified as a covariate is a moderator in all equations.".
  end if.
  do if (errcode(i,1)=52 and mcerpt=0).
compute mcerpt=1.
    print/title = "ERROR: A variable specified as multicategorical must have at least three categories".
  end if.
  do if (errcode(i,1) = 47). 
    do if (iterrmod=1).
      print/title = "ERROR: Iteration for Y model didn't converge to a solution. Interpret results with caution.".
      print/title = "       Try increasing the number of iterations, though this will slow down computation."/space=0. 
    end if.
    do if (bootiter=1).
      print/title = "ERROR: Nonconvergence during bootstrapping. Interpet bootstrap results with caution.".
    end if.
  end if.
  do if (errcode(i,1) = 61).
      print/title = "ERROR: Variable names should be no more than eight characters in length.".       
  end if.
  do if (errcode(i,1) = 62).
      print/title = "ERROR: After listwise deletion of cases with missing data, too few cases remain.".       
  end if.
  do if (errcode(i,1) = 63).
      print/title = "ERROR: The XMINT option is available only for model 4.".       
  end if.
   do if (errcode(i,1) = 64).
      print/title = "ERROR: Incorrect number of values specified in CDEVAL option.".       
  end if.
  do if (errcode(i,1) = 65).
      print/title = "ERROR: Only indicator or sequential coding of X is allowed with the XMINT option.".       
  end if.
  do if (errcode(i,1) = 66).
    print/title = "ERROR: A reference value of X is required for this model.".       
  end if.
  do if (errcode(i,1) = 67).
    print/title = "ERROR: Too many elements provided in XREFVAL option.".       
  end if.
  do if (errcode(i,1)=68).
    print/title = "ERROR: Covariate assignment is not allowed with XMINT option.".
  end if.
  do if (errcode(i,1) = 69).
      print/title = "ERROR: Incorrect number of values specified in COVAL option.".       
  end if.
  do if (errcode(i,1) = 70).
      print/title = "ERROR: Incorrect value(s) in XREFVAL for this dichotomous X variable.".       
  end if.
  do if (errcode(i,1) = 71).
      print/title = "ERROR: The CENTER option is not available when using the XMINT option.".       
  end if.
  do if (errcode(i,1) = 72).
      print/title = "ERROR: The XMINT option is not available for models with a dichotomous Y.".       
  end if.
  do if errcode(i,1)=73.
    print/title = "ERROR: The SPLINE option allows no more than 8 joints.".
  end if.
  do if errcode(i,1)=74.
    print/title = "ERROR: Spline joint locations must be listed in ascending order with no ties.".
  end if.
  do if errcode(i,1)=75.
    print/title = "ERROR: All spline joint locations must be larger than the minimum observed value.".
  end if.
  do if errcode(i,1)=76.
    print/title = "ERROR: Each spline segment must contain at least two cases.".
  end if.
  do if errcode(i,1)=77.
    print/title = "ERROR: The SPLINE option is not available for this model.".
  end if.  
  do if errcode(i,1)=78.
    print/title = "ERROR: CLUSTER requires use of the ROBUSTSE or CLUSBOOT options.".
  end if.  
  do if errcode(i,1)=79.
    print/title = "ERROR: ROBUSTSE and CLUSBOOT options require a clustering variable.".
  end if.  
  do if errcode(i,1)=80.
    print/title = "ERROR: A variable in your PROCESS command does not exist in the data.".
  end if.  
  do if errcode(i,1)=81.
    print/title = "ERROR: Cluster robust standard errors not available with a dichotomous Y. ".
  end if.  
  do if errcode(i,1)=82.
    print/title = "ERROR: Errors-in-variables estimation not allowed for this model.".
  end if.  
  do if errcode(i,1)=83.
    print/title = "ERROR: Incorrect number of reliabilities provided.".
  end if.  
  do if errcode(i,1)=84.
    print/title = "ERROR: Reliabilites must be larger than 0 and no larger than 1.".
  end if.  
  do if errcode(i,1)=85.
    print/title = "ERROR: One of the assumed reliabilities is too small to estimate the model.".
  end if.  
  do if errcode(i,1)=86.
    print/title = "ERROR: Errors-in-variables estimation not available with dichotomous Y.".
  end if. 
  do if errcode(i,1)=87.
    print/title = "ERROR: The ROBUSTSE and HC options are not available with errors-in".
    print/title="       variables estimation. Use the eiv option or let it default to eiv=3."/space=0.
  end if. 
  do if errcode(i,1)=88.
    print/title = "ERROR: Invalid argument for EIV option.".
  end if. 
end loop.
END MATRIX.
set printback=on.
restore.
!ENDDEFINE.
process activate=1.

