/* MODPROBE for SPSS version 2.0 */.
/* This SPSS macro is used to estimate and probe two-way interactions */.
/* in OLS and logistic regression models */.
/* Written by Andrew F. Hayes */.
/* Department of Psychology */.
/* The Ohio State University */.
/* hayes.338@osu.edu */. 
/* Copyright 2015 */.
/* Features for estimating and probing three-way interactions */.
/* are available in PROCESS.  See http://www.guilford.com/p/hayes3 */.

set printback=off.

/* This code should not be posted online except through afhayes.com */.
/* without written permission. Commercial distribution is not authorized */.


DEFINE LOGITER (pt1lp=!charend ('/')/xlp=!charend('/')/ylp=!charend('/')/iterate=!charend('/')/converge=!charend('/')).
compute pt1lp=!pt1lp.
compute xlp=!xlp.
compute ylp=!ylp.
loop jjj = 1 to !iterate.
  compute vt1 = mdiag(pt1lp&*(1-pt1lp)).
  compute b = bt1+inv(t(xlp)*vt1*xlp)*t(xlp)*(ylp-pt1lp).
  compute pt1lp = 1/(1+exp(-(xlp*b))).
  compute itprob = csum((pt1lp < .00000000000001) or (pt1lp > .99999999999999)).
  do if (itprob > 0).
    loop kkk = 1 to nrow(pt1lp).
      do if (pt1lp(kkk,1) = 1).
        compute pt1lp(kkk,1) = .99999999999999.
      end if.
      do if (pt1lp(kkk,1) = 0).
        compute pt1lp(kkk,1) = .00000000000001.
      end if.
    end loop.
    compute itprob = 0.
  end if.
  do if (itprob = 0).
    compute LL = ylp&*ln(pt1lp)+(1-ylp)&*ln(1-pt1lp).
    compute LL2 = -2*csum(ll).
  end if.
  do if (abs(LL1-LL2) < !converge).
    compute vt1 = mdiag(pt1lp&*(1-pt1lp)).
    compute varb = inv(t(xlp)*vt1*xlp).
    compute seb = sqrt(diag(varb)).
    break.
  end if.
  compute bt1 = b.
  compute LL1 = LL2.
end loop.
!ENDDEFINE.

DEFINE MODPROBE (y = !charend ('/')/x = !charend ('/')/modval = !charend ('/') !default (9999)/jn = !charend ('/') !default (0)/alpha = !charend ('/') !default (.05)
  /center = !charend('/') !default (0)/est = !charend ('/') !default (0)/iterate = !charend ('/') !default(10000)
  /converge = !charend ('/') !default(.0000001)/change = !charend ('/') !default(1)/mcfoc=!charend('/') !default(0)/mcmod = !charend('/') !default(0)
  /nomod=!charend('/') !default(0)/save=!charend('/') !default(0)/decimals = !charend('/') !default(F10.4)/ptiles=!charend('/') !default(0)/hc3=!charend('/') !default(0)).
PRESERVE.
set mxloop = 1000000.
set printback = off.
matrix.
get dd/variables = !y !x/names = nm/MISSING = OMIT.
get x/variables = !x/names=tpx/missing=omit.
get tempy/variables = !y/names=tpy/missing=omit.
compute n = nrow(dd).
compute nx=1.
compute criterr=0.
compute errs=0.
compute errsm=make(10,1,0).
compute nyv=ncol(tpy).
compute nxv=ncol(x).
compute nomod=(trunc(!nomod) =1).
compute savplot=(!save = 1).
compute ptiles=(!ptiles=1).
compute hc3=(!hc3=1).
do if (nyv <> 1).
  compute criterr=1.
  compute errs=errs+1.
  compute errsm(errs,1)=1.
end if.
do if (nxv < 2).
  compute criterr=1.
  compute errs=errs+1.
  compute errsm(errs,1)=2.
end if.
compute desc1 = (nrow(dd)*sscp(dd))-(t(csum(dd))*(csum(dd))).
compute desc1 = desc1/(nrow(dd)*(nrow(dd)-1)).
compute desc2=csum(diag(desc1) = 0).
do if (desc2 > 0).
  compute criterr=1.
  compute errs=errs+1.
  compute errsm(errs,1)=5.
end if.
compute ddd = {"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"}.
compute ddd1 = {"int_1", "int_2", "int_3", "int_4", "int_5", "int_6", "int_7", "int_8", "int_9"}.
compute dumok = 0.
compute mcfoc=trunc(!mcfoc).
compute mcmod=trunc(!mcmod).
do if (mcfoc < 0).
  compute mcfoc=0.
end if.
do if (mcfoc > 4).
  compute mcfoc=0.
end if.
do if (mcmod < 0).
  compute mcmod=0.
end if.
do if (mcmod > 4).
  compute mcmod=0.
end if.
do if (mcmod > 0 and mcfoc>0).
  compute errs=errs+1.
  compute errsm(errs,1)=3.
  compute criterr=1.
end if.
compute mcloc=(mcfoc>0).
compute con = make(n,1,1).
compute ncovs=ncol(x)-2.
do if (mcfoc > 0 or mcmod > 0).
  compute temp = dd.
  compute temp(GRADE(dd(:,(ncol(dd)-mcloc))),:) = dd.
  compute dd = temp.
  compute dummy = design(dd(:,ncol(dd)-mcloc)).
  compute nvls = ncol(dummy).
  compute nnvls = csum(dummy).
  compute toosmall=rsum(nnvls < 2).
  compute mnvls = cmin(t(nnvls)).
  do if (rsum(nnvls < 2)) > 0).
    compute criterr = 1.
    compute errs=errs+1.
    compute errsm(errs,1)=6.
  end if.
  do if (nvls > 10).
    compute criterr = 1.
    compute errs=errs+1.
    compute errsm(errs,1)=4.
  end if.
  do if (criterr=0).
    compute dumok = 1.
    compute nnvls=make(nvls,1,0).
    compute nnvls(1,1)=dd(1,ncol(dd)-mcloc).
    compute temp = 2.
    loop i = 2 to n.
      do if (dd(i,ncol(dd)-mcloc) <> nnvls((temp-1),1)).
        compute nnvls(temp,1)=dd(i,ncol(dd)-mcloc).
        compute temp = temp+1.
      end if.
    end loop.
    compute dummy = dummy(:,2:ncol(dummy)).
    do if (mcfoc=4 or mcmod=4).
      compute minus1=make(1,ncol(dummy),-1).
      loop k = 1 to n.
        do if (rsum(dummy(k,:)) = 0).
          compute dummy(k,:)=minus1.
        end if.
      end loop.
    end if.
    do if (mcfoc = 2 or mcfoc = 3) or (mcmod = 2 or mcmod = 3)).
      loop k = 1 to n.
        do if (rsum(dummy(k,:)) > 0).
          loop i = 1 to ncol(dummy).
            do if (dummy(k,i) = 0).
              compute dummy(k,i) = 1.
            else.
              break.
            end if.
          end loop.
        end if.
      end loop.
      do if (mcfoc = 3 or mcmod=3).
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
          compute dummy(k,:)=conmat1((rsum(dummy(k,:))+1),:). 
        end loop.
      end if.
    end if.
    compute nx = ncol(dummy).
    compute xname = ddd(1,1:nx).
    compute xdata = dummy.
    compute xname = ddd(1,1:nx).
    compute xname2={nm(1,(ncol(nm)-mcloc)), xname}.   
    compute indlbs = t(xname).
    compute dummat = make((nx+1),nx,0).
    compute dummat((2:nrow(dummat)),:)=ident(nx).
    do if (mcfoc = 2 or mcmod=2).
      loop i = 2 to nrow(dummat).
        loop j = 1 to (i-1).
          compute dummat(i,j) = 1.
        end loop.
      end loop.
    end if.
    do if (mcfoc = 3 or mcmod=3).
      compute dummat=conmat1.
    end if.
    do if (mcfoc = 4 or mcmod = 4).
       compute dummat(1,:)=minus1.
    end if.
    compute dummat={nnvls, dummat}.
  end if.
end if.

print/title = "**************** MODPROBE Procedure for SPSS Release 2.00 *****************"/space=newpage.
print/title = "          Written by Andrew F. Hayes, Ph.D.       www.afhayes.com".
print/title = "***************************************************************************".
do if (criterr=0).
  compute y = dd(:,1).
  compute ovals = ncol(design(dd(:,1))).
  compute jnerr = 0.
  compute itprob = 0.
  do if (ovals = 2).
    compute omx = cmax(y(:,1)).
    compute omn = cmin(y(:,1)).
    compute y(:,1) = (y(:,1) = omx).
    compute rcd = {omn, 0; omx, 1}.
  end if.
  compute jn = !jn.
  do if (ovals = 2 and jn > 1).
    compute jn = 1.
    compute jnerr = 1.
  end if.
  compute xd = 1.959963984540054.
  compute cilm = .05.
  compute alperr = 0.
  compute sstotal = csum((y-(csum(y)/n))&**2).
  compute outv = t(nm(1,1)).
  compute mdtr = nm(1,ncol(dd)).
  compute fciv = nm(1,(ncol(dd)-1)).
  compute centerv={" "}.
  do if (!center = 1 and mcmod = 0).
    compute dd(:,ncol(dd)) = dd(:,ncol(dd))-(csum(dd(:,ncol(dd)))/n).
    compute centerv={centerv,mdtr}.
  end if.
  do if (!center = 1 and mcfoc=0).
    compute dd(:,(ncol(dd)-1)) = dd(:,(ncol(dd)-1))-(csum(dd(:,(ncol(dd)-1)))/n).
    compute centerv={centerv, fciv}.
  end if.
  do if (mcfoc=0 and mcmod=0).
    compute inter = dd(:,(ncol(dd)-1))&*dd(:,ncol(dd)).
    compute x = {con,dd(:,2:ncol(dd))}.
    do if (nomod = 0).
      compute x={x,inter}.
    end if.
    compute nms = t(nm(1,2:ncol(dd))).
  end if.
  do if (mcfoc <> 0 or mcmod <> 0).
    compute inter=make(n,nx,0).
    loop i = 1 to nx.
      compute inter(:,i)=dummy(:,i)&*dd(:,ncol(dd)-(1-mcloc)).
    end loop.
    do if (ncovs = 0).
      compute x={con,dd(:,ncol(dd)-(1-mcloc)),dummy}.
      do if (nomod = 0).
        compute x={x,inter}.
      end if.
      compute nms={nm(1,ncol(dd)-(1-mcloc)); t(ddd(1,1:ncol(dummy)));t(ddd1(1,1:ncol(dummy)))}.
    end if.
    do if (ncovs > 0).
      compute x={con,dd(:,2:(1+ncovs)),dd(:,ncol(dd)-(1-mcloc)),dummy}.
      do if (nomod = 0).
        compute x={x,inter}.
      end if.   
      compute nms={t(nm(1,2:(1+ncovs))); nm(1,ncol(dd)-(1-mcloc)); t(ddd(1,1:ncol(dummy)));t(ddd1(1,1:ncol(dummy)))}.
    end if.
  end if.
  compute xmns = csum(x)&/n.
  compute focvals = {1,0,0}.
  compute highwarn = 0.
  compute lowwarn = 0.
  do if (ncol(x) > (2*nx+2)).
    compute covmns = {1, xmns(1,2:(ncol(x)-1-(2*nx)))}.
    compute focvals = {covmns,0,0}.
  end if.
  compute dfres = n-ncol(x).
  compute tval = 5.990731+18.568800*(1/dfres).
  do if (jn = 2 and dfres < 50).
    compute jn = 0.
    print/title = "Simultaneous inference is not available at this sample size.".
  end if.
  do if (ovals = 2).
    compute pt2 = make(n,1,(csum(y)/n)).
    compute LL3 = y&*ln(pt2)+(1-y)&*ln(1-pt2).
    compute LL3 = -2*csum(LL3).
    compute pt1 = make(n,1,0.5).
    compute bt1 = make(ncol(x),1,0).
    compute LL1 = 0.
    LOGITER pt1lp=pt1/xlp=x/ylp=y/iterate=!iterate/converge=!converge.
    do if (jjj > !iterate).
      compute itprob = 2.
    end if.
    do if (!alpha = .10).
      compute tval = 2.70554345409541221538.
      compute cilm = .10.
    else if (!alpha = .01).
      compute tval = 6.634896601021211218135.
      compute cilm = .01.
    else if (!alpha = .05).
      compute tval = 3.8414588206941250351219.
      compute cilm = .05.
    else.
      compute alperr = 1.
    end if.
  else if (ovals <> 2).
    compute invXtX=inv(t(x)*x).
    compute b = invXtX*t(x)*y.
    compute resid = y-(x*b).
    compute ssresid=csum(resid&*resid).
    compute msresid=ssresid/dfres.
    compute varb = msresid*invXtX.
    compute seint=sqrt(varb(nrow(varb),nrow(varb))).
    compute k3 = nrow(b).
    do if (hc3 = 1).
      compute xhc=x.
      compute h = xhc(:,1).
      loop i3=1 to n.
        compute h(i3,1)= xhc(i3,:)*invXtX*t(xhc(i3,:)). 
      end loop.
      loop i3=1 to k3.
        compute xhc(:,i3) = (resid(:,ncol(resid))&/(1-h))&*xhc(:,i3).
      end loop.
      compute varb=(invXtX*t(xhc)*xhc*invXtX).
    end if.
    compute seb = sqrt(diag(varb)).
    compute r2 = 1-(ssresid/sstotal).
    do if (mcfoc <> 0 or mcmod <> 0).
      compute ytp=y.
      compute xtp=x(:,1:(ncol(x)-nx)).
      compute btp = inv(t(xtp)*xtp)*t(xtp)*ytp.
      compute residt = ytp-(xtp*btp).
      compute ssresidt=csum(residt&*residt).
      compute r2noint = 1-(ssresidt/sstotal).
    end if.
    compute pr = ncol(x)-1.
    compute lmat = ident(nrow(b)).
    compute lmat = lmat(:,2:ncol(lmat)).
    compute f = (t(t(lmat)*b)*inv(t(lmat)*varb*lmat)*((t(lmat)*b)))/pr).
    compute pf = 1-fcdf(f,pr,dfres).
    compute pf = {r2,f,pr,dfres,pf, n}.
    do if (!alpha = .10).
      compute xd = 1.644853626951472.
      compute tval = 4.604839+10.883656*(1/dfres).
      compute cilm = .10.
    else if (!alpha = .01).
      compute xd = 2.5758293035489.
      compute tval = 9.207688+44.698748*(1/dfres).
      compute cilm = .01.
    else if (!alpha = .05).
      compute xd = 1.959963984540054.
      compute tval = 5.990731+18.568800*(1/dfres).
      compute cilm = .05.
    else.
      compute alperr = 1.
    end if.
  end if.
  do if (itprob = 0).
    compute tstat = b&/seb.
    do if (ovals <> 2).
      compute p = 2*(1-tcdf(abs(tstat), dfres)).
    else if (ovals = 2).
      compute p = 2*(1-cdfnorm(abs(tstat))).
    end if.
    compute outp = {b,seb,tstat,p}.
    compute nms = {"constant"; nms; "interact"}.
    do if (jn < 2 and ovals <> 2).
      compute tval =  (dfres* (exp((dfres-(5/6))*((xd/(dfres-(2/3)+(.11/dfres)))*(xd/(dfres-(2/3)+(.11/dfres)))))-1)).
      compute outp={outp,(b-sqrt(tval)&*seb),(b+sqrt(tval)&*seb)}.
    end if.
    compute bb = tval.
    print outv/title = "Outcome Variable"/format a8.
    print fciv/title = "Focal Predictor Variable"/format a8.
    print mdtr/title = "Moderator Variable"/format a8.
    do if (dumok = 1).
      do if (mcfoc > 0).
        print dummat/title = "Coding of categorical focal predictor for analysis:"/cnames = xname2/format = F5.2.
      else if (mcmod > 0).
        print dummat/title = "Coding of categorical moderator variable for analysis:"/cnames = xname2/format = F5.2.
      end if.
    end if.
    print/title = "***************************************************************************".
    do if (ovals = 2).
      compute nmsd = {outv, "Analysis"}.
      print rcd/title = "Coding of binary Y for analysis:"/cnames = nmsd/format = !decimals.
      compute outp = {outp,(outp(:,1)-sqrt(tval)&*outp(:,2)),(outp(:,1)+sqrt(tval)&*outp(:,2))}.
      compute LLdiff = LL3-LL2.
      compute pvalue=1-chicdf(LLdiff,nrow(b)).
      compute LL4 = LL2.
      compute mcF = LLdiff/LL3.
      compute cox = 1-exp(-LLdiff/n).
      compute nagel = cox/(1-exp(-(LL3)/n)).
      compute pf = {LL2, LLdiff, pvalue, mcF, cox, nagel, n}.
      print pf/title = "Logistic Regression Summary"/clabels = "-2LL" "Model LL" "p-value" "McFadden" "CoxSnell" "Nagelkrk" "n"/format !decimals.
      print outp/title = "Logistic Regression Model"/rnames = nms/clabels "Coeff" "se" "Z" "p" "LLCI" "ULCI"/format !decimals.
      compute varbtmp=varb.
      compute btmp=b.
      compute LL2f=LL2.
      compute pt1lp = make(n,1,0.5).
      compute ylp=y.
      compute xlp=x(:,1:(ncol(x)-nx)).
      compute bt1 = make(ncol(xlp),1,0).
      LOGITER pt1lp=pt1lp/xlp=xlp/ylp=ylp/iterate=!iterate/converge=!converge.
      compute LLdiff=LL2-LL2f.
      compute b=btmp.
      compute varb=varbtmp.
      compute pchi=1-chicdf(lldiff,nx).
      compute rcha={LLdiff,nx,pchi}.
    end if.
    do if (ovals <> 2).
      print pf/title = "Complete Model Regression Summary"/clabels = "R-sq" "F" "df1" "df2" "p" "n"/format !decimals.
      compute lmat=make(nrow(b),1,0).
      compute lmat(nrow(lmat),1)=1.
      compute fcha = (t(t(lmat)*b)*inv(t(lmat)*varb*lmat)*((t(lmat)*b)))/1).
      compute rcha=((b(nrow(b),1)/seint)*(b(nrow(b),1)/seint))*(1-r2)/dfres.
      compute rcha = {rcha, fcha, 1, dfres, outp((pr+1),4)}.
      print outp/title = "Regression Model"/rnames = nms/clabels "Coeff" "se" "t" "p" "LLCI" "ULCI"/format !decimals.
      do if (mcfoc > 0 or mcmod > 0).
        compute rcha=r2-r2noint.
        compute lmat=make((nrow(b)-nx),nx,0).
        compute lmat2=ident(nx).
        compute lmat={lmat;lmat2}.
        compute fcha = (t(t(lmat)*b)*inv(t(lmat)*varb*lmat)*((t(lmat)*b)))/nx).
        compute pvalr2c=1-fcdf(fcha,nx,dfres).
        compute rcha = {rcha, fcha, nx, dfres,pvalr2c}.
      end if.
    end if.
    do if (mcfoc = 0 and mcmod=0 and nomod=0).
      compute cprod = {nms((ncol(dd)-1),1), "X",  nms(ncol(dd),1)}.
      print cprod/title = "Interact is defined as:"/format = a8.
    end if.
    do if ((mcfoc > 0 or mcmod > 0) and nomod=0).
      compute intkey = {"a", "b", "c", "d", "e"}.
      loop i = 1 to nx.
        compute intkey={intkey; ddd1(1,i), " : ", ddd(1,i), " X ", nm(1,ncol(dd)-(1-mcloc))}.
      end loop.
      compute intkey=intkey(2:nrow(intkey),:).
      print intkey/title="Product terms key:"/format = a8.
    end if.
    do if (!change <> 0 and nomod = 0).
      do if (ovals <> 2).
        print rcha/title = "R-square increase due to interaction:"/clabels "R2-chng" "F" "df1" "df2" "p"/format !decimals.
      else if (ovals = 2).
        print rcha/title = "Likelihood ratio test for interaction:"/clabels "Chi-sq" "df" "p"/format !decimals.
      end if.
    end if.
    print/title = "***************************************************************************".
    do if (nomod=0).
      compute mdvar = x(:,(ncol(x)-1)).
      do if (mcfoc = 0 and mcmod = 0).
        compute g1 = b((ncol(x)-2),1).
        compute g3 = b(ncol(x),1).
        compute vg1 = varb((ncol(x)-2),(ncol(x)-2)).
        compute vg3 = varb(ncol(x),ncol(x)).
        compute covg1g3 = varb((ncol(x)-2), ncol(x)).
      end if.
      do if (mcfoc > 0 or mcmod > 0).
        compute mdvar = x(:,(ncol(x)-(2*nx))).
      end if.
      compute mdmin = cmin(mdvar).
      compute mdmax = cmax(mdvar).
      compute fvar = x(:,(ncol(x)-2)).
      compute nval = ncol(design(mdvar)).
      compute fvmin = cmin(fvar).
      compute fvmax = cmax(fvar).
      do if (!modval = 9999 and jn < 1).
        compute mnmd = csum(mdvar)/n.
        compute tmp = make(n,1,mnmd).
        compute sdmd = sqrt(csum(((mdvar-tmp)&**2))/(n-1)).
        compute probeval = {mnmd-sdmd; mnmd; mnmd+sdmd}.
        do if (ptiles = 1).
          compute tmp = mdvar.
          compute tmp(GRADE(mdvar(:,1)),:) = mdvar.
          compute mdvar = tmp.
          compute probeval={mdvar(trunc(0.25*n),1);mdvar(trunc(0.5*n),1);mdvar(trunc(0.75*n),1)}.
        end if.
      end if.
      do if (nval = 2).
        compute probeval = make(2,1,0).
        /* compute probeval(1,1)=mdvar(1,1) */.
        compute probeval(1,1)=cmin(mdvar).
        compute jn = 0.
        loop i = 1 to n.
          do if (mdvar(i,1) <> probeval(1,1)).
            compute probeval(2,1) = mdvar(i,1).
            BREAK.
          end if.
        end loop.
      end if.
      do if (!modval <> 9999 and jn < 1).
        compute probeval = !modval.
      end if.
      do if (jn < 1).
        compute outp = make(nrow(probeval),7,0).
        do if (mcfoc = 0 and mcmod=0).
          loop i = 1 to nrow(probeval).
            compute x2 = probeval(i,1).
            compute w1 = g1+g3*x2.
            compute varw1 = vg1+(2*x2*covg1g3)+((x2*x2)*vg3).
            compute sew1 = sqrt(varw1).
            compute t1 = w1/sew1.
            compute LLCI = (w1-sqrt(tval)&*sew1).
            compute ULCI = (w1+sqrt(tval)&*sew1).
            do if (ovals <> 2).
              compute p = 2*(1-tcdf(abs(t1), dfres)).
            else if (ovals = 2).
              compute p = 2*(1-cdfnorm(abs(t1))).
            end if.
            loop j = 0 to 20.
              compute temp = (fvmin+j*((fvmax-fvmin)/20)).
              do if (ncol(x) > 4).
                compute focvals = {focvals; covmns, temp, probeval(i,1)}.
              else.
                compute focvals = {focvals; 1, temp, probeval(i,1)}.
              end if.
            end loop.
            compute outp(i,:) = {x2, w1, sew1, t1, p, LLCI, ULCI}.
          end loop.
          compute focvals = focvals(2:nrow(focvals),:).
          compute inter2 = focvals(:,(ncol(focvals)-1))&*focvals(:,ncol(focvals)).
          compute focvals = {focvals, inter2}.
          compute yhat = focvals*b.
          compute focvals = {focvals(:,(ncol(focvals)-2):(ncol(focvals)-1)), yhat}.
        end if.
        do if (mcfoc > 0 or mcmod > 0).
          compute focvals=make(1,ncol(x)+1,1).
          do if (mcfoc > 0).
            print/title = "Conditional Effect of Focal Predictor at Values of the Moderator Variable".
            print/title = " "/space=0.
            compute rnn2=mdtr.
            compute matt=make(nx,6,0).
          end if.
          loop jj=1 to nrow(probeval).
            do if (mcfoc > 0).
              loop ii=1 to nx.
                compute g1=b((ncol(x)-(2*nx)+ii),1).
                compute g3=b((ncol(x)-nx+ii),1).
                compute vg1=varb((ncol(x)-(2*nx)+ii),(ncol(x)-(2*nx)+ii)).
                compute vg3=varb((ncol(x)-nx+ii),(ncol(x)-nx+ii)).
                compute covg1g3=varb((ncol(x)-(2*nx)+ii),(ncol(x)-nx+ii)).
                compute x2 = probeval(jj,1).
                compute w1 = g1+g3*x2.
                compute varw1 = vg1+(2*x2*covg1g3)+((x2*x2)*vg3).
                compute sew1 = sqrt(varw1).
                compute t1 = w1/sew1.
                do if (ovals <> 2).
                  compute LLCI = (w1-sqrt(tval)&*sew1).
                  compute ULCI = (w1+sqrt(tval)&*sew1).
                  compute p = 2*(1-tcdf(abs(t1), dfres)).
                  compute cnms = {"Coeff", "se", "t", "p", "LLCI", "ULCI"}.
                  compute matt(ii,:)={w1,sew1,t1,p,llci,ulci}.
                else if (ovals = 2).
                  compute LLCI = (w1-sqrt(tval)&*sew1).
                  compute ULCI = (w1+sqrt(tval)&*sew1).
                  compute p = 2*(1-cdfnorm(abs(t1))).        
                  compute cnms = {"Coeff", "se", "Z", "p", "LLCI", "ULCI"}.
                  compute matt(ii,:)={w1,sew1,t1,p,llci,ulci}.
                end if.
              end loop.
              compute rnms=t(xname).
              compute mdvalpr=probeval(jj,1).
              print mdvalpr/title = "Moderator value:"/rnames=rnn2/format=!decimals/space=0.
              print matt/title=" "/cnames=cnms/rnames=rnms/format=!decimals/space=0.
              compute xprob=x(:,(ncol(x)-(2*nx)))-mdvalpr.
              loop kk = 1 to nx.
                compute xprob={xprob, (xprob(:,1)&*x(:,(ncol(x)-(2*nx)+kk)))}.          
              end loop.
              compute xprob={x(:,1:((ncol(x)-(2*nx))-1)),xprob}.
              do if (ovals <> 2).
                compute bmultc = inv(t(xprob)*xprob)*t(xprob)*y.
                compute residc = y-(xprob*bmultc).
                compute ssresidc=csum(residc&*residc).
                compute r2c = r2-(1-(ssresidc/sstotal)).
                compute fcha2=(dfres*r2c)/(nx*(1-r2)).
                do if (hc3 = 1).
                  loop kk = 1 to nx.
                    compute xprob={xprob, x(:,(ncol(x)-(2*nx)+kk))}.          
                  end loop.
                  compute bmultc = inv(t(xprob)*xprob)*t(xprob)*y.
                  compute k3 = nrow(bmultc).
                  compute xhc=xprob.
                  compute h = xhc(:,1).
                  loop i3=1 to n.
                    compute h(i3,1)= xhc(i3,:)*inv(t(xprob)*xprob)*t(xhc(i3,:)). 
                  end loop.
                  loop i3=1 to k3.
                    compute xhc(:,i3) = (resid(:,ncol(resid))&/(1-h))&*xhc(:,i3).
                  end loop.
                  compute varbc=(inv(t(xprob)*xprob)*t(xhc)*xhc*inv(t(xprob)*xprob)).
                  compute lmat=make((nrow(bmultc)-nx),nx,0).
                  compute lmat2=ident(nx).
                  compute lmat={lmat;lmat2}.
                  compute fcha2 = (t(t(lmat)*bmultc)*inv(t(lmat)*varbc*lmat)*((t(lmat)*bmultc)))/nx).
                end if.
                compute pvalr2cc=1-fcdf(fcha2,nx,dfres).
                compute rcha2 = {r2c, fcha2, nx, dfres,pvalr2cc}.
                print rcha2/title = "Test of equality of conditional means at this value of the moderator"/clabels "R2-chng" "F" "df1" "df2" "p"/format !decimals.
              else if (ovals = 2).
                compute btmp=b.
                compute pt1lp = make(n,1,0.5).
                compute ylp=y.
                compute xlp=xprob.
                compute bt1 = make(ncol(xlp),1,0).
                LOGITER pt1lp=pt1lp/xlp=xlp/ylp=ylp/iterate=!iterate/converge=!converge.
                compute LLdiff=LL2-LL4.
                compute b=btmp.
                compute varb=varbtmp.
                compute pchi=1-chicdf(lldiff,nx).
                compute rcha2 = {LLdiff, nx,pchi}.
                print rcha2/title = "Test of equality of log odds conditioned on this moderator value"/clabels "Chi-sq" "df" "p"/format !decimals.
              end if.
            end if.
            /* end of do if mcfoc > 0 */.
            compute ttttt=make(nrow(dummat),1,probeval(jj,1)).
            compute ttttt={ttttt,dummat(:,2:ncol(dummat))}.
            loop kkk=1 to nx.
              compute ttttt={ttttt,ttttt(:,1)&*ttttt(:,1+kkk)}.
            end loop.
            compute ones=make(nrow(dummat),1,1).
            do if (ncol(x) > (2*nx+2)).
              compute covmnmat=make(nrow(ttttt),ncol(covmns),0).
              loop kkk=1 to nrow(ttttt).
                compute covmnmat(kkk,:)=covmns.  
              end loop.
              compute ttttt={covmnmat,ttttt}.
            else.
              compute ttttt={ones,ttttt}.
            end if.
            compute focvals={focvals;dummat(:,1),ttttt}.  
            do if (mcfoc > 0).
              compute yhat={dummat(:,1),(ttttt*b)}.
              compute cmnms={fciv, "yhat"}.
              do if (ovals <> 2).
                print yhat/title = "Estimated conditional means at this value of the moderator"/cnames=cmnms/format !decimals.
              end if.
              do if (ovals = 2).
                print yhat/title = "Estimated conditional log odds at this value of the moderator"/cnames=cmnms/format !decimals.
              end if.
              do if (jj <> nrow(probeval)).
                print/title="-------------"/space=0.
              end if.
            end if.
          end loop.
          compute focvals=focvals(2:nrow(focvals),:).
          compute yhat=focvals(:,2:ncol(focvals))*b.
          compute focvals={focvals(:,1),focvals(:,(ncol(focvals)-(2*nx))),yhat}.
          compute cnms={fciv,mdtr,"yhat"}.
          do if (mcmod > 0).
            compute cnms={mdtr,fciv,"yhat"}.
          end if.
        end if.
        do if (mcmod > 0).
          compute outp=make((nx+1),7,0).
          compute bcatm={b((ncol(x)-(2*nx)),1);b((ncol(x)-nx+1):ncol(x),1)}.
          compute outp(:,1)=dummat(:,1).
          compute bcatcov=varb((ncol(x)-nx):ncol(x),(ncol(x)-nx):ncol(x)).
          compute bcatcov(1,1)=varb((ncol(x)-(2*nx)),(ncol(x)-(2*nx))).
          compute bcatcov(2:nrow(bcatcov),1)=varb((ncol(x)-nx+1):ncol(x),(ncol(x)-2*nx)).
          compute bcatcov(1,2:nrow(bcatcov))=t(varb((ncol(x)-nx+1):ncol(x),(ncol(x)-2*nx))).
          loop i = 1 to nrow(dummat).
            compute catmval={1,dummat(i,2:ncol(dummat))}.
            compute condeff=catmval*bcatm.
            compute condse=sqrt(catmval*bcatcov*t(catmval)).
            compute outp(i,2:3)={condeff,condse}.
          end loop.
          compute outp(:,4)=outp(:,2)&/outp(:,3).
          compute outp(:,5) = 2*(1-tcdf(abs(outp(:,4)), dfres)).
          compute outp(:,6) = (outp(:,2)-sqrt(tval)&*outp(:,3)).
          compute outp(:,7) = (outp(:,2)+sqrt(tval)&*outp(:,3)).
          do if (ovals <> 2).
            compute cnmms = {xname2(1,1), "Coeff", "se", "t", "p", "LLCI", "ULCI"}.
          end if.
          do if (ovals = 2).
            compute outp(:,5) = 2*(1-cdfnorm(abs(outp(:,4)))).
            compute cnmms = {xname2(1,1), "Coeff", "se", "Z", "p", "LLCI", "ULCI"}.
          end if.
          print outp/title = "Conditional Effect of Focal Predictor in Groups Defined by the Moderator Variable:"/cnames=cnmms/format = !decimals.
        end if.
        do if (mcfoc = 0 and mcmod = 0).
          do if (ovals <> 2).
            compute cnms = {nms((ncol(x)-1),1), "Coeff", "se", "t", "p", "LLCI", "ULCI"}.
          else if (ovals = 2).
            compute cnms = {nms((ncol(x)-1),1), "Coeff", "se", "Z", "p", "LLCI", "ULCI"}.
          end if.
          print outp/title = "Conditional Effect of Focal Predictor at Values of the Moderator Variable"/cnames = cnms/format = !decimals.
        end if.
        do if (probeval(1,1) < mdmin).
          compute lowwarn = 1.
        end if.
        do if (probeval(nrow(probeval),1) > mdmax).
          compute highwarn = 1.
        end if.
        do if (nval > 2 and (!modval = 9999) and mcmod = 0).
          do if (ptiles <> 1).
            print/title = "Moderator values are the sample mean and plus/minus one SD from mean".
          else.
            print/title = "Moderator values are 25th, 50th, and 75th percentiles of the moderator distribution".
          end if.
          do if (highwarn = 1 and ptiles = 0).
            print/title = "Warning: One SD above the mean is beyond the available data".
          end if.
          do if (lowwarn = 1 and ptiles = 0).
            print/title = "Warning: One SD below the mean is beyond the available data".
          end if.
        end if.
        do if (nval = 2 and (!modval = 9999)).
          print/title = "The moderator variable is dichotomous".
        end if.
      end if.
      do if (jn > 0 and nval > 2 and (mcmod=0 and mcfoc=0)).
        compute ajn =(bb*vg3)-(g3*g3).
        compute bjn = 2*((bb*covg1g3)-(g1*g3)).
        compute cjn = (bb*vg1)-(g1*g1).
        compute radarg = (bjn*bjn)-(4*ajn*cjn).
        compute den = 2*ajn.
        compute nrts = 0.
        do if (radarg >= 0 and den <> 0 and nval <> 2).
          compute x21 = (-bjn+sqrt(radarg))/den.
          compute x22 = (-bjn-sqrt(radarg))/den.
          compute roots = 0.
          do if (x21 >= mdmin and x21 <= mdmax).
            compute nrts = 1.
            compute roots = {roots; x21}.
          end if.
          do if (x22 >= mdmin and x22 <= mdmax).
            compute nrts = nrts + 1.
            compute roots = {roots; x22}.
          end if.
          do if (nrts > 0).
            compute roots = roots(2:nrow(roots),1).
            compute cuts=make(nrow(roots),2,0).
            loop j=1 to nrow(roots).
              compute cuts(j,1)=csum((dd(:,ncol(dd)) < roots(j,1)))/(.01*n).
              compute cuts(j,2)=csum((dd(:,ncol(dd)) > roots(j,1)))/(.01*n).
            end loop.
            compute roots={roots,cuts}.
            do if (jn = 1).
              print roots/title = "Moderator Value(s) Defining Nonsimultaneous Johnson-Neyman Significance Region(s)"/
                 clabels "Value" "% below" "% above"/format F10.4.
            else if (jn = 2).
              print roots/title = "Moderator Value(s) Defining Simultaneous Johnson-Neyman Significance Region(s)"/
                  clabels "Value" "% below" "% above"/format F10.4.
            end if.
          end if.
          do if (nrts = 0).
            print/title = "There are no regions of significance for the focal predictor within the observed range of the moderator".
          end if.
        end if.

/* Here is the JN stuff */.
        compute probeval = 0.
        loop j = 0 to 20.
          compute temp = (mdmin+j*((mdmax-mdmin)/20)).
          compute probeval = {probeval; temp}.
        end loop.
        compute probeval = {probeval; (mdmax+5)}.
        do if (nrts > 0).
          loop i = 1 to nrts.
            loop j = 1 to (nrow(probeval)-1).
              do if ((roots(i,1) > probeval(j,1)) and (roots(i,1) < probeval((j+1),1))).
                compute probeva2 = {probeval(1:j,1); roots(i,1); probeval((j+1):nrow(probeval),1)}.
              end if.
            end loop.
            compute probeval = probeva2.
          end loop.
        end if.
        compute probeval = probeval(2:(nrow(probeval)-1),1).
        compute outp = make(nrow(probeval),7,0).
        loop i = 1 to nrow(probeval).
          compute x2 = probeval(i,1).
          compute w1 = g1+g3*x2.
          compute varw1 = vg1+(2*x2*covg1g3)+((x2*x2)*vg3).
          compute sew1 = sqrt(varw1).
          compute t1 = w1/sew1.
          compute LLCI = (w1-sqrt(tval)&*sew1).
          compute ULCI = (w1+sqrt(tval)&*sew1).
          do if (ovals <> 2).
            compute p = 2*(1-tcdf(abs(t1), dfres)).
            compute cnms = {nms((ncol(x)-1),1), "Coeff", "se", "t", "p", "LLCI", "ULCI"}.
            compute outp(i,:) = {x2, w1, sew1, t1, p, LLCI, ULCI}.
          else if (ovals = 2).
            compute p = 2*(1-cdfnorm(abs(t1))).
            compute cnms = {nms((ncol(x)-1),1), "Coeff", "se", "Z", "p", "LLCI", "ULCI"}.
            compute outp(i,:) = {x2, w1, sew1, t1, p, LLCI, ULCI}.
          end if.
        end loop.
        do if (jn = 2).
          compute outp = {outp(:,1), outp(:,6:7)}.
          compute cnms = {nms((ncol(x)-1),1), "LLCI", "ULCI"}.
        end if.
        print outp/title = "Conditional Effect of Focal Predictor at Values of Moderator Variable"/cnames = cnms/format = !decimals.
        print cilm/title = "Alpha level used for Johnson-Neyman method and confidence intervals:"/format = F4.2.
      end if.
      do if (!modval <> 9999 and mcmod = 0 and ((!modval < mdmin) or (!modval > mdmax))).
        print/title = "Warning: MODVAL is outside of the range of the data".
      end if.
      do if (mcfoc = 0 and mcmod=0).
        do if (jn < 1).
          compute fvdes = ncol(design(fvar)).
          do if (fvdes = 2).
            compute fv1 = cmin(fvar).
            compute fv2 = cmax(fvar).
            compute r = 1.
            loop j = 1 to nrow(focvals).
              do if ((focvals(j,1) = fv1) or (focvals(j,1) = fv2)).
                compute focvals(r,:)=focvals(j,:).
                compute r = r + 1.
              end if.
            end loop.
            compute focvals = focvals(1:(r-1),:).
          end if.
        end if.
      end if.
      do if (!est = 1 and jn < 1).
        do if (mcfoc = 0 and mcmod=0).
          compute cnms = {t(nms((ncol(x)-2):(ncol(x)-1),1)), "yhat"}.
        end if.
        do if (ovals = 2).
          compute prob = exp(focvals(:,3))&/(1+exp(focvals(:,3))).
          compute focvals = {focvals, prob}.
          compute cnms = {cnms, "prob"}.
        end if.
        print/title = "***************************************************************************".
        print focvals/title = "Data for Visualizing Conditional Effect of Focal Predictor"/cnames = cnms/format = !decimals.
        do if (savplot= 1).
          save focvals/outfile = */names = cnms.
        end if.
        do if (ncovs > 0).
          print/title = "NOTE: For data above, covariates are set to their sample means.".
        end if.
      end if.
    end if.
    print/title = "********************* ANALYSIS NOTES AND WARNINGS *************************".
    print cilm/title = "Alpha level used for confidence intervals:"/format = F4.2.
    do if (alperr = 1 and jn > 0).
      print/title = "ERROR: Inappropriate alpha level requested.  Alpha set to 0.05.".
    end if.
    do if (jnerr = 1).
      print/title = "NOTE: Simultaneous inference not available for logistic regression.  Nonsimultaneous results are printed.".
    end if.
    do if (hc3 = 1).
      print/title = "NOTE: The HC3 standard error estimator was used.".
    end if.
  end if.
  do if (itprob = 1).
    print/title = "ERROR: There was a problem during iteration.".
  end if.
  do if (itprob = 2).
    print/title = "ERROR: The convergence criterion was not attained.".
  end if.
  do if (ncol(centerv) > 1).
    compute centerv=centerv(:,2:ncol(centerv)).
    print centerv/title = "NOTE: The following variables were mean centered prior to analysis:"/format = A8.
  end if.
end if. 
print.
loop i= 1 to errs.
  do if errsm(i,1)=1.
    print/title = "ERROR: Only one variable can be specifed as the outcome Y."/space=0.
  end if.
  do if errsm(i,1)=2.
    print/title = "ERROR: You must specify at least two variables in the X= list."/space=0.
  end if.
  do if errsm(i,1)=3.
    print/title = "ERROR: Focal predictor and moderator cannot both be specified as multicategorical."/space=0.
  end if.
  do if errsm(i,1)=4.
    print/title = "ERROR: Categorical variables cannot have more than 10 categories."/space=0.
  end if.
  do if errsm(i,1)=5.
    print/title = "ERROR: One of the variables in the model exhibits no variation."/space=0.
  end if.
  do if errsm(i,1)=6.
    print/title = "ERROR: Each group must have at least two cases."/space=0.
  end if.
end loop.
do if (!modval <> 9999 and mcmod <> 0).
  print/title = "NOTE: MODVAL option not available for use with a multicategorical moderator".
end if.
end matrix.
set printback = on.
RESTORE.
!ENDDEFINE.




