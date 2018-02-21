/* MEDIATE for SPSS */.
/* Written by Andrew F. Hayes */.
/* The Ohio State University */.
/* http://www.afhayes.com */.
/* Version 200715 */.
set printback=off.
define bcboot2 (databcbt = !charend ('/')/estmte = !charend ('/') !default(9999)).
compute temp = !databcbt.
compute temp(GRADE(!databcbt)) = !databcbt.
compute badlo = 0.
compute badhi = 0.
do if (!estmte <> 9999).
  compute pv=csum(temp < !estmte)/reps.
  compute ppv = pv. 
  do if (pv > .5).
    compute ppv = 1-pv.
  end if.
  compute y5=sqrt(-2*ln(ppv)).
  compute xp=y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0).
  do if (pv <= .5).
    compute xp = -xp.
  end if.
  compute cilow=trunc(reps*(cdfnorm(2*xp+xp2))).
  compute cihigh=trunc(reps*(cdfnorm(2*xp+(-xp2))))+1.
  do if (cilow < 1).
    compute cilow = 1.
     compute booterr=1.
     compute badlo = 1.
  end if.
  do if (cihigh > reps).
    compute cihigh = reps.
    compute booterr=1.
     compute badhi = 1.
  end if.
  compute llcit=temp(cilow,1).
  compute ulcit=temp(cihigh,1).
end if.
do if (!estmte = 9999).
   compute llcit=temp(cilow,1).
   compute ulcit=temp(cihigh,1).
end if.
!enddefine.

DEFINE MEDIATE (y = !charend('/')/x = !charend('/')/m = !charend('/')/c=!charend('/') !default(xxxxx)/samples =!charend('/') !default(1000)
/ciconf = !charend('/') !default(95)/cimethod = !charend('/') !default(1)/omnibus = !charend('/') !default(0)/total = !charend('/') !default(0)
/catx = !charend('/') !default(0)/seed = !charend('/') !default(random)/percent = !charend('/') !default(0)).
SET PRINTBACK = OFF.
SET LENGTH = NONE.
SET MXLOOPS = 10000001.
SET SEED = !seed.
SET MXCELLS = 100000000.
MATRIX.
get dat/file = */variables = !y !x !m/names = vnames/missing = 9999.
compute ninit = nrow(dat). 
get dd/variables = !y !x !m/names = name/MISSING = OMIT.
compute method = !cimethod.
compute bconoff=(!percent <> 1).
compute percent = (!percent <> 0).
compute samples = !samples.
compute omnibus = !omnibus.
compute booterr = 0.
compute stratify = 0.
compute conf = !ciconf.
compute catx = trunc(!catx).
do if (catx < 0 or catx > 4).
   compute catx = 1.
end if.
compute ovals = ncol(design(dd(:,1))).
compute nte = 1.
compute criterr = 0.
compute note = make(50,1,0).
compute adjust = 0.
compute xskip = 0.
compute dumok = 0.
compute badend = 0.
compute p0=-.322232431088.
compute p1 = -1.
compute p2 = -.342242088547.
compute p3 = -.0204231210245.
compute p4 = -.0000453642210148.
compute q0 = .0993484626060.
compute q1 = .588581570495.
compute q2 = .531103462366.
compute q3 = .103537752850.
compute q4 = .0038560700634.
compute alpha2 = (1-(conf/100))/2.
compute y5=sqrt(-2*ln(alpha2)).
compute xp2=-(y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0)).
compute priorlo = -9999999.
compute priorhi = 9999999.
compute ddd = {"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"}.
do if (trunc(conf) ge 100 or (trunc(conf) le 50)).
   compute conf = 95.
end if.
compute total = !total.
compute reps = samples.
compute errs = 0.
compute runerrs = make(50,1,0).
do if (ovals = 2).
  compute errs= errs+1.
  compute runerrs(errs,1)=6.
  compute criterr = 1.
end if.
do if (method = 2).
  do if (omnibus = 1).
    compute note(nte,1) = 5.
    compute nte=nte+1.
  end if.
  compute samples = 0.
  do if (reps < 1000).
     compute reps = 1000.
  end if.
end if.
loop.
  compute cilow = trunc(reps*(1-(conf/100))/2).
  compute cihigh = trunc((reps*(conf/100)+(reps*(1-(conf/100))/2)))+1.
  do if (cilow < 1 or cihigh > reps).
     compute reps=trunc((reps+1000)/1000)*1000.
   print reps.
     compute adjust = 1.
  end if.
end loop if (cilow gt 0 and cihigh le reps).
do if (method = 1).
   compute samples  = reps.
end if.
do if (adjust = 1).
  compute note(nte,1)=1.
  compute nte = nte+1.
end if.
compute temp = ncol(dd).
compute cc = !quote(!c).
do if (cc = "xxxxx").
  compute nc = 0.  
  else.
    get dd/variables = !y !x !m !c/names = name/MISSING = OMIT.
    compute nc = ncol(dd)-temp.
    compute cname = name(1,(ncol(name)-nc+1):ncol(name)).
    get dat/file = */variables = !y !x !m/names = vnames/missing = 9999.
    compute ninit = nrow(dat). 
end if.
compute n = nrow(dd).
do if (n < ninit).
  compute nmiss = ninit-n.
  compute note(nte,1) = 8.
  compute nte = nte + 1.
end if.
compute ones = make(n,1,1).
get temp/variables = !x/MISSING = OMIT.
compute nxd = ncol(temp).
compute nx = nxd.
do if (nx > 1 and catx > 0).
  compute errs= errs+1.
  compute runerrs(errs,1)=2.
  compute criterr = 1.
end if.
do if (nx = 1 and catx > 0).
  do if (method = 3).
     compute stratify = 1.
  end if.
  compute temp = dd.
  compute temp(GRADE(dd(:,2)),:) = dd.
  compute dd = temp.
  compute dummy = design(dd(:,2)).
  compute nvls = ncol(dummy).
  compute nnvls = csum(dummy).
  compute mnvls = cmin(t(nnvls)).
  do if (mnvls < 2).
    compute errs = errs + 1.
    compute runerrs(errs,1) = 3.
    compute criterr = 1.
  end if.
  do if (nvls > 9).
    compute errs = errs+1.
    compute runerrs(errs,1) = 4.
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
    do if (catx > 0).
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
      do if (stratify <> 0).
        compute note(nte,1) = 6.
        compute nte=nte+1.
      end if.
      do if (catx = 1).
        compute note(nte,1)=2.
        compute nte = nte+1.
      end if.
      do if (catx = 2).
        compute note(nte,1)=3.
        compute nte = nte+1.
        loop k = 1 to n.
          do if (rsum(x(k,:)) = 0).
            compute x(k,:) = minus1.
          end if.
        end loop.
      end if.
      do if (catx = 3 or catx = 4).
         compute note(nte,1)=4.
         do if (catx=4).
           compute note(nte,1)=7.
         end if.
         compute nte = nte+1.
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
         do if (catx = 4).
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
      compute xdata = x.
      compute xskip = 1.
      compute xname = ddd(1,1:nx).
      compute indlbs = t(xname).
      compute dummat = make((nx+1),nx,0).
      compute dummat((2:nrow(dummat)),:)=ident(nx).
      do if (catx = 2).
         compute dummat(1,:) = minus1.
      end if.
      do if (catx = 3).
         loop i = 2 to nrow(dummat).
            loop j = 1 to (i-1).
               compute dummat(i,j) = 1.
            end loop.
          end loop.
      end if.
      do if (catx = 4).
        compute dummat=conmat1.
      end if.
      compute dummat={nnvls, dummat}.
      compute xname2={name(1,2), xname}.   
    end if.
  end if.
end if.
get temp/variables = !m/MISSING = OMIT.
compute nm = ncol(temp).
do if (nm > 15).
    compute errs = errs+1.
    compute runerrs(errs,1) = 5.
    compute criterr = 1.
end if.
compute r2b = make(nm,1,-999).
compute indeff = make((nx*nm),1,-999).
compute cov = make((nx*nm+nm),(nx*nm+nm),0).
compute mns=make(nrow(cov),1,0).
compute ny = ncol(dd)-nm-nxd-nc.
compute y = dd(:,1).
compute ydata = y.
compute yname = name(1,1).
do if (xskip = 0).
  compute x = dd(:,2:(1+nxd)).
  compute xdata = x.
  compute xname = name(1,2:(1+nxd)).
  compute indlbs = t(xname).
end if.

compute m = dd(:,(2+nxd):(1+nxd+nm)).
compute mdata = m.
compute mcheck = 0.
loop i = 1 to nm.
  compute dichm = (ncol(design(m(:,i)))).
  do if (dichm = 2 and mcheck = 0).
    compute errs= errs+1.
    compute runerrs(errs,1) = 1.
    compute mcheck = 1.
  end if.
end loop.
compute mlab = {"   M1 = "; "   M2 ="; "   M3 ="; "   M4 ="; "   M5 ="; "   M6 ="; "   M7 = "; "   M8 = "; "   M9 = "; "  M10 = "}.
compute mlab = {mlab; " M11 = "; " M12 = "; " M13 = "; " M14 = "; " M15 = "}.
compute mname = name(1,(2+nxd):(1+nxd+nm)).
do if (nc > 0).
   compute c = dd(:,(2+nxd+nm):(1+nxd+nm+nc)).
   compute x = {x, c}.
   compute xdata = x.
   compute cdata = c.
   compute xname = {xname, cname}.
end if.
compute nms={yname; t(mname); t(xname(1,1:(ncol(xname)-nc)))}.
do if (dumok = 1).
  compute nms={yname; t(mname); name(1,2)}.
end if.
do if (nm < 16).
  print/title = "*************** MEDIATE Procedure for SPSS Release 050213 ****************".
  print/title = "        Written by Andrew F. Hayes, Ph.D.   http://www.afhayes.com".
  print /title = "**************************************************************************".
  compute cnms={"    Y = "; mlab(1:ncol(mname),1); "    X = "}.
  print nms/title = "VARIABLES IN THE FULL MODEL:"/format = A8/rnames = cnms.
  do if (nc > 0).
    print cname/title = "COVARIATES:"/format = a8.
  end if.
  do if (dumok = 1).
    print dummat/title = "CODING OF CATEGORICAL X FOR ANALYSIS:"/cnames = xname2/format = F7.4.
  end if.
end if.
compute boots=make((samples+1),((nx+1)*nm),-999).
do if (errs = 0 and criterr = 0).
  loop k = 1 to (samples+1).
    do if (k > 1).
      compute u=trunc(uniform(n,1)*n)+1.
      do if (catx > 0 and stratify <> 0).
         loop i = 1 to nrow(xdes).
            compute u((xdes(i,2):xdes(i,3)),1)=trunc(uniform(xdes(i,4),1)*xdes(i,4))+xdes(i,2).
         end loop.
      end if.
      compute xdata = x(u,:).
      compute mdata = m(u,:).
      compute ydata = y(u,:).
      do if (nc > 0).
       compute cdata = c(u,:).
      end if.
    end if.

    /* estimate coefficients and covariance matrices */.

    do if (k = 1 and total = 1).
      print/title = "**************************************************************************".
      print yname/title = "OUTCOME VARIABLE:"/format = A10.
      compute xmat = {ones, xdata}.
      compute ymat = ydata.
      compute cmat = inv(t(xmat)*xmat)*t(xmat)*ymat.
      compute sstotal = t(ymat-(csum(ymat)/n))*(ymat-(csum(ymat)/n)).
      compute ssresid = csum((ymat-xmat*cmat)&**2). 
      compute r2 = (sstotal-ssresid)/sstotal.
      compute adjr2 = 1-((1-r2)*(n-1)/(n-ncol(xmat))).
      compute fratio = ((n-ncol(xmat))*r2)/((1-r2)*(ncol(xmat)-1)).
      compute pfr = 1-fcdf(fratio,(ncol(xmat)-1),(n-ncol(xmat))).
      compute op = {sqrt(r2), r2, adjr2, fratio, (ncol(xmat)-1), (n-ncol(xmat)), pfr}. 
      print op/title = "MODEL SUMMARY (TOTAL  EFFECTS MODEL)"/format = F10.4/clabels = "R", "R-sq", "Adj R-sq", "F", "df1", "df2", "p".
      compute covmat = (ssresid/(n-ncol(xmat)))*inv(t(xmat)*xmat).
      compute tratio = cmat&/sqrt(diag(covmat)).
      compute op = {cmat, sqrt(diag(covmat)), tratio, (2*(1-tcdf(abs(tratio), n-ncol(xmat))))}.
      compute rnms = {"Constant"; t(xname)}.
      print op/title = "MODEL COEFFICIENTS (TOTAL EFFECTS MODEL)"/clabels = "Coeff.", "s.e.", "t", "p"/rnames = rnms/format = F10.4.
      compute ssresid2 = ssresid.
      compute r2total = r2.
      do if (nc > 0).
        compute xmat = {ones, cdata}.
        compute ctmat = inv(t(xmat)*xmat)*t(xmat)*ymat.
        compute ssresid2 = csum((ymat-xmat*ctmat)&**2).
        compute r2total = r2-((sstotal-ssresid2)/sstotal).
      end if.
      do if (omnibus = 1).
        compute fratio = ((n-1-nx-nc)*r2total)/((1-r2)*nx). 
        compute pfr = 1-fcdf(fratio,(nx),(n-1-nx-nc)).
        compute r2direct = {r2total, fratio, nx, (n-1-nx-nc), pfr}.
        print r2direct/title = "OMNIBUS TEST OF TOTAL EFFECT"/clabels = "R-sq", "F", "df1", "df2", "p"/format = F10.4.
      end if.
    end if.


    loop i = 1 to nm.
      compute start = 1+(i-1)*nx.
      compute xmat = {ones, xdata}.
      compute ymat = mdata(:,i).
      compute amat = inv(t(xmat)*xmat)*t(xmat)*ymat.
      compute sstotal = t(ymat-(csum(ymat)/n))*(ymat-(csum(ymat)/n)).
      compute ssresid = csum((ymat-xmat*amat)&**2). 
      compute r2 = (sstotal-ssresid)/sstotal.
      compute adjr2 = 1-((1-r2)*(n-1)/(n-ncol(xmat))).
      compute r2cha = r2.
      compute adjr2cha = adjr2.
      do if (nc > 0).
        compute xmat = {ones, cdata}.
        compute atmat = inv(t(xmat)*xmat)*t(xmat)*ymat.
        compute ssresid2 = csum((ymat-xmat*atmat)&**2).
        compute r2step1 = (sstotal-ssresid2)/sstotal.
        compute r2cha = r2-r2step1.
        compute adjr2cha = adjr2- (1-((1-r2step1)*(n-1)/(n-ncol(xmat)))).
        compute xmat = {ones, xdata}.
      end if.
      compute r2b(i,1)=adjr2cha.
      /* compute r2b(i,1)=r2cha */.
      /* this below is pertinent to printing obtained output */.
      do if (k = 1).
        print/title = "**************************************************************************".
        compute covmat = (ssresid/(n-ncol(xmat)))*inv(t(xmat)*xmat).
        compute tratio = amat&/sqrt(diag(covmat)).
        compute rnms = mname(1,i).
        print rnms/title = "OUTCOME VARIABLE:"/format = A10.
        compute fratio = ((n-ncol(xmat))*r2)/((1-r2)*(ncol(xmat)-1)).
        compute pfr = 1-fcdf(fratio,(ncol(xmat)-1),(n-ncol(xmat))).
        compute op = {sqrt(r2), r2, adjr2, fratio, (ncol(xmat)-1), (n-ncol(xmat)), pfr}. 
        print op/title = "MODEL SUMMARY"/format = F10.4/clabels = "R", "R-sq", "Adj R-sq", "F", "df1", "df2", "p".
        compute op = {amat, sqrt(diag(covmat)), tratio, (2*(1-tcdf(abs(tratio), n-ncol(xmat))))}.
        compute rnms = {"Constant"; t(xname)}.
        print op/title = "MODEL COEFFICIENTS"/clabels = "Coeff.", "s.e.", "t", "p"/rnames = rnms/format = F10.4.
        compute cov(start:(start+(nx-1)),start:(start+(nx-1)))=covmat(2:(ncol(covmat)-nc),2:(ncol(covmat)-nc)).
      end if.
      compute mns(start:(start+(nx-1)),1)=amat(2:(1+nx)).
    end loop.

    compute xmat = {ones,mdata,xdata}.
    compute ymat = ydata.
    compute bmat =  inv(t(xmat)*xmat)*t(xmat)*ymat.
    compute sstotal = t(ymat-(csum(ymat)/n))*(ymat-(csum(ymat)/n)).
    compute ssresid = csum((ymat-xmat*bmat)&**2).

    do if (k = 1).
      print/title = "**************************************************************************".
      compute covmat = (ssresid/(n-ncol(xmat)))*inv(t(xmat)*xmat).
      compute r2 = (sstotal-ssresid)/sstotal.
      compute adjr2 = 1-((1-r2)*(n-1)/(n-ncol(xmat))).
      print yname/title = "OUTCOME VARIABLE:"/format = A10.
      compute fratio = ((n-ncol(xmat))*r2)/((1-r2)*(ncol(xmat)-1)).
      compute pfr = 1-fcdf(fratio,(ncol(xmat)-1),(n-ncol(xmat))).
      compute df2=(n-ncol(xmat)).
      compute op = {sqrt(r2), r2, adjr2, fratio, (ncol(xmat)-1), df2, pfr}. 
      print op/title = "MODEL SUMMARY"/format = F10.4/clabels = "R", "R-sq", "adj R-sq", "F", "df1", "df2", "p".
      compute tratio = bmat&/sqrt(diag(covmat)).
      compute rnms = {"Constant"; t(mname); t(xname)}.
      compute op = {bmat, sqrt(diag(covmat)), tratio, (2*(1-tcdf(abs(tratio), n-ncol(xmat)))) }.
      print op/title = "MODEL COEFFICIENTS"/clabels = "Coeff.", "s.e.", "t", "p"/rnames = rnms/format = F10.4.
      compute cov((ncol(cov)-nm+1):ncol(cov),(ncol(cov)-nm+1):ncol(cov))=covmat(2:(1+nm),2:(1+nm)).
      
      /* Here we test for homogeneity of regression */.
      compute homgen=make((nm+1),5,0).
      compute xinttemp = make(n,(nx*nm),0).
      compute temp=1.
      loop v1=1 to nm.
        loop v2=1 to nx.
          compute xinttemp(:,temp)=xmat(:,(1+v1))&*xmat(:,(nm+1+v2)).
          compute temp=temp+1.
        end loop.  
      end loop.
      loop v1=1 to (nm+1).
        do if (v1 = (nm+1)).
          compute temp={xmat,xinttemp}.
          compute df1temp=ncol(xinttemp).
        else.
          compute temp={xmat,xinttemp(:,(1+((v1-1)*nx)):((1+((v1-1)*nx))+(nx-1)))}.
          compute df1temp=nx.
        end if.
        compute btemp =  inv(t(temp)*temp)*t(temp)*ymat.
        compute restemp = csum((ymat-temp*btemp)&**2).
        compute r2temp=(sstotal-restemp)/sstotal.
        compute homgen(v1,1)=r2temp-r2.
        compute homgen(v1,3)=df1temp.
        compute homgen(v1,4)=df2-df1temp.
        compute homgen(v1,2)=(homgen(v1,4)*homgen(v1,1))/((1-r2temp)*homgen(v1,3)).
        compute homgen(v1,5)=1-fcdf(homgen(v1,2),homgen(v1,3),homgen(v1,4)).
      end loop.
     do if (nm = 1).
       compute mnamet = mname(1,1).
       compute homgen=homgen(1,:).
     end if.
     do if (nm > 1).
       compute mnamet={mname,"OMNIBUS"}.
     end if.
     print homgen/title = "TEST OF HOMOGENEITY OF REGRESSION (X*M INTERACTION)"/rnames = mnamet/clabels="R-sq", "F", "df1", "df2", "p"/format = F10.4.
     /* end of homogeneity test */.
  
      compute xmat = {ones,mdata}.
      do if (nc > 0).
        compute xmat = {xmat, cdata}.
      end if.
      compute btmat = inv(t(xmat)*xmat)*t(xmat)*ymat.
      compute ssresid2 = csum((ymat-xmat*btmat)&**2).
      compute r2direct = r2-((sstotal-ssresid2)/sstotal).
      do if (omnibus = 1).
        compute fratio = ((n-1-nx-nm-nc)*r2direct)/((1-r2)*nx).
        compute pfr = 1-fcdf(fratio,(nx),(n-1-nx-nm-nc)).
        compute r2direct = {r2direct, fratio, nx, (n-1-nx-nm-nc), pfr}.
        print r2direct/title = "OMNIBUS TEST OF DIRECT EFFECT"/clabels = "R-sq", "F", "df1", "df2", "p"/format = F10.4.
      end if.
    end if.
    compute mns((nrow(mns)-nm+1):nrow(mns),1)=bmat(2:(1+nm),1).

    /* calculate relative indirect effects */.
    compute tmp = 1.
    loop i = 1 to (nm*nx).
      compute boots(k,i)=mns(i,1)&*mns(((nm*nx)+tmp),1).
      do if ((i/nx)=trunc(i/nx)).
        compute boots(k,((nm*nx)+tmp))=r2b(tmp,1)*mns(((nm*nx)+tmp),1).
        compute tmp = tmp+1.
      end if.
    end loop.
  end loop.

  do if (method = 1 or method = 3).
    compute llci=make(1,ncol(boots),-999).
    compute ulci=make(1,ncol(boots),-999).
    compute indeff = t(boots(1,1:(nm*nx))).
    compute indeff22=t(boots(1,:)).
    compute omni = t(boots(1,((nm*nx)+1) : ncol(boots))).
    compute boots = boots(2:nrow(boots),:).
    compute bootse = t(sqrt(((reps*cssq(boots))-(csum(boots)&**2))/((reps-1)*reps))).
    loop i = 1 to ncol(boots).
      bcboot2 databcbt = boots(:,i)/estmte=(indeff22(i,1)*bconoff)+(9999*(1-bconoff)).
      compute llci(1,i)=llcit.
      compute ulci(1,i)=ulcit.
      do if (badlo = 1 and llcit <> priorlo).
        compute badend={badend, llcit}.
        compute priorlo = llcit.
      end if.
      do if (badhi = 1 and ulcit <> priorhi).
        compute badend={badend, ulcit}.
        compute priorhi = ulcit.
      end if.
    end loop.
    compute llci = t(llci).
    compute ulci= t(ulci).
    compute llci={llci,ulci}.
  end if.

  /* generate monte carlo samples */.
  do if (method = 2).
    compute tmp=1.
    loop i = 1 to (nm*nx).
      compute indeff(i,1)=mns(i,1)&*mns(((nm*nx)+tmp),1).
      do if ((i/nx)=trunc(i/nx)).
        compute tmp = tmp+1.
      end if.
    end loop.
    compute boots = sqrt(-2*ln(uniform(reps,nrow(cov))))&*cos((2*3.14159265358979)*uniform(reps,nrow(cov))).
    compute boots=boots*chol(cov).
    loop i = 1 to ncol(boots).
      compute boots(:,i)=boots(:,i)+mns(i,1).
    end loop.
    compute tmp=1.
    loop i = 1 to (nm*nx).
      compute boots(:,i)=boots(:,i)&*boots(:,((nm*nx)+tmp)).
      compute temp = boots(:,i).
      compute temp(GRADE(boots(:,i))) = boots(:,i).
      compute boots(:,i) = temp.
      do if ((i/nx)=trunc(i/nx)).
        compute tmp = tmp+1.
      end if.
    end loop.
    compute boots = boots(:,1:(nm*nx)).
    compute bootse = t(sqrt(((reps*cssq(boots))-(csum(boots)&**2))/((reps-1)*reps))).
    compute llci = t({(boots(cilow,:));(boots(cihigh,:))}).
  end if.



  print/title = "**************************************************************************".
  do if (omnibus = 1 and (method = 1 or method = 3)).
    compute indlbs = {indlbs; "OMNIBUS"}.
  end if.
  loop i = 1 to nm.
    compute op = mname(1,i).
    print op/title = "INDIRECT EFFECT(S) THROUGH:"/format = A8.
    do if (method = 2).
      compute clbs = {"Effect", "SE(mc)", "LLCI", "ULCI"}. 
      compute op = {indeff, bootse(1:(nm*nx)), llci}.
    end if.
    do if (method = 1 or method = 3).
      compute clbs = {"Effect", "SE(boot)", "LLCI", "ULCI"}.
      compute op = {indeff, bootse(1:(nm*nx),:), llci(1:(nm*nx),:)}.
    end if.
    compute temp = op((((i-1)*nx)+1):(nx*i),:).
    do if ((method = 1 or method = 3) and omnibus = 1).
      compute temp = {temp; omni(i,1), bootse(((nm*nx)+i),:), llci(((nm*nx)+i),:)}.
    end if.
    print temp/title = " "/rnames = indlbs/cnames = clbs/format = F10.4.
    print/title = "----------".
  end loop.
end if.
print/title = "********************** ANALYSIS NOTES AND WARNINGS ***************************".
loop i = 1 to errs.
   do if (runerrs(i,1) = 1).
     print/title = "CRITICAL ERROR: One of your declared mediators is dichotomous.  This procedure can't be used.".
   end if.
   do if (runerrs(i,1) = 2).
     print/title = "CRITICAL ERROR: Use of categorical option assumes only a single X is listed".
   end if.
   do if (runerrs(i,1) = 3).
     print/title = "CRITICAL ERROR: One of your categories has only 1 case".
   end if.
   do if (runerrs(i,1) = 4).
     print/title = "CRITICAL ERROR: A categorical variable must have fewer than 10 categories".
   end if.
   do if (runerrs(i,1) = 5).
     print/title = "CRITICAL ERROR: The procedure allows a maximum of 15 proposed mediators".
   end if.
   do if (runerrs(i,1) = 6).
     print/title = "CRITICAL ERROR: Your outcome variable is dichotomous.  This procedure can't be used".
   end if.
end loop.
loop i = 1 to nte.
  do if (note(i,1) = 1 and criterr = 0).
    print/title = "NOTE: The number of samples was adjusted upward given your desired confidence".
  end if.
  do if (note(i,1) = 2 and criterr = 0).
    print/title = "NOTE: Indicator coding is used for categorical X".
  end if.
  do if (note(i,1) = 3 and criterr = 0).
    print/title = "NOTE: Effect coding is used for categorical X".
  end if.
  do if (note(i,1) = 4 and criterr = 0).
    print/title = "NOTE: Sequential coding is used for categorical X".
  end if.
  do if (note(i,1) = 7 and criterr = 0).
    print/title = "NOTE:Helmert coding is used for categorical X".
  end if.
  do if (note(i,1) = 5 and criterr = 0).
    print/title = "NOTE: Omnibus test for indirect effect is not available with the Monte Carlo option".
  end if.
  do if (note(i,1) = 6 and criterr = 0).
    print/title = "NOTE: Stratified bootstrap sampling was used".
  end if.
  do if (note(i,1) = 8) and criterr = 0).
    print nmiss/title = "NOTE: Some cases were deleted due to missing data.  The number of such cases was:".
  end if.
end loop.
do if (errs = 0 and criterr = 0).
  print reps/title = "Number of samples used for indirect effect confidence intervals:"/format = F7.0.
  print conf/title = "Level of confidence for confidence intervals:"/format = F8.4.
  do if ((method = 1 or method = 3) and percent = 1).
    print/title = "Percentile bootstrap confidence intervals for indirect effects are printed in output".
  end if.
  do if ((method = 1 or method = 3) and percent = 0).
    print/title = "Bias corrected bootstrap confidence intervals for indirect effects are printed in output".
  end if.
  do if (method = 2).
    print/title = "Monte Carlo confidence intervals for indirect effects are printed in output".
  end if.
  do if (booterr = 1).
    compute badend = badend(1,2:ncol(badend)).
    print badend/title = "WARNING: Bootstrap CI endpoints below not trustworthy.  Decrease confidence or increase bootstraps"/format = F10.4.
  end if.
end if.
END MATRIX.
SET PRINTBACK = ON.
!ENDDEFINE.


